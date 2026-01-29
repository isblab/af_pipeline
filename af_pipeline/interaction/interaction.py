import os
import warnings
import numpy as np
import pandas as pd
from typing import Dict
from itertools import product
from scipy.spatial import distance_matrix
from af_pipeline.utils.misc_utils import (
    get_key_from_res_range,
    generate_cmap,
)
from af_pipeline._initialize import _Initialize
from af_pipeline.tools.matrix_patches import MatrixPatches

class Interaction(_Initialize):
    """ Class to handle interaction data for the predicted structure. \n
    One can obtain:
    1. Interaction map: A binary contact map or distance map.
    2. Restraints: Contacts as pairwise residues in a DataFrame format.
    3. Interacting patches: contiguous regions in the interaction map obtained in (1).
    """

    def __init__(
        self,
        data_file_path: str,
        structure_file_path: str,
        af_offset: dict = {},
        idr_chains: list = [],
        rep_atom_dict: dict = {},
        average_token_pae: bool = False,
        average_token_plddt: bool = False,
        metric_level: str = "per_token",
    ):

        super().__init__(
            data_file_path=data_file_path,
            structure_file_path=structure_file_path,
            af_offset=af_offset,
            rep_atom_dict=rep_atom_dict,
            average_token_pae=average_token_pae,
            average_token_plddt=average_token_plddt,
            metric_level=metric_level,
        )

        self.contact_threshold = 8.0  # Distance threshold in (Angstorm) to define a contact between residue pairs.
        self.plddt_cutoff = 70.0  # pLDDT cutoff to consider a confident prediction.
        self.idr_plddt_cutoff = 50.0  # pLDDT cutoff for IDR chains.
        self.pae_cutoff = 5.0 # PAE cutoff to consider a confident prediction.
        self.idr_chains = idr_chains # List of chains that are disordered

        self.save_plot = False
        self.save_table = False

    @staticmethod
    def get_contact_map(
        coords1: np.ndarray,
        coords2: np.ndarray,
        contact_threshold: float = 8.0,
    ) -> np.ndarray:
        """ Get the contact map between two arrays of coordinates.

        Args:

            coords1 (np.ndarray):
                Coordinates of the chain1 residues.

            coords2 (np.ndarray):
                Coordinates of the chain2 residues.

            contact_threshold (float, optional):
                Threshold for defining a contact between residues in 
                Angstroms. Defaults to 8.0.

        Returns:

            contact_map (np.ndarray):
                Binary contact map where 1 indicates a contact and 0 indicates
                no contact.
        """

        distance_map = distance_matrix(coords1, coords2)

        contact_map = np.where(distance_map < contact_threshold, 1, 0)

        return contact_map

    def create_regions_of_interest(self):
        """
        Create regions of interest for all possible chain pairs.

        Returns:
            regions_of_interest (list): list of regions of interest
        """

        regions_of_interest = []
        token_chain_ids = self.token_chain_ids
        chain_pairs = set()

        for chain1 in set(token_chain_ids):
            for chain2 in set(token_chain_ids):
                if chain1 != chain2:
                    pair = tuple(sorted([chain1, chain2]))
                    chain_pairs.add(pair)

        chain_pairs = list(chain_pairs)

        for chain1, chain2 in chain_pairs:

            ch1_start = self.renumber.renumber_chain_res_num(
                chain_res_num=1,
                chain_id=chain1
            )
            ch1_end = self.renumber.renumber_chain_res_num(
                chain_res_num=self.lengths_dict[chain1],
                chain_id=chain1
            )
            ch2_start = self.renumber.renumber_chain_res_num(
                chain_res_num=1,
                chain_id=chain2
            )
            ch2_end = self.renumber.renumber_chain_res_num(
                chain_res_num=self.lengths_dict[chain2],
                chain_id=chain2
            )

            region_of_interest = {
                chain1: [ch1_start, ch1_end],
                chain2: [ch2_start, ch2_end],
            }

            regions_of_interest.append(region_of_interest)

        return regions_of_interest

    def get_interaction_data(self, region_of_interest: Dict):
        """
        Get the interaction amp, pLDDT, and PAE for the region of interest.

        Args:
            region_of_interest (Dict): Dictionary containing the chain IDs and the residue numbers for the region of interest.

        Returns:
            interaction_map (np.array): binary contact map or distance map
            plddt1 (np.array): plddt values for chain 1
            plddt2 (np.array): plddt values for chain 2
            pae (np.array): PAE matrix for the region of interest
        """

        chain1, chain2 = list(region_of_interest.keys())
        p1_region, p2_region = (
            region_of_interest[chain1],
            region_of_interest[chain2],
        )

        token_rep_chain_ids = self.structure_parser.get_token_chain_ids(
            structure=self.structure,
            rep_atom_dict=self.rep_atom_dict,
            only_representative=True,
        )
        token_rep_res_ids = self.structure_parser.get_token_res_ids(
            structure=self.structure,
            rep_atom_dict=self.rep_atom_dict,
            only_representative=True,
        )

        c1_res_nums = list(range(p1_region[0], p1_region[1] + 1))
        c2_res_nums = list(range(p2_region[0], p2_region[1] + 1))

        c1_res_idxs = [
            idx for idx, (ch_id, res_id)
            in enumerate(zip(token_rep_chain_ids, token_rep_res_ids))
            if (
                ch_id == chain1 and 
                self.renumber.renumber_chain_res_num(res_id, chain1) in c1_res_nums
            )
        ]
        c2_res_idxs = [
            idx for idx, (ch_id, res_id)
            in enumerate(zip(token_rep_chain_ids, token_rep_res_ids))
            if (
                ch_id == chain2 and 
                self.renumber.renumber_chain_res_num(res_id, chain2) in c2_res_nums
            )
        ]

        avg_pae = self.avg_pae[np.ix_(c1_res_idxs, c2_res_idxs)]

        coords1 = np.array([self.token_coords[idx] for idx in c1_res_idxs])
        coords2 = np.array([self.token_coords[idx] for idx in c2_res_idxs])

        plddt1 = {
            chain1: np.array([self.token_plddts[idx] for idx in c1_res_idxs])
        }
        plddt2 = {
            chain2: np.array([self.token_plddts[idx] for idx in c2_res_idxs])
        }

        # Create a contact map or distance map as specified.
        interaction_map = self.get_contact_map(
            coords1=coords1.reshape(-1, 3),
            coords2=coords2.reshape(-1, 3),
            contact_threshold=self.contact_threshold,
        )

        return interaction_map, plddt1, plddt2, avg_pae

    def apply_confidence_cutoffs(
        self,
        plddt1: dict,
        plddt2: dict,
        avg_pae: np.array
    ):
        """
        mask low-confidence interactions.

        Args:
            plddt1 (dict): pLDDT values for chain 1
            plddt2 (dict): pLDDT values for chain 2

        Returns:
            plddt_matrix (np.array): binary matrix for plddt values >= plddt_cutoff
            avg_pae (np.array): binary matrix for avg_pae values <= pae_cutoff
        """

        chain1, chain2 = next(iter(plddt1)), next(iter(plddt2))
        plddt1, plddt2 = plddt1[chain1], plddt2[chain2]

        assert (
            isinstance(plddt1, np.ndarray) and isinstance(plddt2, np.ndarray)
        ), "pLDDT values should be numpy arrays."

        plddt1, plddt2 = plddt1.reshape(-1, 1), plddt2.reshape(-1,1)

        ch1_cutoff = ch2_cutoff = self.plddt_cutoff

        if chain1 in self.idr_chains:
            ch1_cutoff = self.idr_plddt_cutoff

        if chain2 in self.idr_chains:
            ch2_cutoff = self.idr_plddt_cutoff

        plddt1 = np.where(plddt1 >= ch1_cutoff, 1, 0)
        plddt2 = np.where(plddt2 >= ch2_cutoff, 1, 0)

        plddt_matrix = plddt1 * plddt2.T

        avg_pae = np.where(avg_pae <= self.pae_cutoff, 1, 0)

        return plddt_matrix, avg_pae

    def get_confident_interaction_map(self, region_of_interest: Dict):
        """
        For the specified regions in the predicted structure, obtain all
        confident interacting residue pairs.

        Returns:
            confident_interactions (np.array): binary map of confident
            interacting residues
        """

        interaction_map, plddt1, plddt2, avg_pae = self.get_interaction_data(
            region_of_interest=region_of_interest
        )

        plddt_matrix, pae_matrix = self.apply_confidence_cutoffs(
            plddt1=plddt1, plddt2=plddt2, avg_pae=avg_pae
        )

        confident_interactions = interaction_map * plddt_matrix * pae_matrix

        return confident_interactions

    def get_interacting_patches(
        self,
        contact_map: np.array,
        region_of_interest: dict,
    ):
        """This is a dirty implementation to get the interacting patches. \n
        This is a temporary solution until we find a better way to get interacting
        patches for the given contact map.

        Args:
            contact_map (np.array): binary contact map.
            region_of_interest (dict): region of interest for the protein pair.

        Returns:
            patches (dict): interacting patches for the given region of interest of the protein pair.
        """

        patches = {}

        chain1, chain2 = region_of_interest.keys()
        p1_region, p2_region = region_of_interest[chain1], region_of_interest[chain2]

        if np.unique(contact_map).tolist() == [0]: # No interactions found.
            warnings.warn(
                f"No interacting patches found for {chain1}:{p1_region} and {chain2}:{p2_region}."
            )
            return patches

        matrix_patches = MatrixPatches(
            matrix=contact_map,
            row_obj=chain1,
            col_obj=chain2,
        )

        patches_df = matrix_patches.get_patches_from_matrix()

        for patch_idx, patch in patches_df.iterrows():

            ch1_patch = patch[chain1]
            ch2_patch = patch[chain2]

            ch1_patch = sorted([int(x) for x in ch1_patch])
            ch2_patch = sorted([int(x) for x in ch2_patch])

            ch1_patch = np.array(ch1_patch) + region_of_interest[chain1][0]
            ch2_patch = np.array(ch2_patch) + region_of_interest[chain2][0]

            patches[patch_idx] = {
                chain1: np.array(ch1_patch),
                chain2: np.array(ch2_patch),
            }

        return patches

    def save_ppair_interaction(
        self,
        region_of_interest: dict,
        output_dir: str,
        save_plot: bool = False,
        plot_type: str = "static",
        p1_name: str | None = None,
        p2_name: str | None = None,
        concat_residues: bool = True,
        contact_probability: bool = True,
    ):
        """Save the interacting patches for the given region of interest of the protein pair.

        Args:
            region_of_interest (Dict): Dictionary containing the chain IDs and the residue indices for the region of interest.
            save_plot (bool, optional): Outputs the plot if True. Defaults to False.
            plot_type (str, optional): Type of plot to be saved. Defaults to "static"; options: ["static", "interactive", "both"].
            p1_name (str, optional): Name of the first protein. Defaults to None.
            p2_name (str, optional): Name of the second protein. Defaults to None.
            concat_residues (bool, optional): Whether to concatenate the residues into residue ranges. Defaults to True.
            contact_probability (bool, optional): Whether to add contact probability column to the output. Defaults to True.
        """

        os.makedirs(output_dir, exist_ok=True)

        chain1, chain2 = list(region_of_interest.keys())
        p1_region, p2_region = (
            region_of_interest[chain1],
            region_of_interest[chain2],
        )

        contact_map = self.get_confident_interaction_map(
            region_of_interest=region_of_interest
        )

        interacting_patches = self.get_interacting_patches(
            contact_map=contact_map,
            region_of_interest=region_of_interest,
        )

        if p1_name and p2_name:
            p_names = {
                chain1: p1_name,
                chain2: p2_name,
            }
            dir_name_to_replace = "_".join([p1_name, p2_name])
        else:
            p_names = {
                chain1: chain1,
                chain2: chain2,
            }
            dir_name_to_replace = None

        if len(interacting_patches) > 0:

            file_name = "_".join([
                f"{p_names[k]}_{k}:{v[0]}-{v[1]}"
                for k, v in region_of_interest.items()
            ])

            save_map(
                contact_map=contact_map,
                avg_contact_probs_mat=self.contact_probs,
                patches=interacting_patches,
                chain1=chain1,
                chain2=chain2,
                p1_name=p_names[chain1],
                p2_name=p_names[chain2],
                p1_region=p1_region,
                p2_region=p2_region,
                out_file=os.path.join(output_dir, f"patches_{file_name}.html"),
                save_plot=save_plot,
                plot_type=plot_type,
                concat_residues=concat_residues,
                contact_probability=contact_probability,
                num_to_idx=self.rep_num_to_idx,
            )


def save_map(
    contact_map: np.ndarray,
    avg_contact_probs_mat: np.ndarray | None,
    patches: dict,
    chain1: str,
    chain2: str,
    p1_region: tuple,
    p2_region: tuple,
    out_file: str,
    save_plot=False,
    plot_type="static",
    p1_name: str | None = None,
    p2_name: str | None = None,
    concat_residues: bool = True,
    contact_probability: bool = False,
    num_to_idx: dict = None,
):
    """Save the interacting patches and the contact map to a file.

    Args:
        contact_map (np.ndarray): binary contact map or contact map
        avg_contact_probs_mat (np.ndarray): average contact_probs_mat map
        patches (dict): interacting patches from the map
        interacting_region (dict): interacting region specified by the user
        out_file (str): path to save the output file
        save_plot (bool, optional): save the plot. Defaults to False.
    """

    if contact_probability:
        assert avg_contact_probs_mat is not None; "avg_contact_probs_mat must be provided if contact_probability is True"

    out_dir = os.path.dirname(out_file)
    file_name = os.path.basename(out_file).split(".")[0]

    csv_outfile = os.path.join(out_dir, f"{file_name}.csv")

    print(f"Writing interacting patches to {csv_outfile}")

    save_interaction_patches_df(
        avg_contact_probs_mat,
        patches,
        chain1,
        chain2,
        p1_region,
        p2_region,
        p1_name,
        p2_name,
        concat_residues,
        contact_probability,
        num_to_idx,
        csv_outfile
    )

    if save_plot:
        if plot_type == "interactive" or plot_type == "both":

            fig = plot_map(
                contact_map=contact_map,
                chain1=chain1 if p1_name is None else f"{p1_name}_{chain1}",
                chain2=chain2 if p2_name is None else f"{p2_name}_{chain2}",
                p1_region=p1_region,
                p2_region=p2_region,
                plot_type="interactive",
            )
            out_file = os.path.join(out_dir, f"{file_name}.html")
            fig.write_html(
                out_file,
                full_html=False,
            )

            if plot_type == "both":
                plot_type = "static"

        if plot_type == "static":

            fig = plot_map(
                contact_map=contact_map,
                chain1=chain1 if p1_name is None else f"{p1_name}_{chain1}",
                chain2=chain2 if p2_name is None else f"{p2_name}_{chain2}",
                p1_region=p1_region,
                p2_region=p2_region,
                plot_type="static",
            )

            out_file = os.path.join(out_dir, f"{file_name}.png")
            fig.figure.savefig(out_file)


def save_interaction_patches_df(
    avg_contact_probs_mat,
    patches,
    chain1,
    chain2,
    p1_region,
    p2_region,
    p1_name,
    p2_name,
    concat_residues,
    contact_probability,
    num_to_idx,
    csv_outfile
):
    import pandas as pd
    df_rows = []
    for _, patch in patches.items():
        ch1_res_range = patch[chain1].tolist()
        ch2_res_range = patch[chain2].tolist()

        if concat_residues:
            if contact_probability:
                res1_idxs = np.array([num_to_idx[chain1][res_num] for res_num in ch1_res_range])
                res2_idxs = np.array([num_to_idx[chain2][res_num] for res_num in ch2_res_range])
                contact_probs_mat_res_range = avg_contact_probs_mat[np.ix_(res1_idxs, res2_idxs)]
                avg_contact_prob = np.round(np.mean(contact_probs_mat_res_range), 2)

            ch1_res_range = get_key_from_res_range(ch1_res_range)
            ch2_res_range = get_key_from_res_range(ch2_res_range)

            if contact_probability:
                df_rows.append([ch1_res_range, ch2_res_range, avg_contact_prob])

            else:
                df_rows.append([ch1_res_range, ch2_res_range])

        else:
            for res1, res2 in product(ch1_res_range, ch2_res_range):
                if contact_probability:
                    res1_idx = res1 - p1_region[0]
                    res2_idx = res2 - p2_region[0] + p1_region[1] - p1_region[0] + 1

                    contact_probs_mat_res_range = avg_contact_probs_mat[res1_idx, res2_idx]
                    avg_contact_prob = np.mean(contact_probs_mat_res_range)

                    df_rows.append([res1, res2, avg_contact_prob])

                else:
                    df_rows.append([res1, res2])

    column_names = [f"{chain1}", f"{chain2}"]

    if contact_probability:
        column_names.append("avg_contact_probability")

    if p1_name and p2_name:
        column_names[0] = f"{p1_name}_{chain1}"
        column_names[1] = f"{p2_name}_{chain2}"

    df = pd.DataFrame(df_rows, columns=column_names)

    df.to_csv(csv_outfile, index=False, header=True, sep=",")

def plot_map(
    contact_map: np.ndarray,
    chain1: str,
    chain2: str,
    p1_region: tuple,
    p2_region: tuple,
    plot_type: str
):
    """Plot the contact map

    Args:
        contact_map (np.ndarray): binary contact map or segmented map with labels
    """

    xtick_vals = np.arange(
        0, p2_region[1] - p2_region[0] + 1
    )
    xtick_labels = [str(x+p2_region[0]) for x in xtick_vals]

    ytick_vals = np.arange(
        0, p1_region[1] - p1_region[0] + 1
    )
    ytick_labels = [str(x+p1_region[0]) for x in ytick_vals]

    num_unique_patches = len(np.unique(contact_map))

    colorscale = generate_cmap(
        n=num_unique_patches,
        scheme="binary" if num_unique_patches == 2 else "soft-warm",
    )

    if plot_type == "interactive":

        import plotly.graph_objects as go

        fig = go.Figure(
            data=go.Heatmap(
                z=contact_map,
                colorscale=colorscale,
                xgap=0.2,
                ygap=0.2,
            )
        )

        fig.update_layout(
            title="Contact Map",
            yaxis_title=f"Residue number of {chain1}",
            xaxis_title=f"Residue number of {chain2}",
            xaxis=dict(
                tickmode="array",
                tickformat=".0f",
                tickvals=xtick_vals,
                ticktext=xtick_labels,
            ),
            yaxis=dict(
                tickmode="array",
                tickformat=".0f",
                tickvals=ytick_vals,
                ticktext=ytick_labels,
            ),
        )

    elif plot_type == "static":

            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors

            cmap = mcolors.ListedColormap(colorscale)
            fig, ax = plt.subplots()

            fig = ax.imshow(
                contact_map,
                cmap=cmap,
                interpolation="nearest",
            )

            xtick_labels = [xtick_labels[0], xtick_labels[-1]]
            ytick_labels = [ytick_labels[0], ytick_labels[-1]]
            xtick_vals = [xtick_vals[0], xtick_vals[-1]]
            ytick_vals = [ytick_vals[0], ytick_vals[-1]]

            ax.set_xticks(xtick_vals)
            ax.set_xticklabels(xtick_labels)
            ax.set_yticks(ytick_vals)
            ax.set_yticklabels(ytick_labels)

            ax.set_xlabel(f"Residue number of {chain2}")
            ax.set_ylabel(f"Residue number of {chain1}")
            ax.set_title("Contact Map")

    return fig