"""
Interaction module
==============================
Interaction class to handle interaction data for the predicted structure. \n
One can obtain:
1. Interaction map: A binary contact map or distance map.
2. Restraints: Contacts as pairwise residues in a DataFrame format.
3. Interacting patches: contiguous regions in the interaction map obtained in (1).
"""
import os
from textwrap import dedent
import warnings
import numpy as np
import pandas as pd
from itertools import product
from typing import Dict, Optional
from scipy.spatial import distance_matrix
from af_pipeline.utils.misc_utils import (
    get_key_from_res_range,
    generate_cmap,
)
from af_pipeline.parser.initialize import Initialize
from af_pipeline.tools.matrix_patches import MatrixPatches
from af_pipeline.constants.af_constants import InteractionConstants as IntCons
from af_pipeline.constants.af_constants import (
    ColorMapScheme,
    MetricLevel,
    KeywordArg,
    FileFormat,
)
import plotly.graph_objects as go
import matplotlib.pyplot as plt

_error_not_set_up = """
The RigidBodies instance is not set up yet. Please set up the instance
by calling the `set_attributes_from` method with an instance of the
Initialize class.

Alternatively, if you know what attributes to set, you can set the attributes of
the RigidBodies instance directly without using the `set_attributes_from` method.
and set the `is_set_up` attribute to True.

The following attributes need to be set for the RigidBodies instance to work properly:
- structure: Structure
- structure_parser: StructureParser
- token_chain_ids: list
- token_coords: np.ndarray
- token_plddts: np.ndarray
- lengths_dict: dict
- rep_atom_dict: dict
- pae: np.ndarray
- avg_pae: np.ndarray
- contact_probs: np.ndarray
- renumber: af_pipeline.tools.structure_tools.RenumberResidues
- rep_num_to_idx: dict


See specific methods for more details on the required attributes for each method.
"""

class Interaction:
    """ Class to handle interaction data for the predicted structure. """

    contact_threshold: float = IntCons.contact_threshold
    """ Threshold for defining a contact between residues in Angstroms. """

    plddt_cutoff: float = IntCons.plddt_cutoff
    """ Threshold for defining a confident residue based on pLDDT values. """

    pae_cutoff: float = IntCons.pae_cutoff
    """ Threshold for defining a confident residue pair based on PAE values. """

    plddt_cutoff_idr: Optional[float] = IntCons.plddt_cutoff_idr
    """ Threshold for defining a confident residue based on pLDDT values for
    intrinsically disordered regions.
    Only takes effect if `idr_chains` is a non-empty list.
    """

    idr_chains: Optional[list] = []
    """ List of chains that are intrinsically disordered. """

    save_plot: Optional[bool] = IntCons.save_plot
    """ Whether to save the plot of the interaction map. """

    save_table: Optional[bool] = IntCons.save_table
    """ Whether to save the table of interacting residues. """

    def __init__(
        self,
        contact_threshold: float = IntCons.contact_threshold,
        plddt_cutoff: float = IntCons.plddt_cutoff,
        pae_cutoff: float = IntCons.pae_cutoff,
        **kwargs,
    ):

        self.contact_threshold = contact_threshold
        self.plddt_cutoff = plddt_cutoff
        self.plddt_cutoff_idr = kwargs.get(
            KeywordArg.PLDDT_CUTOFF_IDR, IntCons.plddt_cutoff_idr
        )
        self.pae_cutoff = pae_cutoff
        self.idr_chains = kwargs.get(KeywordArg.IDR_CHAINS, []) # List of chains that are disordered

        self.save_plot = kwargs.get(KeywordArg.SAVE_PLOT, IntCons.save_plot)
        self.save_table = kwargs.get(KeywordArg.SAVE_TABLE, IntCons.save_table)

        self._is_set_up = False
        setup_instance = kwargs.get(KeywordArg.SETUP_INSTANCE, None)
        if isinstance(setup_instance, Initialize):
            self.set_attributes_from(instance=setup_instance)

    def check_is_set_up(self):
        """ Check if the Interaction instance is set up. """

        if not self._is_set_up:
            raise ValueError(_error_not_set_up)

    def set_attributes_from(self, instance: Initialize):
        """ Set the attributes of Interaction instance from the Initialize instance.
        This method can be used to set the attributes of Interaction instance

        ## Arguments:

        - **instance (Initialize)**:<br />
            An instance of the Initialize class.
        """

        if instance.metric_level != MetricLevel.REPRESENTATIVE_TOKEN:
            raise NotImplementedError(dedent(f"""
                Currently, Interaction class only considers interactions at the
                residue-level and not at atomic-level.
                Hence, Initialize instance should be initialized with:
                metric_level = "representative_token" """)
            )

        self.structure = instance.structure
        self.structure_parser = instance.structure_parser

        self.token_chain_ids = instance.token_chain_ids
        self.token_coords = instance.token_coords
        self.token_plddts = instance.token_plddts

        self.lengths_dict = instance.lengths_dict

        self.rep_atom_dict = instance.rep_atom_dict

        self.pae = instance.pae
        self.avg_pae = instance.avg_pae

        self.contact_probs = instance.contact_probs

        self.renumber = instance.renumber
        self.rep_num_to_idx = instance.rep_num_to_idx

        self._is_set_up = True

    @staticmethod
    def get_contact_map(
        coords1: np.ndarray,
        coords2: np.ndarray,
        contact_threshold: float = IntCons.contact_threshold,
    ) -> np.ndarray:
        """ Get the contact map between two arrays of coordinates.

        ## Arguments:

        - **coords1 (np.ndarray)**:<br />
            Coordinates of the chain1 residues.

        - **coords2 (np.ndarray)**:<br />
            Coordinates of the chain2 residues.

        - **contact_threshold (float, optional):**:<br />
            Threshold for defining a contact between residues in
            Angstroms. Defaults to 8.0.

        ## Returns:

        - **np.ndarray**:<br />
            Binary contact map where 1 indicates a contact and 0 indicates
            no contact.
        """

        distance_map = distance_matrix(coords1, coords2)

        contact_map = np.where(distance_map < contact_threshold, 1, 0)

        return contact_map

    def create_regions_of_interest(self) -> list:
        """
        Create regions of interest for all possible chain pairs.

        ## Returns:

        - **list**:<br />
            A list of dictionaries, where each dictionary contains the chain IDs
            and the start and end residue numbers for the region of interest
            for a pair of chains.
        """

        self.check_is_set_up()
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

    def get_interaction_data(self, region_of_interest: Dict) -> tuple:
        """ Get the interaction amp, pLDDT, and PAE for the region of interest.

        ## Arguments:

        - **region_of_interest (Dict)**:<br />
            Dictionary containing the chain IDs and the residue numbers for the
            region of interest. The keys should be the chain IDs and the values
            should be lists of the form [start_res_num, end_res_num].

        ## Returns:

        - **tuple**:<br />
            A tuple containing the following elements:
            - **interaction_map (np.array)**:<br />
                Binary contact map or distance map for the region of interest.
            - **plddt1 (dict)**:<br />
                pLDDT values for the residues in chain 1 of the region of interest.
            - **plddt2 (dict)**:<br />
                pLDDT values for the residues in chain 2 of the region of interest.
            - **avg_pae (np.array)**:<br />
                Average PAE values for the region of interest.
        """

        self.check_is_set_up()
        chains = list(region_of_interest.keys())

        assert len(chains) == 2, f"Region of interest is for a chain-pair. Got {len(chains)} chains"

        chain1, chain2 = chains[0], chains[1]
        p1_region = region_of_interest[chain1]
        p2_region = region_of_interest[chain2]

        c1_res_nums = list(range(p1_region[0], p1_region[1] + 1))
        c2_res_nums = list(range(p2_region[0], p2_region[1] + 1))

        c1_res_idxs = [self.rep_num_to_idx[chain1][token_num] for token_num in c1_res_nums]
        c2_res_idxs = [self.rep_num_to_idx[chain2][token_num] for token_num in c2_res_nums]

        avg_pae = self.avg_pae[np.ix_(c1_res_idxs, c2_res_idxs)]

        coords1 = np.array([self.token_coords[idx] for idx in c1_res_idxs])
        coords2 = np.array([self.token_coords[idx] for idx in c2_res_idxs])

        plddt1 = {chain1: np.array([self.token_plddts[idx] for idx in c1_res_idxs])}
        plddt2 = {chain2: np.array([self.token_plddts[idx] for idx in c2_res_idxs])}

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
        avg_pae: np.ndarray
    ) -> tuple:
        """ Mask low-confidence interactions based on pLDDT and PAE cutoffs.

        ## Arguments:

        - **plddt1 (dict)**:<br />
            pLDDT values for chain 1.

        - **plddt2 (dict)**:<br />
            pLDDT values for chain 2.

        - **avg_pae (np.ndarray)**:<br />
            Average PAE values matrix.

        ## Returns:

        - **tuple**:<br />
            Binary matrices for pLDDT and PAE values that are above the respective cutoffs.
        """

        chain1, chain2 = next(iter(plddt1)), next(iter(plddt2))
        plddt1, plddt2 = plddt1[chain1], plddt2[chain2]

        assert (
            isinstance(plddt1, np.ndarray) and isinstance(plddt2, np.ndarray)
        ), "pLDDT values should be numpy arrays."

        plddt1, plddt2 = plddt1.reshape(-1, 1), plddt2.reshape(-1,1)

        ch1_cutoff = ch2_cutoff = self.plddt_cutoff

        if chain1 in self.idr_chains:
            ch1_cutoff = self.plddt_cutoff_idr

        if chain2 in self.idr_chains:
            ch2_cutoff = self.plddt_cutoff_idr

        plddt1 = np.where(plddt1 >= ch1_cutoff, 1, 0)
        plddt2 = np.where(plddt2 >= ch2_cutoff, 1, 0)

        plddt_matrix = plddt1 * plddt2.T

        avg_pae = np.where(avg_pae <= self.pae_cutoff, 1, 0)

        return plddt_matrix, avg_pae

    def get_confident_interaction_map(self, region_of_interest: Dict) -> np.ndarray:
        """ For the specified regions in the predicted structure, obtain all
        confident interacting residue pairs.

        ## Arguments:

        - **region_of_interest (Dict)**:<br />
            Dictionary containing the chain IDs and the residue numbers for the
            region of interest. The keys should be the chain IDs and the values
            should be lists of the form [start_res_num, end_res_num].

        ## Returns:

        - **np.ndarray**:<br />
            Binary map of confident interacting residues, where 1 indicates a
            confident interaction and 0 indicates no confident interaction.
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
        contact_map: np.ndarray,
        region_of_interest: dict,
    ) -> dict:
        """ This is a dirty implementation to get the interacting patches. \n
        This is a temporary solution until we find a better way to get interacting
        patches for the given contact map.

        ## Arguments:

        - **contact_map (np.ndarray)**:<br />
            Binary contact map for the region of interest, where 1 indicates a
            contact and 0 indicates no contact.

        - **region_of_interest (dict)**:<br />
            Dictionary containing the chain IDs and the residue numbers for the
            region of interest. The keys should be the chain IDs and the values
            should be lists of the form [start_res_num, end_res_num].

        ## Returns:

        - **dict**:<br />
            Dictionary containing the interacting patches for the given region of
            interest of the protein pair. The keys are the patch indices and the
            values are dictionaries containing the chain IDs and the residue numbers
            for the interacting patches.
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
        save_plot: bool = IntCons.save_plot,
        plot_type: str = IntCons.plot_type,
        p1_name: str | None = None,
        p2_name: str | None = None,
        concat_residues: bool = True,
        contact_probability: bool = True,
    ):
        """ Save the interacting patches for the given region of interest of the protein pair.

        ## Arguments:

        - **region_of_interest (dict)**:<br />
            Dictionary containing the chain IDs and the residue numbers for the
            region of interest. The keys should be the chain IDs and the values
            should be lists of the form [start_res_num, end_res_num].

        - **output_dir (str)**:<br />
            Directory to save the output files.

        - **save_plot (bool, optional):**:<br />
            Whether to save the plot of the interaction map. Defaults to False.

        - **plot_type (str, optional):**:<br />
            Type of plot to be saved. Defaults to "static". Valid options are:
            "static", "interactive", and "both".

        - **p1_name (str | None, optional):**:<br />
            Name of the first protein. Defaults to None.

        - **p2_name (str | None, optional):**:<br />
            Name of the second protein. Defaults to None.

        - **concat_residues (bool, optional):**:<br />
            Whether to concatenate the residues into residue ranges. Defaults to True.

        - **contact_probability (bool, optional):**:<br />
            Whether to add contact probability column to the output. Defaults to True.
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
                out_file=os.path.join(output_dir, f"patches_{file_name}.{FileFormat.HTML}"),
                save_plot=save_plot,
                plot_type=plot_type,
                concat_residues=concat_residues,
                contact_probability=contact_probability,
                num_to_idx=self.rep_num_to_idx,
            )


    """

    Args:
        contact_map (np.ndarray): binary contact map or contact map
        avg_contact_probs_mat (np.ndarray): average contact_probs_mat map
        patches (dict): interacting patches from the map
        interacting_region (dict): interacting region specified by the user
        out_file (str): path to save the output file
        save_plot (bool, optional): save the plot. Defaults to False.
    """
def save_map(
    contact_map: np.ndarray,
    avg_contact_probs_mat: np.ndarray | None,
    patches: dict,
    chain1: str,
    chain2: str,
    p1_region: tuple,
    p2_region: tuple,
    out_file: str,
    save_plot=IntCons.save_plot,
    plot_type=IntCons.plot_type,
    p1_name: str | None = None,
    p2_name: str | None = None,
    concat_residues: bool = True,
    contact_probability: bool = False,
    num_to_idx: dict = None,
):
    """ Save the interacting patches and the contact map to a file.

    ## Arguments:

    - **contact_map (np.ndarray)**:<br />
        Binary contact map or contact map for the region of interest, where 1
        indicates a contact and 0 indicates no contact.

    - **avg_contact_probs_mat (np.ndarray | None)**:<br />
        Average contact probabilities matrix for the region of interest.
        This is used to calculate the average contact probability for each patch
        if contact_probability is True.

    - **patches (dict)**:<br />
        Dictionary containing the interacting patches for the given region of
        interest of the protein pair. The keys are the patch indices and the
        values are dictionaries containing the chain IDs and the residue numbers
        for the interacting patches.

    - **chain1 (str)**:<br />
        Chain ID for the first chain in the region of interest.

    - **chain2 (str)**:<br />
        Chain ID for the second chain in the region of interest.

    - **p1_region (tuple)**:<br />
        Tuple containing the start and end residue numbers for the first chain in
        the region of interest.

    - **p2_region (tuple)**:<br />
        Tuple containing the start and end residue numbers for the second chain in
        the region of interest.

    - **out_file (str)**:<br />
        Path to save the output file. The file will be saved in CSV format and
        the plot will be saved in the format specified by plot_type.

    - **save_plot (_type_, optional):**:<br />
        Whether to save the plot of the interaction map. Defaults to False.

    - **plot_type (_type_, optional):**:<br />
        Type of plot to be saved. Defaults to "static". Valid options are:
        "static", "interactive", and "both".

    - **p1_name (str | None, optional):**:<br />
        Name of the first chain in the region of interest.

    - **p2_name (str | None, optional):**:<br />
        Name of the second chain in the region of interest.

    - **concat_residues (bool, optional):**:<br />
        Whether to concatenate the residues into residue ranges. Defaults to True.

    - **contact_probability (bool, optional):**:<br />
        Whether to add contact probability column to the output. Defaults to False.

    - **num_to_idx (dict, optional):**:<br />
        Dictionary mapping residue numbers to indices in the contact probabilities
         matrix. This is required if contact_probability is True. Defaults to None.
    """

    if contact_probability:
        assert avg_contact_probs_mat is not None; (
            "avg_contact_probs_mat must be provided if contact_probability is True"
        )

    out_dir = os.path.dirname(out_file)
    file_name = os.path.basename(out_file).split(".")[0]

    csv_outfile = os.path.join(out_dir, f"{file_name}.{FileFormat.CSV}")

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
        if plot_type in IntCons.valid_plot_types[1:3]: # interactive or both

            fig = plot_map(
                contact_map=contact_map,
                chain1=chain1 if p1_name is None else f"{p1_name}_{chain1}",
                chain2=chain2 if p2_name is None else f"{p2_name}_{chain2}",
                p1_region=p1_region,
                p2_region=p2_region,
                plot_type=IntCons.valid_plot_types[1],
            )
            out_file = os.path.join(out_dir, f"{file_name}.{FileFormat.HTML}")
            fig.write_html(
                out_file,
                full_html=False,
            )

            if plot_type == IntCons.valid_plot_types[2]:
                plot_type = IntCons.valid_plot_types[0] # static

        if plot_type == IntCons.valid_plot_types[0]:

            fig = plot_map(
                contact_map=contact_map,
                chain1=chain1 if p1_name is None else f"{p1_name}_{chain1}",
                chain2=chain2 if p2_name is None else f"{p2_name}_{chain2}",
                p1_region=p1_region,
                p2_region=p2_region,
                plot_type=IntCons.valid_plot_types[0],
            )

            out_file = os.path.join(out_dir, f"{file_name}.{FileFormat.PNG}")
            fig.figure.savefig(out_file)


def save_interaction_patches_df(
    avg_contact_probs_mat: np.ndarray | None,
    patches: dict,
    chain1: str,
    chain2: str,
    p1_region: tuple,
    p2_region: tuple,
    p1_name: str | None,
    p2_name: str | None,
    concat_residues: bool,
    contact_probability: bool,
    num_to_idx: dict | None,
    csv_outfile: str,
):
    """ Save the interacting patches df as csv file.

    ## Arguments:

    - **avg_contact_probs_mat (np.ndarray | None)**:<br />
        Average contact probabilities matrix for the region of interest. This is
        used to calculate the average contact probability for each patch if
        contact_probability is True.

    - **patches (dict)**:<br />
        Dictionary containing the interacting patches for the given region of
        interest of the protein pair. The keys are the patch indices and the
        values are dictionaries containing the chain IDs and the residue numbers
        for the interacting patches.

    - **chain1 (str)**:<br />
        Chain ID for the first chain in the region of interest.

    - **chain2 (str)**:<br />
        Chain ID for the second chain in the region of interest.

    - **p1_region (tuple)**:<br />
        Tuple containing the start and end residue numbers for the first chain in
        the region of interest.

    - **p2_region (tuple)**:<br />
        Tuple containing the start and end residue numbers for the second chain in
        the region of interest.

    - **p1_name (str | None)**:<br />
        Name of the first protein in the region of interest.

    - **p2_name (str | None)**:<br />
        Name of the second protein in the region of interest.

    - **concat_residues (bool)**:<br />
        Whether to concatenate the residues into residue ranges.

    - **contact_probability (bool)**:<br />
        Whether to add contact probability column to the output.

    - **num_to_idx (dict | None)**:<br />
        Dictionary mapping residue numbers to indices in the contact probabilities

    - **csv_outfile (str)**:<br />
        Path to save the output CSV file.
    """

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
) -> go.Figure | plt.Figure:
    """ Plot the contact map.

    ## Arguments:

    - **contact_map (np.ndarray)**:<br />
        Binary contact map

    - **chain1 (str)**:<br />
        Chain ID for the first chain in the region of interest.

    - **chain2 (str)**:<br />
        Chain ID for the second chain in the region of interest.

    - **p1_region (tuple)**:<br />
        Tuple containing the start and end residue numbers for the first chain in
        the region of interest.

    - **p2_region (tuple)**:<br />
        Tuple containing the start and end residue numbers for the second chain in
        the region of interest.

    - **plot_type (str)**:<br />
        Type of plot to generate. Options are "static" or "interactive".

    ## Returns:

    - ****:<br />
        The generated plot object. The type of the plot object depends on the
        plot_type argument.
        - For "static", it returns a matplotlib figure object.
        - For "interactive", it returns a plotly figure object.
    """
    xtick_vals = np.arange(0, p2_region[1] - p2_region[0] + 1)
    xtick_labels = [str(x+p2_region[0]) for x in xtick_vals]

    ytick_vals = np.arange(0, p1_region[1] - p1_region[0] + 1)
    ytick_labels = [str(x+p1_region[0]) for x in ytick_vals]

    num_unique_patches = len(np.unique(contact_map))

    colorscale = generate_cmap(
        n=num_unique_patches,
        scheme=ColorMapScheme.BINARY if num_unique_patches == 2 else ColorMapScheme.SOFT_WARM,
    )

    if plot_type not in IntCons.valid_plot_types:
        raise ValueError(
            f"Invalid plot type: {plot_type}. Valid options are: {IntCons.valid_plot_types}"
        )

    if plot_type == IntCons.valid_plot_types[1]: # interactive

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

    elif plot_type == IntCons.valid_plot_types[0]: # static

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