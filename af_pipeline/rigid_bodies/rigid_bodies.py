"""
Rigid Bodies extraction module
==============================
RigidBodies class with methods to extract rigid bodies from AlphaFold predictions.
"""
import os
import time
import copy
import json
import warnings
import numpy as np
from itertools import product, combinations_with_replacement
from typing import Dict, List
from collections import defaultdict
from Bio.PDB.Structure import Structure
from af_pipeline._initialize import _Initialize
import matplotlib.patches
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from af_pipeline.tools.structure_tools import (
    save_structure_obj,
    ResidueSelect,
    get_interaction_map,
)
from af_pipeline.pae_to_domains.pae_to_domains import (
    domains_from_pae_matrix_igraph,
    domains_from_pae_matrix_networkx,
    domains_from_pae_matrix_label_propagation
)
# from af_pipeline.rigid_bodies.output_rigid_bodies import (
#     save_rigid_bodies_txt,
#     save_rigid_bodies_json
# )
from af_pipeline.utils.misc_utils import (
    fill_up_the_blanks,
    extract_protein_chain_mapping,
    get_key_from_res_range,
)
from af_pipeline.tools.structure_tools import has_modifications
from af_pipeline.rigid_bodies.rigid_body_assessment import RigidBodyAssessment

class RigidBodies(_Initialize):
    """ Class to extract rigid bodies from AlphaFold prediction."""

    data_file_path: str
    """ Path to AF2/AF3 data file (e.g. .json, .pkl)"""

    structure_file_path: str
    """ Path to AF2/AF3 structure file (e.g. .pdb, .cif)"""

    af_offset: dict | None
    """ Dictionary containing the offset for AF2/AF3 numbering."""

    idr_chains: list
    """ List of chain IDs that represent IDR protein chains."""

    rep_atom_dict: dict
    """ Dictionary of representative atoms for a residue/token."""

    average_token_pae: bool
    """ Whether to average the PAE in case of per atom tokens."""

    average_token_plddt: bool
    """ Whether to average the pLDDT in case of per atom tokens."""

    metric_level: str
    """ Metric level of the instance. Can be "per_token" or "representative_token"."""

    library: str
    """ Library to use for graph-based community detection.
    ('igraph' or 'networkx' or 'label_propagation')"""

    pae_power: int
    """ Exponent to raise the PAE matrix to."""

    pae_cutoff: float
    """ PAE cutoff to consider an edge between two tokens."""

    resolution: float
    """ Resolution parameter for graph-based community detection."""

    plddt_cutoff: float
    """ pLDDT cutoff to filter residues in rigid bodies."""

    plddt_cutoff_idr: float
    """ pLDDT cutoff to filter residues in rigid bodies for IDR chains."""

    random_seed: int
    """ Random seed for label propagation method."""

    def __init__(
        self,
        data_file_path: str,
        structure_file_path: str,
        af_offset: dict = {},
        idr_chains: list = [],
        rep_atom_dict: dict = {},
        average_token_pae: bool = True,
        average_token_plddt: bool = True,
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

        self.library = "networkx"
        self.pae_power = 1
        self.pae_cutoff = 12
        self.resolution = 0.5
        self.plddt_cutoff = 70
        self.plddt_cutoff_idr = 50
        self.random_seed = 47
        self.idr_chains = idr_chains

    def extract_rigid_bodies(
        self,
        num_res: int = 1,
        num_proteins: int = 1,
        plddt_filter: bool = True
    ) -> list[dict[str, list[tuple[str, int]]]]:
        """Extract Rigid bodies from a PAE file.

        Three implementations are available:
        ```python
        - igraph # (community detection using Leiden algorithm)
        - networkx # (community detection using Clauset-Newman-Moore greedy modularity maximization)
        - label_propagation # (community detection using fast label propagation algorithm)
        ```

        Based on the PAE matrix, a graph is constructed where the nodes are
        the residues/tokens and the edges are formed based on the PAE cutoff.
        Communities are detected using the specified implementation.
        Each community is considered as a pseudo-domain.

        Arguments:

        - **num_res (int)**:<br />
            Minimum number of residues in a rigid body.

        - **num_proteins (int)**:<br />
            Minimum number of proteins in a rigid body.

        - **plddt_filter (bool)**:<br />
            Filter the residues based on the pLDDT cutoff.

        Returns:

        - **rigid_bodies (list)**:<br />
            List of extracted rigid bodies.
        """

        print("Extracting rigid bodies...")

        pae_matrix = self.pae

        if self.library == "igraph":
            pseudo_domains = domains_from_pae_matrix_igraph(
                pae_matrix,
                pae_power=self.pae_power,
                pae_cutoff=self.pae_cutoff,
                graph_resolution=self.resolution,
            )

        elif self.library == "networkx":
            pseudo_domains = domains_from_pae_matrix_networkx(
                pae_matrix,
                pae_power=self.pae_power,
                pae_cutoff=self.pae_cutoff,
                graph_resolution=self.resolution,
            )

        elif self.library == "label_propagation":
            pseudo_domains = domains_from_pae_matrix_label_propagation(
                pae_matrix,
                pae_power=self.pae_power,
                pae_cutoff=self.pae_cutoff,
                random_seed=self.random_seed,
            )

        else:
            raise ValueError(
                """
                Invalid library specified.
                Use 'igraph' or 'networkx' or 'label_propagation'.
                """
            )

        # domains is a list of lists
        # each list contains token indices in a domain
        rigid_bodies = []

        for pseudo_domain in pseudo_domains:

            # domain_dict is a dictionary of rigid bodies
            # each rigid body is represented as a dictionary with chain_id as
            # the key and a list of residue numbers as the value
            domain_dict = self._convert_domain_to_dict(pseudo_domain=pseudo_domain)

            # removing residues with pLDDT score below the cutoff
            # different cutoffs can be used for IDR and non-IDR chains
            if plddt_filter:
                rb_dict = self._filter_by_plddt(domain_dict=domain_dict)
            else:
                rb_dict = domain_dict

            # Remove domains with number of proteins less than certain size
            # The size is determind by `num_res` or `num_proteins`
            rb_dict = RigidBodies._filter_by_domain_size(
                rb_dict=rb_dict,
                num_res=num_res,
                num_proteins=num_proteins,
            )

            if len(rb_dict) > 0:
                rigid_bodies.append(rb_dict)

        return rigid_bodies

    def _convert_domain_to_dict(
        self,
        pseudo_domain: list
    ) -> dict[str, list[tuple[str, int]]]:
        """Convert the pseudo-domain list to a dictionary format.

        Example:
        ```python
        pseudo_domain = [0, 1, 5, 6, 8] # these represent token indices
        ```
        will be converted to:
        ```python
        domain_dict = {
            "A": [("CA", 1), ("CB", 2)], # for token indices 0, 1
            "B": [("CA", 1), ("CB", 2), ("CB", 4)] # for token indices 5, 6, 8
        }
        ```
        Where, the tuple is `(atom_name, token_num)`.\n
        The `atom_name` corresponds to the representative atom for the token.\n
        The `token_num` corresponds to the token number within the chain,
        In case of protein, it corresponds to the residue number as per the
        UniProt sequence (provided the `offset`) in case of protein.

        *(See :py:class:`af_pipeline.tools.structure_tools.RenumberResidues` to
        check how the residue numbering is done based on the `offset`.)*

        Arguments:

        - **pseudo_domain (list)**:<br />
            Token indices in a rigid body.

        Returns:

        - **domain_dict (dict)**:<br />
            `{chain_id: [(atom_name, token_num), ...]}`.
        """

        domain_dict = defaultdict(list)

        for token_idx in pseudo_domain:

            token_num = self.idx_to_num[token_idx].get("token_num")
            chain_id = self.idx_to_num[token_idx].get("chain_id")
            atom_name = self.idx_to_num[token_idx].get("atom_name")

            if chain_id not in domain_dict:
                domain_dict[chain_id] = [(atom_name, token_num)]
            else:
                domain_dict[chain_id].append((atom_name, token_num))

        return domain_dict

    def _filter_by_plddt(
        self,
        domain_dict: dict,
    ) -> dict[str, list[tuple[str, int]]]:
        """Filter the residues in the pseudo-domains based on the pLDDT cutoff.

        Only keep the residues with pLDDT >= cutoff in the `domain_dict`.
        Different cutoffs can be used for IDR and non-IDR chains.\n
        *(See :py:attr:`plddt_cutoff_idr` and :py:attr:`plddt_cutoff`.)*

        Arguments:

        - **domain_dict (dict)**:<br />
            Dictionary of pseudo-domains.

        Returns:

        - **rb_dict (dict)**:<br />
            pLDDT filtered dictionary of rigid bodies.
        """

        rb_dict = {}

        # Filter each chain in the pseudo-domain based on the pLDDT cutoff
        for chain_id, atom_name_token_list in domain_dict.items():

            plddt_cutoff = self.plddt_cutoff

            # pLDDT is indexed by token indices
            chain_token_idxs = np.array([
                self.num_to_idx[chain_id][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list
            ])

            if chain_id in self.idr_chains:
                plddt_cutoff = self.plddt_cutoff_idr

            chain_plddt_mask = np.array(self.token_plddts)[
                chain_token_idxs
            ] >= plddt_cutoff

            confident_tokens = [
                (
                    self.idx_to_num[token_idx].get("atom_name"),
                    self.idx_to_num[token_idx].get("token_num")
                )
                for token_idx in chain_token_idxs[chain_plddt_mask]
            ]

            # Update the rigid body dictionary with the confident residues
            domain_dict[chain_id] = confident_tokens

        for chain_id, confident_tokens in domain_dict.items():

            # Remove chains which have no confident residues
            if len(confident_tokens) == 0:
                continue

            rb_dict[chain_id] = confident_tokens

        return rb_dict

    @staticmethod
    def _filter_by_domain_size(
        rb_dict: dict,
        num_res: int,
        num_proteins: int,
    ) -> dict:
        """Filter the domain based on the size of the domain.

        Only keep the domain if it exceeds certain size.
        The size is determined by `num_res` or `num_proteins`.
        ```python
        - num_res # Minimum number of residues in a rigid body.
        - num_proteins # Minimum number of proteins in a rigid body.
        ```

        Arguments:

        - **rb_dict (dict)**:<br />
            Dictionary of a rigid body.

        - **num_res (int)**:<br />
            Minimum number of residues in a rigid body.

        - **num_proteins (int)**:<br />
            Minimum number of proteins in a rigid body.

        Returns:

        - **rb_dict (dict)**:<br />
            Filtered dictionary of a rigid body.
        """

        if len(rb_dict) < num_proteins:
            return {}

        total_residues = 0
        for chain_id, chain_atom_res_list in rb_dict.items():
            chain_res_set = set()
            for atom_name, res_num in chain_atom_res_list:
                chain_res_set.add(res_num)
            total_residues += len(chain_res_set)

        if total_residues < num_res:
            return {}

        return rb_dict

    @staticmethod
    def _keep_residue_numbers_only(
        rigid_bodies: list[dict[str, list[tuple[str, int]]]]
    ) -> list[dict[str, list[int]]]:
        """ Convert the rigid bodies to a list of residue numbers only.

        By default, the rigid body is in the following format.

            {"A": [("CA", 1), ("CB", 2)]}

        where the key is the chain ID and the value is a list of tuples
        containing the `atom_name` and `res_num`.
        This function converts the rigid body to the following format.

            {"A": [1, 2]}

        Arguments:

        - **rigid_bodies (list)**:<br />
            List of rigid bodies, where each rigid body is a dictionary
            with chain IDs as keys and a list of tuples containing
            atom names and residue numbers as values.

        Returns:

        - **rigid_bodies (list)**:<br />
            List of rigid bodies, where each rigid body is a dictionary
            with chain IDs as keys and a list of residue numbers as values.
        """

        if len(rigid_bodies) == 0:

            return rigid_bodies

        elif isinstance(
            list(rigid_bodies[0].values())[0][0], int
        ):
            return rigid_bodies

        rigid_bodies = copy.deepcopy(rigid_bodies)

        for idx, rb_dict in enumerate(rigid_bodies):
            # Convert the rigid body dictionary to a list of residue numbers

            for chain_id, chain_res_num_list in rb_dict.items():

                only_res_num_list = []
                for atom_name, res_num in chain_res_num_list:
                    only_res_num_list.append(res_num)

                only_res_num_list.sort()

                rb_dict[chain_id] = list(set(only_res_num_list))
            # Update the rigid body dictionary with the sorted list
            rigid_bodies[idx] = rb_dict

        return rigid_bodies

    def save_rigid_bodies(
        self,
        domains: list,
        output_dir: str,
        output_format: str = "txt",
        save_structure: bool = True,
        structure_file_type: str = "cif",
        no_plddt_filter_for_structure: bool = False,
        pae_plot: bool = False,
        rb_assessment: dict | None = None,
        protein_chain_map: dict | None = None,
    ):
        """ Save the rigid bodies to a file and/or save the structure of the
        rigid bodies and assess the rigid bodies.

        Output options:
        - The rigid bodies are saved in a plain text format with the chain IDs
          and residue numbers.
        - The structure of the rigid bodies can be saved in PDB or CIF format.
          For rigid bodies with modifications, it is recommended to use PDB format.
        - The PAE plot can be saved to visualize the rigid bodies in the PAE matrix.
        - The rigid bodies can be assessed based on the interface residues,
          number of contacts, interface PAE and pLDDT, average PAE and plDDT and minimum PAE.
        - The assessment is saved in an Excel file.

        parameters for rigid body assessment:\n
        - `as_average`:
        whether to report only the average of assessment metric to the output file. \n
        - `symmetric_pae`:
        whether to report a single average PAE value or assymetric PAE value for PAE assessment metrics. \n

        Arguments:
            domains (list): list of rigid bodies, where each rigid body is a dictionary with chain IDs as keys and residue numbers as values.
            output_dir (str): Directory to save the output files.
            output_format (str, optional): Defaults to "txt". ("txt" or "csv")
            save_structure (bool, optional): Whether to save the structure of the rigid bodies. Defaults to True.
            structure_file_type (str, optional): File type to save the structure. Defaults to "pdb". ("pdb" or "cif")
            no_plddt_filter_for_structure (bool, optional): Whether to save the structure without filtering based on pLDDT. Defaults to False.
            pae_plot (bool, optional): Whether to save the PAE plot for the rigid bodies. Defaults to False.
            rb_assessment (dict | None, optional): Dictionary containing parameters for rigid body assessment.
            protein_chain_map (dict | None, optional): Protein-to-chain mapping dictionary.
        """

        dir_name = os.path.basename(self.structure_file_path).split(".")[0]
        output_dir = os.path.join(output_dir, dir_name)
        os.makedirs(output_dir, exist_ok=True)

        file_name = f"{dir_name}_rigid_bodies"

        protein_chain_map = extract_protein_chain_mapping(
            protein_chain_mapping=protein_chain_map
        )

        ##################################################

        # Save the PAE plot for the rigid bodies
        # the region of the PAE matrix corresponding to the rigid bodies will
        # be highlighted
        if pae_plot:

            for rb_idx, rb_dict in enumerate(domains):

                pae_patches = PAEPatches(
                    num_to_idx=self.num_to_idx,
                    pae=self.pae,
                    lengths_dict=self.lengths_dict,
                    af_offset=self.af_offset,
                    rb_idx=rb_idx,
                )

                # patches are the highlighted rectangles in the PAE matrix
                patches = pae_patches.extract_pae_patches(rb_dict=rb_dict)
                pae_patches.plot_pae_patches(
                    patches=patches,
                    output_dir=output_dir,
                )

        if rb_assessment is not None:

            assessment_file_name = (
                os.path.basename(self.structure_file_path).split(".")[0]
                + "_rb_assessment.xlsx"
            )
            save_path = os.path.join(output_dir, assessment_file_name)

            coords = np.array(self.token_coords)

            contact_map = get_interaction_map(
                coords1=coords,
                coords2=coords,
                contact_threshold=8,
                map_type="contact",
            )

            for rb_idx, rb_dict in enumerate(domains):

                rb_save_path = save_path.replace(
                    ".xlsx", f"_rb_{rb_idx}.xlsx"
                )

                rb_assess = RigidBodyAssessment(
                    rb_dict=rb_dict,
                    num_to_idx=self.num_to_idx,
                    idx_to_num=self.idx_to_num,
                    contact_map=contact_map,
                    plddt_list=self.token_plddts,
                    pae=self.pae,
                    lengths_dict=self.lengths_dict,
                    save_path=rb_save_path,
                    symmetric_pae=rb_assessment.get("symmetric_pae", True),
                    as_average=rb_assessment.get("as_average", False),
                    idr_chains=self.idr_chains,
                    protein_chain_map=protein_chain_map,
                )

                rb_assess.save_rb_assessment()

        # txt or json output format
        if output_format == "txt":

            domains = self._keep_residue_numbers_only(domains)

            save_rigid_bodies_txt(
                output_dir=output_dir,
                domains=domains,
                protein_chain_map=protein_chain_map,
                file_name=file_name,
            )

        elif output_format == "json":

            save_rigid_bodies_json(
                output_dir=output_dir,
                domains=domains,
                protein_chain_map=protein_chain_map,
                file_name=file_name,
            )

        if save_structure:

            domains = self._keep_residue_numbers_only(domains)

            if structure_file_type == "cif" and has_modifications(self.structure):
                warnings.warn(
                    """

                    Protein or nucleotide modifications are stored as HETATM
                    for which sequence connectivity is lost in CIF format due
                    to a bug in Biopython MMCIFParser.
                    Please use PDB format to save the structure with
                    modifications.
                    """
                )

            # Renumber the structure to match the actual sequence numbering
            # if af_offset is provided
            structure = self.renumber.renumber_structure(
                structure=self.structure
            )

            for idx, rb_dict in enumerate(domains):

                # In the following case, the txt or json ouput will have pLDDT
                # filtered residues but, the structure file will ignore this
                # filter use this flag when you don't want missing residues in
                # the structure file
                if no_plddt_filter_for_structure:
                    for chain_id, res_list in rb_dict.items():
                        if len(res_list) > 0:
                            res_list = fill_up_the_blanks(res_list)
                            rb_dict[chain_id] = res_list

                output_path = os.path.join(
                    output_dir, f"rigid_body_{idx}.{structure_file_type}"
                )

                save_structure_obj(
                    structure=structure,
                    out_file=output_path,
                    res_select_obj=ResidueSelect(rb_dict),
                    save_type=structure_file_type,
                    preserve_header_footer=True,
                )


class PAEPatches:
    """ Class to extract and plot PAE patches for rigid bodies."""

    num_to_idx: dict
    """ Dictionary mapping residues to token indices."""

    pae: np.ndarray
    """ PAE matrix."""

    lengths_dict: dict
    """ Dictionary of lengths of each chain."""

    af_offset: dict | None
    """ Offset describing start and end residue number for each chain in
    the predicted structure.\n
    example: `{'A': [1, 100], 'B': [101, 200]}`."""

    rb_idx: int
    """ Rigid body index."""

    def __init__(
        self,
        num_to_idx,
        pae,
        lengths_dict,
        af_offset,
        rb_idx: int,
    ):

        self.num_to_idx = num_to_idx
        self.pae = pae
        self.lengths_dict = lengths_dict
        self.af_offset = af_offset
        self.rb_idx = rb_idx

    def extract_pae_patches(self, rb_dict: dict) -> list[list]:
        """ Extract PAE patches for the rigid body.

        Arguments:

        - **rb_dict (dict)**:<br />
            Dictionary of a rigid body with following `key`:`value` pair:
            `chain_id`: `(atom_name, res_num)`.

        Returns:

        - **patches (list)**:<br />
            List of patches in the format:
            `[[xy (tuple), height (int), width (int)], ...]`.
        """

        patches = []

        rb_list = [(k, v) for k, v in rb_dict.items()]

        chain_combinations = product(rb_list, repeat=2)

        for (ch1, atom_tokens1), (ch2, atom_tokens2) in chain_combinations:
            res_idxs_1 = [
                self.num_to_idx[ch1][res_num][atom_name]
                for (atom_name, res_num) in atom_tokens1
            ]
            res_idxs_2 = [
                self.num_to_idx[ch2][res_num][atom_name]
                for (atom_name, res_num) in atom_tokens2
            ]

            res_idx_range_1 = get_key_from_res_range(
                res_range=res_idxs_1, as_list=True
            )
            res_idx_range_2 = get_key_from_res_range(
                res_range=res_idxs_2, as_list=True
            )

            res_pair_combinations = product(res_idx_range_1, res_idx_range_2)

            for res_idx_1, res_idx_2 in res_pair_combinations:

                if "-" not in res_idx_1 or "-" not in res_idx_2:
                    continue
                # if "-" in res_idx_1 and "-" in res_idx_2:
                res1_y0 = int(res_idx_1.split("-")[0])
                res1_y1 = int(res_idx_1.split("-")[1])
                res2_x0 = int(res_idx_2.split("-")[0])
                res2_x1 = int(res_idx_2.split("-")[1])

                xy_ = (res2_x0, res1_y0) # xy (0,0) coordinates for rectangle
                h_ = res1_y1 - res1_y0 + 1 # patch height
                w_ = res2_x1 - res2_x0 + 1 # patch width

                if h_ > 0 and w_ > 0:
                    patches.append([xy_, h_, w_])

        return patches

    def plot_pae_patches(
        self,
        patches: list,
        output_dir: str,
    ):
        """ Plot the PAE patches for the rigid body.

        Arguments:

        - **patches (list)**:<br />
            List of patches in the format:
            `[[xy (tuple), height (int), width (int)], ...]`.

        - **output_dir (str)**:<br />
            Directory to save the output plot.
        """

        fig = plt.figure(figsize=(20, 20))
        plt.rcParams['font.size'] = 16
        plt.rcParams['axes.titlesize'] = 28
        plt.rcParams['axes.labelsize'] = 22
        plt.rcParams['xtick.labelsize'] = 13
        plt.rcParams['ytick.labelsize'] = 13
        plt.imshow(
            self.pae,
            # cmap="Greens_r",
            cmap="Greys_r",
            vmax=31.75,
            vmin=0,
            interpolation="nearest",
            )

        for xy, h, w in patches:
            rect = matplotlib.patches.Rectangle(
                xy,
                w,
                h,
                linewidth=0,
                # edgecolor="green",
                facecolor="lime",
                alpha=0.5,
            )
            plt.gca().add_patch(rect)

        cumu_len = 0
        ticks = []
        ticks_labels = []

        for chain_id, p_length in self.lengths_dict.items():
            if chain_id == "total":
                continue

            cumu_len += p_length

            if cumu_len != self.pae.shape[1]:
                plt.axhline(
                    y=cumu_len,
                    color='red',
                    linestyle='--',
                    linewidth=0.75,
                )
                plt.axvline(
                    x=cumu_len,
                    color='red',
                    linestyle='--',
                    linewidth=0.75,
                )

            if self.af_offset is not None:
                ticks_labels.extend([
                    f"\n{self.af_offset.get(chain_id, [1, 1])[0]}" ,
                    f"{self.af_offset.get(chain_id, [1, 1])[1]}\n"
                ])
            else:
                ticks_labels.extend([
                    "\n1",
                    f"{self.lengths_dict[chain_id]}\n"
                ])

            if cumu_len-p_length not in ticks:
                ticks.extend([cumu_len-p_length, cumu_len])
            else:
                ticks.extend([cumu_len-p_length+1, cumu_len])

        plt.xlim(0, self.pae.shape[0])
        plt.ylim(0, self.pae.shape[1])

        plt.gca().invert_yaxis()
        plt.yticks(ticks, ticks_labels)

        plt.xticks(ticks, ticks_labels, rotation=90, ha='center')
        plt.title(f"Predicted aligned error (PAE)", pad=20)
        plt.xlabel("Scored residue")
        plt.ylabel("Aligned residue")

        ax = plt.gca()

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=1.2)
        plt.colorbar(
            label="Predicted Alignment Error (PAE)",
            orientation="horizontal",
            cax=cax,
        )

        plt.savefig(
            os.path.join(output_dir, f"rigid_body_{self.rb_idx}.png"),
            transparent=True
        )
        plt.close(fig)

def save_rigid_bodies_txt(
    output_dir: str,
    domains: list,
    protein_chain_map: dict,
    file_name: str = "rigid_bodies",
):
    """ Save rigid bodies to a text file.

    This function writes the rigid bodies information to a text file in a human-readable format.

    The output file will contain the rigid body index, chain ID, protein name (if available), and the residue range.

    Args:

        output_dir (str):
            Directory where the output file will be saved.

        domains (list):
            List of dictionaries, where each dictionary represents a rigid body.

        protein_chain_map (dict):
            A mapping of chain IDs to protein names.

        file_name (str, optional):
            Name of the output file without extension.
            Defaults to "rigid_bodies".
    """

    file_name += ".txt"
    output_path = os.path.join(output_dir, file_name)

    with open(output_path, "w") as f:

        for idx, rb_dict in enumerate(domains):
            f.write(f"Rigid Body {idx}\n")

            for chain_id, res_list in rb_dict.items():

                protein_name = protein_chain_map.get(chain_id, None)

                if len(res_list) > 0:
                    if protein_name:
                        f.write(
                            f"{protein_name}_{chain_id}: {get_key_from_res_range(res_range=res_list)}\n"
                        )
                    else:
                        f.write(
                            f"{chain_id}:{get_key_from_res_range(res_range=res_list)}\n"
                        )

            f.write("\n")

def save_rigid_bodies_json(
    output_dir: str,
    domains: list,
    protein_chain_map: dict,
    file_name: str = "rigid_bodies",
):
    """ Save rigid bodies to a JSON file.

    This function writes the rigid bodies information to a JSON file.

    For per-atom tokens JSON format is recommended over text format.

    Args:

        output_dir (str):
            Directory where the output file will be saved.

        domains (list):
            List of dictionaries, where each dictionary represents a rigid body.

        protein_chain_map (dict):
            A mapping of chain IDs to protein names.

        file_name (str, optional):
            Name of the output file without extension.
            Defaults to "rigid_bodies".
    """

    file_name += ".json"
    output_path = os.path.join(output_dir, file_name)

    rigid_bodies = []
    for idx, rb_dict in enumerate(domains):
        ch_dict = {}
        for chain_id, res_num_list in rb_dict.items():
            protein_name = protein_chain_map.get(chain_id, None)
            if not protein_name:
                protein_name = "Unknown"
            ch_dict[chain_id] = {
                "protein": protein_name,
                "residues": []
            }
            for atom_name, res_num in res_num_list:
                ch_dict[chain_id]["residues"].append(
                    (atom_name, res_num)
                )
        rigid_bodies.append(ch_dict)

    with open(output_path, "w") as f:
        json.dump(rigid_bodies, f, indent=4)