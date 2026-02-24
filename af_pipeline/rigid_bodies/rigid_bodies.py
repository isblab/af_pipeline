"""
Rigid Bodies extraction module
==============================
RigidBodies class with methods to extract rigid bodies from AlphaFold predictions.
"""
import os
import copy
import json
import warnings
import numpy as np
from itertools import product
from typing import Dict, List
from collections import defaultdict
from Bio.PDB.Structure import Structure
from af_pipeline.initialize import Initialize
import matplotlib.patches
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from af_pipeline.tools.structure_tools import (
    save_structure_obj,
    ResidueSelect,
)
from af_pipeline.pae_to_domains.pae_to_domains import (
    domains_from_pae_matrix_igraph,
    domains_from_pae_matrix_networkx,
    domains_from_pae_matrix_label_propagation
)
from af_pipeline.utils.misc_utils import (
    fill_up_the_blanks,
    extract_protein_chain_mapping,
    get_key_from_res_range,
)
from af_pipeline.tools.structure_tools import has_modifications
from af_pipeline.rigid_bodies.rigid_body_assessment import RigidBodyAssessment
from af_pipeline.constants.af_constants import RigidBodiesConstants as RBCons
from af_pipeline.constants.af_constants import (
    RES_SEPARATOR,
    FileFormat,
    KeywordArg,
    MiscStrEnum,
    CommunityDetectionLibrary,
    ResidueMapKeys,
)

_error_not_set_up = """
The RigidBodies instance is not set up yet. Please set up the instance
by calling the `set_from_initializer` method with an instance of the
Initialize class.

Alternatively, if you know what attributes to set, you can set the attributes of
the RigidBodies instance directly without using the `set_from_initializer` method.
and set the `is_set_up` attribute to True.

The following attributes need to be set for the RigidBodies instance to work properly:
- structure_file_path: str
- structure: Structure
- af_offset: dict | None
- token_coords: list
- token_plddts: list
- pae: np.ndarray
- avg_pae: np.ndarray
- lengths_dict: dict
- renumber: af_pipeline.tools.structure_tools.RenumberResidues
- idx_to_num: dict (can be obtained from renumber object)
- num_to_idx: dict (can be obtained from renumber object)

See specific methods for more details on the required attributes for each method.
"""

class RigidBodies:
    """ Class to extract rigid bodies from AlphaFold prediction."""

    library: str
    """ Library to use for graph-based community detection.
    ('igraph' or 'networkx' or 'label_propagation')"""

    pae_cutoff: float
    """ PAE cutoff to consider an edge between two tokens."""

    pae_power: int
    """ Exponent to raise the PAE matrix to."""

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
        library: str = RBCons.library,
        pae_cutoff: float = RBCons.pae_cutoff,
        pae_power: int = RBCons.pae_power,
        resolution: float = RBCons.resolution,
        plddt_cutoff: float = RBCons.plddt_cutoff,
        **kwargs,
    ):

        self.library = library
        self.pae_cutoff = pae_cutoff
        self.pae_power = pae_power
        self.resolution = resolution
        self.plddt_cutoff = plddt_cutoff
        self.idr_chains = kwargs.get(KeywordArg.IDR_CHAINS, [])
        self.plddt_cutoff_idr = kwargs.get(
            KeywordArg.PLDDT_CUTOFF_IDR, RBCons.plddt_cutoff_idr
        )
        self.random_seed = kwargs.get(
            KeywordArg.RANDOM_SEED, RBCons.random_seed
        )

        self._is_set_up = False

    def check_is_set_up(self):
        """ Check if the RigidBodies instance is set up. """

        if not self._is_set_up:
            raise ValueError(_error_not_set_up)

    def set_attributes_from(self, instance: Initialize):
        """ Set the attributes of RigidBodies instance from the initializer instance.

        This method can be used to set the attributes of RigidBodies instance
        from the initializer instance after the RigidBodies instance is created.
        This is useful when the initializer instance is created after the
        RigidBodies instance is created.

        ## Arguments:

        - **instance (Initialize)**:<br />
            An instance of the Initialize class.
        """

        assert isinstance(instance, Initialize), "The instance should be of type Initialize."

        self.structure = instance.structure
        self.structure_file_path = instance.structure_file_path
        self.af_offset = instance.af_offset

        self.token_coords = instance.token_coords
        self.token_plddts = instance.token_plddts

        self.pae = instance.pae
        self.avg_pae = instance.avg_pae

        self.lengths_dict = instance.lengths_dict

        self.renumber = instance.renumber
        self.idx_to_num = instance.idx_to_num
        self.num_to_idx = instance.num_to_idx

        self._is_set_up = True

    def extract_rigid_bodies(
        self,
        pae_matrix: np.ndarray,
        min_res: int = RBCons.min_res,
        min_proteins: int = RBCons.min_proteins,
        plddt_filter: bool = RBCons.plddt_filter,
    ) -> list[dict[str, list[tuple[str, int]]]]:
        """Extract Rigid bodies from a PAE file.

        Three implementations for community detection are available:
        ```python
        - igraph # (Leiden algorithm)
        - networkx # (Clauset-Newman-Moore greedy modularity maximization)
        - label_propagation # (fast label propagation algorithm)
        ```

        Based on the PAE matrix, a graph is constructed where the nodes are
        the residues/tokens and the edges are formed based on the PAE cutoff.
        Communities are detected using the specified implementation.
        Each community is considered as a pseudo-domain.

        ## Arguments:

        - **min_res (int)**:<br />
            Minimum number of residues in a rigid body.

        - **min_proteins (int)**:<br />
            Minimum number of proteins in a rigid body.

        - **plddt_filter (bool)**:<br />
            Filter the residues based on the pLDDT cutoff.

        ## Returns:

        - **rigid_bodies (list)**:<br />
            List of extracted rigid bodies.
        """

        self.check_is_set_up()
        print("Extracting rigid bodies...")
        if self.library not in RBCons.valid_libraries:
            raise ValueError(
                "Invalid library specified."\
                f"Use one of {RBCons.valid_libraries}."
            )

        if self.library == CommunityDetectionLibrary.IGRAPH:
            pseudo_domains = domains_from_pae_matrix_igraph(
                pae_matrix,
                pae_power=self.pae_power,
                pae_cutoff=self.pae_cutoff,
                graph_resolution=self.resolution,
            )

        elif self.library == CommunityDetectionLibrary.NETWORKX:
            pseudo_domains = domains_from_pae_matrix_networkx(
                pae_matrix,
                pae_power=self.pae_power,
                pae_cutoff=self.pae_cutoff,
                graph_resolution=self.resolution,
            )

        elif self.library == CommunityDetectionLibrary.LABEL_PROPAGATION:
            pseudo_domains = domains_from_pae_matrix_label_propagation(
                pae_matrix,
                pae_power=self.pae_power,
                pae_cutoff=self.pae_cutoff,
                random_seed=self.random_seed,
            )

        # domains is a list of lists
        # each list contains token indices in a domain
        rigid_bodies = []

        for pseudo_domain in pseudo_domains:

            # domain_dict is a dictionary of rigid bodies
            # each rigid body is represented as a dictionary with chain_id as
            # the key and a list of residue numbers as the value
            domain_dict = self._convert_domain_to_dict(pseudo_domain)

            # removing residues with pLDDT score below the cutoff
            # different cutoffs can be used for IDR and non-IDR chains
            if plddt_filter:
                rb_dict = self._filter_by_plddt(
                    domain_dict=domain_dict,
                    token_plddts=self.token_plddts,
                )
            else:
                rb_dict = domain_dict

            # Remove domains with number of proteins less than certain size
            # The size is determind by `min_res` or `min_proteins`
            rb_dict = RigidBodies._filter_by_domain_size(
                rb_dict=rb_dict,
                min_res=min_res,
                min_proteins=min_proteins,
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

        ## Arguments:

        - **pseudo_domain (list)**:<br />
            Token indices in a rigid body.

        ## Returns:

        - **domain_dict (dict)**:<br />
            `{chain_id: [(atom_name, token_num), ...]}`.
        """

        self.check_is_set_up()
        domain_dict = defaultdict(list)

        for token_idx in pseudo_domain:

            token_num = self.idx_to_num[token_idx].get(ResidueMapKeys.TOKEN_NUM)
            chain_id = self.idx_to_num[token_idx].get(ResidueMapKeys.CHAIN_ID)
            atom_name = self.idx_to_num[token_idx].get(ResidueMapKeys.ATOM_NAME)

            if chain_id not in domain_dict:
                domain_dict[chain_id] = [(atom_name, token_num)]
            else:
                domain_dict[chain_id].append((atom_name, token_num))

        return domain_dict

    def _filter_by_plddt(
        self,
        domain_dict: dict,
        token_plddts: list,
    ) -> dict[str, list[tuple[str, int]]]:
        """Filter the residues in the pseudo-domains based on the pLDDT cutoff.

        Only keep the residues with pLDDT >= cutoff in the `domain_dict`.
        Different cutoffs can be used for IDR and non-IDR chains.\n
        *(See :py:attr:`plddt_cutoff_idr` and :py:attr:`plddt_cutoff`.)*

        ## Arguments:

        - **domain_dict (dict)**:<br />
            Dictionary of pseudo-domains.

        ## Returns:

        - **rb_dict (dict)**:<br />
            pLDDT filtered dictionary of rigid bodies.
        """

        self.check_is_set_up()
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

            chain_plddt_mask = np.array(token_plddts)[
                chain_token_idxs
            ] >= plddt_cutoff

            confident_tokens = [
                (
                    self.idx_to_num[token_idx].get(ResidueMapKeys.ATOM_NAME),
                    self.idx_to_num[token_idx].get(ResidueMapKeys.TOKEN_NUM)
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
        min_res: int,
        min_proteins: int,
    ) -> dict:
        """Filter the domain based on the size of the domain.

        Only keep the domain if it exceeds certain size.
        The size is determined by `min_res` or `min_proteins`.
        ```python
        - min_res # Minimum number of residues in a rigid body.
        - min_proteins # Minimum number of proteins in a rigid body.
        ```

        ## Arguments:

        - **rb_dict (dict)**:<br />
            Dictionary of a rigid body.

        - **min_res (int)**:<br />
            Minimum number of residues in a rigid body.

        - **min_proteins (int)**:<br />
            Minimum number of proteins in a rigid body.

        ## Returns:

        - **rb_dict (dict)**:<br />
            Filtered dictionary of a rigid body.
        """

        if len(rb_dict) < min_proteins:
            return {}

        total_residues = 0
        for _chain_id, chain_atom_res_list in rb_dict.items():
            chain_res_set = set()
            for _atom_name, res_num in chain_atom_res_list:
                chain_res_set.add(res_num)
            total_residues += len(chain_res_set)

        if total_residues < min_res:
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

        ## Arguments:

        - **rigid_bodies (list)**:<br />
            List of rigid bodies, where each rigid body is a dictionary
            with chain IDs as keys and a list of tuples containing
            atom names and residue numbers as values.

        ## Returns:

        - **rigid_bodies (list)**:<br />
            List of rigid bodies, where each rigid body is a dictionary
            with chain IDs as keys and a list of residue numbers as values.
        """

        if len(rigid_bodies) == 0:

            return rigid_bodies

        elif isinstance(list(rigid_bodies[0].values())[0][0], int): #! this is a weak check

            return rigid_bodies

        rigid_bodies = copy.deepcopy(rigid_bodies)

        # Convert the rigid body dictionary to a list of residue numbers
        for idx, rb_dict in enumerate(rigid_bodies):

            for chain_id, chain_res_num_list in rb_dict.items():

                only_res_num_list = []
                for _atom_name, res_num in chain_res_num_list:
                    only_res_num_list.append(res_num)

                only_res_num_list.sort()

                rb_dict[chain_id] = list(set(only_res_num_list))

            rigid_bodies[idx] = rb_dict

        return rigid_bodies

    def save_rigid_bodies(
        self,
        domains: list,
        output_dir: str,
        rb_out_fmt: str = RBCons.rb_out_fmt,
        save_structure: bool = RBCons.save_structure,
        rb_struct_fmt: str = RBCons.rb_struct_fmt,
        filter_struct_by_plddt: bool = RBCons.filter_struct_by_plddt,
        protein_chain_map: dict = {},
    ):
        """ Save the rigid bodies to a file and/or save the structure of the
        rigid bodies and assess the rigid bodies.

        Output options:
        - The rigid bodies are saved in a plain text or JSON format with the
          chain IDs and residue numbers.
        - The structure of the rigid bodies can be saved in PDB or CIF format.
          For rigid bodies with modifications, it is recommended to use PDB
          format.
        - The PAE plot can be saved to visualize the rigid bodies in the PAE
          matrix. Set `pae_plot=True` to save the PAE plot.
        - The rigid bodies can be assessed based on the interface residues,
          number of contacts, interface PAE and pLDDT, average PAE and plDDT
          and minimum PAE. Set `rb_assessment` dictionary with the following
          keys:
            - `as_average`:<br />
                Whether to report only the average of assessment metric to the
                output file.<br />
            - `symmetric_pae`:<br />
                Whether to report a single average PAE value or assymetric PAE
                value for PAE assessment metrics.
        - The assessment is saved in an Excel file.

        ## Arguments:

        - **domains (list)**:<br />
            list of rigid bodies, where each rigid body is a dictionary with
            chain IDs as keys and residue numbers as values.

        - **output_dir (str)**:<br />
            Directory to save the output files.

        - **rb_out_fmt (str, optional)**:<br />
            Defaults to "txt". ("txt" or "csv")

        - **save_structure (bool, optional)**:<br />
            Whether to save the structure of the rigid bodies.

        - **rb_struct_fmt (str, optional)**:<br />
            File type to save the structure ("pdb" or "cif").

        - **filter_struct_by_plddt (bool, optional)**:<br />
            Whether to save the structure without filtering based on pLDDT.

        - **pae_plot (bool, optional)**:<br />
            Whether to save the PAE plot for the rigid bodies.

        - **rb_assessment (dict | None, optional)**:<br />
            Dictionary containing parameters for rigid body assessment.<br />
            Parameters for rigid body assessment:<br />
            - `as_average`:<br />
            Whether to report only the average of assessment metric to the
            output file.<br />
            - `symmetric_pae`:<br />
            Whether to report a single average PAE value or assymetric PAE
            value for PAE assessment metrics. \n

        - **protein_chain_map (dict | None, optional)**:<br />
            Protein-to-chain mapping dictionary.
        """

        self.check_is_set_up()
        os.makedirs(output_dir, exist_ok=True)

        dir_name = os.path.basename(
            os.path.dirname(os.path.dirname(self.structure_file_path))
        )

        file_name = f"{dir_name}_rigid_bodies"

        protein_chain_map = extract_protein_chain_mapping(
            protein_chain_mapping=protein_chain_map
        )

        ##################################################

        if rb_out_fmt not in RBCons.valid_rb_out_fmts:
            raise ValueError(
                "Invalid rigid body output format specified."\
                f"Use one of {RBCons.valid_rb_out_fmts}."
            )

        # txt or json output format
        if rb_out_fmt == FileFormat.TXT:

            domains = self._keep_residue_numbers_only(domains)

            save_rigid_bodies_txt(
                output_dir=output_dir,
                domains=domains,
                protein_chain_map=protein_chain_map,
                file_name=file_name,
            )

        elif rb_out_fmt == FileFormat.JSON:

            save_rigid_bodies_json(
                output_dir=output_dir,
                domains=domains,
                protein_chain_map=protein_chain_map,
                file_name=file_name,
            )

        if save_structure:

            domains = self._keep_residue_numbers_only(domains)

            if rb_struct_fmt not in RBCons.valid_rb_struct_fmts:
                raise ValueError(
                    "Invalid rigid body structure output format specified."\
                    f"Use one of {RBCons.valid_rb_struct_fmts}."
                )

            if (
                rb_struct_fmt == FileFormat.CIF and
                has_modifications(self.structure)
            ):
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
                if filter_struct_by_plddt is False:
                    for chain_id, res_list in rb_dict.items():
                        if len(res_list) > 0:
                            res_list = fill_up_the_blanks(res_list)
                            rb_dict[chain_id] = res_list

                output_path = os.path.join(
                    output_dir,
                    f"{RBCons.rb_name_template.substitute(rb_idx=idx)}.{rb_struct_fmt}"
                )

                save_structure_obj(
                    structure=structure,
                    out_file=output_path,
                    res_select_obj=ResidueSelect(rb_dict),
                    save_type=rb_struct_fmt,
                    preserve_header_footer=True,
                )

    def show_rigid_bodies_on_pae_matrix(
        self,
        domains: list[dict],
        output_dir: str,
    ):
        for rb_idx, rb_dict in enumerate(domains):

            pae_patches = PAEPatches(
                num_to_idx=self.num_to_idx,
                pae=self.pae,
                lengths_dict=self.lengths_dict,
                rb_idx=rb_idx,
                af_offset=self.af_offset,
            )

            # patches are the highlighted rectangles in the PAE matrix
            patches = pae_patches.extract_pae_patches(rb_dict=rb_dict)
            pae_patches.plot_pae_patches(
                patches=patches,
                output_dir=output_dir,
            )

    def assess_rigid_bodies(
        self,
        domains: list[dict],
        output_dir: str,
        protein_chain_map: dict = {},
        symmetric_pae: bool = True,
        as_average: bool = True,
    ):

        for rb_idx, rb_dict in enumerate(domains):

            rb_save_path = os.path.join(
                output_dir,
                f"{RBCons.rb_assessment_name_template.substitute(rb_idx=rb_idx)}.{FileFormat.XLSX}"
            )

            rb_assess = RigidBodyAssessment(
                rb_dict=rb_dict,
                symmetric_pae=symmetric_pae,
                as_average=as_average,
                idr_chains=self.idr_chains,
                protein_chain_map=protein_chain_map,
            )

            rb_assess.set_attributes_from(instance=self)
            rb_assess.perform_assessment()
            rb_assess.save_rb_assessment(
                rb_c_assess=rb_assess.rb_c_assess,
                rb_cp_assess=rb_assess.rb_cp_assess,
                overall_assessment=rb_assess.overall_assessment,
                save_path=rb_save_path,
            )


class PAEPatches:
    """ Class to extract and plot PAE patches for rigid bodies."""

    num_to_idx: dict
    """ Dictionary mapping token numbers to token indices."""

    pae: np.ndarray
    """ PAE matrix."""

    lengths_dict: dict
    """ Dictionary containing the chain lengths and total length."""

    rb_idx: int
    """ Rigid body index."""

    af_offset: dict | None
    """ Offset describing start and end residue number for each chain in
    the predicted structure.\n
    example: `{'A': [1, 100], 'B': [101, 200]}`."""

    def __init__(
        self,
        num_to_idx: dict,
        pae: np.ndarray,
        lengths_dict: dict,
        rb_idx: int,
        af_offset: dict | None = None,
    ):

        self.num_to_idx = num_to_idx
        self.pae = pae
        self.lengths_dict = lengths_dict
        self.af_offset = af_offset
        self.rb_idx = rb_idx

    def extract_pae_patches(self, rb_dict: dict) -> list[list]:
        """ Extract PAE patches for the rigid body.

        ## Arguments:

        - **rb_dict (dict)**:<br />
            Dictionary of a rigid body with following `key`:`value` pair:
            `chain_id`: `(atom_name, res_num)`.

        ## Returns:

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

                if (
                    f"{RES_SEPARATOR}" not in res_idx_1 or
                    f"{RES_SEPARATOR}" not in res_idx_2
                ):
                    continue
                res1_y0 = int(res_idx_1.split(RES_SEPARATOR)[0])
                res1_y1 = int(res_idx_1.split(RES_SEPARATOR)[1])
                res2_x0 = int(res_idx_2.split(RES_SEPARATOR)[0])
                res2_x1 = int(res_idx_2.split(RES_SEPARATOR)[1])

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
        """ Show the PAE patches for the rigid body on the PAE matrix plot.

        ## Arguments:

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
            if chain_id == MiscStrEnum.TOTAL:
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

    This function writes the rigid bodies information to a text file in a
    human-readable format.<br />
    The output file will contain the rigid body index, chain ID, protein name
    (if available), and the residue range.

    ## Arguments:

    - **output_dir (str)**:<br />
        Directory where the output file will be saved.

    - **domains (list)**:<br />
        List of dictionaries, where each dictionary represents a rigid body.

    - **protein_chain_map (dict)**:<br />
        A mapping of chain IDs to protein names.

    - **file_name (str, optional)**:<br />
        Name of the output file without extension.
    """

    file_name += ".txt"
    output_path = os.path.join(output_dir, file_name)

    with open(output_path, "w") as f:

        for idx, rb_dict in enumerate(domains):
            f.write(f"Rigid Body {idx}\n")

            for chain_id, res_list in rb_dict.items():

                protein_name = protein_chain_map.get(chain_id, None)
                res_range = get_key_from_res_range(res_range=res_list)

                if len(res_list) > 0:
                    if protein_name:
                        f.write(f"{protein_name}_{chain_id}: {res_range}\n")
                    else:
                        f.write(f"{chain_id}:{res_range}\n")

            f.write("\n")

def save_rigid_bodies_json(
    output_dir: str,
    domains: list,
    protein_chain_map: dict,
    file_name: str = "rigid_bodies",
):
    """ Save rigid bodies to a JSON file.

    This function writes the rigid bodies information to a JSON file.<br />
    For per-atom tokens JSON format is recommended over text format.

    ## Arguments:

    - **output_dir (str)**:<br />
        Directory where the output file will be saved.

    - **domains (list)**:<br />
        List of dictionaries, where each dictionary represents a rigid body.

    - **protein_chain_map (dict)**:<br />
        A mapping of chain IDs to protein names.

    - **file_name (str, optional)**:<br />
        Name of the output file without extension.
    """

    file_name += FileFormat.JSON
    output_path = os.path.join(output_dir, file_name)

    rigid_bodies = []
    for idx, rb_dict in enumerate(domains):
        ch_dict = {}
        for chain_id, res_num_list in rb_dict.items():
            protein_name = protein_chain_map.get(chain_id, None)
            if not protein_name:
                protein_name = MiscStrEnum.UNKNOWN
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