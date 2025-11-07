"""
Rigid Bodies extraction module
==============================
RigidBodies class with methods to extract rigid bodies from AlphaFold predictions.
"""
import os
import copy
import warnings
import time
import numpy as np
from collections import defaultdict
from Bio.PDB.Structure import Structure
from af_pipeline._initialize import _Initialize
from af_pipeline.tools.structure_tools import (
    save_structure_obj,
    ResidueSelect
)
from af_pipeline.pae_to_domains.pae_to_domains import (
    domains_from_pae_matrix_igraph,
    domains_from_pae_matrix_networkx,
    domains_from_pae_matrix_label_propagation
)
from af_pipeline.rigid_bodies.output_rigid_bodies import (
    save_rigid_bodies_txt,
    save_rigid_bodies_json
)
from af_pipeline.tools.misc_tools import (
    fill_up_the_blanks,
)

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

    state: str
    """ State of the instance. Can be "per_token" or "per_residue"."""

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
        state: str = "per_token",
    ):

        super().__init__(
            data_file_path=data_file_path,
            structure_file_path=structure_file_path,
            af_offset=af_offset,
            rep_atom_dict=rep_atom_dict,
            average_token_pae=average_token_pae,
            average_token_plddt=average_token_plddt,
            state=state,
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
            1. igraph based
            2. networkx based
            3. label_propagation based

        Args:
            num_res (int): Minimum number of residues in a rigid body
            num_proteins (int): Minimum number of proteins in a rigid body
            plddt_filter (bool): Filter the residues based on the pLDDT cutoff

        Returns:
            `rigid_bodies (list)`: List of extracted rigid bodies
        """

        print("Extracting rigid bodies...")
        start_time = time.time()

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

        for rb in pseudo_domains:

            # rb_dict is a dictionary of rigid bodies
            # each rigid body is represented as a dictionary with chain_id as
            # the key and a list of residue numbers as the value
            rb_dict = self.rb_to_rb_dict(rb=rb)

            # removing residues with pLDDT score below the cutoff
            if plddt_filter:
                rb_dict = self.filter_by_plddt(rb_dict=rb_dict)

            rigid_bodies.append(rb_dict)

        # Remove domains with number of proteins less than `num_proteins`
        rigid_bodies = [
            rb_dict
            for rb_dict in rigid_bodies
            if len(rb_dict) >= num_proteins
        ]

        # Remove domains with number of residues less than `num_res`
        rigid_bodies = [
            rb_dict
            for rb_dict in rigid_bodies
            if sum([len(res_list) for res_list in rb_dict.values()]) >= num_res
        ]

        end_time = time.time()

        print(
            f"Done extracting rigid bodies in \
            {end_time - start_time:.2f} seconds"
        )

        return rigid_bodies

    def rb_to_rb_dict(
        self,
        rb: list
    ) -> dict[str, list[tuple[str, int]]]:
        """Convert the domain list to a dictionary of rigid bodies.

        Args:
            rb (list): List of token indices in a rigid body.

        Returns:
            `rb_dict (dict)`:
                Dictionary of rigid bodies with chain IDs as keys and
                a list of tuples containing atom names and residue numbers
                as values.
        """

        rb_dict = defaultdict(list)

        for token_idx in rb:

            res_num = self.idx_to_num[token_idx].get("res_num")
            chain_id = self.idx_to_num[token_idx].get("chain_id")
            atom_name = self.idx_to_num[token_idx].get("atom_name")

            if chain_id not in rb_dict:
                rb_dict[chain_id] = [(atom_name, res_num)]
            else:
                rb_dict[chain_id].append((atom_name, res_num))

        return rb_dict

    def filter_by_plddt(
        self,
        rb_dict: dict,
    ) -> dict[str, list[tuple[str, int]]]:
        """Filter the residues in the rigid bodies based on the pLDDT cutoff.

        If the pLDDT score of a residue is less than the cutoff, it is removed
        from the rigid body.

        Args:
            rb_dict (dict): Dictionary of rigid bodies

        Returns:
        - **rb_dict (dict)**: pLDDT filtered dictionary of rigid bodies.
        """

        # Filter the residues in each chain in the rigid body based on the pLDDT cutoff
        for rb_ch_id, rb_ch_res_num_list in rb_dict.items():

            confident_residues = []

            plddt_res_num_arr = np.array([
                self.num_to_idx[rb_ch_id][res_num][atom_name]
                for atom_name, res_num in rb_ch_res_num_list
            ])

            if rb_ch_id in self.idr_chains:
                tf_plddt_filtered = np.array(self.token_plddts)[
                    plddt_res_num_arr
                ] >= self.plddt_cutoff_idr
            else:
                tf_plddt_filtered = np.array(self.token_plddts)[
                    plddt_res_num_arr
                ] >= self.plddt_cutoff

            confident_residues = plddt_res_num_arr[tf_plddt_filtered]
            confident_residues = [
                (self.idx_to_num[token_idx].get("atom_name"),
                 self.idx_to_num[token_idx].get("res_num"))
                for token_idx in confident_residues
            ]

            # Update the rigid body dictionary with the confident residues
            rb_dict[rb_ch_id] = confident_residues

        # Remove chains which have no confident residues
        empty_chains = []

        for chain_id, confident_residues in rb_dict.items():
            if len(confident_residues) == 0:
                empty_chains.append(chain_id)

        for chain_id in empty_chains:
            del rb_dict[chain_id]

        return rb_dict

    @staticmethod
    def keep_residue_numbers_only(
        rigid_bodies: list[dict[str, list[tuple[str, int]]]] | list
    ) -> list[dict[str, list[int]]] | list:
        """ Convert the rigid bodies to a list of residue numbers only.

        By default, the rigid body is in the following format:
        `{"A": [("CA", 1), ("CB", 2)]}`
        where the key is the chain ID and the value is a list of tuples
        containing the atom name and residue number.

        This function converts the rigid body to a list of residue numbers only:
        `{"A": [1, 2]}`

        Args:
            rigid_bodies (list[dict[str, list[tuple[str, int]]]] | list):
                List of rigid bodies, where each rigid body is a dictionary
                with chain IDs as keys and a list of tuples containing
                atom names and residue numbers as values.

        Returns:
        - **rigid_bodies (list[dict[str, list[int]]] | list)**:
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

            for chain_id, res_num_list in rb_dict.items():

                rb_list = []
                for atom_name, res_num in res_num_list:
                    rb_list.append(res_num)

                # Sort the list of residue numbers
                rb_list.sort()

                rb_dict[chain_id] = list(set(rb_list))
            # Update the rigid body dictionary with the sorted list
            rigid_bodies[idx] = rb_dict

        return rigid_bodies

    def extract_protein_chain_mapping(
        self,
        protein_chain_mapping: dict | None = None
    ) -> dict[str, str]:
        """ Extract the protein chain mapping from the provided dictionary.

        For e.g., if the user provides the following mapping:
        ```python
        {
            "ProteinA": "A,B",
            "ProteinB": "C"
        }
        ```
        The function will return the following dictionary:
        ```python
        {
            "A": "ProteinA",
            "B": "ProteinA",
            "C": "ProteinB"
        }
        ```

        Args:
            protein_chain_mapping (dict): Protein-to-chain map.

        Returns:
        - **protein_chain_map (dict)**:
            Dictionary with chain IDs as keys and protein names as values.
        """

        protein_chain_map = {}

        if protein_chain_mapping is None:
            return protein_chain_map

        for p_c_maps in protein_chain_mapping:
            protein_name, chain_ids = p_c_maps.split(":")
            chain_ids = chain_ids.split(",")
            for chain_id in chain_ids:
                if chain_id not in protein_chain_map:
                    protein_chain_map[chain_id] = protein_name

        return protein_chain_map

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

        Args:
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

        file_name = (
            os.path.basename(self.structure_file_path).split(".")[0] + "_rigid_bodies"
        )

        protein_chain_map = self.extract_protein_chain_mapping(
            protein_chain_mapping=protein_chain_map
        )

        ##################################################

        # txt or json output format
        if output_format == "txt":

            domains = self.keep_residue_numbers_only(domains)

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

            domains = self.keep_residue_numbers_only(domains)

            if structure_file_type == "cif":
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
            if not isinstance(self.structure, Structure):
                raise TypeError(
                    f"""
                    The structure should be a Bio.PDB.Structure.Structure
                    object. Got {type(self.structure)} instead.
                    """
                )


            structure = self.renumber.renumber_structure(
                structure=self.structure
            )

            for idx, rb_dict in enumerate(domains):

                # In the following case, the txt or json ouput will have pLDDT filtered residues
                # but, the structure file will ignore this filter
                # use this flag when you don't want missing residues in the structure file

                if no_plddt_filter_for_structure:
                    for chain_id, res_list in rb_dict.items():
                        if len(res_list) > 0:
                            res_list = fill_up_the_blanks(res_list)
                            rb_dict[chain_id] = res_list

                output_path = os.path.join(output_dir, f"rigid_body_{idx}.{structure_file_type}")

                save_structure_obj(
                    structure=structure,
                    out_file=output_path,
                    res_select_obj=ResidueSelect(rb_dict),
                    save_type=structure_file_type,
                    preserve_header_footer=True,
                )

        # Save the PAE plot for the rigid bodies
        # the region of the PAE matrix corresponding to the rigid bodies will be highlighted
        if pae_plot:

            raise NotImplementedError(
                "PAE plot for rigid bodies is not implemented yet."
            )

        if rb_assessment is not None:

            raise NotImplementedError(
                "Rigid body assessment is not implemented yet."
            )