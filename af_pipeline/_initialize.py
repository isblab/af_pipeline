from typing import Dict
from af_pipeline.parser1.structure_parser1 import StructureParser
from af_pipeline.parser1.data_parser import DataParser
import Bio.PDB
import Bio.PDB.Structure
import numpy as np
from af_pipeline.utils.obj_helpers import (
    get_duplicate_indices,
    update_matrix_row_col,
    symmetrize_matrix,
    create_mask,
)
from af_pipeline.tools.structure_tools import RenumberResidues

# rep_atom_dict is "res_name": "rep_atom_id"
# average_atom_pae will only take effect if state is "per_residue"
# average_atom_pae supersedes rep_atom_dict

class _Initialize:

    def __init__(
        self,
        data_file_path: str,
        structure_file_path: str | None = None,
        af_offset: dict | None = None,
        rep_atom_dict: dict | None = None,
        average_atom_pae: bool = True,
        state: str = "per_token",
    ):

        self.structure_file_path = structure_file_path
        self.data_file_path = data_file_path
        self.af_offset = af_offset
        self.rep_atom_dict = rep_atom_dict
        self.average_atom_pae = average_atom_pae
        self.state = state

        self.structure_parser = StructureParser(
            struct_file_path=self.structure_file_path,
            preserve_header_footer=False,
            which_parser="biopython",
        )

        self.data_parser = DataParser(
            data_file_path=self.data_file_path,
        )

        self.renumber = RenumberResidues(
            af_offset=self.af_offset
        )

        self.get_attributes(state=state)

    def get_attributes(self, state: str):

        data = self.data_parser.get_data_dict()
        structure = self.structure_parser.get_structure(
            parser=self.structure_parser.get_parser()
        )

        # get pae matrix
        self.pae = data["pae"]

        self.contact_probs = data["contact_probs"]

        if state == "per_token":

            self.token_chain_ids = self.structure_parser.get_token_chain_ids(
                structure=structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=False,
            )
            self.token_res_ids = self.structure_parser.get_token_res_ids(
                structure=structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=False,
            )

            self.token_plddts = self.structure_parser.get_plddt(
                structure=structure,
                per_atom=False,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=False,
            )

            self.token_coords = self.structure_parser.get_coordinates(
                structure=structure,
                per_atom=False,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=False,
           )

        elif state == "per_residue":

            self.token_chain_ids = self.structure_parser.get_token_chain_ids(
                structure=structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=True,
            )

            self.token_res_ids = self.structure_parser.get_token_res_ids(
                structure=structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=True,
            )

            self.token_plddts = self.structure_parser.get_plddt(
                structure=structure,
                per_atom=False,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=True,
            )

            self.token_coords = self.structure_parser.get_coordinates(
                structure=structure,
                per_atom=False,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=True,
            )

            self.idxs_to_keep = self.get_idxs_to_keep(
                structure=structure,
                rep_atom_dict=self.rep_atom_dict,
            )

            self.pae = self.update_pae(
                pae=self.pae,
                token_res_ids=self.token_res_ids,
                token_chain_ids=self.token_chain_ids,
                average_atom_pae=self.average_atom_pae,
                idxs_to_keep=self.idxs_to_keep,
            )

            self.contact_probs = self.update_contact_probs(
                contact_probs_mat=self.contact_probs,
                token_chain_ids=self.token_chain_ids,
                token_res_ids=self.token_res_ids,
                idxs_to_keep=self.idxs_to_keep,
            )

        else:
            raise Exception(
                "State should be either 'per_token' or 'per_residue'."
            )

        self.lengths_dict = self.get_chain_lengths(
            token_chain_ids=self.token_chain_ids
        )

        self.idx_to_num, self.num_to_idx = self.renumber.residue_map(
            token_chain_ids=self.token_chain_ids,
            token_res_ids=self.token_res_ids,
        )

    @staticmethod
    def get_chain_lengths(token_chain_ids: list) -> Dict[str, int]:
        """Get the chain lengths.

        lengths_dict is a dictionary containing the chain lengths. \n
        {chain_id: length} \n
        "total" is the total length of the system. \n
        For example, if the system has 2 chains A and B, \n
        lengths_dict = {"A": 100, "B": 50, "total": 150} \n

        Args:

            token_chain_ids (list):
                tokenized chain IDs.

        Returns:

            lengths_dict (Dict):
                dict containing the chain lengths.

        Example:

            >>> token_chain_ids = ['A', 'A', 'B', 'B', 'B']
            >>> lengths_dict = AfParser.get_chain_lengths(token_chain_ids)
            >>> print(lengths_dict)
            {'total': 5, 'A': 2, 'B': 3}
        """

        lengths_dict = {}
        lengths_dict["total"] = 0

        for chain_id in token_chain_ids:
            if chain_id not in lengths_dict:
                lengths_dict[chain_id] = 1
            else:
                lengths_dict[chain_id] += 1
            lengths_dict["total"] += 1

        return lengths_dict

    @staticmethod
    def get_idxs_to_keep(
        structure: Bio.PDB.Structure.Structure,
        rep_atom_dict: dict = {},
    ) -> dict:
        """Get the indices to keep for the PAE matrix.

        Args:

            token_chain_ids (list):
                tokenized chain IDs.

            token_res_ids (list):
                tokenized residue IDs.

            rep_atom_dict (dict, optional):
                dictionary with representative atom IDs.
                Defaults to {}.

        Returns:

            idxs_to_keep (dict):
                dictionary with indices to keep.
        """

        idxs_to_keep = {}

        for residue, ch_id in StructureParser.get_residues(structure):
            rep_atom = rep_atom_dict.get(
                residue.get_resname(),
                StructureParser.get_rep_atom(residue=residue)
            )
            quants = StructureParser.extract_perresidue_quantities(
                residue=residue,
                quantities=["res_pos", "chain_id", "res_name", "rep_atom", "rep_atom_local_idx"],
                rep_atom=rep_atom,
            )
            key = (quants["chain_id"], quants["res_pos"])

            if key not in idxs_to_keep:
                idxs_to_keep[key] = quants["rep_atom_local_idx"]

        return idxs_to_keep

    @staticmethod
    def update_pae(
        pae: np.ndarray,
        token_res_ids: list,
        token_chain_ids: list,
        average_atom_pae: bool = False,
        idxs_to_keep: dict = {},
    ) -> np.ndarray:
        """Update the PAE matrix based on the keyword.

        If average_atom_pae is set to True, the repeated residue
        IDs are removed. \n
        PAE values for the repeated residue IDs are replaced with
        the mean of the PAE values. \n

        Args:

            pae (np.ndarray):
                PAE matrix.

            token_res_ids (list):
                tokenized residue IDs.

            token_chain_ids (list):
                tokenized chain IDs.

            **average_atom_pae (bool, optional):
                If True, the repeated residue IDs are removed. \n
                Defaults to False.

        Returns:

            pae (np.ndarray):
                updated PAE matrix.
        """

        token_ids = list(zip(token_chain_ids, token_res_ids))

        dup_token_indices = get_duplicate_indices(
            my_li=token_ids,
            return_type="dict",
            keep_which=None,  # Keep all duplicates
        )

        pae = update_matrix_row_col(
            matrix=pae,
            idxs_to_update=dup_token_indices,
            replace_with_avg=average_atom_pae,
            idxs_to_keep=idxs_to_keep,
        )
        return pae

    @staticmethod
    def update_contact_probs(
        contact_probs_mat: np.ndarray,
        token_chain_ids: list,
        token_res_ids: list,
        idxs_to_keep: dict = {},
    ) -> np.ndarray:
        """Update the contact probabilities matrix based on the keyword.

        If average_atom_pae is set to True, the repeated residue
        IDs are removed. \n
        Contact probabilities for the repeated residue IDs are replaced with
        the mean of the contact probabilities. \n

        Args:

            contact_probs_mat (np.ndarray):
                Contact probabilities matrix.

            token_chain_ids (list):
                tokenized chain IDs.

            token_res_ids (list):
                tokenized residue IDs.

            **average_atom_pae (bool, optional):
                If True, the repeated residue IDs are removed. \n
                Defaults to False.

        Returns:

            contact_probs_mat (np.ndarray):
                updated contact probabilities matrix.
        """

        token_ids = list(zip(token_chain_ids, token_res_ids))

        dup_token_indices = get_duplicate_indices(
            my_li=token_ids,
            return_type="dict",
            keep_which=None,  # Keep all duplicates
        )

        contact_probs_mat = update_matrix_row_col(
            matrix=contact_probs_mat,
            idxs_to_update=dup_token_indices,
            replace_with_avg=False,
            idxs_to_keep=idxs_to_keep,
        )
        return contact_probs_mat

    @staticmethod
    def get_min_pae(
        pae_matrix: np.ndarray,
        lengths_dict: Dict,
        along_axis: int | None = None,
        hide_interactions: str = "intrachain",
        return_type: str = "array",
    ) -> np.ndarray | Dict[str, list] | list:
        """ Per-residue minimum PAE values.

        Given the PAE matrix, obtain minimum PAE values for each residue. \n

        If `hide_interactions=="intrachain"`, only the interchain interactions
        are considered. \n
        If `hide_interactions=="interchain"`, only the intrachain interactions
        are considered.

        If `return_type==dict`, a dictionary containing the min PAE values for
        each chain is returned. \n
        If `return_type=="array"` or `return_type=="list"`, a numpy array or a
        list containing the min PAE values for all residues is returned
        respectively. \n

        Args:

            pae_matrix (np.ndarray):
                Average PAE matrix.

            lengths_dict (Dict):
                Dictionary containing the chain lengths.

            along_axis (int | None):
                Axis along which to get the min PAE values.
                If None, average PAE is calculated and along_axis is set to 1.

            hide_interactions (str):
                Hide intrachain or interchain interactions.

            return_type (str):
                Whether to return min_pae as dict or list or array.

        Returns:

            min_pae_dict (Dict):
                Dictionary containing the min PAE values for each chain.

            min_pae (np.ndarray):
                Minimum PAE values for all residues in a 2D numpy array.

        Examples:

            >>> pae_matrix = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
            >>> lengths_dict = {"A": 2, "B": 1, "total": 3}
            >>> min_pae = AfParser.get_min_pae(pae_matrix, lengths_dict)
            >>> print(min_pae)
            [0 1 2]
            >>>
            >>> min_pae = AfParser.get_min_pae(
            ... pae_matrix, lengths_dict, along_axis=0
            ... )
            >>> print(min_pae)
            [0 1 0]
            >>>
            >>> min_pae = AfParser.get_min_pae(
            ... pae_matrix, lengths_dict, return_type="dict"
            ... )
            >>> print(min_pae)
            {'A': [0, 1], 'B': [2]}
        """

        if along_axis is None:
            pae_matrix = symmetrize_matrix(matrix=pae_matrix)
            along_axis = 1

        hide_keys = {
            "intrachain": "intra_part",
            "interchain": "inter_part",
        }

        interchain_mask = create_mask(
            partition_dict=lengths_dict,
            hide_interactions=hide_keys[hide_interactions],
            masked_value=1,
            unmasked_value=0,
        )

        masked_pae_matrix = np.ma.masked_array(
            pae_matrix, mask=interchain_mask
        )

        min_pae = np.min(masked_pae_matrix, axis=along_axis)

        if return_type == "array":
            return min_pae.data

        elif return_type == "list":
            return min_pae.tolist()

        elif return_type == "dict":

            min_pae_dict = {}
            start = 0
            min_pae = min_pae.tolist()

            for chain_id in lengths_dict:
                if chain_id != "total":

                    end = start + lengths_dict[chain_id]
                    min_pae_dict[chain_id] = min_pae[start:end]
                    start = end
            return min_pae_dict

        else:
            raise Exception(
                "return_type should be either 'array', 'list' or 'dict'."
            )
