import numpy as np
from typing import Dict
import Bio.PDB.Structure
from Bio.PDB.Structure import Structure
from af_pipeline.parser1.structure_parser1 import StructureParser
from af_pipeline.parser1.data_parser import DataParser
from af_pipeline.tools.structure_tools import RenumberResidues
from af_pipeline.tools.misc_tools import (
    get_duplicate_indices,
    update_matrix_row_col,
    symmetrize_matrix,
    create_mask,
)

#TODO: Put this at appropriate place
# rep_atom_dict is "res_name": "rep_atom_id"
# average_token_pae will only take effect if state is "per_residue"
# average_token_pae and average_token_plddt supersedes rep_atom_dict

class _Initialize:
    """

    Attributes:

        structure_file_path (str):
            Path to the structure file (PDB or CIF).

        data_file_path (str):
            Path to the data file (json or pkl).

        af_offset (dict):
            Dictionary containing the offset for the prediction
            Defaults to {}.

        rep_atom_dict (dict):
            Dictionary containing the representative atoms for residues.
            Defaults to {}.

        average_token_pae (bool):
            If True, average PAE values for residues with per-atom tokens.
            Defaults to True.

        average_token_plddt (bool):
            If True, average pLDDT values for residues with per-atom tokens.
            Defaults to True.

        state (str):
            State of the parser, either "per_token" or "per_residue".
            Defaults to "per_token".
    """

    def __init__(
        self,
        data_file_path: str,
        structure_file_path: str,
        af_offset: dict = {},
        rep_atom_dict: dict = {},
        average_token_pae: bool = True,
        average_token_plddt: bool = True,
        state: str = "per_token",
    ):
        """ Initialize the _Initialize class.

        Args:

            data_file_path (str):
                Path to the data file (json or pkl).

            structure_file_path (str):
                Path to the structure file (PDB or CIF).

            af_offset (dict, optional):
                Dictionary containing the offset for the prediction.
                Defaults to {}.

        rep_atom_dict (dict):
            Dictionary containing the representative atoms for residues.
            Defaults to {}.

        average_token_pae (bool):
            If True, average PAE values for residues with per-atom tokens.
            Defaults to True.

        average_token_plddt (bool):
            If True, average pLDDT values for residues with per-atom tokens.
            Defaults to True.

        state (str):
            State of the parser, either "per_token" or "per_residue".
            Defaults to "per_token".
        """

        self.structure_file_path = structure_file_path
        self.data_file_path = data_file_path
        self.af_offset = af_offset
        self.rep_atom_dict = rep_atom_dict
        self.average_token_pae = average_token_pae
        self.average_token_plddt = average_token_plddt
        self.state = state
        self.structure = None

        self.structure_parser = StructureParser(
            structure_file_path=self.structure_file_path,
            preserve_header_footer=False,
        )

        self.data_parser = DataParser(
            data_file_path=self.data_file_path,
        )

        self.get_attributes(state=self.state)

        self.avg_pae = symmetrize_matrix(matrix=self.pae)

        self.lengths_dict = self.get_chain_lengths(
            token_chain_ids=self.token_chain_ids
        )

        self.renumber = RenumberResidues(
            af_offset=self.af_offset
        )

        self.idx_to_num, self.num_to_idx = self.renumber.residue_map(
            token_chain_ids=self.token_chain_ids,
            token_res_ids=self.token_res_ids,
            token_atom_names=self.token_atom_names,
        )

    def get_attributes(self, state: str) -> None:
        """ Get the attributes of the class based on the state.

        ! Add description of the state here.

        Args:
            state (str):
                State of the parser, either "per_token" or "per_residue".
        """

        data = self.data_parser.get_data_dict()
        self.structure = self.structure_parser.get_structure(
            parser=self.structure_parser.get_parser()
        )

        self.pae = self.data_parser.get_pae(data)

        self.contact_probs = self.data_parser.get_contact_probs_mat(data)

        if not isinstance(self.structure, Structure):
            raise TypeError(
                f"""

                Structure should be a Bio.PDB.Structure.Structure object.
                Got {type(self.structure)} instead. \n
                """
            )

        if state == "per_token":

            self.only_representative = False

            self.token_chain_ids = self.structure_parser.get_token_chain_ids(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=self.only_representative,
            )

            self.token_res_ids = self.structure_parser.get_token_res_ids(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=self.only_representative,
            )

            self.token_plddts = self.structure_parser.get_plddt(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
                average_token_plddt=False,
                only_representative=self.only_representative,
            )

            self.token_coords = self.structure_parser.get_coordinates(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=self.only_representative,
           )

            self.token_atom_names = self.structure_parser.get_token_atom_names(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=self.only_representative,
            )

        elif state == "per_residue":

            self.only_representative = True

            self.token_chain_ids = self.structure_parser.get_token_chain_ids(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=self.only_representative,
            )

            self.token_res_ids = self.structure_parser.get_token_res_ids(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=self.only_representative,
            )

            self.token_plddts = self.structure_parser.get_plddt(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
                average_token_plddt=self.average_token_plddt,
                only_representative=self.only_representative,
            )

            self.token_coords = self.structure_parser.get_coordinates(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=self.only_representative,
            )

            self.token_atom_names = self.structure_parser.get_token_atom_names(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
                only_representative=self.only_representative,
            )

            self.idxs_to_keep = self.get_idxs_to_keep(
                structure=self.structure,
                rep_atom_dict=self.rep_atom_dict,
            )

            self.pae = self.update_pae(
                token_res_ids=self.data_parser.get_token_res_ids(data),
                token_chain_ids=self.data_parser.get_token_chain_ids(data),
            )

            self.contact_probs = self.update_contact_probs(
                token_chain_ids= self.data_parser.get_token_chain_ids(data),
                token_res_ids=self.data_parser.get_token_res_ids(data),
            )

        else:
            raise Exception(
                f"""

                State should be either 'per_token' or 'per_residue'.
                Got '{state}' instead. \n
                """
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
                Token chain IDs.

        Returns:

            lengths_dict (Dict):
                Dictionary containing the chain lengths and total length.

        Example:

            >>> token_chain_ids = ['A', 'A', 'B', 'B', 'B']
            >>> lengths_dict = _Initialize.get_chain_lengths(token_chain_ids)
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
        structure: Structure,
        rep_atom_dict: dict = {},
    ) -> Dict[tuple, int]:
        """Get the indices to keep in the PAE matrix.

        Args:

            token_chain_ids (list):
                Token chain IDs.

            token_res_ids (list):
                Token residue IDs.

            rep_atom_dict (dict, optional):
                Dictionary with residue names as keys and representative
                atoms as values. \n
                If only_representative is True, this dictionary is used to get
                the representative atom for the specified residue.

        Returns:

            idxs_to_keep (dict):
                Dictionary with indices to keep.
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

    def update_pae(
        self,
        token_res_ids: list | None,
        token_chain_ids: list | None,
    ) -> np.ndarray:
        """Update the PAE matrix based on the keyword.

        If average_token_pae is set to True, the repeated residue
        IDs are removed. \n
        PAE values for the repeated residue IDs are replaced with
        the mean of the PAE values. \n

        Args:

            pae (np.ndarray):
                PAE matrix.

            token_res_ids (list):
                Token residue IDs.

            token_chain_ids (list):
                Token chain IDs.

        Returns:

            pae (np.ndarray):
                Updated PAE matrix.
        """

        if token_chain_ids is None or token_res_ids is None:

            return self.pae

        token_ids = list(zip(token_chain_ids, token_res_ids))

        dup_token_indices = get_duplicate_indices(
            my_li=token_ids,
            return_type="dict",
            keep_which=None,  # Keep all duplicates
        )

        if not isinstance(dup_token_indices, dict):
            raise TypeError(
                f"""

                Duplicate indices should be a dictionary.
                Got {type(dup_token_indices)} instead. \n
                """
            )

        pae = update_matrix_row_col(
            matrix=self.pae,
            idxs_to_update=dup_token_indices,
            replace_with_avg=self.average_token_pae,
            idxs_to_keep=self.idxs_to_keep,
        )
        return pae

    def update_contact_probs(
        self,
        token_chain_ids: list | None,
        token_res_ids: list | None,
    ) -> np.ndarray | None:
        """Update the contact probabilities matrix based on the keyword.

        If average_token_pae is set to True, the repeated residue
        IDs are removed. \n
        Contact probabilities for the repeated residue IDs are replaced with
        the mean of the contact probabilities. \n

        Args:

            contact_probs_mat (np.ndarray):
                Contact probabilities matrix.

            token_chain_ids (list):
                Token chain IDs.

            token_res_ids (list):
                Token residue IDs.

        Returns:

            contact_probs_mat (np.ndarray):
                Updated contact probabilities matrix.
        """

        if token_chain_ids is None or token_res_ids is None:

            return self.contact_probs

        token_ids = list(zip(token_chain_ids, token_res_ids))

        dup_token_indices = get_duplicate_indices(
            my_li=token_ids,
            return_type="dict",
            keep_which=None,  # Keep all duplicates
        )

        if not isinstance(dup_token_indices, dict):
            raise TypeError(
                f"""

                Duplicate indices should be a dictionary.
                Got {type(dup_token_indices)} instead. \n
                """
            )

        contact_probs_mat = update_matrix_row_col(
            matrix=self.contact_probs,
            idxs_to_update=dup_token_indices,
            replace_with_avg=False,
            idxs_to_keep=self.idxs_to_keep,
        )
        return contact_probs_mat

    @staticmethod
    def get_min_pae(
        pae_matrix: np.ndarray,
        lengths_dict: Dict,
        along_axis: int | None = None,
        return_type: str = "array",
    ) -> np.ndarray | Dict[str, list] | list:
        """ Per-residue minimum PAE values.

        Given the PAE matrix, obtain minimum PAE values for each residue. \n

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

            return_type (str):
                Whether to return min_pae as dict or list or array.

        Returns:

            min_pae_dict (Dict):
                Dictionary containing the min PAE values for each chain.

            min_pae (np.ndarray):
                Minimum PAE values for all residues in a 2D numpy array.

        Examples:

            >>> pae_matrix = np.array([
            ... [0, 1, 1],
            ... [1, 0, 3],
            ... [2, 3, 0]
            ... ])
            >>> lengths_dict = {"A": 2, "B": 1, "total": 3}
            >>> _Initialize.get_min_pae(
            ...     pae_matrix,
            ...     lengths_dict
            ... )
            array([1.5, 3. , 1.5])

            >>> _Initialize.get_min_pae(
            ... pae_matrix,
            ... lengths_dict,
            ... along_axis=0,
            ... return_type="list"
            ... )
            [2, 3, 1]

            >>> _Initialize.get_min_pae(
            ... pae_matrix,
            ... lengths_dict,
            ... along_axis=1,
            ... return_type="dict"
            ... )
            {'A': [1, 3], 'B': [2]}
        """

        if along_axis is None:
            pae_matrix = symmetrize_matrix(matrix=pae_matrix)
            along_axis = 1

        interchain_mask = create_mask(
            partition_dict=lengths_dict,
            hide_interactions="intra_part",
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
                f"""

                return_type should be either 'array', 'list' or 'dict'.
                Got '{return_type}' instead.
                """
            )

if __name__ == "__main__":
    import doctest
    doctest.testmod()