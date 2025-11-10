"""
Data Parser module
==================
DataParser class with methods to parse the JSON or PKL data files

Currently supports the following:
- AlphaFold2 data files (.pkl or .json)
- AlphaFold3 data files (.json)
- ColabFold data files (.json)

To be implemented:
- Boltz data files (.npy)
"""
import warnings
import numpy as np
import os
import json
import pickle as pkl
from typing import Dict

class DataParser:
    """Class with methods to parse the prediction data files."""

    data_file_path: str
    """ Path to the data file."""

    def __init__(
        self,
        data_file_path: str
    ):
        self.data_file_path = data_file_path

    def get_data_dict(self) -> Dict:
        """Get the data from the data file.

        Args:
            data_file_path (str): Path to the data file.

        Returns:
            `data (Dict)`: Data dictionary from the data file.
        """

        ext = os.path.splitext(self.data_file_path)[1]

        # AF2 data file
        if "pkl" in ext:
            with open(self.data_file_path, "rb") as f:
                data = pkl.load(f)

        # AF3 or ColabFold data file
        elif "json" in ext:
            with open(self.data_file_path, "r") as f:
                data = json.load(f)

            if isinstance(data, list):
                data = data[0]

        elif "npy" in ext:
            raise NotImplementedError(
                """

                Boltz predictions are not supported yet.
                """
            )

        else:
            raise Exception(
                """

                Incorrect file format.. Suported .pkl/.json only.
                """
            )

        return data

    @staticmethod
    def get_token_chain_ids(data: Dict) -> list | None:
        """Get the token chain IDs from the data dictionary.

        This is specific to AF3: `"token_chain_ids"` key.
        For others, `None` is returned.

        In general, each token is a residue/nucleotide/ion in the chain.
        In case of modified residues or nucleotides or ions or glycan chains,
        each atom is a token.

        Args:
            data (Dict): Data dictionary from the data file.

        Returns:
            `token_chain_ids (list)`: Token chain IDs.
        """

        if "token_chain_ids" in data:
            token_chain_ids = data["token_chain_ids"]

        else:
            warnings.warn(
                """

                Chain IDs not found, data file might be AF2.
                Structure file is required for AF2.
                """
            )
            token_chain_ids = None

        return token_chain_ids

    @staticmethod
    def get_token_res_ids(data: Dict) -> list | None:
        """Get the token residue IDs from the data dictionary.

        This is specific to AF3: `"token_res_ids"` key.
        For others, `None` is returned.

        Atom-level tokens have the same token residue IDs.

        Args:
            data (Dict): Data dictionary from the data file.

        Returns:
            `token_res_ids (list)`: Token residue IDs.
        """

        if "token_res_ids" in data:
            token_res_ids = data["token_res_ids"]

        else:
            warnings.warn(
                """

                Residue IDs not found, data file might be AF2.
                Structure file is required for AF2.
                """
            )
            token_res_ids = None

        return token_res_ids

    @staticmethod
    def get_pae(data: Dict) -> np.ndarray:
        """Return the PAE matrix from the data dictionary.

        Size of the PAE matrix is NxN, where N is the number of tokens.

        Args:
            data (Dict): Data dictionary from the data file.

        Returns:
            `pae (np.ndarray)`: PAE matrix.
        """

        # For AF2.
        if "predicted_aligned_error" in data:
            pae = np.array(data["predicted_aligned_error"])

        # For AF3.
        elif "pae" in data:
            pae = np.array(data["pae"])

        else:
            raise Exception("PAE matrix not found...")

        return pae

    @staticmethod
    def get_contact_probs_mat(data: Dict) -> np.ndarray | None:
        """Get the contact probabilities from the data dictionary.

        This is specific to AF3: `"contact_probs"` key.
        For others, `None` is returned.

        Args:
            data (Dict): Data dictionary from the data file.

        Returns:
            `contact_probs_mat (np.ndarray)`:
                Contact probabilities matrix from AlphaFold3 output.
        """

        if "contact_probs" in data:
            contact_probs_mat = np.array(data["contact_probs"])

        else:
            warnings.warn(
                """

                Contact probabilities not found, data file might not be AF3.
                """
            )
            contact_probs_mat = None

        return contact_probs_mat

    @staticmethod
    def get_atom_chain_ids(data: Dict) -> list | None:
        """Get per atom chain IDs from the data dictionary.

        This is specific to AF3: `"atom_chain_ids"` key.
        For others, `None` is returned.

        Args:
            data (Dict): Data dictionary from the data file.

        Returns:
            `atom_chain_ids (list)`: Per atom chain IDs.
        """

        if "atom_chain_ids" in data:
            atom_chain_ids = data["atom_chain_ids"]

        else:
            warnings.warn(
                """
                Per atom chain IDs not found, data file might be AF2.
                Structure file is required for AF2.
                """
            )
            atom_chain_ids = None

        return atom_chain_ids

    @staticmethod
    def get_atom_plddts(data: Dict) -> np.ndarray | None:
        """Get per atom pLDDT scores from the data dictionary.

        This is specific to AF3: `"atom_plddts"` key.
        If data is AF2, `None` is returned.

        Args:
            data (Dict): Data dictionary from the data file.

        Returns:
            `atom_plddts (np.ndarray)`: Per atom pLDDT scores.
        """

        if "atom_plddts" in data:
            atom_plddts = data["atom_plddts"]

        else:
            warnings.warn(
                """
                Per atom pLDDT scores not found, data file might be AF2.
                Structure file is required for AF2.
                """
            )
            atom_plddts = None

        return atom_plddts