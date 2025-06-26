import warnings
import numpy as np
import os
import json
import pickle as pkl
from typing import Dict

class DataParser:
    """Class containing methods to parse the AF2/3 data file.

    Attributes:

        data_file_path (str):
            Path to the AF2/3 data file.
    """

    def __init__(
        self,
        data_file_path: str
    ):
        """Initialize the DataParser class.

        Args:
            data_file_path (str):
                Path to the AF2/3 data file.
        """

        self.data_file_path = data_file_path

    def get_data_dict(self) -> Dict:
        """Parse the AF2/3 data file.

        Args:

            data_file_path (str):
                path to the data file.

        Returns:

            data (Dict):
                data dict from the data file.
        """

        ext = os.path.splitext(self.data_file_path)[1]

        # AF2 data file
        if "pkl" in ext:
            with open(self.data_file_path, "rb") as f:
                data = pkl.load(f)

        # AF3 data file
        elif "json" in ext:
            with open(self.data_file_path, "r") as f:
                data = json.load(f)

            if isinstance(data, list):
                data = data[0]

        else:
            raise Exception(
                "Incorrect file format.. Suported .pkl/.json only."
            )

        return data

    @staticmethod
    def get_token_chain_ids(data: Dict) -> list:
        """Get the token chain IDs from the data dict.

        This is specific to AF3: "token_chain_ids" key. \n
        If data is AF2, None is returned. \n

        In general, each token is a residue/nucleotide/ion in the chain.
        In case of modified residues or nucleotides or ions or glycan chains,
        each atom is a token.

        Args:

            data (Dict):
                data dictionary from the data file.

        Returns:

            token_chain_ids (list):
                token chain IDs.
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
    def get_token_res_ids(data: Dict) -> list:
        """Get the token residue IDs from the data dict.

        This is specific to AF3: "token_res_ids" key. \n
        If data is AF2, None is returned. \n

        Args:

            data (Dict):
                data dictionary from the data file.

        Returns:

            token_res_ids (list):
                token residue IDs.
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
    def get_pae(data: Dict):
        """Return the PAE matrix from the data dict.

        Size of the PAE matrix is NxN, where N is the number of tokens.

        Args:

            data (Dict):
                data dictionary from the data file.

        Returns:

            pae (np.array):
                PAE matrix.
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
    def get_contact_probs_mat(data: Dict):
        """Get the contact probabilities from the data dict.

        Args:

            data (Dict):
                data dictionary from the data file.

        Returns:

            contact_probs_mat (np.array):
                contact probabilities matrix from AlphaFold3 output.
        """

        if "contact_probs" in data:
            contact_probs_mat = np.array(data["contact_probs"])

        else:
            warnings.warn(
                "Contact probabilities not found, data file might not be AF3."
            )
            contact_probs_mat = None

        return contact_probs_mat

    @staticmethod
    def get_atom_chain_ids(data: Dict) -> list:
        """Get per atom chain IDs from the data dict.

        This is specific to AF3: "atom_chain_ids" key. \n
        If data is AF2, None is returned. \n

        Args:

            data (Dict):
                data dict from the data file.

        Returns:

            atom_chain_ids (list):
                per atom chain IDs.
        """

        if "atom_chain_ids" in data:
            atom_chain_ids = data["atom_chain_ids"]

        else:
            warnings.warn(
                """
                Atomized chain IDs not found, data file might be AF2.
                Structure file is required for AF2.
                """
            )
            atom_chain_ids = None

        return atom_chain_ids

    @staticmethod
    def get_atom_plddts(data: Dict) -> np.array:
        """Get per atom pLDDT scores from the data dict.

        This is specific to AF3: "atom_plddts" key. \n
        If data is AF2, None is returned. \n
        However, similar information can be obtained from the structure file. \n
        see :py:meth:`Parser.StructureParser.get_atom_plddts`.

        Args:

            data (Dict):
                data dict from the data file.

        Returns:

            atom_plddts (np.array):
                per atom pLDDT scores.
        """

        if "atom_plddts" in data:
            atom_plddts = np.array(data["atom_plddts"])

        else:
            warnings.warn(
                """
                Atomized pLDDT scores not found, data file might be AF2.
                Structure file is required for AF2.
                """
            )
            atom_plddts = None

        return atom_plddts