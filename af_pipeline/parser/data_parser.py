"""
Data Parser module
==================
Methods to parse the `JSON` or `PKL` data files[^data_file] provided with the structure
predictions from AlphaFold2, AlphaFold3, and ColabFold.

Currently supports the following:
- AlphaFold2 data files (.pkl or .json)
- AlphaFold3 data files (.json)
- ColabFold data files (.json)

[^data_file]: *The "data file" in this context refers to the file that contains \
  prediction metrics such as PAE matrices, pLDDT scores, (token chain IDs, \
  token residue IDs, and contact probabilities for AF3), among others.*
"""

import warnings
import numpy as np
import os
from typing import Dict, Callable, List

from af_pipeline.constants.af_constants import (
    ALLOWED_DATA_FORMATS,
    AVAILABLE_DATA_READERS,
)

class DataParser:
    """Class with methods to parse the prediction data files."""

    data_file_path: str
    """ Path to the data file provided with the prediction. """

    def __init__(self, data_file_path: str):

        self.data_file_path = data_file_path

    @property
    def data_type(self) -> str:
        """ Data file type based on the file extension.

        See af_pipeline.constants.af_constants.ALLOWED_DATA_FORMATS for
        supported formats.
        """

        ext = os.path.splitext(self.data_file_path)[1].replace(".", "")

        if ext not in ALLOWED_DATA_FORMATS:
            raise Exception(
                f"""

                Incorrect file format: {ext}.
                Suported formats are {ALLOWED_DATA_FORMATS}.
                """
            )

        return ext

    @property
    def parser(self) -> Callable[[str], Dict | List]:
        """Parser function based on the file extension.

        See af_pipeline.constants.af_constants.AVAILABLE_DATA_READERS for
        supported parsers.
        """

        ext = self.data_type

        parser = AVAILABLE_DATA_READERS[ext]
        assert callable(parser), f"Parser for '{ext}' is not a function."

        return parser

    def get_data_dict(self) -> Dict:
        """Get the data from the data file.
        Currently supports `.pkl` and `.json` data files.

        Returns:

        - **data (Dict)**:<br />
            Data the data file in the dictionary format.
        """

        data = self.parser(self.data_file_path)

        if isinstance(data, list):
            data = data[0]

        return data

    @staticmethod
    def get_token_chain_ids(data: Dict) -> list | None:
        """Get the token chain IDs from the data dictionary.

        This is specific to AF3: `"token_chain_ids"` key.
        For others, `None` is returned.

        Each token can be a residue/nucleotide/ion in the chain.<br />
        In case of modified residues or nucleotides or glycans,
        each atom is a token.

        Arguments:

        - **data (Dict)**:<br />
            Data dictionary from the data file.

        Returns:

        - **token_chain_ids (list)**:<br />
            Token chain IDs.
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

        > [!NOTE]
        > Atom-level tokens of a given modified residue or nucleotide or glycans
        > have the same token residue IDs.

        Arguments:

        - **data (Dict)**:<br />
            Data dictionary from the data file.

        Returns:

        - **token_res_ids (list)**:<br />
            Token residue IDs provided in the AF3 `JSON` file.
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

        Size of the PAE matrix is `NxN`, where `N` is the number of tokens.

        Arguments:

        - **data (Dict)**:<br />
            Data dictionary from the data file.

        Returns:

        - **pae (np.ndarray)**:<br />
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
    def get_contact_probs_mat(data: Dict) -> np.ndarray | None:
        """Get the contact probabilities from the data dictionary.

        This is specific to AF3: `"contact_probs"` key.
        For others, `None` is returned.

        Arguments:

        - **data (Dict)**:<br />
            Data dictionary from the data file.

        Returns:

        - **contact_probs_mat (np.ndarray)**:<br />
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

        Arguments:

        - **data (Dict)**:<br />
            Data dictionary from the data file.

        Returns:

        - **atom_chain_ids (list)**: <br />
            Atom chain IDs provided in the AF3 JSON file.
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

        Arguments:

        - **data (Dict)**:<br />
            Data dictionary from the data file.

        Returns:

        - **atom_plddts (np.ndarray)**:<br />
            Per atom pLDDT scores provided in the AF3 JSON file.
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