from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from itertools import combinations, product
from collections import defaultdict
from tqdm import tqdm
import warnings
import copy
from af_pipeline.constants.af_constants import (
    CHAIN_PAIRWISE_ASSESSMENT_COLUMNS,
    CHAINWISE_ASSESSMENT_COLUMNS
)
from af_pipeline.utils.misc_utils import time_it

class RigidBodyAssessment:
    """ Class to assess rigid bodies extracted from AlphaFold predictions."""

    rb_dict: Dict[str, List[Tuple[str, int]]]
    """ Dictionary of rigid bodies, where each rigid body is a dictionary
    with chain IDs as keys and residue numbers as values."""

    num_to_idx: Dict[str, Dict[int, Dict[str, int]]]
    """ Residue number to index mapping."""

    idx_to_num: Dict[int, Dict[str, str|int]]
    """ Index to residue number mapping."""

    contact_map: np.ndarray
    """ Contact map of the structure."""

    plddt_list: np.ndarray
    """ List of pLDDT scores for all residues in the structure."""

    pae: np.ndarray
    """ Predicted aligned error (PAE) matrix."""

    lengths_dict: Dict[str, int]
    """ Dictionary of lengths of chains in the structure."""

    symmetric_pae: bool
    """ Whether to treat PAE symmetrically between chain pairs."""

    as_average: bool
    """ Whether to report only the average of assessment metrics."""

    idr_chains: List[str]
    """ List of chain IDs that are considered disordered (IDRs)."""

    protein_chain_map: Dict[str, str]
    """ Mapping of chain IDs to protein names."""

    unique_chains: List[str]
    """ List of unique chain IDs in the rigid body."""

    chain_pairs: List[Tuple[str, str]]
    """ List of unique chain pairs in the rigid body."""

    rb_res_binary_map: np.ndarray
    """ Binary map of residues in the rigid body."""

    rb_res_pairs: Dict[Tuple[str, str], List[Tuple[int, int]]]
    """ Dictionary of residue pairs in the rigid body."""

    per_chain_plddt: Dict[str, np.ndarray]
    """ Dictionary of per-chain pLDDT scores."""

    # per_chain_avg_plddt: Dict[str, float]
    # """ Dictionary of per-chain average pLDDT scores."""

    pairwise_pae: Dict[Tuple[str, str], Dict[str, np.ndarray]]
    """ Dictionary of pairwise PAE values for each chain pair."""

    interface_res_pairs: Dict[Tuple[str, str], List[Tuple[int, int]]]
    """ Dictionary of interface residue pairs for each chain pair."""

    per_chain_interface_res: Dict[str, List[int]]
    """ Dictionary of interface residues for each chain."""

    num_contacts: Dict[Tuple[str, str], int]
    """ Dictionary of number of contacts for each chain pair."""

    num_interface_residues: Dict[Tuple[str, str], int]
    """ Dictionary of interface residues for each chain pair."""

    pairwise_ipae: Dict[Tuple[str, str], Dict[str, Dict[Tuple[int, int], float]]]
    """ Dictionary of pairwise iPAE values for each chain pair."""

    per_chain_iplddt: Dict[str, Dict[int, float]]
    """ Dictionary of per-chain ipLDDT scores."""

    # per_chain_avg_iplddt: Dict[str, float]
    # """ Dictionary of per-chain average ipLDDT scores."""

    # pairwise_avg_iplddt: Dict[Tuple[str, str], Dict[str, float]]
    # """ Dictionary of pairwise average ipLDDT scores for each chain pair."""

    pairwise_min_pae: Dict[Tuple[str, str], Dict[str, float]]
    """ Dictionary of pairwise minimum PAE values for each chain pair."""

    pairwise_avg_pae: Dict[Tuple[str, str], Dict[str, float]]
    """ Dictionary of pairwise average PAE values for each chain pair."""

    pairwise_avg_ipae: Dict[Tuple[str, str], Dict[str, float]]
    """ Dictionary of pairwise average iPAE values for each chain pair."""

    overall_assessment: dict
    """ Dictionary of overall assessment metrics for the rigid body."""

    save_path: str
    """ Path to save the assessment results."""

    def __init__(
        self,
        rb_dict: dict,
        num_to_idx: dict,
        idx_to_num: dict,
        contact_map: np.ndarray,
        plddt_list: np.ndarray,
        pae: np.ndarray,
        lengths_dict: dict,
        save_path: str,
        **kwargs,
    ):

        self.rb_dict = rb_dict
        self.num_to_idx = num_to_idx
        self.idx_to_num = idx_to_num
        self.symmetric_pae = kwargs.get("symmetric_pae", True)
        self.as_average = kwargs.get("as_average", False)
        self.idr_chains = kwargs.get("idr_chains", [])
        self.protein_chain_map = kwargs.get("protein_chain_map", {})

        self.unique_chains = self.get_unique_chains()
        self.chain_pairs = self.get_chain_pairs()
        self.rb_res_binary_map = self.get_rb_res_binary_map(lengths_dict)
        self.rb_res_pairs = self.get_rb_res_pairs()

        self.per_chain_plddt = self.get_per_chain_plddt(plddt_list)
        self.pairwise_pae = self.get_pairwise_pae(pae=pae)

        self.interface_res_pairs = self.get_interface_res_pairs(contact_map)
        self.per_chain_interface_res = self.get_per_chain_interface_residues()
        self.num_interface_residues = self.get_num_interface_residues()
        self.num_contacts = self.get_num_contacts()

        self.pairwise_ipae = self.get_pairwise_ipae(pae=pae)
        self.per_chain_iplddt = self.get_per_chain_iplddt(plddt_list)

        # self.per_chain_avg_plddt = self.get_per_chain_avg_plddt()
        # self.per_chain_avg_iplddt = self.get_per_chain_average_iplddt()
        # self.pairwise_avg_iplddt = self.get_pairwise_avg_iplddt()

        self.pairwise_min_pae = self.get_pairwise_min_pae(self.symmetric_pae)
        self.pairwise_avg_pae = self.get_pairwise_avg_pae(self.symmetric_pae)
        self.pairwise_avg_ipae = self.get_pairwise_avg_ipae(self.symmetric_pae)

        self.overall_assessment = self.get_overall_assessment()

        self.save_path = save_path

    def _get_protein_name(self, chain_id: str) -> str:
        """ Get the protein name for a given chain ID.

        Arguments:

        - **chain_id (str)**:<br />
            Chain ID.

        Returns:

        - **protein_name (str)**:<br />
            Protein name corresponding to the chain ID.
            If not found, returns "Unknown".
        """

        protein_name = self.protein_chain_map.get(chain_id, "Unknown")

        return protein_name

    def _get_chain_type(self, chain_id: str) -> str:
        """ Get the chain type (IDR or R) for a given chain ID.

        Arguments:

        - **chain_id (str)**:
            Chain ID.

        Returns:

        - **chain_type (str)**:<br />
            Chain type: "IDR" if the chain is in the list of IDR chains,
            otherwise "R".
        """

        chain_type = "IDR" if chain_id in self.idr_chains else "R"

        return chain_type

    def _get_ipae(
        self,
        res_pair: Tuple[int, int],
        chain_pair: Tuple[str, str],
        ipae_type: str
    ):
        """ Get the iPAE value for a given residue pair.

        There are three possible outcomes:
        - **ij**: iPAE value of residue 1 wrt 2.
        - **ji**: iPAE value of residue 2 wrt 1.
        - **average**: Average of ij and ji iPAE values.

        Arguments:

        - **res_pair (tuple)**:<br />
            Tuple of residue indices (res1, res2).

        - **chain_pair (tuple)**:<br />
            Tuple of chain IDs (chain1, chain2).

        - **ipae_type (str)**:<br />
            Type of iPAE value to retrieve. One of "ij", "ji", or "average".

        Returns:

        - **(float)**:<br />
            The requested iPAE value.
        """
        res1, res2 = res_pair

        ipae_ij = self.pairwise_ipae.get(chain_pair, {}).get("ij", {}).get(
            (res1, res2), np.nan
        )
        ipae_ji = self.pairwise_ipae.get(chain_pair, {}).get("ji", {}).get(
            (res2, res1), np.nan
        )

        ipae_dict = {
            "ij": ipae_ij,
            "ji": ipae_ji,
            "average": (ipae_ij + ipae_ji) / 2,
        }

        if ipae_type not in ipae_dict:
            raise ValueError(
                f"iPAE type '{ipae_type}' not recognized for residue pair "
                f"{res_pair} in chain pair {chain_pair}."
            )

        return ipae_dict[ipae_type]

    def _get_avg_ipae(
        self,
        chain_pair: Tuple[str, str],
        ipae_type: str
    ) -> float:
        """ Get the average iPAE value for a given chain pair.

        Arguments:

        - **chain_pair (tuple)**:<br />
            Tuple of chain IDs (chain1, chain2).

        - **ipae_type (str)**:<br />
            Type of average iPAE value to retrieve.
            One of "ij", "ji", or "average".

        Returns:

        - **(float)**:<br />
            The requested average iPAE value.
        """

        # ipae_dict = self.pairwise_ipae.get(chain_pair, {})
        # if len(ipae_dict) == 0:
        #     return np.nan

        # ipae_ij = list(ipae_dict["ij"].values())
        # ipae_ji = list(ipae_dict["ji"].values())

        # if self.symmetric_pae:
        #     ipae_ij = (ipae_ij + ipae_ji) / 2
        #     ipae_ji = ipae_ij

        # attrs_ = {
        #     "ij": np.mean(ipae_ij),
        #     "ji": np.mean(ipae_ji),
        #     "average": np.mean(ipae_ij + ipae_ji) / 2,
        # }

        avg_ipae_ij = self.pairwise_avg_ipae.get(chain_pair, {}).get(
            "ij", np.nan
        )
        avg_ipae_ji = self.pairwise_avg_ipae.get(chain_pair, {}).get(
            "ji", np.nan
        )

        avg_ipae_dict = {
            "ij": avg_ipae_ij,
            "ji": avg_ipae_ji,
            "average": (avg_ipae_ij + avg_ipae_ji) / 2,
        }

        if ipae_type not in avg_ipae_dict:
            raise ValueError(
                f"Average iPAE type '{ipae_type}' not recognized for chain pair"
                f" {chain_pair}."
            )

        return avg_ipae_dict[ipae_type]

    # @time_it
    def get_pairwise_avg_ipae(self, symmetric_pae: bool = True) -> Dict[tuple, Dict[str, float]]:
        """ Get the average iPAE for each chain pair.

        Arguments:

        - **symmetric_pae (bool, optional)**:<br />
            If True, calculates the average iPAE symmetrically for both
            directions (ij and ji).

        Returns:

        - **pairwise_avg_ipae (defaultdict)**:<br />
            A dictionary where keys are chain pairs (tuples) and values are the
            average iPAE values.
        """

        pairwise_avg_ipae = defaultdict(dict)

        for chain_pair, ipae_dict in self.pairwise_ipae.items():

            if symmetric_pae:
                _avg_ipae = np.mean(
                    list(ipae_dict["ij"].values()) +
                    list(ipae_dict["ji"].values())
                ) / 2
                pairwise_avg_ipae[chain_pair]["ij"] = _avg_ipae
                pairwise_avg_ipae[chain_pair]["ji"] = _avg_ipae
            else:
                pairwise_avg_ipae[chain_pair]["ij"] = np.mean(list(ipae_dict["ij"].values()))
                pairwise_avg_ipae[chain_pair]["ji"] = np.mean(list(ipae_dict["ji"].values()))

        pairwise_avg_ipae = dict(pairwise_avg_ipae)

        return pairwise_avg_ipae

    def _get_avg_pae(
        self,
        chain_pair: Tuple[str, str],
        pae_type: str,
    ) -> float:
        """ Get the average PAE value for a given chain pair.

        Arguments:

        - **chain_pair (tuple)**:<br />
            Tuple of chain IDs (chain1, chain2).

        - **pae_type (str)**:<br />
            Type of average PAE value to retrieve.
            One of "ij", "ji", or "average".

        Returns:

        - **(float)**:<br />
            The requested average PAE value.
        """

        # pae_ij = self.pairwise_pae[chain_pair]["ij"]
        # pae_ji = self.pairwise_pae[chain_pair]["ji"]

        # if self.symmetric_pae:
        #     pae_ij = (pae_ij + pae_ji) / 2
        #     pae_ji = pae_ij

        # attrs_ = {
        #     "ij": np.mean(pae_ij),
        #     "ji": np.mean(pae_ji),
        #     "average": np.mean(pae_ij + pae_ji) / 2,
        # }

        avg_pae_ij = self.pairwise_avg_pae.get(chain_pair, {}).get(
            "ij", np.nan
        )
        avg_pae_ji = self.pairwise_avg_pae.get(chain_pair, {}).get(
            "ji", np.nan
        )

        avg_pae_dict = {
            "ij": avg_pae_ij,
            "ji": avg_pae_ji,
            "average": (avg_pae_ij + avg_pae_ji) / 2,
        }

        if pae_type not in avg_pae_dict:
            raise ValueError(
                f"Average PAE type '{pae_type}' not recognized for chain pair"
                f" {chain_pair}."
            )

        return avg_pae_dict[pae_type]

    # @time_it
    def get_pairwise_avg_pae(self, symmetric_pae: bool = True) -> Dict[tuple, Dict[str, float]]:
        """ Get the average PAE for each chain pair.

        Arguments:

        - **symmetric_pae (bool, optional)**:<br />
            If True, calculates the average PAE symmetrically for both
            directions (ij and ji).

        Returns:

        - **pairwise_avg_pae (dict)**:
            A dictionary where keys are chain pairs (tuples) and values are the
            average PAE values.
        """

        pairwise_avg_pae = defaultdict(dict)

        for chain_pair in self.chain_pairs:
            if symmetric_pae:
                _avg_pae = np.mean(
                    self.pairwise_pae[chain_pair]["ij"] +
                    self.pairwise_pae[chain_pair]["ji"]
                ) / 2
                pairwise_avg_pae[chain_pair]["ij"] = _avg_pae
                pairwise_avg_pae[chain_pair]["ji"] = _avg_pae
            else:
                pairwise_avg_pae[chain_pair]["ij"] = np.mean(
                    self.pairwise_pae[chain_pair]["ij"]
                )
                pairwise_avg_pae[chain_pair]["ji"] = np.mean(
                    self.pairwise_pae[chain_pair]["ji"]
                )

        pairwise_avg_pae = dict(pairwise_avg_pae)

        return pairwise_avg_pae

    def _get_min_pae(
        self,
        chain_pair: Tuple[str, str],
        pae_type: str,
    ):
        """ Get the minimum PAE value for a given chain pair.

        Arguments:

        - **chain_pair (tuple)**:<br />
            Tuple of chain IDs (chain1, chain2).

        - **pae_type (str)**:<br />
            Type of minimum PAE value to retrieve.
            One of "ij", "ji", or "average".

        Returns:

        - **(float)**:<br />
            The requested minimum PAE value.
        """

        pae_ij = self.pairwise_pae[chain_pair]["ij"]
        pae_ji = self.pairwise_pae[chain_pair]["ji"]

        if self.symmetric_pae:
            pae_ij = np.minimum(pae_ij, pae_ji)
            pae_ji = pae_ij

        attrs_ = {
            "ij": np.min(pae_ij),
            "ji": np.min(pae_ji),
            "average": np.min(np.concatenate([pae_ij, pae_ji])),
        }

        if pae_type not in attrs_:
            raise ValueError(
                f"Minimum PAE type '{pae_type}' not recognized for chain pair"
                f" {chain_pair}."
            )

        return attrs_[pae_type]

    # # @time_it
    def get_pairwise_min_pae(self, symmetric_pae: bool = True) -> Dict[tuple, Dict[str, float]]:
        """Get the minimum PAE for each chain pair.

        Arguments:

        - **symmetric_pae (bool, optional)**:
            If True, calculates the minimum PAE symmetrically for both
            directions (ij and ji).

        Returns:

        - **pairwise_min_pae (defaultdict)**:<br />
            A dictionary where keys are chain pairs (tuples) and values are the
            minimum PAE values.
            If `symmetric_pae` is `True`, the minimum PAE is calculated as the
            minimum of both directions (ij and ji).
            If `symmetric_pae` is `False`, the minimum PAE is calculated
            separately for each direction.
        """

        pairwise_min_pae = defaultdict(dict)

        for chain_pair, pae_dict in self.pairwise_pae.items():
            if symmetric_pae:
                _min_pae = np.min(
                    [np.min(pae_dict["ij"]), np.min(pae_dict["ji"])]
                )
                pairwise_min_pae[chain_pair]["ij"] = _min_pae
                pairwise_min_pae[chain_pair]["ji"] = _min_pae
            else:
                pairwise_min_pae[chain_pair]["ij"] = np.min(pae_dict["ij"])
                pairwise_min_pae[chain_pair]["ji"] = np.min(pae_dict["ji"])

        return pairwise_min_pae

    def _get_avg_plddt(
        self,
        chain_id: str,
    ) -> float:
        """ Get the average pLDDT value for a given chain ID.

        Arguments:

        - **chain_id (str)**:<br />
            Chain ID.

        Returns:

        - **(float)**:<br />
            The average pLDDT value for the given chain ID.
        """

        return np.nanmean(self.per_chain_plddt[chain_id])

    def _get_avg_iplddt(
        self,
        chain_pair: Tuple[str, str] | None,
        chain_id: str | None = None,
    ) -> float:
        """ Get the average ipLDDT value for a given chain in the context of
        the chain pair or the chain pair itself.

        Arguments:

        - **chain_pair (Tuple[str, str])**:<br />
            Tuple of chain IDs (chain1, chain2).

        - **chain_id (str | None, optional):<br />
            Chain ID for which average ipLDDT is to be calculated.
            If `None`, average ipLDDT for the entire interface is calculated.

        Returns:

        - **(float)**:<br />
            The requested average ipLDDT value.
        """

        assert chain_pair is not None or chain_id is not None, \
            "Either chain_pair or chain_id must be provided."

        if chain_pair is None:
            return np.nanmean(
                list(self.per_chain_iplddt.get(chain_id, {}).values())
            )

        chain1, chain2 = chain_pair
        iplddt_chain1 = self.per_chain_iplddt.get(chain1, {})
        iplddt_chain2 = self.per_chain_iplddt.get(chain2, {})
        interacting_res_pairs = self.interface_res_pairs.get(chain_pair, [])
        iplddt_vals_chain1 = [
            iplddt_chain1.get(res1, np.nan)
            for res1, _ in interacting_res_pairs
        ]
        iplddt_vals_chain2 = [
            iplddt_chain2.get(res2, np.nan)
            for _, res2 in interacting_res_pairs
        ]
        iplddt_vals_interface = iplddt_vals_chain1 + iplddt_vals_chain2

        attrs_ = {
            chain1: iplddt_vals_chain1 if iplddt_vals_chain1 else [np.nan],
            chain2: iplddt_vals_chain2 if iplddt_vals_chain2 else [np.nan],
            None: iplddt_vals_interface if iplddt_vals_interface else [np.nan],
        }

        assert chain_id in attrs_, "Invalid chain_id provided."

        return np.nanmean(attrs_[chain_id])

    def _get_chain_attr(
        self,
        chain_id: str,
        attr_name: str
    ) -> float | str | int:
        """ Get the attribute value for a given chain ID.

        Arguments:

        - **chain_id (str)**:<br />
            Chain ID.

        - **attr_name (str)**:<br />
            Name of the attribute to retrieve.
            Valid attributes are defined in `CHAINWISE_ASSESSMENT_COLUMNS`.

        Returns:

        - **(float | str | int)**:<br />
            The requested attribute value.
        """

        if attr_name not in CHAINWISE_ASSESSMENT_COLUMNS:
            raise ValueError(
                f"Attribute '{attr_name}' not recognized for chain "
                f"'{chain_id}'."
            )

        attrs_ = {
            "Chain ID": chain_id,
            "Protein Name": self._get_protein_name(chain_id),
            "Average pLDDT": self._get_avg_plddt(chain_id),
            "Average ipLDDT": self._get_avg_iplddt(None, chain_id),
            "Interface Residues": len(
                self.per_chain_interface_res.get(chain_id, [])
            ),
            "Chain Type": self._get_chain_type(chain_id),
        }

        return attrs_[attr_name]

    def _get_res_attr(
        self,
        chain_id: str,
        res_num: int,
        attr_name: str
    ) -> float:
        """ Get the attribute value for a given residue in a chain.

        Arguments:

        - **chain_id (str)**:<br />
            Chain ID.

        - **res_num (int)**:<br />
            Residue number.

        - **attr_name (str)**:<br />
            Name of the attribute to retrieve.

        Returns:

        - **(float)**:<br />
            The requested attribute value.
        """

        attrs_ = {
            "iplddt": self.per_chain_iplddt.get(chain_id, {}).get(res_num, np.nan),
        }
        if attr_name not in attrs_:
            raise ValueError(
                f"Attribute '{attr_name}' not recognized for residue "
                f"{res_num} in chain '{chain_id}'."
            )

        return attrs_[attr_name]

    # @time_it
    def _get_cp_attr(
        self,
        chain_pair: tuple,
        attr_name: str,
    ) -> float | str | int:

        chain1, chain2 = chain_pair

        attrs_ = {
            "Chain ID 1": chain1,
            "Chain ID 2": chain2,
            "Protein Name 1": self._get_protein_name(chain1),
            "Protein Name 2": self._get_protein_name(chain2),
            "Chain Type 1": self._get_chain_type(chain1),
            "Chain Type 2": self._get_chain_type(chain2),
            "Interface Residues": self.num_interface_residues.get(chain_pair, 0),
            "Number of contacts": self.num_contacts.get(chain_pair, 0),
            "Average ipLDDT": self._get_avg_iplddt(chain_pair, None),
            "Average ipLDDT chain1": self._get_avg_iplddt(chain_pair, chain1),
            "Average ipLDDT chain2": self._get_avg_iplddt(chain_pair, chain2),
            "Average PAE": self._get_avg_pae(chain_pair, "average"),
            "Average PAE ij": self._get_avg_pae(chain_pair, "ij"),
            "Average PAE ji": self._get_avg_pae(chain_pair, "ji"),
            "Average iPAE": self._get_avg_ipae(chain_pair, "average"),
            "Average iPAE ij": self._get_avg_ipae(chain_pair, "ij"),
            "Average iPAE ji": self._get_avg_ipae(chain_pair, "ji"),
            "Minimum PAE": self._get_min_pae(chain_pair, "average"),
            "Minimum PAE ij": self._get_min_pae(chain_pair, "ij"),
            "Minimum PAE ji": self._get_min_pae(chain_pair, "ji"),
        }

        if attr_name not in attrs_:
            raise ValueError(
                f"Attribute '{attr_name}' not recognized for chain pair "
                f"{chain_pair}."
            )

        return attrs_[attr_name]

    def _get_rp_attr(self, res_pair: tuple, chain_pair: str, attr_name: str):

        res1, res2 = res_pair
        chain1, chain2 = chain_pair

        attrs_ = {
            "Chain ID 1": chain1,
            "Chain ID 2": chain2,
            "Protein Name 1": self._get_protein_name(chain1),
            "Protein Name 2": self._get_protein_name(chain2),
            "Chain Type 1": self._get_chain_type(chain1),
            "Chain Type 2": self._get_chain_type(chain2),
            "Residue 1": res1,
            "Residue 2": res2,
            "ipLDDT res1": self._get_res_attr(chain1, res1, "iplddt"),
            "ipLDDT res2": self._get_res_attr(chain2, res2, "iplddt"),
            "iPAE": self._get_ipae(res_pair, chain_pair, "average"),
            "iPAE ij": self._get_ipae(res_pair, chain_pair, "ij"),
            "iPAE ji": self._get_ipae(res_pair, chain_pair, "ji"),
        }
        if attr_name not in attrs_:
            raise ValueError(
                f"Attribute '{attr_name}' not recognized for residue pair "
                f"{res_pair} in chain pair {chain_pair}."
            )
        return attrs_[attr_name]

    # @time_it
    def save_rb_assessment(self):
        """ Save the assessment of the rigid bodies to an Excel file.

        The assessment includes:
        - **Per chain assessment**:<br />
            Average pLDDT, Average iLDDT, interface residues,
            Chain type (IDR or R).
        - **Per chain pair assessment**:<br />
            interface residues, Number of contacts, Average PAE,
            Average iPAE, Minimum PAE, Average iLDDT for each chain,
            Chain type (IDR or R) for each chain.
        - **Overall assessment**:<br />
            Average pLDDT, Average iLDDT, interface residues,
            Chain type (IDR or R).

        The assessment is saved in an Excel file with three sheets:
        - "Chain Wise Assessment": Contains per chain assessment data.
        - "Chain Pairwise Assessment": Contains per chain pair assessment data.
        - "Overall Assessment": Contains overall assessment data.
        """

        chain_wise_assessment_rows = []
        chain_pairwise_assessment_rows = []
        overall_assessment_rows = []

        for chain_id in self.unique_chains:
            chain_wise_assessment_rows.append({
                k: self._get_chain_attr(chain_id, k)
                for k in CHAINWISE_ASSESSMENT_COLUMNS
            })

        attr_state = (self.as_average, self.symmetric_pae)

        for chain_pair in self.chain_pairs:

            if self.as_average is True:
                chain_pairwise_assessment_rows.append({
                    k: self._get_cp_attr(chain_pair, k)
                    for k in CHAIN_PAIRWISE_ASSESSMENT_COLUMNS[attr_state]
                })

            elif self.as_average is False:
                for res_pair in self.interface_res_pairs[chain_pair]:
                    chain_pairwise_assessment_rows.append({
                        k: self._get_rp_attr(res_pair, chain_pair, k)
                        for k in CHAIN_PAIRWISE_ASSESSMENT_COLUMNS[attr_state]
                    })

        overall_assessment_keys = {
            "Number of Chains": "num_chains",
            "Number of Interacting Chain Pairs": "num_interacting_chain_pairs",
            "Interface Residues": "num_interface_residues",
            "Number of Contacts": "num_contacts",
            "Average ipLDDT": "avg_iplddt",
            "Average IDR ipLDDT": "avg_idr_iplddt",
            "Average iPAE ij": "avg_ipae_ij",
            "Average iPAE ji": "avg_ipae_ji",
        }

        for col_head, key in overall_assessment_keys.items():
            if self.overall_assessment.get(key, np.nan) is not np.nan:
                overall_assessment_rows.append({
                    "Key": col_head,
                    "Value": self.overall_assessment.get(key)
                })

        chain_pairwise_assessment_df = pd.DataFrame(chain_pairwise_assessment_rows)
        chainwise_assessment_df = pd.DataFrame(chain_wise_assessment_rows)
        overall_assessment_df = pd.DataFrame(overall_assessment_rows)

        df_dict = {
            "chain_pairwise_assessment": chain_pairwise_assessment_df,
            "chainwise_assessment": chainwise_assessment_df,
            "overall_assessment": overall_assessment_df,
        }

        for k, df_ in df_dict.items():
            df_dict[k] = df_.fillna(np.nan)
            df_dict[k] = df_.map(lambda x: round(x, 2) if isinstance(x, (int, float)) else x)

        with pd.ExcelWriter(self.save_path, engine='openpyxl', mode='w') as writer:
            for sheet_name, df in df_dict.items():

                if df.empty:
                    warnings.warn(f"Skipping empty DataFrame for sheet: {sheet_name}")
                    continue

                df.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    index=False,
                )

    # @time_it
    def get_unique_chains(self) -> List[str]:
        """Get unique chains in the rigid body.

        Returns:

        - **unique_chains (list)**:<br />
            List of unique chain IDs in the rigid body.
        """

        unique_chains = [
            chain_id
            for chain_id in self.rb_dict.keys()
            if len(self.rb_dict[chain_id]) > 0
        ]

        return unique_chains

    # @time_it
    def get_chain_pairs(self) -> List[Tuple[str, str]]:
        """Get all unique chain pairs in the rigid body.

        Returns:

        - **chain_pairs (list)**:<br />
            List of tuples containing unique chain pairs.
            Each tuple contains two chain IDs.
        """

        chain_pairs = list(combinations(self.unique_chains, 2))

        return [tuple(pair) for pair in chain_pairs]

    # @time_it
    def get_rb_res_binary_map(self, lengths_dict: Dict[str, int]) -> np.ndarray:
        """Get a binary map of residues in the rigid body.

        Arguments:

        - **lengths_dict (dict)**:<br />
            Dictionary containing lengths of chains in the structure.

        Returns:

        - **rb_res_binary_map (np.ndarray)**:<br />
            A binary map of residues in the rigid body.
            The shape is (`total_length`, `total_length`) where `total_length`
            is the sum of lengths of all chains.
            The value is 1 if the residue is part of the rigid body, 0 otherwise.
        """

        total_len = lengths_dict.get("total", 0)
        rb_res_binary_map = np.zeros((total_len, total_len), dtype=int)
        all_rb_interface_res_idxs = []

        for chain_id, atom_name_token_list in self.rb_dict.items():

            res_idxs = [
                self.num_to_idx[chain_id][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list
            ]
            all_rb_interface_res_idxs.extend(res_idxs)

        all_rb_interface_res_idxs = np.unique(all_rb_interface_res_idxs)

        rb_res_binary_map[
            np.ix_(all_rb_interface_res_idxs, all_rb_interface_res_idxs)
        ] = 1

        return rb_res_binary_map

    # @time_it
    def get_rb_res_pairs(self) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
        """Get all unique residue pairs in the rigid body.

        Returns:

        - **rb_res_pairs (defaultdict)**:<br />
            A dictionary where keys are chain pairs (tuples) and values are
            lists of residue index pairs.
            Each residue index pair is a tuple of indices from the two chains
            in the rigid body.
        """

        rb_res_pairs = defaultdict(list)

        for chain_pair in self.chain_pairs:

            chain1, chain2 = chain_pair

            atom_name_token_list1 = self.rb_dict[chain1]
            atom_name_token_list2 = self.rb_dict[chain2]

            res1_idxs = [
                self.num_to_idx[chain1][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list1
            ]
            res2_idxs = [
                self.num_to_idx[chain2][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list2
            ]

            # Create pairs of residues from the two chains
            pairs = list(product(res1_idxs, res2_idxs))
            rb_res_pairs[chain_pair].extend(pairs)

        return rb_res_pairs

    # @time_it
    def get_interface_res_pairs(
        self,
        contact_map: np.ndarray,
    ) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
        """ Get interface residue pairs from the contact map.

        Arguments:

        - **contact_map (np.ndarray)**:<br />
            A binary contact map where 1 indicates a contact between residues
            and 0 indicates no contact.

        Returns:

        - **interface_res_pairs (defaultdict)**:<br />
            A dictionary where keys are chain pairs (tuples) and values are
            lists of residue index pairs.
        """

        interface_res_pairs = defaultdict(list)
        contacting_res_indices = np.argwhere(contact_map == 1)

        for chain1, chain2 in tqdm(self.chain_pairs):

            atom_name_token_list1 = self.rb_dict[chain1]
            atom_name_token_list2 = self.rb_dict[chain2]

            res1_idxs = [
                self.num_to_idx[chain1][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list1
            ]
            res2_idxs = [
                self.num_to_idx[chain2][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list2
            ]

            mask1 = np.isin(contacting_res_indices[:, 0], res1_idxs)
            mask2 = np.isin(contacting_res_indices[:, 1], res2_idxs)

            mask = mask1 & mask2
            contacting_res_pairs = set(
                map(tuple, contacting_res_indices[mask])
            )

            if len(contacting_res_pairs) > 0:
                interface_res_pairs[(chain1, chain2)].extend(
                    list(contacting_res_pairs)
                )

        return interface_res_pairs

    # @time_it
    def get_per_chain_interface_residues(self) -> Dict[str, List[int]]:
        """Get interface residues for each chain.

        Returns:

        - **per_chain_interface_residues (defaultdict)**:<br />
            A dictionary where keys are chain IDs and values are lists of
            residue indices.
            Each list contains the indices of residues that are part of any of
            the interfaces that the chain is involved in.
        """

        per_chain_interface_res = defaultdict(list)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            chain1, chain2 = chain_pair

            for res1_idx, res2_idx in interacting_res_pairs:

                per_chain_interface_res[chain1].append(
                    res1_idx
                ) if res1_idx not in per_chain_interface_res[chain1] else None

                per_chain_interface_res[chain2].append(
                    res2_idx
                ) if res2_idx not in per_chain_interface_res[chain2] else None

        per_chain_interface_res = dict(per_chain_interface_res)

        return per_chain_interface_res

    # @time_it
    def get_num_interface_residues(self) -> Dict[tuple, int]:
        """Get the interface residues for each chain pair.

        Returns:

        - **num_interface_residues (dict)**:<br />
            A dictionary where keys are chain pairs (tuples) and values are the
            number of unique interface residues.
            Each key is a tuple of two chain IDs, and the value is the count of
            unique residues that interact between those chains.
        """

        num_interface_residues = {}

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            unique_interface_residues = np.unique(
                np.array(interacting_res_pairs).flatten()
            )
            num_interface_residues[chain_pair] = len(unique_interface_residues)

        return num_interface_residues

    # @time_it
    def get_num_contacts(self) -> Dict[tuple, int]:
        """Get the number of contacts for each chain pair.

        Returns:

        - **num_contacts (dict)**:<br />
            A dictionary where keys are chain pairs (tuples) and values are the
            number of contacts.
            Each key is a tuple of two chain IDs, and the value is the count of
            contacts between those chains.
        """

        num_contacts = {}

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            num_contacts[chain_pair] = len(interacting_res_pairs)

        return num_contacts

    # @time_it
    def get_per_chain_plddt(self, plddt_list: list) -> Dict[str, np.ndarray]:
        """ Get per-chain pLDDT scores from a list of pLDDT scores.

        Arguments:

        - **plddt_list (list)**:<br />
            A list of pLDDT scores for all residues in the structure.

        Returns:

        - **per_chain_plddt (dict)**:<br />
            A dictionary where keys are chain IDs and values are numpy arrays
            of pLDDT scores for residues in that chain.
        """

        per_chain_plddt = {}

        for chain_id, atom_name_token_list in self.rb_dict.items():

            res_idxs = [
                self.num_to_idx[chain_id][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list
            ]

            plddt_scores = np.array(plddt_list)[res_idxs]
            per_chain_plddt[chain_id] = plddt_scores

        return per_chain_plddt

    # @time_it
    def get_per_chain_iplddt(self, plddt_list: list) -> Dict[str, Dict[int, float]]:
        """ Get per-chain ipLDDT scores from a list of pLDDT scores.

        Arguments:

        - **plddt_list (list)**:<br />
            A list of pLDDT scores for all residues in the structure.

        Returns:

        - **per_chain_iplddt (dict)**:<br />
            A dictionary where keys are chain IDs and values are dictionaries
            mapping residue indices to their pLDDT scores.
        """

        per_chain_iplddt = defaultdict(dict)

        for chain_id, interface_res_idxs in self.per_chain_interface_res.items():

            for res_idx in interface_res_idxs:

                per_chain_iplddt[chain_id][res_idx] = plddt_list[res_idx]

        per_chain_iplddt = dict(per_chain_iplddt)

        return per_chain_iplddt

    # @time_it
    def get_pairwise_pae(
        self,
        pae: np.ndarray
    ) -> Dict[Tuple[str, str], Dict[str, np.ndarray]]:
        """ Get pairwise PAE values for each chain pair.

        Arguments:

        - **pae (np.ndarray)**:<br />
            A 2D numpy array representing the PAE matrix.

        Returns:

        - **pairwise_pae (defaultdict)**:
            A dictionary where keys are chain pairs (tuples) and values are
            dictionaries containing PAE values for residue pairs.
        """

        pairwise_pae = defaultdict(np.ndarray)

        for chain_pair in self.chain_pairs:

            rb_chain_pair_res = self.rb_res_pairs[chain_pair]

            rb_pae_vals_ij = [
                pae[res1_idx, res2_idx]
                for res1_idx, res2_idx in rb_chain_pair_res
            ]

            rb_pae_vals_ji = [
                pae[res2_idx, res1_idx]
                for res1_idx, res2_idx in rb_chain_pair_res
            ]

            if len(rb_pae_vals_ij) > 0:
                pairwise_pae[chain_pair] = {
                    "ij": rb_pae_vals_ij,
                    "ji": rb_pae_vals_ji,
                }

        return pairwise_pae

    # @time_it
    def get_pairwise_ipae(self, pae) -> Dict[Tuple[str, str], Dict[str, Dict[Tuple[int, int], float]]]:
        """ Get pairwise iPAE values for each chain pair.

        Arguments:

        - **pae (np.ndarray)**:<br />
            A 2D numpy array representing the PAE matrix.

        Returns:

        - **pairwise_ipae (defaultdict)**:<br />
            A dictionary where keys are chain pairs (tuples) and values are
            dictionaries containing iPAE values for residue pairs.
        """

        pairwise_ipae = {}

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            pairwise_ipae[chain_pair] = {
                "ij" : {
                    (res1_idx, res2_idx): pae[res1_idx, res2_idx]
                    for res1_idx, res2_idx in interacting_res_pairs
                },
                "ji" : {
                    (res2_idx, res1_idx): pae[res2_idx, res1_idx]
                    for res1_idx, res2_idx in interacting_res_pairs
                }
            }

        return pairwise_ipae

    # @time_it
    def get_overall_assessment(self):
        """ Get overall assessment of the rigid body.

        Returns:

        - **overall_assessment (dict)**:<br />
            A dictionary containing overall statistics about the rigid body.
            It includes the number of chains, number of interacting chain pairs,
            interface residues, number of contacts, average ipLDDT,
            average IDR ipLDDT, average iPAE ij, and average iPAE ji.
        """

        overall_assessment = {}

        overall_assessment["num_chains"] = len(self.unique_chains)

        overall_assessment["num_interacting_chain_pairs"] = len(self.interface_res_pairs)

        overall_assessment["num_interface_residues"] = sum(
            len(res_list)
            for res_list in self.per_chain_interface_res.values()
        )

        overall_assessment["num_contacts"] = sum(
            len(contact_pairs)
            for contact_pairs in self.interface_res_pairs.values()
        )

        global_iplddt_scores = [
            iplddt
            for iplddt_scores in self.per_chain_iplddt.values()
            for iplddt in iplddt_scores.values()
        ]

        overall_assessment["avg_iplddt"] = (
            np.mean(global_iplddt_scores) if global_iplddt_scores else np.nan
        )

        global_idr_iplddt_scores = [
            iplddt
            for chain_id, iplddt_scores in self.per_chain_iplddt.items()
            for iplddt in iplddt_scores.values()
            if chain_id in self.idr_chains
        ]

        overall_assessment["avg_idr_iplddt"] = (
            np.mean(global_idr_iplddt_scores) if global_idr_iplddt_scores else np.nan
        )

        global_ipae_ij_scores = [
            ipae
            for ipae_dict in self.pairwise_ipae.values()
            for ipae in ipae_dict["ij"].values()
        ]

        global_ipae_ji_scores = [
            ipae
            for ipae_dict in self.pairwise_ipae.values()
            for ipae in ipae_dict["ji"].values()
        ]

        overall_assessment["avg_ipae_ij"] = (
            np.mean(global_ipae_ij_scores) if global_ipae_ij_scores else np.nan
        )

        overall_assessment["avg_ipae_ji"] = (
            np.mean(global_ipae_ji_scores) if global_ipae_ji_scores else np.nan
        )

        return overall_assessment

    # @time_it
    # def get_per_chain_avg_plddt(self) -> Dict[str, float]:
    #     """ Get the average pLDDT score for each chain.

    #     Returns:

    #     - **per_chain_avg_plddt (dict)**:<br />
    #         A dictionary where keys are chain IDs and values are the average
    #         pLDDT scores for that chain.
    #     """

    #     return {
    #         chain_id: np.mean(plddt_scores)
    #         for chain_id, plddt_scores in self.per_chain_plddt.items()
    #     }

    # @time_it
    # def get_per_chain_average_iplddt(self) -> Dict[str, float]:
    #     """ Get the average ipLDDT score for each chain.

    #     Returns:

    #     - **per_chain_avg_iplddt (dict)**:<br />
    #         A dictionary where keys are chain IDs and values are the average
    #         ipLDDT scores for that chain.
    #     """

    #     return {
    #         chain_id: np.mean(list(iplddt_scores.values()))
    #         for chain_id, iplddt_scores in self.per_chain_iplddt.items()
    #     }

    # @time_it
    # def get_pairwise_avg_iplddt(self) -> Dict[tuple, Dict[str, float]]:
    #     """ Get the average ipLDDT for each chain pair.

    #     Returns:

    #     - **pairwise_avg_iplddt (dict)**:<br />
    #         A dictionary where keys are chain pairs (tuples) and values are
    #         dictionaries containing average ipLDDT values for each chain in the
    #         pair.
    #     """

    #     pairwise_avg_iplddt = defaultdict(dict)

    #     for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

    #         chain1, chain2 = chain_pair

    #         iplddt1_values = [
    #             self.per_chain_iplddt[chain1].get(res1_idx, np.nan)
    #             for res1_idx, res2_idx in interacting_res_pairs
    #         ]

    #         iplddt2_values = [
    #             self.per_chain_iplddt[chain2].get(res2_idx, np.nan)
    #             for res1_idx, res2_idx in interacting_res_pairs
    #         ]

    #         pairwise_avg_iplddt[chain_pair][chain1] = np.mean(iplddt1_values)
    #         pairwise_avg_iplddt[chain_pair][chain2] = np.mean(iplddt2_values)

    #     pairwise_avg_iplddt = dict(pairwise_avg_iplddt)

    #     return pairwise_avg_iplddt
