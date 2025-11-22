from typing import Dict, List, Tuple
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from itertools import combinations, combinations_with_replacement, product
from collections import defaultdict
from tqdm import tqdm
import warnings
import copy
from af_pipeline.constants.af_constants import (
    CHAIN_PAIRWISE_ASSESSMENT_COLUMNS,
    CHAINWISE_ASSESSMENT_COLUMNS
)
from af_pipeline.utils.misc_utils import time_it
from af_pipeline.utils.misc_utils import create_mask

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

    rb_mask: np.ndarray
    """ Binary map of residues in the rigid body."""

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

        assert contact_map.shape == pae.shape, (
            "Contact map and PAE matrix must have the same shape."
        )

        self.rb_dict = rb_dict
        self.num_to_idx = num_to_idx
        self.idx_to_num = idx_to_num
        self.lengths_dict = lengths_dict
        self.save_path = save_path
        self.interchain_mask = create_mask(
            partition_dict=lengths_dict,
            hide_interactions="intra_part",
            masked_value=1,
            unmasked_value=0,
        )
        self.contact_map_mask_2d = np.ma.masked_array(
            contact_map,
            mask=self.interchain_mask,
        )
        # plt.imshow(self.contact_map_mask_2d)
        # plt.show()
        # exit()
        self.contact_map_mask_1d = self.contact_map_mask_2d.any(axis=0)
        # plt.scatter(np.arange(len(self.contact_map_mask_1d)), self.contact_map_mask_1d)
        # plt.show()
        # exit()
        # print(np.unique(self.contact_map_mask_2d, return_counts=True))
        self.plddt_list = np.array(plddt_list)
        self.rb_dict_idxs = self.transform_rb_dict_to_idxs(rb_dict)

        self.symmetric_pae = kwargs.get("symmetric_pae", True)
        self.as_average = kwargs.get("as_average", False)

        self.idr_chains = kwargs.get("idr_chains", [])
        self.protein_chain_map = kwargs.get("protein_chain_map", {})

        self.unique_chains = self.get_unique_chains()
        self.chain_pairs = self.get_chain_pairs()
        self.rb_mask_2d = self.get_rb_mask(lengths_dict, 2) # use for PAE
        self.rb_mask_1d = self.get_rb_mask(lengths_dict, 1) # use for pLDDT

        self.chain_mask_stack_1d = self.get_chain_mask_stack(1)
        self.chain_mask_stack_2d = self.get_chain_mask_stack(2)

        self.chain_pair_mask_stack_1d = self.get_chain_pair_mask_stack(1)
        self.chain_pair_mask_stack_2d = self.get_chain_pair_mask_stack(2)

        # self.chain_contact_masks = self.get_chain_contact_masks()
        # self.chain_pair_contact_masks = self.get_chain_pair_contact_masks()

        self.rb_ch_assess = RigidBodyChainAssessment(
            unique_chains=self.unique_chains,
            idx_to_num=self.idx_to_num,
            as_average=self.as_average,
            chain_mask_stack_1d=self.chain_mask_stack_1d,
            rb_mask_1d=self.rb_mask_1d,
            contact_map_mask_1d=self.contact_map_mask_1d,
            contact_map_mask_2d=self.contact_map_mask_2d,
            plddt_list=self.plddt_list,
            idr_chains=self.idr_chains,
            protein_chain_map=self.protein_chain_map,
        )


    def get_unique_chains(self) -> List[str]:
        """Get unique chains in the rigid body.

        Returns:

        - **unique_chains (list)**:<br />
            List of unique chain IDs in the rigid body.
        """

        unique_chains = sorted([
            chain_id
            for chain_id in self.rb_dict.keys()
            if len(self.rb_dict[chain_id]) > 0
        ])

        return unique_chains

    def get_chain_pairs(self) -> List[Tuple[str, str]]:
        """Get all unique chain pairs in the rigid body.

        Returns:

        - **chain_pairs (list)**:<br />
            List of tuples containing unique chain pairs.
            Each tuple contains two chain IDs.
        """

        chain_pairs = list(combinations(self.unique_chains, 2))

        return sorted([tuple(pair) for pair in chain_pairs])

    def transform_rb_dict_to_idxs(
        self,
        rb_dict: Dict[str, List[Tuple[str, int]]],
    ) -> Dict[str, List[Tuple[str, int]]]:
        """Transform rigid body dictionary from residue numbers to indices.

        Arguments:

        - **rb_dict (dict)**:<br />
            Dictionary of rigid bodies with residue numbers.

        Returns:

        - **rb_dict_idxs (dict)**:<br />
            Dictionary of rigid bodies with residue indices.
        """

        rb_dict_idxs = {
            chain_id: [
                self.num_to_idx[chain_id][token_num][atom_name]
                for atom_name, token_num in rb_dict[chain_id]
            ]
            for chain_id in rb_dict.keys()
        }

        return rb_dict_idxs

    def get_rb_mask(
        self,
        lengths_dict: Dict[str, int],
        dimensions: int = 2,
    ) -> np.ndarray:
        """Get a binary map of residues in the rigid body.

        Arguments:

        - **lengths_dict (dict)**:<br />
            Dictionary containing lengths of chains in the structure.

        Returns:

        - **rb_mask (np.ndarray)**:<br />
            A binary map of residues in the rigid body.
            The shape is (`total_length`, `total_length`) where `total_length`
            is the sum of lengths of all chains.
            The value is 1 if the residue is part of the rigid body, 0 otherwise.
        """

        total_len = lengths_dict.get("total", 0)

        if dimensions == 1:
            rb_mask = np.zeros((total_len,), dtype=int)
        elif dimensions == 2:
            rb_mask = np.zeros((total_len, total_len), dtype=int)

        rb_res_idxs = [
            token_idx
            for chain_id in self.rb_dict_idxs.keys()
            for token_idx in self.rb_dict_idxs[chain_id]
        ]

        if dimensions == 1:
            rb_mask[rb_res_idxs] = 1
        elif dimensions == 2:
            rb_mask[np.ix_(rb_res_idxs, rb_res_idxs)] = 1

        rb_mask = np.ma.make_mask(rb_mask)

        return rb_mask

    def get_chain_mask(
        self,
        chain_id: str,
        lengths_dict: Dict[str, int],
        dimensions: int = 2,
    ) -> np.ndarray:
        """Get a binary map of residues in a chain.

        Arguments:

        - **chain_id (str)**:<br />
            Chain ID for which the mask is to be generated.

        - **lengths_dict (dict)**:<br />
            Dictionary containing lengths of chains in the structure.

        Returns:

        - **chain_mask (np.ndarray)**:<br />
            A binary map of residues in the chain.
            The shape is (`total_length`, `total_length`) where `total_length`
            is the sum of lengths of all chains.
            The value is 1 if the residue is part of the chain, 0 otherwise.
        """

        total_len = lengths_dict.get("total", 0)

        if dimensions == 1:
            chain_mask = np.zeros((total_len,), dtype=int)
        elif dimensions == 2:
            chain_mask = np.zeros((total_len, total_len), dtype=int)

        chain_res_idxs = self.rb_dict_idxs.get(chain_id, [])

        if dimensions == 1:
            chain_mask[chain_res_idxs] = 1
        elif dimensions == 2:
            chain_mask[np.ix_(chain_res_idxs, chain_res_idxs)] = 1

        chain_mask = np.ma.make_mask(chain_mask)

        return chain_mask

    def get_chain_pair_mask(
        self,
        chain_id_1: str,
        chain_id_2: str,
        lengths_dict: Dict[str, int],
        dimensions: int = 2,
    ) -> np.ndarray:
        """Get a binary map of residues in a chain pair.

        Arguments:

        - **chain_id_1 (str)**:<br />
            First chain ID of the chain pair.

        - **chain_id_2 (str)**:<br />
            Second chain ID of the chain pair.

        - **lengths_dict (dict)**:<br />
            Dictionary containing lengths of chains in the structure.

        Returns:

        - **chain_pair_mask (np.ndarray)**:<br />
            A binary map of residues in the chain pair.
            The shape is (`total_length`, `total_length`) where `total_length`
            is the sum of lengths of all chains.
            The value is 1 if the residue is part of either chain, 0 otherwise.
        """

        total_len = lengths_dict.get("total", 0)

        if dimensions == 1:
            chain_pair_mask = np.zeros((total_len,), dtype=int)
        elif dimensions == 2:
            chain_pair_mask = np.zeros((total_len, total_len), dtype=int)

        chain_1_res_idxs = self.rb_dict_idxs.get(chain_id_1, [])
        chain_2_res_idxs = self.rb_dict_idxs.get(chain_id_2, [])

        if dimensions == 1:
            chain_pair_mask[chain_1_res_idxs] = 1
            chain_pair_mask[chain_2_res_idxs] = 1
        elif dimensions == 2:
            chain_pair_mask[np.ix_(chain_1_res_idxs, chain_2_res_idxs)] = 1
            chain_pair_mask[np.ix_(chain_2_res_idxs, chain_1_res_idxs)] = 1

        chain_pair_mask = np.ma.make_mask(chain_pair_mask)

        return chain_pair_mask

    def get_chain_mask_stack(self, dimensions: int = 2):
        """Get a stack of binary maps for all chains in the rigid body.

        Returns:

        - **chain_mask_stack (np.ndarray)**:<br />
            A stack of binary maps for all chains in the rigid body.
            The shape is (`num_chains`, `total_length`, `total_length`)
            where `num_chains` is the number of unique chains in the rigid body
            and `total_length` is the sum of lengths of all chains.
            The value is 1 if the residue is part of the chain, 0 otherwise.
        """

        chain_mask_stack = np.array([
            self.get_chain_mask(
                chain_id,
                self.lengths_dict,
                dimensions=dimensions,
            )
            for chain_id in self.unique_chains
        ])

        return chain_mask_stack

    def get_chain_pair_mask_stack(self, dimensions: int = 2):
        """Get a stack of binary maps for all chain pairs in the rigid body.

        Returns:

        - **chain_pair_mask_stack (np.ndarray)**:<br />
            A stack of binary maps for all chain pairs in the rigid body.
            The shape is (`num_chain_pairs`, `total_length`, `total_length`)
            where `num_chain_pairs` is the number of unique chain pairs in the
            rigid body and `total_length` is the sum of lengths of all chains.
            The value is 1 if the residue is part of either chain, 0 otherwise.
        """

        chain_pair_mask_stack = np.array([
            self.get_chain_pair_mask(
                chain_id_1,
                chain_id_2,
                self.lengths_dict,
                dimensions=dimensions,
            )
            for chain_id_1, chain_id_2 in self.chain_pairs
        ])

        return chain_pair_mask_stack

    def save_rb_assessment(self):
        """ Save the assessment of the rigid bodies to an Excel file.

        The assessment includes:
        - **Per chain assessment**:<br />
            Average pLDDT, Average iLDDT, interface residues count,
            Chain type (IDR or R).
        - **Per chain pair assessment**:<br />
            interface residues count, Number of contacts, Average PAE,
            Average iPAE, Minimum PAE, Average iLDDT for each chain,
            Chain type (IDR or R) for each chain.
        - **Overall assessment**:<br />
            Average pLDDT, Average iLDDT, interface residues count,
            Chain type (IDR or R).

        The assessment is saved in an Excel file with three sheets:
        - "Chain Wise Assessment": Contains per chain assessment data.
        - "Chain Pairwise Assessment": Contains per chain pair assessment data.
        - "Overall Assessment": Contains overall assessment data.
        """

        chain_wise_assessment_rows = []
        chain_pairwise_assessment_rows = []
        overall_assessment_rows = []

        if self.as_average:

            iterators = [(chain_id,) for chain_id in self.unique_chains]
            func = self.rb_ch_assess.get_chain_attr

        else:

            iterators = [
                (ch_id, res_idx)
                for ch_id in self.unique_chains
                for res_idx in self.rb_ch_assess.per_chain_interface_res[ch_id]
            ]
            func = self.rb_ch_assess.get_res_attr

        chain_wise_assessment_rows.append({
            k: func(*it, k)
            for it in iterators
            for k in CHAINWISE_ASSESSMENT_COLUMNS[self.as_average]
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


    # def get_chain_contact_masks(self, dimensions: int = 2) -> np.ndarray:
    #     """Get contact masks for all chains in the rigid body.

    #     Returns:

    #     - **chain_contact_masks (np.ndarray)**:<br />
    #         A stack of contact masks for all chains in the rigid body.
    #         The shape is (`num_chains`, `total_length`, `total_length`)
    #         where `num_chains` is the number of unique chains in the rigid body
    #         and `total_length` is the sum of lengths of all chains.
    #         The value is 1 if the residue pair is in contact, 0 otherwise.
    #     """

    #     if dimensions == 1:
    #         chain_mask_stack = self.chain_mask_stack_1d
    #         contact_map_mask = self.contact_map_mask_1d

    #     elif dimensions == 2:
    #         chain_mask_stack = self.chain_mask_stack_2d
    #         contact_map_mask = self.contact_map_mask_2d

    #     chain_contact_masks = np.ma.array(
    #         chain_mask_stack,
    #         mask=np.broadcast_to(contact_map_mask, chain_mask_stack.shape)
    #     )

    #     return chain_contact_masks

    # def get_chain_pair_contact_masks(self, dimensions: int = 2) -> np.ndarray:
    #     """Get contact masks for all chain pairs in the rigid body.

    #     Returns:

    #     - **chain_pair_contact_masks (np.ndarray)**:<br />
    #         A stack of contact masks for all chain pairs in the rigid body.
    #         The shape is (`num_chain_pairs`, `total_length`, `total_length`)
    #         where `num_chain_pairs` is the number of unique chain pairs in the
    #         rigid body and `total_length` is the sum of lengths of all chains.
    #         The value is 1 if the residue pair is in contact, 0 otherwise.
    #     """

    #     if dimensions == 1:
    #         cp_mask_stack = self.chain_pair_mask_stack_1d
    #         contact_map_mask = self.contact_map_mask_1d

    #     elif dimensions == 2:
    #         cp_mask_stack = self.chain_pair_mask_stack_2d
    #         contact_map_mask = self.contact_map_mask_2d

    #     chain_pair_contact_masks = np.ma.array(
    #         cp_mask_stack,
    #         mask=np.broadcast_to(contact_map_mask, cp_mask_stack.shape)
    #     )

    #     return chain_pair_contact_masks

class RigidBodyChainAssessment:

    def __init__(
        self,
        unique_chains: List[str],
        idx_to_num: Dict[int, Dict[str, str|int]],
        as_average: bool,
        chain_mask_stack_1d: np.ndarray,
        rb_mask_1d: np.ndarray,
        contact_map_mask_1d: np.ndarray,
        contact_map_mask_2d: np.ndarray,
        plddt_list: np.ndarray,
        idr_chains: List[str] = [],
        protein_chain_map: Dict[str, str] = {},
    ):

        self.unique_chains = unique_chains
        self.idx_to_num = idx_to_num
        self.as_average = as_average
        self.chain_mask_stack_1d = chain_mask_stack_1d
        self.rb_mask_1d = rb_mask_1d
        self.contact_map_mask_1d = contact_map_mask_1d
        self.contact_map_mask_2d = contact_map_mask_2d
        self.plddt_list = plddt_list
        self.idr_chains = idr_chains
        self.protein_chain_map = protein_chain_map

        self.per_chain_plddt = self.get_per_chain_plddt(self.as_average)
        self.per_chain_iplddt = self.get_per_chain_iplddt(self.as_average)

        self.per_chain_interface_res = self.get_per_chain_interface_residues(
            self.as_average
        )

    def get_per_chain_plddt(self, only_avg: bool = True) -> Dict[str, float]:
        """Get average pLDDT score per chain in the rigid body.

        Returns:

        - **per_chain_plddt (dict)**:<br />
            Dictionary with chain IDs as keys and average pLDDT scores
            as values.
        """

        per_chain_plddt = {}

        for i, ch_id in enumerate(self.unique_chains):
            chain_mask = self.chain_mask_stack_1d[i, :]
            plddt_mask = chain_mask * self.rb_mask_1d
            if only_avg:
                per_chain_plddt[ch_id] = np.mean(self.plddt_list[plddt_mask])
                continue
            per_chain_plddt[ch_id] = self.plddt_list[plddt_mask].tolist()

        # print(per_chain_plddt)

        return per_chain_plddt

    def get_per_chain_iplddt(self, only_avg: bool = True) -> Dict[str, float]:
        """Get average inter-chain pLDDT score per chain in the rigid body.

        Returns:

        - **per_chain_iplddt (dict)**:<br />
            Dictionary with chain IDs as keys and average inter-chain
            pLDDT scores as values.
        """

        per_chain_iplddt = {}

        for i, ch_id in enumerate(self.unique_chains):
            chain_mask = self.chain_mask_stack_1d[i, :]
            iplddt_mask = chain_mask * self.rb_mask_1d * self.contact_map_mask_1d
            if only_avg:
                per_chain_iplddt[ch_id] = np.mean(self.plddt_list[iplddt_mask])
                continue
            per_chain_iplddt[ch_id] = self.plddt_list[iplddt_mask].tolist()

        # print(per_chain_iplddt)

        return per_chain_iplddt

    def get_per_chain_interface_residues(
        self,
        only_count: bool = True
    ) -> Dict[str, List[int]]:
        """Get interface residues per chain in the rigid body.

        Returns:

        - **per_chain_interface_residues (dict)**:<br />
            Dictionary with chain IDs as keys and list of residue numbers
            at the interface as values.
        """

        per_chain_interface_residues = {}

        for i, ch_id in enumerate(self.unique_chains):
            chain_mask = self.chain_mask_stack_1d[i, :]
            interface_mask = (
                chain_mask
                * self.rb_mask_1d
                * self.contact_map_mask_1d
            )

            if only_count:
                res_count = np.sum(interface_mask)
                per_chain_interface_residues[ch_id] = res_count
                continue

            res_idxs = np.where(interface_mask)[0]
            per_chain_interface_residues[ch_id] = res_idxs

        return per_chain_interface_residues

    def get_chain_type(self, chain_id: str) -> str:
        """Get the type of a chain (IDR or ordered).

        Arguments:

        - **chain_id (str)**:<br />
            Chain ID for which the type is to be determined.

        Returns:

        - **chain_type (str)**:<br />
            "IDR" if the chain is disordered, "ordered" otherwise.
        """

        chain_type = "IDR" if chain_id in self.idr_chains else "R"

        return chain_type

    def get_protein_name(self, chain_id: str) -> str:
        """Get the protein name for a chain.

        Arguments:

        - **chain_id (str)**:<br />
            Chain ID for which the protein name is to be determined.

        Returns:

        - **protein_name (str)**:<br />
            Protein name corresponding to the chain ID.
            If not found, returns the chain ID itself.
        """

        protein_name = self.protein_chain_map.get(chain_id, chain_id)

        return protein_name

    def get_chain_attr(
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

        if attr_name not in CHAINWISE_ASSESSMENT_COLUMNS[self.as_average]:
            raise ValueError(
                f"Attribute '{attr_name}' not recognized for chain "
                f"'{chain_id}'."
            )

        attrs_ = {
            "Chain ID": chain_id,
            "Protein Name": self.get_protein_name(chain_id),
            "Average pLDDT": self.per_chain_plddt[chain_id],
            "Average ipLDDT": self.per_chain_iplddt[chain_id],
            "Interface Residues": self.per_chain_interface_res[chain_id],
            "Chain Type": self.get_chain_type(chain_id),
        }

        return attrs_[attr_name]

    def get_res_attr(
        self,
        chain_id: str,
        res_idx: int,
        attr_name: str
    ):
        """ Get the attribute value for a given residue in a chain.

        Arguments:

        - **chain_id (str)**:<br />
            Chain ID.

        - **res_num (int)**:<br />
            Residue number.

        Returns:

        - **(float | str | int)**:<br />
            The requested attribute value.
        """

        res_num = self.idx_to_num[res_idx]["token_num"]

        if attr_name not in CHAINWISE_ASSESSMENT_COLUMNS[self.as_average]:
            raise ValueError(
                f"Attribute '{attr_name}' not recognized for residue "
                f"'{res_num}' in chain '{chain_id}'."
            )

        local_idx = list(self.per_chain_interface_res[chain_id]).index(res_idx)

        attrs_ = {
            "Chain ID": chain_id,
            "Residue Number": res_num,
            "Protein Name": self.get_protein_name(chain_id),
            "Average pLDDT": self.per_chain_plddt[chain_id][local_idx],
            "Average ipLDDT": self.per_chain_iplddt[chain_id][local_idx],
            "Chain Type": self.get_chain_type(chain_id),
        }

        return attrs_[attr_name]