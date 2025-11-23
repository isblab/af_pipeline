from typing import Dict, List, Tuple
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from itertools import combinations
from collections import defaultdict
from tqdm import tqdm
import warnings
import copy
from af_pipeline.constants.af_constants import (
    CHAIN_PAIRWISE_ASSESSMENT_COLUMNS,
    CHAINWISE_ASSESSMENT_COLUMNS,
    OVERALL_ASSESSMENT_COLUMNS,
)
from af_pipeline.utils.misc_utils import (
    time_it,
    create_mask,
    symmetrize_matrix,
)

class _Mask:

    def __init__(
        self,
        rb_dict: dict,
        num_to_idx: dict,
        idx_to_num: dict,
        contact_map: np.ndarray,
        plddt_list: np.ndarray,
        pae: np.ndarray,
        avg_pae: np.ndarray,
        lengths_dict: dict,
        as_average: bool = True,
        symmetric_pae: bool = True,
        idr_chains: List[str] = [],
        protein_chain_map: Dict[str, str] = {},
    ):

        assert contact_map.shape == pae.shape, (
            "Contact map and PAE matrix must have the same shape."
        )

        self.rb_dict = rb_dict
        self.pae = pae
        self.avg_pae = avg_pae
        self.num_to_idx = num_to_idx
        self.idx_to_num = idx_to_num
        self.lengths_dict = lengths_dict

        self.symmetric_pae = symmetric_pae
        self.as_average = as_average

        self.idr_chains = idr_chains
        self.protein_chain_map = protein_chain_map

        self.interchain_mask = create_mask(
            partition_dict=lengths_dict,
            hide_interactions="intra_part",
            masked_value=0,
            unmasked_value=1,
        )
        # self.contact_map_mask_2d = np.ma.masked_array(
        #     contact_map,
        #     mask=self.interchain_mask,
        # )
        self.contact_map_mask_2d = np.ma.make_mask(
            contact_map * self.interchain_mask
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

        self.unique_chains = self.get_unique_chains()
        self.chain_pairs = self.get_chain_pairs()

        self.rb_mask_2d = self.get_rb_mask(lengths_dict, 2) # use for PAE
        self.rb_mask_1d = self.get_rb_mask(lengths_dict, 1) # use for pLDDT

        self.chain_mask_stack_1d = self.get_chain_mask_stack(1)
        self.chain_mask_stack_2d = self.get_chain_mask_stack(2)

        self.chain_pair_mask_stack_1d = self.get_chain_pair_mask_stack(1)
        self.chain_pair_mask_stack_2d = self.get_chain_pair_mask_stack(2)

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

class RigidBodyChainAssessment:

    def __init__(
        self,
        _mask: _Mask,
    ):

        self.unique_chains = _mask.unique_chains
        self.idx_to_num = _mask.idx_to_num
        self.as_average = _mask.as_average
        self.chain_mask_stack_1d = _mask.chain_mask_stack_1d
        self.rb_mask_1d = _mask.rb_mask_1d
        self.contact_map_mask_1d = _mask.contact_map_mask_1d
        self.plddt_list = _mask.plddt_list
        self.idr_chains = _mask.idr_chains
        self.protein_chain_map = _mask.protein_chain_map

        self.per_chain_plddt = self.get_per_chain_plddt(self.as_average, False)
        self.per_chain_iplddt = self.get_per_chain_plddt(self.as_average, True)
        self.per_chain_interface_res = self.get_per_chain_interface_residues(
            self.as_average
        )

    def get_per_chain_plddt(
        self,
        only_avg: bool = True,
        only_interface: bool = False,
    ) -> Dict[str, list | float]:
        """Get average pLDDT score per chain in the rigid body.

        Returns:

        - **per_chain_plddt (dict)**:<br />
            Dictionary with chain IDs as keys and average pLDDT scores
            as values.
        """

        per_chain_plddt = {}

        for i, ch_id in enumerate(self.unique_chains):
            chain_mask_1d = self.chain_mask_stack_1d[i, :]
            plddt_mask_1d = chain_mask_1d * self.rb_mask_1d
            if only_interface:
                plddt_mask_1d = plddt_mask_1d * self.contact_map_mask_1d
            if only_avg:
                per_chain_plddt[ch_id] = np.mean(self.plddt_list[plddt_mask_1d])
                continue
            per_chain_plddt[ch_id] = self.plddt_list[plddt_mask_1d].tolist()

        # print(per_chain_plddt)

        return per_chain_plddt

    def get_per_chain_interface_residues(
        self,
        only_count: bool = True
    ) -> Dict[str, List[int]]:
        """Get interface residues per chain in the rigid body.

        Returns:

        - **per_chain_interface_res (dict)**:<br />
            Dictionary with chain IDs as keys and list of residue numbers
            at the interface as values.
        """

        per_chain_interface_res = {}

        for i, ch_id in enumerate(self.unique_chains):
            chain_mask_1d = self.chain_mask_stack_1d[i, :]
            interface_mask_1d = (
                chain_mask_1d
                * self.rb_mask_1d
                * self.contact_map_mask_1d
            )

            if only_count:
                res_count = np.sum(interface_mask_1d)
                per_chain_interface_res[ch_id] = res_count
                continue

            res_idxs = np.where(interface_mask_1d)[0]
            per_chain_interface_res[ch_id] = res_idxs

        return per_chain_interface_res

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
            "Protein Name": self.protein_chain_map.get(chain_id, chain_id),
            "Average pLDDT": self.per_chain_plddt[chain_id],
            "Average ipLDDT": self.per_chain_iplddt[chain_id],
            "Interface Residues": self.per_chain_interface_res[chain_id],
            "Chain Type": "IDR" if chain_id in self.idr_chains else "R",
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
            "Protein Name": self.protein_chain_map.get(chain_id, chain_id),
            "Average pLDDT": self.per_chain_plddt[chain_id][local_idx],
            "Average ipLDDT": self.per_chain_iplddt[chain_id][local_idx],
            "Chain Type": "IDR" if chain_id in self.idr_chains else "R",
        }

        return attrs_[attr_name]

    def get_chain_assessment(self):
        """Get chain-wise assessment for the rigid body.

        Returns:

        - **chain_wise_assessment (list)**:<br />
            List of dictionaries containing chain-wise assessment.
        """

        chain_wise_assessment_rows = []
        if self.as_average:

            iterators = [(chain_id,) for chain_id in self.unique_chains]
            func = self.get_chain_attr
        else:

            iterators = [
                (ch_id, res_idx)
                for ch_id in self.unique_chains
                for res_idx in self.per_chain_interface_res[ch_id]
            ]
            func = self.get_res_attr

        chain_wise_assessment_rows = [
            {
                k: func(*it, k)
                for k in CHAINWISE_ASSESSMENT_COLUMNS[self.as_average]
            } for it in iterators
        ]

        return pd.DataFrame(chain_wise_assessment_rows)

class RigidBodyChainPairAssessment:

    def __init__(
        self,
        _mask: _Mask,
    ):

        self.unique_chains = _mask.unique_chains
        self.chain_pairs = _mask.chain_pairs
        self.idx_to_num = _mask.idx_to_num
        self.as_average = _mask.as_average
        self.symmetric_pae = _mask.symmetric_pae

        self.chain_mask_stack_1d = _mask.chain_mask_stack_1d
        self.chain_mask_stack_2d = _mask.chain_mask_stack_2d

        self.chain_pair_mask_stack_1d = _mask.chain_pair_mask_stack_1d
        self.chain_pair_mask_stack_2d = _mask.chain_pair_mask_stack_2d

        self.rb_mask_1d = _mask.rb_mask_1d
        self.rb_mask_2d = _mask.rb_mask_2d

        self.contact_map_mask_1d = _mask.contact_map_mask_1d
        self.contact_map_mask_2d = _mask.contact_map_mask_2d

        self.plddt_list = _mask.plddt_list
        self.pae = _mask.pae
        self.avg_pae = _mask.avg_pae
        self.symmetric_pae = _mask.symmetric_pae
        self.idr_chains = _mask.idr_chains
        self.protein_chain_map = _mask.protein_chain_map

        self.chain_pair_interface_res = self.get_chain_pair_interface(
            self.as_average, False,
        )
        self.chain_pair_contacts = self.get_chain_pair_interface(
            self.as_average, True,
        )
        # self.chain_pair_plddt = self.get_chain_pair_plddt(
        #     self.as_average, False,
        # )
        self.chain_pair_iplddt = self.get_chain_pair_plddt(
            self.as_average, True,
        )

        self.chain_pair_pae = self.get_chain_pair_pae(
            self.as_average, False, (True, "")
        )
        self.chain_pair_pae_ij = self.get_chain_pair_pae(
            self.as_average, False, (False, "ij")
        )
        self.chain_pair_pae_ji = self.get_chain_pair_pae(
            self.as_average, False, (False, "ji")
        )

        self.chain_pair_ipae = self.get_chain_pair_pae(
            self.as_average, True, (True, "")
        )
        self.chain_pair_ipae_ij = self.get_chain_pair_pae(
            self.as_average, True, (False, "ij")
        )
        self.chain_pair_ipae_ji = self.get_chain_pair_pae(
            self.as_average, True, (False, "ji")
        )

    def get_chain_pair_attr(
        self,
        chain_pair: Tuple[str, str],
        attr_name: str,
    ):
        # print(chain_pair)
        chain1, chain2 = chain_pair
        attr_state = (self.as_average, self.symmetric_pae)
        if attr_name not in CHAIN_PAIRWISE_ASSESSMENT_COLUMNS[attr_state]:
            raise ValueError(f"Invalid attribute name: {attr_name}")

        attrs_ = {
            "Chain ID": (chain1, chain2),

            "Protein Name": (
                self.protein_chain_map.get(chain1, chain1),
                self.protein_chain_map.get(chain2, chain2),
            ),

            "Chain Type": (
                "IDR" if chain1 in self.idr_chains else "R",
                "IDR" if chain2 in self.idr_chains else "R",
            ),

            "Interface Residues": self.chain_pair_interface_res[chain_pair],
            "Number of contacts": self.chain_pair_contacts[chain_pair],

            # "pLDDT Chain 1": self.chain_pair_plddt[chain_pair][0],
            # "pLDDT Chain 2": self.chain_pair_plddt[chain_pair][1],

            "Average ipLDDT chain1": self.chain_pair_iplddt[chain_pair][0],
            "Average ipLDDT chain2": self.chain_pair_iplddt[chain_pair][1],

            "Average PAE": self.chain_pair_pae[chain_pair],
            "Average iPAE": self.chain_pair_ipae[chain_pair],

            "Average PAE ij": self.chain_pair_pae_ij[chain_pair],
            "Average PAE ji": self.chain_pair_pae_ji[chain_pair],

            "Average iPAE ij": self.chain_pair_ipae_ij[chain_pair],
            "Average iPAE ji": self.chain_pair_ipae_ji[chain_pair],
        }

        return attrs_[attr_name]

    def get_res_pair_attr(
        self,
        chain_pair: Tuple[str, str],
        res_pair: Tuple[int, int],
        attr_name: str,
    ):
        attr_state = (self.as_average, self.symmetric_pae)
        if attr_name not in CHAIN_PAIRWISE_ASSESSMENT_COLUMNS[attr_state]:
            raise ValueError(f"Invalid attribute name: {attr_name}")
        # print(chain_pair, res_pair)
        ch_id_1, ch_id_2 = chain_pair
        res_idx_1, res_idx_2 = res_pair
        attrs_ = {
            "Chain ID 1": ch_id_1,
            "Chain ID 2": ch_id_2,
            "Residue 1": self.idx_to_num[res_idx_1]["token_num"],
            "Residue 2": self.idx_to_num[res_idx_2]["token_num"],
            "pLDDT 1": self.plddt_list[res_idx_1],
            "pLDDT 2": self.plddt_list[res_idx_2],
            "ipLDDT 1": self.plddt_list[res_idx_1],
            "ipLDDT 2": self.plddt_list[res_idx_2],
        }

        return attrs_[attr_name]

    def get_chain_pair_interface(
        self,
        only_count: bool = True,
        only_contacts: bool = False,
    ) -> Dict[str, List[int]]:
        """Get interface residues for the chain pair in the rigid body.

        Returns:

        - **chain_pair_interface_residues (dict)**:<br />
            Dictionary with chain IDs as keys and list of residue numbers
            at the interface as values.
        """

        chain_pair_interface_residues = {}

        for i, (ch_id_1, ch_id_2) in enumerate(self.chain_pairs):
            chain_pair_mask = self.chain_pair_mask_stack_2d[i, :, :]
            interface_mask = (
                chain_pair_mask
                * self.rb_mask_2d
                * self.contact_map_mask_2d
            )
            if only_contacts:
                num_contacts = np.sum(interface_mask)
                chain_pair_interface_residues[(ch_id_1, ch_id_2)] = num_contacts
                continue

            interface_mask_1d = interface_mask.any(axis=0)
            ch1_idx = self.unique_chains.index(ch_id_1)
            chain1_mask = interface_mask_1d * self.chain_mask_stack_1d[ch1_idx, :]
            ch2_idx = self.unique_chains.index(ch_id_2)
            chain2_mask = interface_mask_1d * self.chain_mask_stack_1d[ch2_idx, :]

            if only_count:
                res_count1 = np.sum(chain1_mask)
                res_count2 = np.sum(chain2_mask)
                chain_pair_interface_residues[(ch_id_1, ch_id_2)] = (res_count1, res_count2)
                continue

            res_idx_pairs = np.argwhere(interface_mask)
            chain_pair_interface_residues[(ch_id_1, ch_id_2)] = res_idx_pairs

        return chain_pair_interface_residues

    def get_chain_pair_plddt(
        self,
        only_avg: bool = True,
        only_interface: bool = False,
    ):

        chain_pair_plddt = {}

        for i, (ch_id_1, ch_id_2) in enumerate(self.chain_pairs):
            chain_pair_mask = self.chain_pair_mask_stack_2d[i, :, :]
            plddt_mask = chain_pair_mask * self.rb_mask_2d

            if only_interface:
                plddt_mask = plddt_mask * self.contact_map_mask_2d

            plddt_mask_1d = plddt_mask.any(axis=0)
            ch1_idx = self.unique_chains.index(ch_id_1)
            chain1_mask = plddt_mask_1d * self.chain_mask_stack_1d[ch1_idx, :]
            ch2_idx = self.unique_chains.index(ch_id_2)
            chain2_mask = plddt_mask_1d * self.chain_mask_stack_1d[ch2_idx, :]

            if only_avg:
                chain_pair_plddt[(ch_id_1, ch_id_2)] = (
                    np.mean(self.plddt_list[chain1_mask]),
                    np.mean(self.plddt_list[chain2_mask])
                )
                continue

            res_pair_idxs = np.argwhere(plddt_mask)
            ch1_plddts = self.plddt_list[res_pair_idxs[:, 0]]
            ch2_plddts = self.plddt_list[res_pair_idxs[:, 1]]

            chain_pair_plddt[(ch_id_1, ch_id_2)] = list(zip(ch1_plddts, ch2_plddts))

        return chain_pair_plddt

    # @time_it
    def get_chain_pair_pae(
        self,
        only_avg: bool = True,
        only_interface: bool = False,
        symmetric: Tuple[bool, str] = (True, ""),
    ):
        chain_pair_pae = {}
        sym_bool, sym_type = symmetric

        assert symmetric in [(True, ""), (False, "ij"), (False, "ji")], (
            """`symmetric` must be one of the following.
            (True, ''), (False, 'ij'), (False, 'ji')
            """
        )

        # (sym_bool, sym_type)
        pae_attrs_ = {
            (True, ""): lambda a: self.avg_pae,
            (False, "ij"): np.triu,
            (False, "ji"): np.tril,
        }
        # only_interface
        mask_attrs_ = {
            True: self.rb_mask_2d * self.contact_map_mask_2d,
            False: self.rb_mask_2d,
        }
        # sym_bool
        sym_bool_attrs_ = {
            True: lambda a: a,
            False: lambda a: a + a.T - np.diag(np.diag(a)),
        }

        pae_matrix = pae_attrs_[(sym_bool, sym_type)](self.pae)
        pae_matrix = sym_bool_attrs_[sym_bool](pae_matrix)
        pae_mask_multiplier = mask_attrs_[only_interface]

        for i, (ch_id_1, ch_id_2) in enumerate(self.chain_pairs):

            chain_pair_mask = self.chain_pair_mask_stack_2d[i, :, :]
            pae_mask_2d = chain_pair_mask * pae_mask_multiplier
            pae_idxs = np.where(pae_mask_2d)
            pae_vals = pae_matrix[pae_idxs].flatten()

            if only_avg:
                chain_pair_pae[(ch_id_1, ch_id_2)] = np.mean(pae_vals)
                continue

            chain_pair_pae[(ch_id_1, ch_id_2)] = pae_vals.tolist()

        return chain_pair_pae

    def get_chain_pair_assessment(self):
        """Get chain-pair-wise assessment for the rigid body.

        Returns:
        - **chain_pairwise_assessment (list)**:<br />
            List of dictionaries containing chain-pair-wise assessment.
        """

        chain_pairwise_assessment_rows = []


        if self.as_average:

            iterators_cp = [(chain_pair,) for chain_pair in self.chain_pairs]
            func_cp = self.get_chain_pair_attr
            # print(iterators_cp)

        else:

            iterators_cp = [
                (chain_pair, res_pair)
                for chain_pair in self.chain_pairs
                for res_pair in self.chain_pair_interface_res[chain_pair]
            ]
            func_cp = self.get_res_pair_attr

        chain_pairwise_assessment_rows = [
            {
                k: func_cp(*it, k)
                for k in CHAIN_PAIRWISE_ASSESSMENT_COLUMNS[(self.as_average, self.symmetric_pae)]
            } for it in iterators_cp
        ]
        return pd.DataFrame(chain_pairwise_assessment_rows)

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
        avg_pae: np.ndarray,
        lengths_dict: dict,
        save_path: str,
        **kwargs,
    ):

        self.save_path = save_path

        _mask = _Mask(
            rb_dict=rb_dict,
            num_to_idx=num_to_idx,
            idx_to_num=idx_to_num,
            contact_map=contact_map,
            plddt_list=plddt_list,
            pae=pae,
            avg_pae=avg_pae,
            lengths_dict=lengths_dict,
            as_average=kwargs.get("as_average", False),
            symmetric_pae=kwargs.get("symmetric_pae", True),
            idr_chains=kwargs.get("idr_chains", []),
            protein_chain_map=kwargs.get("protein_chain_map", {}),
        )

        self.rb_ch_assess = RigidBodyChainAssessment(_mask=_mask)
        self.rb_ch_pair_assess = RigidBodyChainPairAssessment(_mask=_mask)
        self.overall_assessment = self.get_overall_assessment()

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

        overall_assessment_rows = []

        for col_head, key in OVERALL_ASSESSMENT_COLUMNS.items():
            if self.overall_assessment.get(key, np.nan) is not np.nan:
                overall_assessment_rows.append({
                    "Key": col_head,
                    "Value": self.overall_assessment.get(key)
                })
        overall_assessment_df = pd.DataFrame(overall_assessment_rows)

        c_assessment_df = self.rb_ch_assess.get_chain_assessment()
        cp_assessment_df = self.rb_ch_pair_assess.get_chain_pair_assessment()

        df_dict = {
            "chain_pairwise_assessment": cp_assessment_df,
            "chainwise_assessment": c_assessment_df,
            "overall_assessment": overall_assessment_df,
        }

        for k, df_ in df_dict.items():
            df_dict[k] = df_.fillna(np.nan)
            df_dict[k] = df_.map(
                lambda x: round(x, 2) if isinstance(x, (int, float)) else x
            )

        with pd.ExcelWriter(
            self.save_path, engine='openpyxl', mode='w'
        ) as writer:

            for sheet_name, df in df_dict.items():

                if df.empty:
                    warnings.warn(
                        f"Skipping empty DataFrame for sheet: {sheet_name}"
                    )
                    continue

                df.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    index=False,
                )


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

        overall_assessment["num_chains"] = len(self.rb_ch_assess.unique_chains)

        overall_assessment["num_interacting_chain_pairs"] = len([
            pair
            for pair in self.rb_ch_pair_assess.chain_pairs
            if len(self.rb_ch_pair_assess.chain_pair_interface_res[pair]) > 0
        ])

        overall_assessment["num_interface_residues"] = sum(
            self.rb_ch_assess.per_chain_interface_res.values()
        )

        overall_assessment["num_contacts"] = sum(
            self.rb_ch_pair_assess.chain_pair_contacts.values()
        )

        # global_iplddt_scores = [
        #     iplddt
        #     for iplddt_scores in self.rb_ch_assess.per_chain_iplddt.values()
        #     if isinstance(iplddt_scores, list)
        #     for iplddt in iplddt_scores
        # ]

        # overall_assessment["avg_iplddt"] = (
        #     np.mean(global_iplddt_scores) if global_iplddt_scores else np.nan
        # )

        # global_idr_iplddt_scores = [
        #     iplddt
        #     for chain_id, iplddt_scores in self.per_chain_iplddt.items()
        #     for iplddt in iplddt_scores.values()
        #     if chain_id in self.idr_chains
        # ]

        # overall_assessment["avg_idr_iplddt"] = (
        #     np.mean(global_idr_iplddt_scores) if global_idr_iplddt_scores else np.nan
        # )

        # global_ipae_ij_scores = [
        #     ipae
        #     for ipae_dict in self.pairwise_ipae.values()
        #     for ipae in ipae_dict["ij"].values()
        # ]

        # global_ipae_ji_scores = [
        #     ipae
        #     for ipae_dict in self.pairwise_ipae.values()
        #     for ipae in ipae_dict["ji"].values()
        # ]

        # overall_assessment["avg_ipae_ij"] = (
        #     np.mean(global_ipae_ij_scores) if global_ipae_ij_scores else np.nan
        # )

        # overall_assessment["avg_ipae_ji"] = (
        #     np.mean(global_ipae_ji_scores) if global_ipae_ji_scores else np.nan
        # )

        return overall_assessment