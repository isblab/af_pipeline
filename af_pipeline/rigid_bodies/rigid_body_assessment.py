from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from itertools import combinations
from matplotlib import pyplot as plt
from typing import Dict, List, Tuple
from af_pipeline.constants.af_constants import (
    CHAIN_PAIRWISE_ASSESSMENT_COLUMNS,
    CHAINWISE_ASSESSMENT_COLUMNS,
    OVERALL_ASSESSMENT_COLUMNS,
)
from af_pipeline.utils.misc_utils import (
    time_it,
    create_mask,
)
from af_pipeline.tools.structure_tools import get_interaction_map
from af_pipeline.initialize import Initialize
from af_pipeline.constants.af_constants import InteractionConstants as IntCons

_error_not_set_up="""
The RigidBodyAssessment instance is not set up yet. Please set up the instance
by calling the `set_fromInitializer` method with an instance of the
Initialize class.

Alternatively, if you know what attributes to set, you can set the attributes of
the RigidBodyAssessment instance directly without using the `set_fromInitializer` method.
and set the `is_set_up` attribute to True.

The following attributes need to be set for the RigidBodyAssessment instance to work properly:
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

# @time_it
class _Mask:
    """ Class to create various masks for rigid body assessment.

    These masks are used to extract a subset (e.g. per-chain or per-chain pair)
    of information from the confidence metrics like PAE, pLDDT, contact map, etc.
    for the rigid body."""

    rb_dict: Dict[str, List[int]]
    """ Dictionary of rigid bodies with residue indices for each chain."""

    num_to_idx: Dict[str, Dict[int, Dict[str, int]]]
    """ Mapping from residue numbers to indices for each chain."""

    idx_to_num: Dict[int, Dict[str, int]]
    """ Mapping from residue indices to numbers."""

    lengths_dict: Dict[str, int]
    """ Dictionary containing lengths of chains in the structure."""

    contact_map: np.ndarray
    """ Contact map of the structure."""

    plddt_list: np.ndarray
    """ pLDDT scores for each residue."""

    symmetric_pae: bool
    """ Whether to consider PAE matrix as symmetric."""

    idr_chains: List[str]
    """ List of chain IDs considered as IDRs."""

    protein_chain_map: Dict[str, str]
    """ Mapping from chain IDs to protein names."""

    pae: np.ndarray
    """ PAE matrix of the structure."""

    avg_pae: np.ndarray
    """ Symmetrized PAE matrix of the structure."""

    contact_map_mask_2d: np.ndarray
    """ Boolean mask of contact map for inter-chain interactions (2D)."""

    contact_map_mask_1d: np.ndarray
    """ Boolean mask of contact map for inter-chain interactions (1D)."""

    rb_mask_2d: np.ndarray
    """ Boolean mask of residues in the rigid body (2D)."""

    rb_mask_1d: np.ndarray
    """ Boolean mask of residues in the rigid body (1D)."""

    chain_mask_stack_1d: np.ndarray
    """ Stack of boolean masks for each chain (1D)."""

    chain_mask_stack_2d: np.ndarray
    """ Stack of boolean masks for each chain (2D)."""

    chain_pair_mask_stack_1d: np.ndarray
    """ Stack of boolean masks for each chain pair (1D)."""

    chain_pair_mask_stack_2d: np.ndarray
    """ Stack of boolean masks for each chain pair (2D)."""

    def __init__(
        self,
        rb_dict: dict,
        pae: np.ndarray,
        avg_pae: np.ndarray,
        plddt_list: np.ndarray,
        contact_map: np.ndarray,
        lengths_dict: dict,
        num_to_idx: dict,
        idx_to_num: dict,
        symmetric_pae: bool = True,
        idr_chains: List[str] = [],
        protein_chain_map: Dict[str, str] = {},
    ):

        self.pae = pae
        self.avg_pae = avg_pae
        self.plddt_list = np.array(plddt_list)

        self.num_to_idx = num_to_idx
        self.idx_to_num = idx_to_num

        self.lengths_dict = lengths_dict

        self.symmetric_pae = symmetric_pae

        self.idr_chains = idr_chains
        self.protein_chain_map = protein_chain_map

        self.rb_dict = self.transform_rb_dict_to_idxs(rb_dict)

        _interchain_mask = create_mask(
            partition_dict=lengths_dict,
            hide_interactions="intra_part",
            masked_value=0,
            unmasked_value=1,
        )

        self.unique_chains = self.get_unique_chains()
        self.chain_pairs = self.get_chain_pairs()

        self.rb_mask_1d = self.get_rb_mask(lengths_dict, 1) # use for pLDDT
        self.rb_mask_2d = self.get_rb_mask(lengths_dict, 2) # use for PAE

        self.contact_map_mask_2d = np.ma.make_mask(
            contact_map * _interchain_mask * self.rb_mask_2d
        )
        self.contact_map_mask_1d = self.contact_map_mask_2d.any(axis=0)

        self.chain_mask_stack_1d = self.get_chain_mask_stack(1)
        self.chain_mask_stack_2d = self.get_chain_mask_stack(2)

        self.chain_pair_mask_stack_1d = self.get_chain_pair_mask_stack(1)
        self.chain_pair_mask_stack_2d = self.get_chain_pair_mask_stack(2)

    def get_unique_chains(self) -> List[str]:
        """Get unique chains in the rigid body.

        ## Returns:

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

        ## Returns:

        - **chain_pairs (list)**:<br />
            List of tuples containing unique chain pairs.
            Each tuple contains two chain IDs.
        """

        chain_pairs = list(combinations(self.unique_chains, 2))

        return sorted([tuple(pair) for pair in chain_pairs])

    # @time_it
    def transform_rb_dict_to_idxs(
        self,
        rb_dict: Dict[str, List[Tuple[str, int]]],
    ) -> Dict[str, List[int]]:
        """Transform rigid body dictionary from residue numbers to indices.

        ## Arguments:

        - **rb_dict (dict)**:<br />
            Dictionary of rigid bodies with residue numbers.

        ## Returns:

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

    # @time_it
    def get_rb_mask(
        self,
        lengths_dict: Dict[str, int],
        dimensions: int = 2,
    ) -> np.ndarray:
        """Get a binary map of residues in the rigid body.

        ## Arguments:

        - **lengths_dict (dict)**:<br />
            Dictionary containing lengths of chains in the structure.

        - **dimensions (int)**:<br />
            Dimensions of the mask to be generated.

        ## Returns:

        - **rb_mask (np.ndarray)**:<br />
            A binary map of residues in the rigid body.
            The shape is (`total_length`, `total_length`) where `total_length`
            is the sum of lengths of all chains.<br />
            `True` if the residue is part of the rigid body, `False` otherwise.
        """

        _Mask.sanity_check_mask_dimensions(dimensions)
        total_len = lengths_dict.get("total", 0)

        rb_res_idxs = [
            token_idx
            for chain_id in self.rb_dict.keys()
            for token_idx in self.rb_dict[chain_id]
        ]

        if dimensions == 1:
            rb_mask = np.zeros((total_len,), dtype=int)
            rb_mask[rb_res_idxs] = 1

        elif dimensions == 2:
            rb_mask = np.zeros((total_len, total_len), dtype=int)
            rb_mask[np.ix_(rb_res_idxs, rb_res_idxs)] = 1

        rb_mask = np.ma.make_mask(rb_mask)

        return rb_mask

    # @time_it
    def get_chain_mask(
        self,
        chain_id: str,
        lengths_dict: Dict[str, int],
        dimensions: int = 2,
    ) -> np.ndarray:
        """Get a binary map of residues in a chain.

        ## Arguments:

        - **chain_id (str)**:<br />
            Chain ID for which the mask is to be generated.

        - **lengths_dict (dict)**:<br />
            Dictionary containing lengths of chains in the structure.

        ## Returns:

        - **chain_mask (np.ndarray)**:<br />
            A binary map of residues in the chain.
            The shape is (`total_length`, `total_length`) where `total_length`
            is the sum of lengths of all chains.<br />
            `True` if the residue is part of the chain, `False` otherwise.
        """

        _Mask.sanity_check_mask_dimensions(dimensions)
        total_len = lengths_dict.get("total", 0)

        chain_res_idxs = self.rb_dict.get(chain_id, [])

        if dimensions == 1:
            chain_mask = np.zeros((total_len,), dtype=int)
            chain_mask[chain_res_idxs] = 1

        elif dimensions == 2:
            chain_mask = np.zeros((total_len, total_len), dtype=int)
            chain_mask[np.ix_(chain_res_idxs, chain_res_idxs)] = 1

        chain_mask = np.ma.make_mask(chain_mask)

        return chain_mask

    # @time_it
    def get_chain_pair_mask(
        self,
        chain_id_1: str,
        chain_id_2: str,
        lengths_dict: Dict[str, int],
        dimensions: int = 2,
    ) -> np.ndarray:
        """Get a binary map of residues in a chain pair.

        ## Arguments:

        - **chain_id_1 (str)**:<br />
            First chain ID of the chain pair.

        - **chain_id_2 (str)**:<br />
            Second chain ID of the chain pair.

        - **lengths_dict (dict)**:<br />
            Dictionary containing lengths of chains in the structure.

        ## Returns:

        - **chain_pair_mask (np.ndarray)**:<br />
            A binary map of residues in the chain pair.
            The shape is (`total_length`, `total_length`) where `total_length`
            is the sum of lengths of all chains.
            `True` if the residue is part of either chain, `False` otherwise.
        """

        _Mask.sanity_check_mask_dimensions(dimensions)
        total_len = lengths_dict.get("total", 0)

        chain_1_res_idxs = self.rb_dict.get(chain_id_1, [])
        chain_2_res_idxs = self.rb_dict.get(chain_id_2, [])

        if dimensions == 1:
            chain_pair_mask = np.zeros((total_len,), dtype=int)
            chain_pair_mask[chain_1_res_idxs] = 1
            chain_pair_mask[chain_2_res_idxs] = 1

        elif dimensions == 2:
            chain_pair_mask = np.zeros((total_len, total_len), dtype=int)
            chain_pair_mask[np.ix_(chain_1_res_idxs, chain_2_res_idxs)] = 1
            chain_pair_mask[np.ix_(chain_2_res_idxs, chain_1_res_idxs)] = 1

        chain_pair_mask = np.ma.make_mask(chain_pair_mask)

        return chain_pair_mask

    # @time_it
    def get_chain_mask_stack(self, dimensions: int = 2) -> np.ndarray:
        """Get a stack of binary maps for all chains in the rigid body.

        ## Arguments:

        - **dimensions (int)**:<br />
            Dimensions of the mask to be generated.

        ## Returns:

        - **chain_mask_stack (np.ndarray)**:<br />
            A stack of masks for all chains in the rigid body.<br />
            See `get_chain_mask` for more details about the mask.
        """

        _Mask.sanity_check_mask_dimensions(dimensions)

        chain_mask_stack = np.array([
            self.get_chain_mask(
                chain_id,
                self.lengths_dict,
                dimensions=dimensions,
            )
            for chain_id in self.unique_chains
        ])

        return chain_mask_stack

    # @time_it
    def get_chain_pair_mask_stack(self, dimensions: int = 2) -> np.ndarray:
        """Get a stack of binary maps for all chain pairs in the rigid body.

        ## Arguments:

        - **dimensions (int)**:<br />
            Dimensions of the mask to be generated.

        ## Returns:

        - **chain_pair_mask_stack (np.ndarray)**:<br />
            A stack of masks for all chain pairs in the rigid body.<br />
            See `get_chain_pair_mask` for more details about the mask.
        """

        _Mask.sanity_check_mask_dimensions(dimensions)
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

    @staticmethod
    def sanity_check_mask_dimensions(dimensions: int):
        """Sanity check for mask dimensions.

        ## Arguments:

        - **dimensions (int)**:<br />
            Dimensions of the mask to be generated.

        Raises:

        - **ValueError**:<br />
            If dimensions is not 1 or 2.
        """

        if dimensions not in [1, 2]:
            raise ValueError(
                f"Invalid dimensions {dimensions}. "
                f"Dimensions must be 1 or 2."
            )

# @time_it
class RigidBodyChainAssessment:
    """ Class to perform chain-wise assessment of a rigid body."""

    as_average: bool
    """ Whether to return average values or residue-level values as lists."""

    unique_chains: List[str]
    """ List of unique chain IDs in the rigid body."""

    idx_to_num: Dict[int, Dict[str, int]]
    """ Mapping from residue indices to numbers."""

    chain_mask_stack_1d: np.ndarray
    """ Stack of boolean masks for each chain (1D)."""

    rb_mask_1d: np.ndarray
    """ Boolean mask of residues in the rigid body (1D)."""

    contact_map_mask_1d: np.ndarray
    """ Boolean mask of contact map for inter-chain interactions (1D)."""

    plddt_list: np.ndarray
    """ pLDDT scores for each token."""

    idr_chains: List[str]
    """ List of chain IDs considered as IDRs."""

    protein_chain_map: Dict[str, str]
    """ Mapping from chain IDs to protein names."""

    per_chain_plddt: Dict[str, list | float]
    """ pLDDT scores per chain."""

    per_chain_iplddt: Dict[str, list | float]
    """ Interface pLDDT scores per chain."""

    per_chain_interface_res: Dict[str, List[int] | int]
    """ Interface residues per chain."""

    def __init__(
        self,
        _mask: _Mask,
        as_average: bool = True,
    ):

        self.as_average = as_average

        self.unique_chains = _mask.unique_chains
        self.idx_to_num = _mask.idx_to_num

        self.chain_mask_stack_1d = _mask.chain_mask_stack_1d
        self.rb_mask_1d = _mask.rb_mask_1d
        self.contact_map_mask_1d = _mask.contact_map_mask_1d

        self.plddt_list = _mask.plddt_list
        self.idr_chains = _mask.idr_chains
        self.protein_chain_map = _mask.protein_chain_map

        self.per_chain_plddt = self.get_per_chain_plddt(
            only_avg=self.as_average,
            only_interface=False,
        )
        self.per_chain_iplddt = self.get_per_chain_plddt(
            only_avg=self.as_average,
            only_interface=True,
        )
        self.per_chain_interface_res = self.get_per_chain_interface_residues(
            only_count=self.as_average
        )

    # @time_it
    def get_per_chain_plddt(
        self,
        only_avg: bool = True,
        only_interface: bool = False,
    ) -> Dict[str, float | List[float]]:
        """Get average pLDDT score per chain in the rigid body.

        ## Arguments:

        - **only_avg (bool)**:<br />
            If True, returns average pLDDT score per chain.

        - **only_interface (bool)**:<br />
            If True, considers only interface residues for pLDDT calculation.

        ## Returns:

        - **per_chain_plddt (dict)**:<br />
            Dictionary with chain IDs as keys and average pLDDT scores
            as values.
        """

        per_chain_plddt = {}

        mask_attrs_ = {
            True: self.rb_mask_1d * self.contact_map_mask_1d,
            False: self.rb_mask_1d,
        }

        plddt_attrs_ = {
            True: lambda x: np.mean(x),
            False: lambda x: x.tolist(),
        }

        for i, ch_id in enumerate(self.unique_chains):

            chain_mask_1d = self.chain_mask_stack_1d[i, :]
            plddt_mask_1d = chain_mask_1d * mask_attrs_[only_interface]
            per_chain_plddt[ch_id] = plddt_attrs_[only_avg](
                self.plddt_list[plddt_mask_1d]
            )

        return per_chain_plddt

    # @time_it
    def get_per_chain_interface_residues(
        self,
        only_count: bool = True
    ) -> Dict[str, int | List[int]]:
        """Get interface residues per chain in the rigid body.

        ## Arguments:

        - **only_count (bool)**:<br />
            If True, returns count of interface residues per chain.

        ## Returns:

        - **per_chain_interface_res (dict)**:<br />
            Dictionary with chain IDs as keys and list of residue numbers
            at the interface as values.
        """

        per_chain_interface_res = {}

        interface_attrs_ = {
            True: lambda x: np.sum(x),
            False: lambda x: np.where(x)[0],
        }

        mask_multiplier = self.rb_mask_1d * self.contact_map_mask_1d

        for i, ch_id in enumerate(self.unique_chains):

            chain_mask_1d = self.chain_mask_stack_1d[i, :]
            interface_mask_1d = chain_mask_1d * mask_multiplier
            per_chain_interface_res[ch_id] = interface_attrs_[only_count](
                interface_mask_1d
            )

        return per_chain_interface_res

    def get_per_chain_residues(
        self,
        only_count: bool = True
    ) -> Dict[str, int | List[int]]:
        """Get residues per chain in the rigid body.

        ## Arguments:

        - **only_count (bool)**:<br />
            If True, returns count of residues per chain.

        ## Returns:
        - **per_chain_residues (dict)**:<br />
            Dictionary with chain IDs as keys and list of residue numbers
            in the rigid body as values.
        """
        per_chain_residues = {}

        residue_attrs_ = {
            True: lambda x: np.sum(x),
            False: lambda x: np.where(x)[0],
        }

        for i, ch_id in enumerate(self.unique_chains):

            chain_mask_1d = self.chain_mask_stack_1d[i, :]
            residue_mask_1d = chain_mask_1d * self.rb_mask_1d
            per_chain_residues[ch_id] = residue_attrs_[only_count](
                residue_mask_1d
            )

        return per_chain_residues

    def get_chain_attr(
        self,
        chain_id: str,
        attr_name: str
    ) -> float | str | int:
        """ Get the attribute value for a given chain ID.

        > [!NOTE]
        > This function should be used when `as_average` is `True`.

        ## Arguments:

        - **chain_id (str)**:<br />
            Chain ID.

        - **attr_name (str)**:<br />
            Name of the attribute to retrieve.
            Valid attributes are defined in `CHAINWISE_ASSESSMENT_COLUMNS`.

        ## Returns:

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
    ) -> float | str | int:
        """ Get the attribute value for a given residue in a chain.

        > [!NOTE]
        > This function should be used when `as_average` is `False`.

        ## Arguments:

        - **chain_id (str)**:<br />
            Chain ID.

        - **res_idx (int)**:<br />
            Token index of the residue in the chain.

        - **attr_name (str)**:<br />
            Name of the attribute to retrieve.
            Valid attributes are defined in `CHAINWISE_ASSESSMENT_COLUMNS`.

        ## Returns:

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

    def get_chain_assessment(self) -> pd.DataFrame:
        """Get chain-wise assessment for the rigid body.

        ## Returns:

        - **(pd.DataFrame)**:<br />
            DataFrame containing chain-wise assessment.
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

        chain_wise_assessment_rows = [{
                k: func(*it, k)
                for k in CHAINWISE_ASSESSMENT_COLUMNS[self.as_average]
            }
            for it in iterators
        ]

        return pd.DataFrame(chain_wise_assessment_rows)

# @time_it
class RigidBodyChainPairAssessment:
    """ Class to perform chain-pair-wise assessment of a rigid body."""

    as_average: bool
    """ Whether to return average values or residue-level values as lists."""

    unique_chains: List[str]
    """ List of unique chain IDs in the rigid body."""

    chain_pairs: List[Tuple[str, str]]
    """ List of unique chain pairs in the rigid body."""

    idx_to_num: Dict[int, Dict[str, int]]
    """ Mapping from residue indices to numbers."""

    symmetric_pae: bool
    """ Whether to consider PAE matrix as symmetric."""

    chain_mask_stack_1d: np.ndarray
    """ Stack of boolean masks for each chain (1D)."""

    chain_mask_stack_2d: np.ndarray
    """ Stack of boolean masks for each chain (2D)."""

    chain_pair_mask_stack_1d: np.ndarray
    """ Stack of boolean masks for each chain pair (1D)."""

    chain_pair_mask_stack_2d: np.ndarray
    """ Stack of boolean masks for each chain pair (2D)."""

    rb_mask_1d: np.ndarray
    """ Boolean mask of residues in the rigid body (1D)."""

    rb_mask_2d: np.ndarray
    """ Boolean mask of residues in the rigid body (2D)."""

    contact_map_mask_1d: np.ndarray
    """ Boolean mask of contact map for inter-chain interactions (1D)."""

    contact_map_mask_2d: np.ndarray
    """ Boolean mask of contact map for inter-chain interactions (2D)."""

    plddt_list: np.ndarray
    """ pLDDT scores for each token."""

    pae: np.ndarray
    """ PAE matrix of the structure."""

    avg_pae: np.ndarray
    """ Symmetrized PAE matrix of the structure."""

    idr_chains: List[str]
    """ List of chain IDs considered as IDRs."""

    protein_chain_map: Dict[str, str]
    """ Mapping from chain IDs to protein names."""

    chain_pair_interface_res: Dict[Tuple[str, str], List[int] | int]
    """ Interface residues per chain pair."""

    chain_pair_contacts: Dict[Tuple[str, str], List[int] | int]
    """ Number of contacts per chain pair."""

    chain_pair_iplddt: Dict[Tuple[str, str], List[float] | float]
    """ Interface pLDDT per chain pair."""

    chain_pair_pae: Dict[Tuple[str, str], List[float] | float]
    """ PAE per chain pair."""

    chain_pair_ipae: Dict[Tuple[str, str], List[float] | float]
    """ iPAE per chain pair."""

    chain_pair_pae_ij: Dict[Tuple[str, str], List[float] | float]
    """ PAE (i->j) per chain pair."""

    chain_pair_pae_ji: Dict[Tuple[str, str], List[float] | float]
    """ PAE (j->i) per chain pair."""

    chain_pair_ipae_ij: Dict[Tuple[str, str], List[float] | float]
    """ iPAE (i->j) per chain pair."""

    chain_pair_ipae_ji: Dict[Tuple[str, str], List[float] | float]
    """ iPAE (j->i) per chain pair."""

    def __init__(
        self,
        _mask: _Mask,
        as_average: bool = True,
    ):

        self.as_average = as_average

        self.unique_chains = _mask.unique_chains
        self.chain_pairs = _mask.chain_pairs
        self.idx_to_num = _mask.idx_to_num
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
        self.chain_pair_pae = {}
        self.chain_pair_pae_ij = {}
        self.chain_pair_pae_ji = {}
        self.chain_pair_ipae = {}
        self.chain_pair_ipae_ij = {}
        self.chain_pair_ipae_ji = {}

        self.chain_pair_interface_res = self.get_chain_pair_interface(
            per_chain=True,
            only_count=self.as_average,
        )

        self.chain_pair_contacts = self.get_chain_pair_interface(
            per_chain=False,
            only_count=True,
        )

        self.chain_pair_iplddt = self.get_chain_pair_plddt(
            only_avg=self.as_average,
            only_interface=True,
        )

        if self.symmetric_pae:
            self.chain_pair_pae = self.get_chain_pair_pae(
                only_avg=self.as_average,
                only_interface=False,
                symmetric=(True, ""),
            )
        else:
            self.chain_pair_pae_ij = self.get_chain_pair_pae(
                only_avg=self.as_average,
                only_interface=False,
                symmetric=(False, "ij"),
            )
            self.chain_pair_pae_ji = self.get_chain_pair_pae(
                only_avg=self.as_average,
                only_interface=False,
                symmetric=(False, "ji"),
            )

        if self.symmetric_pae:
            self.chain_pair_ipae = self.get_chain_pair_pae(
                only_avg=self.as_average,
                only_interface=True,
                symmetric=(True, ""),
            )
        else:
            self.chain_pair_ipae_ij = self.get_chain_pair_pae(
                only_avg=self.as_average,
                only_interface=True,
                symmetric=(False, "ij"),
            )
            self.chain_pair_ipae_ji = self.get_chain_pair_pae(
                only_avg=self.as_average,
                only_interface=True,
                symmetric=(False, "ji"),
            )

    def get_chain_pair_attr(
        self,
        chain_pair: Tuple[str, str],
        attr_name: str,
    ):
        """ Get the attribute value for a given chain pair.

        ## Arguments:

        - **chain_pair (tuple)**:<br />
            Tuple of two chain IDs.

        - **attr_name (str)**:<br />
            Name of the attribute to retrieve.
            Valid attributes are defined in `CHAIN_PAIRWISE_ASSESSMENT_COLUMNS`.

        ## Returns:

        - **(float | str | int)**:<br />
            The requested attribute value.
        """

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

            "Average PAE": (
                self.chain_pair_pae[chain_pair]
                if self.symmetric_pae else None
            ),
            "Average iPAE": (
                self.chain_pair_ipae[chain_pair]
                if self.symmetric_pae else None
            ),

            "Average PAE ij": (
                self.chain_pair_pae_ij[chain_pair]
                if not self.symmetric_pae else None
            ),
            "Average PAE ji": (
                self.chain_pair_pae_ji[chain_pair]
                if not self.symmetric_pae else None
            ),

            "Average iPAE ij": (
                self.chain_pair_ipae_ij[chain_pair]
                if not self.symmetric_pae else None
            ),
            "Average iPAE ji": (
                self.chain_pair_ipae_ji[chain_pair]
                if not self.symmetric_pae else None
            ),
        }

        return attrs_[attr_name]

    def get_res_pair_attr(
        self,
        chain_pair: Tuple[str, str],
        res_pair: Tuple[int, int],
        attr_name: str,
    ):
        """ Get the attribute value for a given residue pair in a chain pair.

        ## Arguments:

        - **chain_pair (tuple)**:<br />
            Tuple of two chain IDs.

        - **res_pair (tuple)**:<br />
            Tuple of two residue token indices.

        - **attr_name (str)**:<br />
            Name of the attribute to retrieve.
            Valid attributes are defined in `CHAIN_PAIRWISE_ASSESSMENT_COLUMNS`.

        ## Returns:

        - **(float | str | int)**:<br />
            The requested attribute value.
        """

        ch_id_1, ch_id_2 = chain_pair
        res_idx_1, res_idx_2 = res_pair
        attr_state = (self.as_average, self.symmetric_pae)

        if attr_name not in CHAIN_PAIRWISE_ASSESSMENT_COLUMNS[attr_state]:
            raise ValueError(f"Invalid attribute name: {attr_name}")

        attrs_ = {
            "Chain ID 1": ch_id_1,
            "Chain ID 2": ch_id_2,

            "Protein Name 1": self.protein_chain_map.get(ch_id_1, ch_id_1),
            "Protein Name 2": self.protein_chain_map.get(ch_id_2, ch_id_2),

            "Chain Type 1": "IDR" if ch_id_1 in self.idr_chains else "R",
            "Chain Type 2": "IDR" if ch_id_2 in self.idr_chains else "R",

            "Residue 1": self.idx_to_num[res_idx_1]["token_num"],
            "Residue 2": self.idx_to_num[res_idx_2]["token_num"],

            "pLDDT 1": self.plddt_list[res_idx_1],
            "pLDDT 2": self.plddt_list[res_idx_2],

            "PAE ij": self.pae[res_idx_1, res_idx_2],
            "PAE ji": self.pae[res_idx_2, res_idx_1],
            "PAE": self.avg_pae[res_idx_1, res_idx_2],
        }

        return attrs_[attr_name]

    # @time_it
    def get_chain_pair_interface(
        self,
        per_chain: bool = True,
        only_count: bool = True,
    ) -> Dict[Tuple[str, str], List[int] | int | Tuple[int, int]]:
        """Get interface residues for the chain pair in the rigid body.

        ## Arguments:

        - **per_chain (bool)**:<br />
            If True, returns count of interface residues per chain pair.

        - **only_count (bool)**:<br />
            If True, returns number of contacts per chain pair.

        ## Returns:

        - **chain_pair_interface_residues (dict)**:<br />
            Dictionary with chain IDs as keys and list of residue numbers
            at the interface as values.
        """

        chain_pair_interface_residues = {}

        attr_state = (per_chain, only_count)
        interface_mask_multiplier = self.rb_mask_2d * self.contact_map_mask_2d

        attrs_ = {
            (True, True): lambda x: (len(x[:, 0]), len(x[:, 1])),
            (False, True): lambda x: len(x),
            (True, False): lambda x: x.tolist(),
        }

        for i, (ch_id_1, ch_id_2) in enumerate(self.chain_pairs):

            chain_pair_mask = self.chain_pair_mask_stack_2d[i, :, :]
            interface_mask = chain_pair_mask * interface_mask_multiplier
            interface_mask = np.triu(interface_mask) # To avoid duplicates
            res_idx_pairs = np.argwhere(interface_mask)
            chain_pair_interface_residues[(ch_id_1, ch_id_2)] = (
                attrs_[attr_state](res_idx_pairs)
            )

        return chain_pair_interface_residues

    def get_chain_pair_plddt(
        self,
        only_avg: bool = True,
        only_interface: bool = False,
    ):
        """ Get pLDDT scores for the chain pair in the rigid body.

        ## Arguments:

        - **only_avg (bool, optional)**:<br />
            If True, returns average pLDDT scores per chain.

        - **only_interface (bool, optional)**:<br />
            If True, considers only interface residues for pLDDT calculation.

        ## Returns:

        - **chain_pair_plddt (dict)**:<br />
            Dictionary with chain pair tuples as keys and pLDDT scores
            as values.
        """

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

            plddt_mask = np.triu(plddt_mask)  # To avoid duplicate pairs
            res_pair_idxs = np.argwhere(plddt_mask)
            ch1_plddts = self.plddt_list[res_pair_idxs[:, 0]]
            ch2_plddts = self.plddt_list[res_pair_idxs[:, 1]]

            chain_pair_plddt[(ch_id_1, ch_id_2)] = list(
                zip(ch1_plddts, ch2_plddts)
            )

        return chain_pair_plddt

    # @time_it
    def get_chain_pair_pae(
        self,
        only_avg: bool = True,
        only_interface: bool = False,
        symmetric: Tuple[bool, str] = (True, ""),
    ):
        """ Get PAE values for the chain pair in the rigid body.

        ## Arguments:

        - **only_avg (bool, optional)**:<br />
            If True, returns average PAE values per chain pair.

        - **only_interface (bool, optional)**:<br />
            If True, considers only interface residues for PAE calculation.

        - **symmetric (Tuple[bool, str], optional)**:<br />
            If True, treats PAE symmetrically between chain pairs.

        ## Returns:

        - **chain_pair_pae (dict)**:<br />
            Dictionary with chain pair tuples as keys and PAE values
            as values.
        """
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

        ## Returns:

        - **(pd.DataFrame)**:<br />
            DataFrame containing chain-pair-wise assessment.
        """

        chain_pairwise_assessment_rows = []

        attr_state = (self.as_average, self.symmetric_pae)

        if self.as_average:

            iterators_cp = [(chain_pair,) for chain_pair in self.chain_pairs]
            func_cp = self.get_chain_pair_attr

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
                for k in CHAIN_PAIRWISE_ASSESSMENT_COLUMNS[attr_state]
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
        as_average: bool = True,
        symmetric_pae: bool = True,
        **kwargs,
    ):

        self.as_average = as_average
        self.symmetric_pae = symmetric_pae
        self.idr_chains = kwargs.get("idr_chains", [])
        self.protein_chain_map = kwargs.get("protein_chain_map", {})
        self.is_set_up = False
        self.rb_dict = rb_dict

    def check_is_set_up(self):
        """ Check if the RigidBodyAssessment instance is set up. """

        if not self.is_set_up:
            raise ValueError(_error_not_set_up)

    def perform_assessment(self):
        """ Perform the assessment of the rigid body."""

        self.check_is_set_up()

        _mask = _Mask(
            rb_dict=self.rb_dict,
            num_to_idx=self.num_to_idx,
            idx_to_num=self.idx_to_num,
            contact_map=self.contact_map,
            plddt_list=self.token_plddts,
            pae=self.pae,
            avg_pae=self.avg_pae,
            lengths_dict=self.lengths_dict,
            symmetric_pae=self.symmetric_pae,
            idr_chains=self.idr_chains,
            protein_chain_map=self.protein_chain_map,
        )

        self.rb_c_assess = RigidBodyChainAssessment(
            _mask=_mask,
            as_average=self.as_average,
        )
        self.rb_cp_assess = RigidBodyChainPairAssessment(
            _mask=_mask,
            as_average=self.as_average,
        )
        self.overall_assessment = self.get_overall_assessment(
            rb_c_assess=self.rb_c_assess,
            rb_cp_assess=self.rb_cp_assess,
        )

    def set_attributes_from(
        self,
        instance: Initialize | RigidBodies
    ):
        """ Set the attributes of the RigidBodyAssessment instance from or based
        on an instance of Initialize or RigidBodies class.

        > [!CAUTION]
        > If instance is of type RigidBodies it should have already been set up.

        ## Arguments:

        - **instance (Initialize | RigidBodies )**:<br />
            An instance of Initialize or RigidBodies class
        """

        self.num_to_idx = instance.num_to_idx
        self.idx_to_num = instance.idx_to_num

        self.token_coords = instance.token_coords
        self.token_plddts = instance.token_plddts

        self.pae = instance.pae
        self.avg_pae = instance.avg_pae

        self.lengths_dict = instance.lengths_dict

        self.contact_map = get_interaction_map(
            coords1=np.array(self.token_coords),
            coords2=np.array(self.token_coords),
            contact_threshold=IntCons.contact_threshold,
            map_type="contact",
        )

        self.is_set_up = True

    def save_rb_assessment(
        self,
        rb_c_assess: RigidBodyChainAssessment,
        rb_cp_assess: RigidBodyChainPairAssessment,
        overall_assessment: dict,
        save_path: str,
    ):
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

        for col, key in OVERALL_ASSESSMENT_COLUMNS[self.symmetric_pae].items():
            if overall_assessment.get(key, np.nan) is not np.nan:
                overall_assessment_rows.append({
                    "Key": col,
                    "Value": overall_assessment.get(key)
                })
        overall_assessment_df = pd.DataFrame(overall_assessment_rows)

        c_assessment_df = rb_c_assess.get_chain_assessment()
        cp_assessment_df = rb_cp_assess.get_chain_pair_assessment()

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

        with pd.ExcelWriter(save_path, engine='openpyxl', mode='w') as writer:

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
    def get_overall_assessment(
        self,
        rb_c_assess: RigidBodyChainAssessment,
        rb_cp_assess: RigidBodyChainPairAssessment
    ):
        """ Get overall assessment of the rigid body.

        ## Returns:

        - **overall_assessment (dict)**:<br />
            A dictionary containing overall statistics about the rigid body.
            It includes the number of chains, number of interacting chain pairs,
            interface residues, number of contacts, average ipLDDT,
            average IDR ipLDDT, average iPAE ij, and average iPAE ji.
        """

        overall_assessment = {}

        # Average pLDDT across all chains in the rigid body
        overall_assessment["avg_plddt"] = np.mean([
            np.mean(plddt_scores)
            for plddt_scores in rb_c_assess.get_per_chain_plddt(
                only_avg=False,
                only_interface=False,
            ).values()
        ]) if rb_c_assess.per_chain_plddt else np.nan

        # Number of chains in the rigid body
        overall_assessment["num_chains"] = len(rb_c_assess.unique_chains)

        # Number of interacting chain pairs in the rigid body
        overall_assessment["num_interacting_chain_pairs"] = len([
            pair
            for pair in rb_cp_assess.chain_pairs
            if len(rb_cp_assess.chain_pair_interface_res[pair]) > 0
        ])

        # Number of interface residues in the rigid body
        per_c_interface = rb_c_assess.get_per_chain_interface_residues(
            only_count=True
        )
        overall_assessment["num_interface_residues"] = sum(
            per_c_interface.values()
        )

        # Total number of residues in the rigid body
        per_c_total_res = rb_c_assess.get_per_chain_residues(
            only_count=True
        )
        overall_assessment["num_total_residues"] = sum(
            per_c_total_res.values()
        )
        overall_assessment["rb_coverage"] = (
            overall_assessment["num_total_residues"] /
            self.lengths_dict["total"]
        )

        # Number of contacts formed in the rigid body
        overall_assessment["num_contacts"] = sum(
            rb_cp_assess.chain_pair_contacts.values()
        )

        # Average ipLDDT scores across all chains in the rigid body
        global_iplddt_scores = np.array([
            np.array([iplddt, chain_id])
            for chain_id, iplddt_scores in rb_c_assess.get_per_chain_plddt(
                only_avg=False,
                only_interface=True,
            ).items()
            for iplddt in iplddt_scores
        ])
        overall_assessment["avg_iplddt"] = (
            np.mean(global_iplddt_scores[:, 0].astype(float))
            if global_iplddt_scores.size > 0 else np.nan
        )

        # Average ipLDDT scores for IDR chains in the rigid body
        if (
            len(rb_c_assess.idr_chains) == 0 or
            global_iplddt_scores.size == 0
        ):
            overall_assessment["avg_idr_iplddt"] = np.nan
        else:
            _idr_chain_mask = np.isin(
                global_iplddt_scores[:, 1],
                rb_c_assess.idr_chains
            )
            global_idr_iplddt_scores = global_iplddt_scores[
                _idr_chain_mask, 0
            ].astype(float)
            overall_assessment["avg_idr_iplddt"] = (
                np.mean(global_idr_iplddt_scores)
                if global_idr_iplddt_scores.size > 0 else np.nan
            )

        # Average iPAE scores across all chain pairs in the rigid body
        attrs_pae = {
            (True, False, "ij"): rb_cp_assess.get_chain_pair_pae,
            (False, False, "ij"): rb_cp_assess.chain_pair_ipae_ij,
            (True, False, "ji"): rb_cp_assess.get_chain_pair_pae,
            (False, False, "ji"): rb_cp_assess.chain_pair_ipae_ji,
            (True, True, ""): rb_cp_assess.get_chain_pair_pae,
            (False, True, ""): rb_cp_assess.chain_pair_ipae,
        }

        col_template = "avg_ipae"
        attr_states_ = {
            True: [
                {
                    "key": (self.as_average, self.symmetric_pae, ""),
                    "col": "",
                },
            ],
            False: [
                {
                    "key": (self.as_average, self.symmetric_pae, "ij"),
                    "col": "_ij",
                },
                {
                    "key": (self.as_average, self.symmetric_pae, "ji"),
                    "col": "_ji",
                }
            ],
        }

        attr_states = attr_states_[self.symmetric_pae]

        for attr_dict in attr_states:

            attr_state = attr_dict["key"]
            col_name = f"{col_template}{attr_dict['col']}"
            pae_direction = attr_dict['col'].lstrip("_")

            if callable(attrs_pae[attr_state]):

                global_pae_scores = [
                    pae
                    for _cp, pae_lst in attrs_pae[attr_state](
                        only_avg=False,
                        only_interface=True,
                        symmetric=(self.symmetric_pae, pae_direction),
                    ).items()
                    for pae in pae_lst
                ]

            else:
                global_pae_scores = [
                    pae
                    for _cp, pae_lst in attrs_pae[attr_state].items()
                    for pae in pae_lst
                ]

            overall_assessment[col_name] = (
                np.mean(global_pae_scores) if global_pae_scores else np.nan
            )

        return overall_assessment