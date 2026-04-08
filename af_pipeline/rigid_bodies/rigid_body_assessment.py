"""
[rigid_body_assessment](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/rigid_bodies/rigid_body_assessment.py)
==============================

RigidBodyAssessment class with methods to assess confidently predicted regions (rigid bodies) extracted from AlphaFold predictions.
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from itertools import combinations
from typing import Dict, List, Tuple
from af_pipeline.utils.misc_utils import create_mask
from af_pipeline.tools.structure_tools import get_interaction_map
from af_pipeline.parser.initialize import Initialize
from af_pipeline.constants.af_constants import InteractionConstants as IntCons, OverallAssessment
from af_pipeline.constants.af_constants import (
    CHAIN_PAIRWISE_ASSESSMENT_COLUMNS,
    CHAINWISE_ASSESSMENT_COLUMNS,
    OVERALL_ASSESSMENT_COLUMNS,
    ChainAssessment,
    ChainPairAssessment,
    ChainType,
    KeywordArg,
    ResidueMapKeys,
    MaskedInteractionType,
    InteractionMapType,
    MiscStrEnum
)

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
            hide_interactions=MaskedInteractionType.INTRA_PART,
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
        total_len = lengths_dict.get(MiscStrEnum.TOTAL, 0)

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
        total_len = lengths_dict.get(MiscStrEnum.TOTAL, 0)

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
        total_len = lengths_dict.get(MiscStrEnum.TOTAL, 0)

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
    """ Whether to return an average value of a metric for the chain or residue-level values as lists."""

    show_interface_residues_only: bool
    """ Whether to show the metrics for the interface residues only in the output
    chain assessment and chain-pair assessment.
    > [!NOTE]
    > This option only takes effect when there are residue-level or
    > residue-pair-level metrics available, i.e. when `as_average` is `False`.
    """

    unique_chains: List[str]
    """ List of unique chain IDs in the rigid body."""

    idx_to_num: Dict[int, Dict[str, int]]
    """ Mapping from residue indices to numbers."""

    _chain_mask_stack_1d: np.ndarray
    """ Stack of boolean masks for each chain (1D)."""

    _rb_mask_1d: np.ndarray
    """ Boolean mask of residues in the rigid body (1D)."""

    _contact_map_mask_1d: np.ndarray
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

    per_chain_residues: Dict[str, List[int] | int]
    """ Residues per chain."""

    per_chain_interface_residues: Dict[str, List[int] | int]
    """ Interface residues per chain."""

    total_residues: int
    """ Total number of residues in the rigid body."""

    total_interface_residues: int
    """ Total number of interface residues in the rigid body."""

    total_idr_residues: int
    """ Total number of residues in the rigid body that within IDRs."""

    total_interface_idr_residues: int
    """ Total number of interface residues in the rigid body that within IDRs."""

    def __init__(
        self,
        _mask: _Mask,
        as_average: bool = True,
        show_interface_residues_only: bool = False,
    ):

        self.as_average = as_average
        self.show_interface_residues_only = show_interface_residues_only

        self.unique_chains = _mask.unique_chains
        self.idx_to_num = _mask.idx_to_num

        self._chain_mask_stack_1d = _mask.chain_mask_stack_1d
        self._rb_mask_1d = _mask.rb_mask_1d
        self._contact_map_mask_1d = _mask.contact_map_mask_1d

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
        self.per_chain_residues = self.get_per_chain_residues(
            only_count=self.as_average,
            only_interface=False,
        )
        self.per_chain_interface_residues = self.get_per_chain_residues(
            only_count=self.as_average,
            only_interface=True,
        )
        self.total_residues = self.get_total_residue_count(
            only_interface=False,
            only_idr=False,
        )
        self.total_interface_residues = self.get_total_residue_count(
            only_interface=True,
            only_idr=False,
        )
        self.total_idr_residues = self.get_total_residue_count(
            only_interface=False,
            only_idr=True,
        )
        self.total_interface_idr_residues = self.get_total_residue_count(
            only_interface=True,
            only_idr=True,
        )
    # @time_it
    def get_per_chain_plddt(
        self,
        only_avg: bool = True,
        only_interface: bool = False,
    ) -> Dict[str, float | List[float]]:
        """Get pLDDT or ipLDDT score per chain in the rigid body.

        There are four possible outcomes depending on the combination of `only_avg`
        and `only_interface` arguments:
        - **only_avg=True** and **only_interface=False**:<br />
            Average pLDDT score per chain for all residues in the chain. (Default)
        - **only_avg=False** and **only_interface=False**:<br />
            List of per-residue pLDDT scores in each chain.
        - **only_avg=True** and **only_interface=True**:<br />
            Average pLDDT score per chain for interface residues only.
        - **only_avg=False** and **only_interface=True**:<br />
            List of per-residue pLDDT scores (only interface residues) in the chain.

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

        # only_interface
        mask_attrs_ = {
            True: self._rb_mask_1d * self._contact_map_mask_1d,
            False: self._rb_mask_1d,
        }

        # only_avg
        plddt_attrs_ = {
            True: np.mean,
            False: lambda x: x.tolist(),
        }

        for i, ch_id in enumerate(self.unique_chains):

            chain_mask_1d = self._chain_mask_stack_1d[i, :]
            plddt_mask_1d = chain_mask_1d * mask_attrs_[only_interface]
            per_chain_plddt[ch_id] = plddt_attrs_[only_avg](
                self.plddt_list[plddt_mask_1d]
            )

        return per_chain_plddt

    # @time_it
    def get_per_chain_residues(
        self,
        only_count: bool = True,
        only_interface: bool = False,
    ) -> Dict[str, int | List[int]]:
        """Get residue indices or count per chain in the rigid body.

        There are four possible outcomes depending on the combination of `only_count`
        and `only_interface` arguments:
        - **only_count=True** and **only_interface=False**:<br />
            Count of residues per chain in the rigid body. (Default)
        - **only_count=False** and **only_interface=False**:<br />
            List of residue indices in each chain within the rigid body.
        - **only_count=True** and **only_interface=True**:<br />
            Count of interface residues per chain in the rigid body.
        - **only_count=False** and **only_interface=True**:<br />
            List of interface residue indices in each chain within the rigid body.

        ## Arguments:

        - **only_count (bool)**:<br />
            If True, returns count of residues per chain.

        - **only_interface (bool)**:<br />
            If True, considers only interface residues for counting.

        ## Returns:
        - **per_chain_residues (dict)**:<br />
            Dictionary with chain IDs as keys and list of residue indices
            in the rigid body as values.
        """

        per_chain_residues = {}

        # only_interface
        mask_attrs_ = {
            True: self._rb_mask_1d * self._contact_map_mask_1d,
            False: self._rb_mask_1d,
        }

        # only_count
        residue_attrs_ = {
            True: np.sum,
            False: lambda x: np.where(x)[0],
        }

        for i, ch_id in enumerate(self.unique_chains):

            chain_mask_1d = self._chain_mask_stack_1d[i, :]
            residue_mask_1d = chain_mask_1d * mask_attrs_[only_interface]
            per_chain_residues[ch_id] = residue_attrs_[only_count](
                residue_mask_1d
            )

        return per_chain_residues

    def get_total_residue_count(
        self,
        only_interface: bool = False,
        only_idr: bool = False
    ) -> int:
        """ Get total residue count in the rigid body.

        There are four possible outcomes depending on the combination of `only_interface`
        and `only_idr` arguments:
        - **only_interface=False** and **only_idr=False**:<br />
            Total count of residues in the rigid body. (Default)
        - **only_interface=True** and **only_idr=False**:<br />
            Total count of interface residues in the rigid body.
        - **only_interface=False** and **only_idr=True**:<br />
            Total count of residues in the rigid body that are within IDRs.
        - **only_interface=True** and **only_idr=True**:<br />
            Total count of interface residues in the rigid body that are within IDRs.

        ## Arguments:

        - **only_interface (bool, optional):**:<br />
            If True, counts only interface residues in the rigid body. Default is False.

        - **only_idr (bool, optional):**:<br />
            If True, counts only residues within IDRs in the rigid body. Default is False.

        ## Returns:

        - **int**:<br />
            Total residue count in the rigid body based on the specified criteria.
        """

        # only_interface, only_idr
        quantity_attrs_ = {
            (True, False): self.per_chain_interface_residues,
            (False, False): self.per_chain_residues,
            (True, True): {
                ch_id: res
                for ch_id, res in self.per_chain_interface_residues.items()
                if ch_id in self.idr_chains
            },
            (False, True): {
                ch_id: res for ch_id, res in self.per_chain_residues.items()
                if ch_id in self.idr_chains
            },
        }

        attr_state = (only_interface, only_idr)
        if self.as_average:
            total_residues = sum(quantity_attrs_[attr_state].values())
        else:
            total_residues = sum([
                len(res_idxs) for res_idxs in quantity_attrs_[attr_state].values()
            ])

        return total_residues

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
            Name of the attribute to retrieve. Valid attributes are defined in
            `af_pipeline.constants.af_constants.CHAINWISE_ASSESSMENT_COLUMNS`.

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
            ChainAssessment.CHAIN_ID: chain_id,
            ChainAssessment.PROTEIN_NAME: self.protein_chain_map.get(chain_id, f"chain_{chain_id}"),
            ChainAssessment.AVERAGE_PLDDT: self.per_chain_plddt[chain_id],
            ChainAssessment.AVERAGE_IPLDDT: self.per_chain_iplddt[chain_id],
            ChainAssessment.INTERFACE_RESIDUES: self.per_chain_interface_residues[chain_id],
            ChainAssessment.CHAIN_TYPE: ChainType.IDR if chain_id in self.idr_chains else ChainType.R,
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
            Name of the attribute to retrieve. Valid attributes are defined in
            `af_pipeline.constants.af_constants.CHAINWISE_ASSESSMENT_COLUMNS`.

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

        local_idx1 = list(self.per_chain_residues[chain_id]).index(res_idx)
        local_idx2 = (
            list(self.per_chain_interface_residues[chain_id]).index(res_idx)
            if res_idx in self.per_chain_interface_residues[chain_id] else None
        )

        attrs_ = {
            ChainAssessment.CHAIN_ID: chain_id,
            ChainAssessment.RESIDUE_NUMBER: res_num,
            ChainAssessment.PROTEIN_NAME: self.protein_chain_map.get(chain_id, f"chain_{chain_id}"),
            ChainAssessment.AVERAGE_PLDDT: self.per_chain_plddt[chain_id][local_idx1],
            ChainAssessment.IS_INTERFACE_RESIDUE: True if local_idx2 is not None else False,
            ChainAssessment.CHAIN_TYPE: ChainType.IDR if chain_id in self.idr_chains else ChainType.R,
        }

        return attrs_[attr_name]

    def get_chain_assessment(self) -> pd.DataFrame:
        """Get chain-wise assessment for the rigid body.

        Chain-wise assessment includes the following:
        - Chain ID
        - Protein name (if available, otherwise defaults to "chain_{chain_id}")
        - pLDDT score
        - Interface pLDDT score
        - Interface residues
        - Chain type (IDR or R)

        ## Returns:

        - **(pd.DataFrame)**:<br />
            DataFrame containing chain-wise assessment.
        """

        chain_wise_assessment_rows = []

        if self.as_average:

            iterators = [(chain_id,) for chain_id in self.unique_chains]
            func = self.get_chain_attr

        else:

            if self.show_interface_residues_only:
                _per_chain_residues = self.per_chain_interface_residues
            else:
                _per_chain_residues = self.per_chain_residues

            iterators = [
                (ch_id, res_idx)
                for ch_id in self.unique_chains
                for res_idx in _per_chain_residues[ch_id]
            ]
            func = self.get_res_attr

        valid_cols = CHAINWISE_ASSESSMENT_COLUMNS[self.as_average]

        if (
            self.show_interface_residues_only and
            ChainAssessment.IS_INTERFACE_RESIDUE in valid_cols
        ):
            valid_cols = [
                col for col in valid_cols
                if col != ChainAssessment.IS_INTERFACE_RESIDUE
            ]

        chain_wise_assessment_rows = [
            {k: func(*it, k) for k in valid_cols}
            for it in iterators
        ]

        return pd.DataFrame(chain_wise_assessment_rows)

# @time_it
class RigidBodyChainPairAssessment:
    """ Class to perform chain-pair-wise assessment of a rigid body."""

    as_average: bool
    """ Whether to return average values or residue-level values as lists."""

    show_interface_residues_only: bool
    """ Whether to show the metrics for the interface residues only in the output
    chain assessment and chain-pair assessment.
    > [!NOTE]
    > This option only takes effect when there are residue-level or
    > residue-pair-level metrics available, i.e. when `as_average` is `False`.
    """

    unique_chains: List[str]
    """ List of unique chain IDs in the rigid body."""

    chain_pairs: List[Tuple[str, str]]
    """ List of unique chain pairs in the rigid body."""

    idx_to_num: Dict[int, Dict[str, int]]
    """ Mapping from residue indices to numbers."""

    symmetric_pae: bool
    """ Whether to consider PAE matrix as symmetric."""

    _chain_mask_stack_1d: np.ndarray
    """ Stack of boolean masks for each chain (1D)."""

    # _chain_mask_stack_2d: np.ndarray
    # """ Stack of boolean masks for each chain (2D)."""

    # _chain_pair_mask_stack_1d: np.ndarray
    # """ Stack of boolean masks for each chain pair (1D)."""

    _chain_pair_mask_stack_2d: np.ndarray
    """ Stack of boolean masks for each chain pair (2D)."""

    # _rb_mask_1d: np.ndarray
    # """ Boolean mask of residues in the rigid body (1D)."""

    _rb_mask_2d: np.ndarray
    """ Boolean mask of residues in the rigid body (2D)."""

    # _contact_map_mask_1d: np.ndarray
    # """ Boolean mask of contact map for inter-chain interactions (1D)."""

    _contact_map_mask_2d: np.ndarray
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

    chain_pair_interface_residues: Dict[Tuple[str, str], List[int] | int]
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

    chain_pair_residues: Dict[Tuple[str, str], List[int] | int] = {}
    """ Residues per chain pair. Only calculated when `show_interface_residues_only` is `False`."""

    chain_pair_residue_counts: Dict[Tuple[str, str], int] = {}
    """ Residue counts per chain pair. Only calculated when `show_interface_residues_only` is `False`."""

    def __init__(
        self,
        _mask: _Mask,
        as_average: bool = True,
        show_interface_residues_only: bool = True,
    ):

        self.as_average = as_average
        self.show_interface_residues_only = show_interface_residues_only

        self.unique_chains = _mask.unique_chains
        self.chain_pairs = _mask.chain_pairs
        self.idx_to_num = _mask.idx_to_num
        self.symmetric_pae = _mask.symmetric_pae

        self._chain_mask_stack_1d = _mask.chain_mask_stack_1d
        # self._chain_mask_stack_2d = _mask.chain_mask_stack_2d

        # self._chain_pair_mask_stack_1d = _mask.chain_pair_mask_stack_1d
        self._chain_pair_mask_stack_2d = _mask.chain_pair_mask_stack_2d

        # self._rb_mask_1d = _mask.rb_mask_1d
        self._rb_mask_2d = _mask.rb_mask_2d

        # self._contact_map_mask_1d = _mask.contact_map_mask_1d
        self._contact_map_mask_2d = _mask.contact_map_mask_2d

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

        self.chain_pair_residues = {}
        self.chain_pair_residue_counts = {}

        if self.show_interface_residues_only is False:

            self.chain_pair_residues = self.get_chain_pair_residues(
                per_chain=True,
                only_count=self.as_average,
                only_interface=False,
            )

            self.chain_pair_residue_counts = self.get_chain_pair_residues(
                per_chain=False,
                only_count=True,
                only_interface=False,
            )

        self.chain_pair_interface_residues = self.get_chain_pair_residues(
            per_chain=True,
            only_count=self.as_average,
            only_interface=True,
        )

        self.chain_pair_contacts = self.get_chain_pair_residues(
            per_chain=False,
            only_count=True,
            only_interface=True,
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
            self.chain_pair_ipae = self.get_chain_pair_pae(
                only_avg=self.as_average,
                only_interface=True,
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

        > [!NOTE]
        > This function should be used when `as_average` is `True`.

        ## Arguments:

        - **chain_pair (tuple)**:<br />
            Tuple of two chain IDs.

        - **attr_name (str)**:<br />
            Name of the attribute to retrieve. Valid attributes are defined in
            `af_pipeline.constants.af_constants.CHAIN_PAIRWISE_ASSESSMENT_COLUMNS`.

        ## Returns:

        - **(float | str | int)**:<br />
            The requested attribute value.
        """

        chain1, chain2 = chain_pair
        attr_state = (self.as_average, self.symmetric_pae)

        if attr_name not in CHAIN_PAIRWISE_ASSESSMENT_COLUMNS[attr_state]:
            raise ValueError(f"Invalid attribute name: {attr_name}")

        attrs_ = {
            ChainPairAssessment.CHAIN_ID: (chain1, chain2),

            ChainPairAssessment.PROTEIN_NAME: (
                self.protein_chain_map.get(chain1, f"chain_{chain1}"),
                self.protein_chain_map.get(chain2, f"chain_{chain2}"),
            ),

            ChainPairAssessment.CHAIN_TYPE: (
                ChainType.IDR if chain1 in self.idr_chains else ChainType.R,
                ChainType.IDR if chain2 in self.idr_chains else ChainType.R,
            ),

            ChainPairAssessment.INTERFACE_RESIDUES: self.chain_pair_interface_residues[chain_pair],
            ChainPairAssessment.NUMBER_OF_CONTACTS: self.chain_pair_contacts[chain_pair],

            ChainPairAssessment.AVERAGE_IPLDDT_CHAIN1: self.chain_pair_iplddt[chain_pair][0],
            ChainPairAssessment.AVERAGE_IPLDDT_CHAIN2: self.chain_pair_iplddt[chain_pair][1],

            ChainPairAssessment.AVERAGE_PAE: (
                self.chain_pair_pae[chain_pair]
                if self.symmetric_pae else None
            ),
            ChainPairAssessment.AVERAGE_IPAE: (
                self.chain_pair_ipae[chain_pair]
                if self.symmetric_pae else None
            ),

            ChainPairAssessment.AVERAGE_PAE_IJ: (
                self.chain_pair_pae_ij[chain_pair]
                if not self.symmetric_pae else None
            ),
            ChainPairAssessment.AVERAGE_PAE_JI: (
                self.chain_pair_pae_ji[chain_pair]
                if not self.symmetric_pae else None
            ),

            ChainPairAssessment.AVERAGE_IPAE_IJ: (
                self.chain_pair_ipae_ij[chain_pair]
                if not self.symmetric_pae else None
            ),
            ChainPairAssessment.AVERAGE_IPAE_JI: (
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

        > [!NOTE]
        > This function should be used when `as_average` is `False`.

        ## Arguments:

        - **chain_pair (tuple)**:<br />
            Tuple of two chain IDs.

        - **res_pair (tuple)**:<br />
            Tuple of two residue token indices.

        - **attr_name (str)**:<br />
            Name of the attribute to retrieve. Valid attributes are defined in
            `af_pipeline.constants.af_constants.CHAIN_PAIRWISE_ASSESSMENT_COLUMNS`.

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
            ChainPairAssessment.CHAIN_ID_1: ch_id_1,
            ChainPairAssessment.CHAIN_ID_2: ch_id_2,

            ChainPairAssessment.PROTEIN_NAME_1: self.protein_chain_map.get(ch_id_1, f"chain_{ch_id_1}"),
            ChainPairAssessment.PROTEIN_NAME_2: self.protein_chain_map.get(ch_id_2, f"chain_{ch_id_2}"),

            ChainPairAssessment.CHAIN_TYPE_1: ChainType.IDR if ch_id_1 in self.idr_chains else ChainType.R,
            ChainPairAssessment.CHAIN_TYPE_2: ChainType.IDR if ch_id_2 in self.idr_chains else ChainType.R,

            ChainPairAssessment.RESIDUE_1: self.idx_to_num[res_idx_1][ResidueMapKeys.TOKEN_NUM],
            ChainPairAssessment.RESIDUE_2: self.idx_to_num[res_idx_2][ResidueMapKeys.TOKEN_NUM],

            ChainPairAssessment.PLDDT_1: self.plddt_list[res_idx_1],
            ChainPairAssessment.PLDDT_2: self.plddt_list[res_idx_2],

            ChainPairAssessment.PAE_IJ: self.pae[res_idx_1, res_idx_2],
            ChainPairAssessment.PAE_JI: self.pae[res_idx_2, res_idx_1],
            ChainPairAssessment.PAE: self.avg_pae[res_idx_1, res_idx_2],
        }

        return attrs_[attr_name]

    def get_chain_pair_residues(
        self,
        per_chain: bool = True,
        only_interface: bool = False,
        only_count: bool = True,
    ) -> Dict[Tuple[str, str], List[int] | int | Tuple[int, int]]:
        """ Get residues in each chain pair within a rigid body.

        There are several possible outcomes depending on the combination of `per_chain`,
        `only_interface`, and `only_count` arguments:
        - **per_chain=True**, **only_interface=False**, and **only_count=True**:<br />
            Count of residues per chain in each chain pair in the rigid body. (Default)
        - **per_chain=True**, **only_interface=False**, and **only_count=False**:<br />
            List of residue indices per chain in each chain pair in the rigid body.
        - **per_chain=True**, **only_interface=True**, and **only_count=True**:<br />
            Count of interface residues per chain in each chain pair in the rigid body.
        - **per_chain=True**, **only_interface=True**, and **only_count=False**:<br />
            List of interface residue indices per chain in each chain pair in the rigid body.
        - **per_chain=False**, **only_interface=False**, and **only_count=True**:<br />
            Count of residue pairs in each chain pair in the rigid body.
        - **per_chain=False**, **only_interface=False**, and **only_count=False**:<br />
            List of residue index pairs as tuples in each chain pair in the rigid body.
        - **per_chain=False**, **only_interface=True**, and **only_count=True**:<br />
            Count of interface residue pairs in each chain pair in the rigid body.
        - **per_chain=False**, **only_interface=True**, and **only_count=False**:<br />
            List of interface residue index pairs as tuples in each chain pair in the rigid body.

        ## Arguments:

        - **per_chain (bool, optional):**:<br />
            If True, returns residues for each chain in the pair separately.
            If False, returns residue pairs as tuples of indices.

        - **only_interface (bool, optional):**:<br />
            If True, considers only interface residues.

        - **only_count (bool, optional):**:<br />
            If True, returns the count of residues instead of their indices.

        ## Returns:

        - **Dict[Tuple[str, str], List[int] | int | Tuple[int, int]]**:<br />
            Dictionary with chain pair tuples as keys and residue indices (or counts)
            as values.
        """

        chain_pair_residues = {}

        # only_interface
        mask_attrs_ = {
            True: self._rb_mask_2d * self._contact_map_mask_2d,
            False: self._rb_mask_2d,
        }

        # per_chain, only_count
        res_pair_attrs_ = {
            (True, True): lambda x: (len(x[:, 0]), len(x[:, 1])),
            (False, True): len,
            (True, False): lambda x: x.tolist(),
            (False, False): lambda x: x.tolist(),
        }

        attr_state = (per_chain, only_count)
        for i, (ch_id_1, ch_id_2) in enumerate(self.chain_pairs):

            chain_pair_mask = self._chain_pair_mask_stack_2d[i, :, :]
            residue_mask = chain_pair_mask * mask_attrs_[only_interface]
            residue_mask = np.triu(residue_mask) # To avoid duplicates
            res_idx_pairs = np.argwhere(residue_mask)
            chain_pair_residues[(ch_id_1, ch_id_2)] = (
                res_pair_attrs_[attr_state](res_idx_pairs)
            )

        return chain_pair_residues

    def get_chain_pair_plddt(
        self,
        only_avg: bool = True,
        only_interface: bool = False,
    ):
        """ Get pLDDT or ipLDDT scores for the chain pair in the rigid body.

        There are four possible outcomes depending on the combination of `only_avg`
        and `only_interface` arguments:
        - **only_avg=True** and **only_interface=False**:<br />
            Average pLDDT score per chain for all residues in the chain pair. (Default)
        - **only_avg=False** and **only_interface=False**:<br />
            List of per-residue pLDDT scores for each residue pair in the chain pair
            in the rigid body.
        - **only_avg=True** and **only_interface=True**:<br />
            Average pLDDT score per chain for interface residues only in the chain pair
            in the rigid body.
        - **only_avg=False** and **only_interface=True**:<br />
            List of per-residue pLDDT scores for each interface residue pair in the chain
            pair in the rigid body.

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
            chain_pair_mask = self._chain_pair_mask_stack_2d[i, :, :]
            plddt_mask = chain_pair_mask * self._rb_mask_2d

            if only_interface:
                plddt_mask = plddt_mask * self._contact_map_mask_2d

            plddt_mask_1d = plddt_mask.any(axis=0)
            ch1_idx = self.unique_chains.index(ch_id_1)
            chain1_mask = plddt_mask_1d * self._chain_mask_stack_1d[ch1_idx, :]
            ch2_idx = self.unique_chains.index(ch_id_2)
            chain2_mask = plddt_mask_1d * self._chain_mask_stack_1d[ch2_idx, :]

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

        There are several possible outcomes depending on the combination of `only_avg`,
        `only_interface`, and `symmetric` arguments:
        - **only_avg=True**, **only_interface=False**, and **symmetric=(True, '')**:<br />
            Average PAE value per chain pair in the rigid body, treating PAE symmetrically. (Default)
        - **only_avg=False**, **only_interface=False**, and **symmetric=(True, '')**:<br />
            List of PAE values for each residue pair in the chain pair in the rigid body, treating PAE symmetrically.
        - **only_avg=True**, **only_interface=False**, and **symmetric=(False, 'ij')**:<br />
            Average PAE value per chain pair in the rigid body, treating PAE asymmetrically (i->j).
        - **only_avg=False**, **only_interface=False**, and **symmetric=(False, 'ij')**:<br />
            List of PAE values for each residue pair in the chain pair in the rigid body, treating PAE asymmetrically (i->j).
        - **only_avg=True**, **only_interface=False**, and **symmetric=(False, 'ji')**:<br />
            Average PAE value per chain pair in the rigid body, treating PAE asymmetrically (j->i).
        - **only_avg=False**, **only_interface=False**, and **symmetric=(False, 'ji')**:<br />
            List of PAE values for each residue pair in the chain pair in the rigid body, treating PAE asymmetrically (j->i).
        - **only_avg=True**, **only_interface=True**, and **symmetric=(True, '')**:<br />
            Average PAE value per chain pair in the rigid body for interface residues only, treating PAE symmetrically.
        - **only_avg=False**, **only_interface=True**, and **symmetric=(True, '')**:<br />
            List of PAE values for each interface residue pair in the chain pair in the rigid body, treating PAE symmetrically.
        - **only_avg=True**, **only_interface=True**, and **symmetric=(False, 'ij')**:<br />
            Average PAE value per chain pair in the rigid body for interface residues only, treating PAE asymmetrically (i->j).
        - **only_avg=False**, **only_interface=True**, and **symmetric=(False, 'ij')**:<br />
            List of PAE values for each interface residue pair in the chain pair in the rigid body, treating PAE asymmetrically (i->j).
        - **only_avg=True**, **only_interface=True**, and **symmetric=(False, 'ji')**:<br />
            Average PAE value per chain pair in the rigid body for interface residues only, treating PAE asymmetrically (j->i).
        - **only_avg=False**, **only_interface=True**, and **symmetric=(False, 'ji')**:<br />
            List of PAE values for each interface residue pair in the chain pair in the rigid body, treating PAE asymmetrically (j->i).

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
            True: self._rb_mask_2d * self._contact_map_mask_2d,
            False: self._rb_mask_2d,
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

            chain_pair_mask = self._chain_pair_mask_stack_2d[i, :, :]
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

        Chain-pair-wise assessment includes the following:
        - Chain IDs
        - Protein names (if available, otherwise defaults to "chain_{chain_id}")
        - Chain types (IDR or R)
        - Number of interface residues
        - Average pLDDT for interface residues
        - Average PAE (if `symmetric_pae` is `True`)
        - Average iPAE (if `symmetric_pae` is `True`)
        - Average PAE (i->j) (if `symmetric_pae` is `False`)
        - Average PAE (j->i) (if `symmetric_pae` is `False`)
        - Average iPAE (i->j) (if `symmetric_pae` is `False`)
        - Average iPAE (j->i) (if `symmetric_pae` is `False`)

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

            if self.show_interface_residues_only is True:
                _chain_pair_residues = self.chain_pair_interface_residues
            else:
                _chain_pair_residues = self.chain_pair_residues

            iterators_cp = [
                (chain_pair, res_pair)
                for chain_pair in self.chain_pairs
                for res_pair in _chain_pair_residues[chain_pair]
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

    as_average: bool
    """ Whether to report only the average of assessment metrics."""

    symmetric_pae: bool
    """ Whether to treat PAE symmetrically between chain pairs."""

    show_interface_residues_only: bool
    """ Whether to show the metrics for the interface residues only in the output
    chain assessment and chain-pair assessment.
    > [!NOTE]
    > This option only takes effect when there are residue-level or
    > residue-pair-level metrics available, i.e. when `as_average` is `False`.
    """

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

    idr_chains: List[str]
    """ List of chain IDs that are considered disordered (IDRs)."""

    protein_chain_map: Dict[str, str]
    """ Mapping of chain IDs to protein names."""

    chain_pairs: List[Tuple[str, str]]
    """ List of unique chain pairs in the rigid body."""

    overall_assessment: dict
    """ Dictionary of overall assessment metrics for the rigid body."""

    save_path: str
    """ Path to save the assessment results."""

    def __init__(
        self,
        rb_dict: dict,
        as_average: bool = True,
        symmetric_pae: bool = True,
        show_interface_residues_only: bool = True,
        **kwargs,
    ):

        self.as_average = as_average
        self.symmetric_pae = symmetric_pae
        self.show_interface_residues_only = show_interface_residues_only
        self.idr_chains = kwargs.get(KeywordArg.IDR_CHAINS, [])
        self.protein_chain_map = kwargs.get(KeywordArg.PROTEIN_CHAIN_MAP, {})
        self._is_set_up = False
        setup_instance = kwargs.get(KeywordArg.SETUP_INSTANCE, None)
        if isinstance(setup_instance, Initialize):
            self.set_attributes_from(instance=setup_instance)
        else:
            from af_pipeline.rigid_bodies.rigid_bodies import RigidBodies
            if isinstance(setup_instance, RigidBodies):
                self.set_attributes_from(instance=setup_instance)
        self.rb_dict = rb_dict

    def check_is_set_up(self):
        """ Check if the RigidBodyAssessment instance is set up. """

        if not self._is_set_up:
            raise ValueError(_error_not_set_up)

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
            map_type=InteractionMapType.CONTACT,
        )

        self._is_set_up = True

    def perform_assessment(self):
        """ Perform the assessment of the rigid body.

        The assessment is performed using the following:
        - af_pipeline.rigid_bodies.rigid_body_assessment.RigidBodyChainAssessment
        - af_pipeline.rigid_bodies.rigid_body_assessment.RigidBodyChainPairAssessment

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
        """

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
            show_interface_residues_only=self.show_interface_residues_only,
        )
        self.rb_cp_assess = RigidBodyChainPairAssessment(
            _mask=_mask,
            as_average=self.as_average,
            show_interface_residues_only=self.show_interface_residues_only,
        )

        self.avg_plddt = self.get_average_plddt(
            only_idr=False,
            only_interface=False,
        )
        self.avg_idr_plddt = self.get_average_plddt(
            only_idr=True,
            only_interface=False,
        )
        self.avg_iplddt = self.get_average_plddt(
            only_idr=False,
            only_interface=True,
        )
        self.avg_idr_iplddt = self.get_average_plddt(
            only_idr=True,
            only_interface=True,
        )
        self.interacting_chains = self.get_interacting_chains()

        self.overall_assessment = self.get_overall_assessment()

    def save_rb_assessment(self, save_path: str):
        """ Save the assessment of the rigid bodies to an Excel file.

        The assessment is saved in an Excel file with three sheets:
        - **"Chain Wise Assessment"**: Contains per chain assessment data.
        - **"Chain Pairwise Assessment"**: Contains per chain pair assessment data.
        - **"Overall Assessment"**: Contains overall assessment data.
        """

        overall_assessment_rows = []

        for col in OVERALL_ASSESSMENT_COLUMNS[self.symmetric_pae]:
            overall_assessment_rows.append({
                "Key": col,
                "Value": self.overall_assessment.get(col, np.nan),
            })
        overall_assessment_df = pd.DataFrame(overall_assessment_rows)

        c_assessment_df = self.rb_c_assess.get_chain_assessment()
        cp_assessment_df = self.rb_cp_assess.get_chain_pair_assessment()

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

    def get_average_plddt(
        self,
        only_idr: bool = False,
        only_interface: bool = False,
    ) -> float:
        """ Get the average pLDDT score across all chains in the rigid body.

        There are four possible outcomes depending on the combination of `only_idr`
        and `only_interface` arguments:
        - **only_idr=False** and **only_interface=False**:<br />
            Average pLDDT score considering all chains and all residues in the rigid body. (Default)
        - **only_idr=True** and **only_interface=False**:<br />
            Average pLDDT score considering only IDR chains and all residues in the rigid body
        - **only_idr=False** and **only_interface=True**:<br />
            Average pLDDT score considering all chains and only interface residues in the rigid body
        - **only_idr=True** and **only_interface=True**:<br />
            Average pLDDT score considering only IDR chains and only interface residues in the rigid body

        ## Arguments:

        - **only_idr (bool, optional)**:<br />
            If True, calculates the average pLDDT score considering only IDR chains.
            If False, considers all chains in the rigid body. (Default)

        - **only_interface (bool, optional)**:<br />
            If True, calculates the average pLDDT score considering only interface residues.
            If False, considers all residues in the rigid body. (Default)

        ## Returns:

        - **avg_plddt (float)**:<br />
            The average pLDDT score across all chains in the rigid body.
        """

        res_selector = {
            True: self.rb_c_assess.per_chain_interface_residues,
            False: self.rb_c_assess.per_chain_residues,
        }
        plddt_selector = {
            True: self.rb_c_assess.per_chain_iplddt,
            False: self.rb_c_assess.per_chain_plddt,
        }

        def get_numerator_from_avg_plddt(only_idr:bool, only_interface:bool) -> float:
            numerator = 0
            for chain_id in self.rb_c_assess.unique_chains:
                if only_idr and chain_id not in self.idr_chains:
                    continue
                weight_factor = res_selector[only_interface].get(chain_id, 0)
                avg_plddt_chain = plddt_selector[only_interface].get(chain_id, 0)
                numerator += weight_factor * avg_plddt_chain
            return numerator

        def get_numerator_from_all_plddt(only_idr:bool, only_interface:bool) -> float:
            return np.sum([
                plddt
                for chain_id, plddt_scores in plddt_selector[only_interface].items()
                if not only_idr or chain_id in self.idr_chains
                for plddt in plddt_scores
            ])

        func_selector = {
            True: get_numerator_from_avg_plddt,
            False: get_numerator_from_all_plddt,
        }

        # only_idr, only_interface
        denominator_selector = {
            (True, False): self.rb_c_assess.total_idr_residues,
            (True, True): self.rb_c_assess.total_interface_idr_residues,
            (False, True): self.rb_c_assess.total_interface_residues,
            (False, False): self.rb_c_assess.total_residues,
        }

        numerator = func_selector[self.as_average](only_idr, only_interface)
        denominator = denominator_selector[(only_idr, only_interface)]
        avg_plddt = numerator / denominator if denominator > 0 else np.nan

        return avg_plddt

    def get_interacting_chains(self) -> List[Tuple[str, str]]:
        """ Get the list of interacting chain pairs in the rigid body.

        ## Returns:

        - **interacting_chains (List[Tuple[str, str]])**:<br />
            A list of tuples, where each tuple contains the IDs of two
            interacting chains in the rigid body.
        """

        interacting_chains = [
            pair
            for pair in self.rb_cp_assess.chain_pairs
            if self.rb_cp_assess.chain_pair_contacts[pair] > 0
        ]
        return interacting_chains

    def get_avg_pae(
        self,
        only_interface: bool = True,
    ) -> Tuple[float, float | None, float | None]:
        """ Get average PAE or iPAE for a rigid body.

        ## Arguments:

        - **only_interface (bool, optional):**:<br />
            If True, calculates the average PAE considering only interface residues.
            If False, considers all residues in the rigid body. (Default)

        ## Returns:

        - **avg_pae (float)**:<br />
            The average PAE or iPAE for the rigid body.
        - **avg_pae_ij (float | None)**:<br />
            The average PAE for i->j across all interacting chain pairs in the rigid body.
        - **avg_pae_ji (float | None)**:<br />
            The average PAE for j->i across all interacting chain pairs in the rigid body.
        """

        avg_pae_ij = avg_pae_ji = None

        # direction, only_interface
        pae_selector = {
            ("sym", True): self.rb_cp_assess.chain_pair_ipae,
            ("sym", False): self.rb_cp_assess.chain_pair_pae,
            ("ij", True): self.rb_cp_assess.chain_pair_ipae_ij,
            ("ij", False): self.rb_cp_assess.chain_pair_pae_ij,
            ("ji", True): self.rb_cp_assess.chain_pair_ipae_ji,
            ("ji", False): self.rb_cp_assess.chain_pair_pae_ji,
        }

        chain_pair_selector = {
            True: [
                pair for pair in self.rb_cp_assess.chain_pairs
                if self.rb_cp_assess.chain_pair_contacts[pair] > 0
            ],
            False: self.rb_cp_assess.chain_pairs,
        }
        residue_selector = {
            True: self.rb_cp_assess.chain_pair_contacts,
            False: self.rb_cp_assess.chain_pair_residue_counts,
        }

        def get_numerator_denominator_from_avg_pae(direction:str, only_interface:bool) -> float:

            numerator = 0
            denominator = 0
            for pair in chain_pair_selector[only_interface]:
                weight_factor = residue_selector[only_interface].get(pair, 0)
                numerator += weight_factor * pae_selector[(direction, only_interface)].get(pair, 0)
                denominator += weight_factor
            return numerator, denominator

        def get_numerator_denominator_from_all_pae(direction:str, only_interface:bool) -> float:

            all_paes = np.array([
                pae
                for pae_vals in pae_selector[(direction, only_interface)].values()
                for pae in pae_vals
            ])

            return np.sum(all_paes), len(all_paes)

        func_selector = {
            True: get_numerator_denominator_from_avg_pae,
            False: get_numerator_denominator_from_all_pae,
        }

        if self.symmetric_pae:
            numerator, denominator = func_selector[self.as_average]("sym", only_interface)
            avg_pae = numerator / denominator if denominator > 0 else np.nan
        else:
            numerator_ij, denominator_ij = func_selector[self.as_average]("ij", only_interface)
            numerator_ji, denominator_ji = func_selector[self.as_average]("ji", only_interface)
            avg_pae_ij = numerator_ij / denominator_ij if denominator_ij > 0 else np.nan
            avg_pae_ji = numerator_ji / denominator_ji if denominator_ji > 0 else np.nan
            if not np.isnan(avg_pae_ij) and not np.isnan(avg_pae_ji):
                avg_pae = (avg_pae_ij + avg_pae_ji) / 2
            else:
                avg_pae = np.nan

        return avg_pae, avg_pae_ij, avg_pae_ji

    # @time_it
    def get_overall_assessment(self):
        """ Get overall assessment of the rigid body.

        Overall assessment includes the following metrics:
        - Average pLDDT across all chains in the rigid body
        - Number of chains in the rigid body
        - Number of interacting chain pairs in the rigid body
        - Number of interface residues in the rigid body
        - Total number of residues in the rigid body
        - Rigid body coverage (total residues in the rigid body / total residues in the prediction)
        - Number of contacts formed in the rigid body
        - Average ipLDDT scores across all chains in the rigid body
        - Average ipLDDT scores for IDR chains in the rigid body
        - Average iPAE scores across all interacting chain pairs in the rigid body
        - Average iPAE (i->j) scores across all interacting chain pairs (if `symmetric_pae` is `False`)
        - Average iPAE (j->i) scores across all interacting chain pairs (if `symmetric_pae` is `False`)

        ## Returns:

        - **overall_assessment (dict)**:<br />
            A dictionary containing overall statistics about the rigid body.
            It includes the number of chains, number of interacting chain pairs,
            interface residues, number of contacts, average ipLDDT,
            average IDR ipLDDT, average iPAE ij, and average iPAE ji.
        """

        overall_assessment = {}

        # Average pLDDT across all chains in the rigid body
        overall_assessment[OverallAssessment.AVERAGE_PLDDT] = self.avg_plddt

        # Number of chains in the rigid body
        overall_assessment[OverallAssessment.NUMBER_OF_CHAINS] = len(self.rb_c_assess.unique_chains)

        # Number of interacting chain pairs in the rigid body
        overall_assessment[OverallAssessment.NUMBER_OF_INTERACTING_CHAIN_PAIRS] = len(self.interacting_chains)

        # Number of interface residues in the rigid body
        overall_assessment[OverallAssessment.INTERFACE_RESIDUES] = self.rb_c_assess.total_interface_residues

        # Total number of residues in the rigid body
        overall_assessment[OverallAssessment.TOTAL_RESIDUES] = self.rb_c_assess.total_residues

        overall_assessment[OverallAssessment.SEQUENCE_COVERAGE] = (
            overall_assessment[OverallAssessment.TOTAL_RESIDUES] /
            self.lengths_dict[MiscStrEnum.TOTAL] # total residues in the prediction
        )

        # Number of contacts formed in the rigid body
        overall_assessment[OverallAssessment.NUMBER_OF_CONTACTS] = sum(self.rb_cp_assess.chain_pair_contacts.values())

        # Average ipLDDT scores across all chains in the rigid body
        overall_assessment[OverallAssessment.AVERAGE_IPLDDT] = self.avg_iplddt

        # Average ipLDDT scores for IDR chains in the rigid body
        overall_assessment[OverallAssessment.AVERAGE_IDR_IPLDDT] = self.avg_idr_iplddt

        # Average iPAE scores across all interacting chain pairs in the rigid body
        avg_ipae, avg_ipae_ij, avg_ipae_ji = self.get_avg_pae(only_interface=True)
        overall_assessment[OverallAssessment.AVERAGE_IPAE] = avg_ipae

        if not self.symmetric_pae:
            overall_assessment[OverallAssessment.AVERAGE_IPAE_IJ] = avg_ipae_ij
            overall_assessment[OverallAssessment.AVERAGE_IPAE_JI] = avg_ipae_ji
        return overall_assessment