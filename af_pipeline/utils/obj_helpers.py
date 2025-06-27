import warnings
import pandas as pd
from collections import Counter
import numpy as np
from typing import Dict
import copy

def create_mask(
    partition_dict: Dict,
    hide_interactions: str = "intra_part",
    masked_value: int = 1,
    unmasked_value: int = 0,
) -> np.ndarray:
    """Create a binary 2D mask.

    Create a binary 2D mask for selecting only inter_part or intra_part
    interactions. \n
    The mask is created by setting the values of the intra_part or
    inter_part interactions to masked values. \n
    if `hide_interactions`=="intra_part"`, intra_part interactions masked.
    if `hide_interactions=="inter_part"`, inter_part interactions masked.

    Args:

        partition_dict (Dict):
            Dictionary containing the chain lengths.

        hide_interactions (str):
            Hide intra_part or inter_part interactions.
            Defaults to "intra_part".

        masked_value (int):
            Value to set for the masked interactions.
            Defaults to 1.

        unmasked_value (int):
            Value to set for the unmasked interactions.
            Defaults to 0.

    Returns:

        mask_ (np.ndarray):
            binary 2D mask for selecting only inter_part interactions.

    Examples:

        >>> partition_dict = {"A": 2, "B": 1, "total": 3}
        >>> mask = AfParser.create_mask(
        ... partition_dict, hide_interactions="intra_part"
        ... )
        >>> print(mask)
        [[1 1 0]
            [1 1 0]
            [0 0 0]]
        >>>
        >>> mask = AfParser.create_mask(
        ... partition_dict, hide_interactions="inter_part"
        ...)
        >>> print(mask)
        [[0 0 1]
            [0 0 1]
            [1 1 0]]
    """

    assert (hide_interactions in [
        "intra_part",
        "inter_part",
    ]), "hide_interactions should be either 'intra_part' or 'inter_part'."

    assert masked_value != unmasked_value, \
        "masked_value and unmasked_value should be different."

    sys_len = partition_dict.get("total", np.sum(list(partition_dict.values())))
    mask_ = np.full((sys_len, sys_len), unmasked_value)

    prev = 0
    for chain in partition_dict:
        if chain == "total":
            continue
        l = partition_dict[chain]
        curr = prev + l
        mask_[prev:curr:, prev:curr] = masked_value
        prev += l

    if hide_interactions == "intra_part":
        return mask_

    elif hide_interactions == "inter_part":
        new_mask_ = np.full((sys_len, sys_len), unmasked_value)
        new_mask_[mask_ == unmasked_value] = masked_value
        return new_mask_

def symmetrize_matrix(matrix: np.ndarray) -> np.ndarray:
    """Symmetrize a matrix by averaging it with its transpose

    Args:
        matrix (np.ndarray): input matrix

    Returns:
        sym_matrix (np.ndarray): symmetrized matrix

    Example:
    matrix = np.array([[1, 2], [3, 4]]) \n
    symmetrize_matrix(matrix) -> np.array([[1, 2.5], [2.5, 4]])
    """

    assert isinstance(matrix, np.ndarray), "Input must be a numpy array"
    assert matrix.ndim == 2, "Input must be a 2D array"

    sym_matrix = (matrix + matrix.T) / 2

    return sym_matrix


def fill_up_the_blanks(li: list) -> list:
    """Fill up the blanks in a list

    Example:
        fill_up_the_blanks([1, 2, 4, 5]) -> [1, 2, 3, 4, 5]

    Args:
        li (list): list with missing numbers

    Returns:
        new_li (list): list with all the missing numbers filled
            up between the minimum and maximum values
    """

    min_li_val = min(li)
    max_li_val = max(li)

    new_li = [x for x in range(min_li_val, max_li_val+1)]

    return new_li


def get_key_from_res_range(
    res_range: list,
    as_list=False
) -> str | list:
    """Returns a residue range string from a list of residue numbers.

    Args:
        res_range (list): List of residue numbers, e.g., [1, 2, 3, 5, 6, 7]

    Returns:
        str: Residue range string, e.g., "1-3,5-7"

    Example:
    get_key_from_res_range([1, 2, 3, 5, 6, 7]) -> "1-3,5-7" \n
    get_key_from_res_range([1, 2, 3, 5, 6, 7], as_list=True) -> ["1-3", "5-7"]
    """

    if not res_range:
        return ""

    res_range = sorted(res_range)
    ranges = []
    start = prev = res_range[0]

    for num in res_range[1:]:
        if num == prev + 1:
            prev = num
        else:
            ranges.append(
                f"{start}-{prev}") if start != prev else ranges.append(str(start)
            )
            start = prev = num

    if start == prev:
        ranges.append(str(start))

    else:
        ranges.append(f"{start}-{prev}")

    if as_list:
        return ranges

    else:
        return ",".join(ranges)


def get_res_range_from_key(
    res_range: str,
    return_type: str = "list"
) -> list | set:
    """Convert a residue range string to a list of residue numbers

    Args:
        res_range (str): residue range string

    Returns:
        res_range_list (list): list of residue numbers

    Example:
    get_res_range_from_key("1-3,5-7") -> [1, 2, 3, 5, 6, 7]
    """

    res_range_list = []

    for res_range in res_range.split(","):

        if "-" in res_range:
            start, end = map(int, res_range.split("-"))
            res_range_list.extend(list(range(start, end+1)))

        else:
            res_range_list.append(int(res_range))

    if return_type == "set":
        return set(res_range_list)

    return res_range_list

def convert_false_to_true(arr: np.ndarray | list, threshold:int=5):
    """ Convert False values in a binary array to True if the patch length is less than or equal to a threshold
        a patch is defined as a sequence of consecutive False values

    Args:
        arr (list): binary array with False values
        threshold (int, optional): _description_. Defaults to 5.

    Returns:
        _type_: _description_
    """

    if isinstance(arr, list):
        arr = np.array(arr)

    where_false = list(np.argwhere(arr == False).flatten())

    false_patches = get_key_from_res_range(where_false, as_list=True)

    for patch in false_patches:

        patch = patch.split("-")
        if len(patch) == 1:
            patch.append(patch[0])

        start = int(patch[0])
        end = int(patch[1])

        if end - start + 1 <= threshold:
            arr[start:end + 1] = True

    return arr

def get_duplicate_indices(
    my_li: list,
    return_type: str = "list",
    keep_which: str | None = "first",
):
    """ Get the indices of duplicate elements in a list.

    Args:

        my_li (list):
            List of elements to check for duplicates.

        return_type (str, optional):
            Output type, either "list" or "dict".

        keep_which (None | str, optional):
            - If "first", keeps the first occurrence of the duplicate element.
            - If "last", keeps the last occurrence of the duplicate element.
            - If None, keeps all occurrences of the duplicate element.
            Defaults to "first".

    Returns:

        duplicate_indices (list | dict):
            - If return_type is "list", returns a list of indices of duplicate 
            elements.
            - If return_type is "dict", returns a dictionary with residue IDs 
            as keys and duplicate indices as values.
            - first or last occurrence of the duplicate element is excluded 
            from the output list or dict based on the keep_which parameter.
    """

    token_counts = Counter(my_li)
    # Get the repeated residue IDs, these are pTMs or small molecules
    repeated_tokens = [
        token_id
        for token_id, count
        in token_counts.items()
        if count > 1
    ]

    duplicate_indices = []

    for token_id in repeated_tokens:
        indices = [i for i, x in enumerate(my_li) if x == token_id]
        if keep_which == "first":
            indices.pop(0)
        elif keep_which == "last":
            indices.pop(-1)
        duplicate_indices.extend(indices)

    if return_type == "list":
        return duplicate_indices

    elif return_type == "dict":

        duplicate_indices = {}

        for token_id in repeated_tokens:
            indices = [i for i, x in enumerate(my_li) if x == token_id]
            if keep_which == "first":
                indices.pop(0)
            elif keep_which == "last":
                indices.pop(-1)
            duplicate_indices[token_id] = indices

        return duplicate_indices

def update_list(
    li: list,
    idxs_to_update: dict,
    replace_with_avg: bool = False,
    idxs_to_keep: dict = {},
):
    """ Update a list by replacing specified indices.

    Args:
        li (list): User-defined list to be updated.
        idxs_to_update (dict): Dictionary of indices to update.
        replace_with_avg (bool, optional): _description_. Defaults to False.
        idxs_to_keep (dict, optional): _description_. Defaults to {}.

    Returns:
        li (list): Updated list with specified indices replaced.
    """

    li = np.array(li, dtype=float)
    li = copy.deepcopy(li)

    idxs_to_update_list = [
        idx for idxs in idxs_to_update.values() for idx in idxs[1:]
    ]

    for token_id, idxs in idxs_to_update.items():
        start, end = idxs[0], idxs[-1]
        end += 1

        to_replace = li[start:end]

        if replace_with_avg is True:
            to_replace = np.mean(to_replace)

        else:
            idx_ = idxs_to_keep.get(token_id, 0)
            to_replace = to_replace[idx_]

        li[start:end] = to_replace

    mask = np.ones(li.shape[0], dtype=bool)
    mask[idxs_to_update_list] = False
    li = li[mask]

    return li.tolist()

def update_matrix_row_col(
    matrix: np.ndarray,
    idxs_to_update: dict,
    replace_with_avg: bool = False,
    idxs_to_keep: dict = {},
):
    """ Update a square matrix by replacing rows and columns.

    This function updates a square matrix by replacing the specified rows and 
    columns.

    idxs_to_update is a dictionary which specifies the indices to update.
    For each key, a list of indices is provided. In the output matrix, the rows 
    and columns corresponding to these indices will be replaced with the 
    average if `replace_with_avg` is True.

    If `replace_with_avg` is False, the rows and columns will be replaced with 
    the values at the specified indices in `idxs_to_keep`.

    If `idxs_to_keep` is not provided, the first index in each list of 
    `idxs_to_update` will be used to replace the columns and rows.

    Args:

        matrix (np.ndarray):
            Square matrix to be updated.

        idxs_to_update (dict):
            Dictionary of indices to update.
            e.g. {token1: [idx1, idx2, ...], token2: [idx3, idx4, ...]}

        replace_with_avg (bool, optional):
            If True, replaces the rows and columns with the average value.
            Defaults to False.

        idxs_to_keep (dict, optional):
            Dictionary of indices to keep.
            If provided, the rows and columns will be replaced with the values 
            at the specified indices in `idxs_to_keep`.
            Defaults to {}.

    Returns:

        matrix (np.ndarray):
            Updated square matrix with specified rows and columns replaced.
    """

    assert matrix.shape[0] == matrix.shape[1], "Matrix must be square (NxN)."

    assert all(isinstance(idxs, list) for idxs in idxs_to_update.values()), \
        "All values in idxs_to_update must be lists of indices."

    if all(len(idxs) <= 1 for idxs in idxs_to_update.values()):
        warnings.warn(
            """
            All indices in idxs_to_update have length <= 1.
            There's nothing to update in the matrix.
            Returning the original matrix.
            """
        )
        return matrix

    assert all(isinstance(idx, int) for idx in idxs_to_keep.values()), \
        "All values in idxs_to_keep must be integers."

    assert set(idxs_to_keep.keys()).issubset(idxs_to_update.keys()), \
        "Keys in idxs_to_keep must be a subset of keys in idxs_to_update."

    for token_id, idx_to_keep in idxs_to_keep.items():
        assert idx_to_keep < len(idxs_to_update[token_id]), \
            f"Index {idx_to_keep} for token {token_id} is out of bounds " \
            f"for the list of indices {idxs_to_update[token_id]}."

    matrix = copy.deepcopy(matrix)
    things_to_update = []

    idxs_to_update_list = [
        idx for idxs in idxs_to_update.values() for idx in idxs[1:]
    ]

    for token_id, idxs in idxs_to_update.items():
        start, end = idxs[0], idxs[-1]
        end += 1

        center_val = matrix[start:end, start:end]
        col_val = matrix[start:end, :]
        row_val = matrix[:, start:end]

        if replace_with_avg is True:

            center_val = np.mean(center_val)
            col_val = np.mean(col_val, axis=0)
            row_val = np.mean(row_val, axis=1)

            for idxs_ in idxs_to_update.values():
                start_, end_ = idxs_[0], idxs_[-1]
                end_ += 1

                row_val[start_:end_] = np.mean(row_val[start_:end_])
                col_val[start_:end_] = np.mean(col_val[start_:end_])

        else:
            idx_ = idxs_to_keep.get(token_id, 0)
            center_val = center_val[idx_][idx_]
            col_val = col_val[idx_, :]
            row_val = row_val[:, idx_]

            for token_id_, idxs_ in idxs_to_update.items():
                idx2_ = idxs_to_keep.get(token_id_, 0)
                start_, end_ = idxs_[0], idxs_[-1]
                end_ += 1

                row_val[start_:end_] = row_val[start_:end_][idx2_]
                col_val[start_:end_] = col_val[start_:end_][idx2_]

        things_to_update.append(
            {
                "pos": start,
                "center_val": center_val,
                "row_val": row_val,
                "col_val": col_val,
            }
        )

    for thing in things_to_update:
        start = thing["pos"]
        matrix[:, start] = thing["row_val"]
        matrix[start, :] = thing["col_val"]
        matrix[start, start] = thing["center_val"]

    mask = np.ones(matrix.shape[0], dtype=bool)
    mask[idxs_to_update_list] = False
    matrix = matrix[mask][:, mask]

    return matrix