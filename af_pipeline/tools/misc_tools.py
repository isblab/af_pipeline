"""
Miscellaneous tools
===================
Utility functions for various tasks.
"""
import copy
import string
import warnings
import pandas as pd
import numpy as np
from typing import Any, Dict
from collections import Counter
from collections import defaultdict

def chain_id_gen():
    """ Generator to sequentially generate 52 alphabets to use as Chain IDs

    TODO: Extend to more than 52 chains if needed

    Yields:
        `i (str)`: Chain ID
    """

    for i in (list(string.ascii_uppercase)):
        yield i
    for i in (list(string.ascii_lowercase)):
        yield i

def create_mask(
    partition_dict: Dict,
    hide_interactions: str = "intra_part",
    masked_value: int = 1,
    unmasked_value: int = 0,
):
    """Create a binary 2D mask.

    Create a binary 2D mask for selecting only inter_part or intra_part interactions. \n
    The mask is created by setting the values of the intra_part or
    inter_part interactions to masked values. \n

    if `hide_interactions=="intra_part"`, intra_part interactions masked.

    if `hide_interactions=="inter_part"`, inter_part interactions masked.

    Args:
        partition_dict (Dict):
            `{"chain_id":chain_length, ..., "total": total_length}`
        hide_interactions (str):
            Hide "intra_part" or "inter_part" interactions.
        masked_value (int):
            Value to set for the masked interactions.
        unmasked_value (int):
            Value to set for the unmasked interactions.

    Returns:
        `new_mask_ (np.ndarray)`:
            Binary 2D mask for selecting only "inter_part" or "intra_part" interactions.

    Examples:

        >>> partition_dict = {"A": 2, "B": 1, "total": 3}
        >>> mask = create_mask(
        ... partition_dict, hide_interactions="intra_part"
        ... )
        >>> print(mask)
        [[1 1 0]
         [1 1 0]
         [0 0 1]]
        >>> mask = create_mask(
        ... partition_dict, hide_interactions="inter_part"
        ... )
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

    assert all(
        isinstance(part_len, int) and part_len > 0
        for part_len in partition_dict.values()
    ), "All chain lengths in partition_dict should be positive integers."

    if "total" not in partition_dict:
        warnings.warn(
            "The partition_dict does not contain 'total'. "
            "Assuming the total length is the sum of all chain lengths.",
            UserWarning,
        )

    sys_len = partition_dict.get("total", sum(partition_dict.values()))

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

def symmetrize_matrix(matrix: np.ndarray | None) -> np.ndarray:
    """Symmetrize a matrix by averaging it with its transpose

    Args:
        matrix (np.ndarray): Input matrix

    Returns:
        `sym_matrix (np.ndarray)`: Symmetrized matrix

    Example:

        >>> matrix = np.array([[1, 2], [3, 4]])
        >>> print(symmetrize_matrix(matrix))
        [[1.  2.5]
         [2.5 4. ]]
    """

    assert isinstance(matrix, np.ndarray), "Input must be a numpy array"
    assert matrix.ndim == 2, "Input must be a 2D array"

    sym_matrix = (matrix + matrix.T) / 2

    return sym_matrix


def fill_up_the_blanks(li: list) -> list:
    """Fill up the blanks in a list

    Args:
        li (list): list with missing numbers

    Returns:
        `new_li (list)`: list with all the missing numbers filled
            up between the minimum and maximum values

    Example:
        >>> fill_up_the_blanks([1, 2, 4, 5])
        [1, 2, 3, 4, 5]
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
        `str`: Residue range string, e.g., "1-3,5-7"

    Example:
        >>> get_key_from_res_range([1, 2, 3, 5, 6, 7])
        '1-3,5-7'
        >>> get_key_from_res_range([1, 2, 3, 5, 6, 7], as_list=True)
        ['1-3', '5-7']
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
                f"{start}-{prev}"
            ) if start != prev else ranges.append(str(start))

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
        `res_range_list (list)`: list of residue numbers

    Example:
        >>> get_res_range_from_key("1-3,5-7")
        [1, 2, 3, 5, 6, 7]
        >>> get_res_range_from_key("1-3,5-7", return_type="set")
        {1, 2, 3, 5, 6, 7}
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

def convert_false_to_true(
    arr: np.ndarray | list,
    threshold:int=2
):
    """
    Convert False values in a binary array to True
    if the patch length is less than or equal to a threshold \n
    A patch is defined as a sequence of consecutive False values

    Args:
        arr (list):
            binary array with False values
        threshold (int, optional):
            Threshold for patch length to convert False to True.
            Defaults to 5.

    Returns:
        `arr (np.ndarray)`:
            binary array with False values converted to True
            if the patch length is less than or equal to threshold

    Example:
        >>> arr = np.array([True, False, False, True, False, False, False])
        >>> convert_false_to_true(arr, threshold=2)
        array([ True,  True,  True,  True, False, False, False])
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
        `duplicate_indices (list | dict)`:\n
            - If return_type is "list", returns a list of indices of duplicate
            elements.\n
            - If return_type is "dict", returns a dictionary with residue IDs
            as keys and duplicate indices as values.\n
            - first or last occurrence of the duplicate element is excluded
            from the output list or dict based on the keep_which parameter.

    Example:
        >>> my_li = [3, 2, 3, 4, 2, 5, 1]
        >>> get_duplicate_indices(my_li)
        [2, 4]
        >>> get_duplicate_indices(my_li, return_type="dict")
        {3: [2], 2: [4]}
        >>> get_duplicate_indices(my_li, keep_which="last")
        [0, 1]
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
    idxs_to_update: dict[Any, list[int]],
    replace_with_avg: bool = False,
    idxs_to_keep: dict = {},
):
    """ Update a list by replacing specified indices.

    Args:
        li (list):
            User-defined list to be updated.
        idxs_to_update (dict):
            Dictionary of indices to update.
        replace_with_avg (bool, optional):
            If `True`, replaces the specified indices with the average value.\n
            If `False`, replaces the specified indices with the values at the
            specified indices in `idxs_to_keep`.\n
            If `idxs_to_keep` is not provided, the first index in each
            list of `idxs_to_update` will be used to replace the indices.
        idxs_to_keep (dict, optional):
            Dictionary of indices to keep.
            If provided, the specified indices will be replaced with the
            values at the specified indices in `idxs_to_keep`.

    Returns:
        `li (list)`: Updated list with specified indices replaced.

    Example:

        >>> li = [1, 2, 3, 4, 5, 6]
        >>> idxs_to_update = {1: [1, 2], 2: [3, 4]}
        >>> idxs_to_keep = {1: 0, 2: 1}
        >>> update_list(
        ... li, idxs_to_update, False, idxs_to_keep
        ... )
        [1.0, 2.0, 5.0, 6.0]
        >>> update_list(
        ... li, idxs_to_update, True, idxs_to_keep
        ... )
        [1.0, 2.5, 4.5, 6.0]
    """

    li1 = np.asarray(copy.deepcopy(li), dtype=float)

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

        li1[start:end] = to_replace

    mask = np.ones(li1.shape[0], dtype=bool)
    mask[idxs_to_update_list] = False
    li1 = li1[mask]

    return li1.tolist()

def update_matrix_row_col(
    matrix: np.ndarray | None,
    idxs_to_update: dict,
    replace_with_avg: bool = False,
    idxs_to_keep: dict = {},
):
    """ Update a square matrix by replacing rows and columns.

    This function updates a square matrix by replacing the specified rows and
    columns.

    `idxs_to_update` is a dictionary which specifies the indices to update.\n
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
            If `True`, replaces the rows and columns with the average value.
        idxs_to_keep (dict, optional):
            Dictionary of indices to keep.
            If provided, the rows and columns will be replaced with the values
            at the specified indices in `idxs_to_keep`.

    Returns:
        `matrix (np.ndarray)`:
            Updated square matrix with specified rows and columns replaced.

    Example:

        >>> matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> idxs_to_update = {1: [0, 1]}
        >>> idxs_to_keep = {1: 0}
        >>> update_matrix_row_col(
        ... matrix, idxs_to_update, False, idxs_to_keep
        ... )
        array([[1, 3],
               [7, 9]])
    """

    if matrix is None:
        raise TypeError(
            "Matrix must be a numpy array, not None."
        )

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

class MatrixPatches:
    """Class to get interacting patches from a binary matrix"""

    matrix: np.ndarray
    """ Binary matrix where rows and columns represent different objects
    (e.g., chains in a protein complex). """

    row_obj: str
    """ Identifier for the rows in the matrix. """

    col_obj: str
    """ Identifier for the columns in the matrix. """

    def __init__(
        self,
        matrix: np.ndarray,
        row_obj: str = "row_obj",
        col_obj: str = "col_obj",
    ):
        self.matrix = matrix
        self.row_obj = row_obj
        self.col_obj = col_obj

    def get_patches_from_matrix(self):
        """Get all interacting patches from a binary matrix

        Args:
            matrix (np.ndarray):
                Binary matrix where rows and columns represent different objects
                (e.g., chains in a protein complex).
            row_obj (str):
                Identifier for the rows in the matrix
            col_obj (str):
                Identifier for the columns in the matrix

        Returns:
            `patches (dict)`: Dictionary of interacting patches

        Example:

            >>> matrix = np.array([
            ... [0, 0, 0, 1],
            ... [0, 1, 1, 1],
            ... [0, 0, 1, 1],
            ... [0, 1, 0, 0],
            ... [0, 1, 0, 1]
            ... ])
            >>> matrix_patches = MatrixPatches(
            ... matrix, row_obj="A", col_obj="B"
            ... )
            >>> matrix_patches.get_patches_from_matrix()
                       A          B
            0  {0, 1, 2}        {3}
            1        {1}  {1, 2, 3}
            2     {1, 2}     {2, 3}
            3     {3, 4}        {1}
            4        {4}        {3}
        """

        assert np.unique(self.matrix).tolist() == [0, 1]
        "Matrix must be binary and non-empty"

        row_sets = self.get_one_sets_from_matrix(self.matrix, axis=0)
        col_sets = self.get_one_sets_from_matrix(self.matrix, axis=1)

        split_row_sets = self.extend_one_sets_by_subsets(row_sets)
        split_col_sets = self.extend_one_sets_by_subsets(col_sets)

        df_row = self.one_sets_to_df(
            split_row_sets, [self.row_obj, self.col_obj]
        )
        df_col = self.one_sets_to_df(
            split_col_sets, [self.col_obj, self.row_obj]
        )

        df_row = self.aggregate_df_rows(df_row, self.col_obj, self.row_obj)
        df_col = self.aggregate_df_rows(df_col, self.row_obj, self.col_obj)

        combined_df = self.combine_dfs(
            df_row, df_col, self.row_obj, self.col_obj
        )

        for col in [self.row_obj, self.col_obj]:
            combined_df[col] = combined_df[col].apply(
                get_res_range_from_key, return_type="set"
            )

        combined_df = self.remove_subset_rows(
            combined_df, self.row_obj, self.col_obj
        )

        return combined_df

    @staticmethod
    def get_one_sets_from_matrix(matrix: np.ndarray, axis: int = 0):
        """Get the indices of 1s in a binary matrix rowwise or columnwise

        Args:
            matrix (np.ndarray):
                Binary matrix where rows and columns represent different objects
                (e.g., chains in a protein complex).
            axis (int, optional):
                0 for rowwise, 1 for columnwise.

        Returns:
            `one_sets (dict)`:
                `{k:v}` where `v` is a set of indices of 1s for key `k`

        Example:

            >>> matrix = np.array([
            ... [1, 0, 1],
            ... [0, 1, 0],
            ... [1, 1, 0]
            ... ])
            >>> MatrixPatches.get_one_sets_from_matrix(matrix=matrix, axis=0)
            {0: {0, 2}, 1: {1}, 2: {0, 1}}
            >>> MatrixPatches.get_one_sets_from_matrix(matrix=matrix, axis=1)
            {0: {0, 2}, 1: {1, 2}, 2: {0}}
        """

        assert np.unique(matrix).tolist() == [0, 1]
        "Matrix must be binary"

        one_sets = {}

        if axis == 0:  # row_sets
            for i in range(matrix.shape[0]):
                one_sets[i] = set(np.where(matrix[i] == 1)[0])

        elif axis == 1:  # col_sets
            for j in range(matrix.shape[1]):
                one_sets[j] = set(np.where(matrix[:, j] == 1)[0])

        return one_sets

    @staticmethod
    def extend_one_sets_by_subsets(one_sets: dict) -> dict:
        """Add the subsets of the sets in list_of_sets to the one_sets

        Args:
            one_sets (dict):
                `{k:v}` where `v` is a set of indices of 1s for key `k`

        Returns:

            new_one_sets (dict):
                `{k:v}` where `v` is a list of sets of indices of 1s for key `k`
                each set is a subset of the original set and is present in
                the values of `one_sets`

        Example:

            >>> one_sets = {
            ... 0: {0, 1, 2, 3, 5, 6},
            ... 1: {1},
            ... 2: {0, 1}
            ... }
            >>> MatrixPatches.extend_one_sets_by_subsets(one_sets)
            {0: [{0, 1, 2, 3}, {5, 6}, {1}, {0, 1}], 1: [{1}], 2: [{1}, {0, 1}]}
        """

        split_sets = MatrixPatches.split_one_sets(one_sets)

        new_one_sets = defaultdict(list)
        list_of_sets = []  # unique sets from split_sets
        list_of_sets = [
            set(x)
            for xs in split_sets.values()
            for x in xs
            if set(x) not in list_of_sets
        ]

        for set1 in list_of_sets:
            for idx, one_set in one_sets.items():
                if set1.issubset(one_set):
                    (
                        new_one_sets[idx].append(set1)
                        if set1 not in new_one_sets[idx]
                        else None
                    )

        return dict(new_one_sets)

    @staticmethod
    def split_one_sets(one_sets: dict) -> dict:
        """Split the sets in `one_sets` into sub-sets such that
        each subset only contains consecutive indices

        Args:
            one_sets (dict):
                `{k:v}` where `v` is a set of indices of 1s for key `k`

        Returns:
            new_one_sets (dict):
                dictionary of lists of lists where each list contains the
                indices of 1s

        Example:

            >>> one_sets = {0: {0, 1, 2, 3, 5, 6}, 1: {1}, 2: {0, 1}}
            >>> MatrixPatches.split_one_sets(one_sets)
            {0: [[0, 1, 2, 3], [5, 6]], 1: [[1]], 2: [[0, 1]]}
        """

        new_one_sets = {}

        for i, one_set in one_sets.items():

            if not isinstance(one_set, set):
                raise TypeError("one_set must be a set")

            sub_sets = MatrixPatches.split_one_set(one_set)
            new_one_sets[i] = sub_sets

        return new_one_sets

    @staticmethod
    def split_one_set(one_set: set | list) -> list:
        """Split a set of indices into sub-sets such that
        each subset only contains consecutive indices

        Args:
            one_set (set | list): Set of indices of 1s

        Returns:
            `sub_sets (list)`:
                List of lists where each list contains the indices of 1s

        Example:

            >>> one_set = {0, 1, 2, 3, 5, 6} \n
            >>> MatrixPatches.split_one_set(one_set)
            [[0, 1, 2, 3], [5, 6]]
        """

        assert isinstance(
            one_set, set | list
        ), "one_set must be a set or a list"

        sub_sets = []

        if isinstance(one_set, list):
            # need to remove duplicates if any
            one_set = set(one_set)

        one_set = sorted(list(one_set))

        for idx, one_pos in enumerate(one_set):

            curr_pos = one_pos
            prev_pos = one_set[idx - 1] if idx > 0 else None

            if idx == 0:
                # If it's the first position, create a new sub-set
                sub_sets.append([curr_pos])

            elif curr_pos - prev_pos == 1:
                # If the current position is consecutive to the previous one
                # add it to the last sub-set
                sub_sets[-1].append(one_pos)

            else:
                # If the current position is not consecutive to the previous one
                # create a new sub-set
                sub_sets.append([curr_pos])

        return sub_sets

    @staticmethod
    def one_sets_to_df(
        one_sets: dict,
        columns: list
    ):
        """Convert a dictionary to a pandas DataFrame

        Args:
            one_sets (dict): Dictionary to convert
            columns (list): Column names

        Returns:
            `df (pd.DataFrame)`:
                DataFrame with the dictionary keys as first column and values
                as second column in columns

        Example:

            >>> one_sets = {
            ... 1: [{1, 2}, {5}],
            ... 2: [{4, 5}, {6}]
            ... }
            >>> columns = ["A", "B"]
            >>> MatrixPatches.one_sets_to_df(one_sets, columns)
               A       B
            0  1  {1, 2}
            1  1     {5}
            2  2  {4, 5}
            3  2     {6}
        """

        if all([isinstance(val, list) for val in one_sets.values()]):

            df_rows = []

            for k, v in one_sets.items():
                for val in v:
                    df_rows.append([str(k), val])

            df = pd.DataFrame(df_rows, columns=columns)

        else:
            raise ValueError("All values in the dictionary must be lists.")

        return df

    @staticmethod
    def aggregate_df_rows(
        df: pd.DataFrame,
        groupby_col: str,
        agg_col: str
    ):
        """Group a DataFrame by a column and aggregate another column

        Args:
            df (pd.DataFrame): DataFrame with groupby_col and agg_col
            groupby_col (str): Column to group by (each value is a set)
            agg_col (str): Column to aggregate (each value is a string)

        Returns:
            `df_group (pd.DataFrame)`:
                Grouped DataFrame with both columns as a set

        Example:

            >>> df = pd.DataFrame({
            ... "A": ["1", "1", "1", "2", "3", "4"],
            ... "B": [{1}, {1,2}, {5}, {4,5}, {1,2}, {1,2}]
            ... })
            >>> MatrixPatches.aggregate_df_rows(df, "B", "A")
                    B          A
            0     {1}        {1}
            1  {1, 2}  {1, 3, 4}
            2  {4, 5}        {2}
            3     {5}        {1}
        """

        df_group = (
            df.groupby(df[groupby_col].map(tuple))[agg_col]
            .apply(",".join)
            .reset_index()
        )
        for idx, row in df_group.iterrows():
            one_set = row[agg_col].split(",")
            one_set = [int(x) for x in one_set]
            one_set = sorted(one_set)
            df_group.at[idx, agg_col] = set(one_set)

        df_group[groupby_col] = df_group[groupby_col].apply(lambda x: set(x))

        return df_group

    @staticmethod
    def combine_dfs(
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        colname_1: str,
        colname_2: str
    ):
        """Combine two DataFrames with columns colname_1 and colname_2
        into a new DataFrame with interacting residues ranges without duplicates

        Args:
            df1 (pd.DataFrame): DataFrame 1
            df2 (pd.DataFrame): DataFrame 2
            colname_1 (str): Column name 1
            colname_2 (str): Column name 2

        Returns:
            `new_df (pd.DataFrame)`:
                Combined DataFrame of interacting residues ranges without
                duplicates

        Example:

            >>> df1 = pd.DataFrame({
            ... "A":[{1, 3, 4}, {1}, {1, 2}, {0, 1, 2, 4}],
            ... "B":[{1}, {1, 2, 3}, {2, 3}, {3}]
            ... })
            >>> df2 = pd.DataFrame({
            ... "A": [{0, 1, 2}, {1}, {1, 2}, {3, 4}, {4}],
            ... "B": [{3}, {1, 2, 3}, {2, 3}, {1}, {1, 3}]
            ... })

            >>> MatrixPatches.combine_dfs(df1, df2, "A", "B")
                 A    B
            0  0-2    3
            1    1  1-3
            2  1-2  2-3
            3  3-4    1
            4    4    1
            5    4    3
            6    1    1
        """

        combined_df = pd.concat([df2, df1], axis=0)
        combined_df.reset_index(drop=True, inplace=True)

        new_df = pd.DataFrame(columns=[colname_1, colname_2])
        df_rows = []

        for _, row in combined_df.iterrows():
            if isinstance(row[colname_1], set) and isinstance(
                row[colname_2], set
            ):
                for res_range1 in get_key_from_res_range(
                    row[colname_1], as_list=True
                ):
                    for res_range2 in get_key_from_res_range(
                        row[colname_2], as_list=True
                    ):
                        df_rows.append([res_range1, res_range2])

        new_df = pd.DataFrame(df_rows, columns=[colname_1, colname_2])
        new_df.drop_duplicates(inplace=True, keep="first")
        new_df.reset_index(drop=True, inplace=True)

        return new_df

    @staticmethod
    def remove_subset_rows(
        df: pd.DataFrame,
        colname_1: str,
        colname_2: str
    ):
        """Remove rows that are subsets of other rows
        (from chatgpt)

        Args:
            df (pd.DataFrame): DataFrame with columns `colname_1` and `colname_2`
            colname_1 (str): column name 1
            colname_2 (str): column name 2

        Returns:
            `filtered_df (pd.DataFrame)`: DataFrame with subset rows removed

        Example:

        >>> df = pd.DataFrame({
        ... "A": [{0, 1, 2}, {1}, {1, 2}, {3, 4}, {4}, {4}, {1}],
        ... "B": [{3}, {1, 2, 3}, {2, 3}, {1}, {1}, {3}, {1}]
        ... })
        >>> MatrixPatches.remove_subset_rows(df, "A", "B")
                   A          B
        0  {0, 1, 2}        {3}
        1        {1}  {1, 2, 3}
        2     {1, 2}     {2, 3}
        3     {3, 4}        {1}
        4        {4}        {3}
        """

        rows_to_keep = []

        for i, row in df.iterrows():

            if not any(
                MatrixPatches.is_subset(
                    row, df.iloc[j], colname_1, colname_2
                )
                for j in range(len(df))
                if i != j
            ):
                rows_to_keep.append(i)

        filtered_df = df.loc[rows_to_keep].reset_index(drop=True)

        return filtered_df

    @staticmethod
    def is_subset(
        row: pd.Series,
        other_row: pd.Series,
        colname_1: str,
        colname_2: str,
    ):
        """Check if row is a subset of other_row for two specified columns

        Args:
            row (pd.Series): Row to check if it is a subset of other_row
            other_row (pd.Series): Row to check against
            colname_1 (str): Column name 1
            colname_2 (str): Column name 2

        Returns:
            `bool`: `True` if row is a subset of `other_row`, `False` otherwise.

        Example:

            >>> row = pd.Series({"A": {0, 1, 2}, "B": {3}})
            >>> other_row = pd.Series({"A": {0, 1, 2, 3}, "B": {3}})
            >>> MatrixPatches.is_subset(row, other_row, "A", "B")
            True
            >>> other_row = pd.Series({"A": {1, 2}, "B": {3}})
            >>> MatrixPatches.is_subset(row, other_row, "A", "B")
            False
        """

        return (
            row[colname_1].issubset(other_row[colname_1])
            and row[colname_2].issubset(other_row[colname_2])
        )

if __name__ == "__main__":

    import doctest
    doctest.testmod()