"""
Miscellaneous utils
===================
Utility functions for various tasks.
"""
import copy
import time
import string
import warnings

import yaml
import numpy as np
from typing import Any, Dict
from collections import Counter
from af_pipeline.constants.af_constants import (
    RES_SEPARATOR,
    RES_SPLITTER,
    SEED_MULTIPLIER,
    RANDOM_SEED,
    MaskedInteractionType,
    MaskedInteractionValue,
    MiscStrEnum,
    ReturnType,
    UpdateConfigMode,
    ColorMapScheme,
    BinaryColorMap,
)
from typing import List
import random
from functools import wraps

def add_attribute(
    config_yaml: dict,
    attribute_name: str,
    attribute_value: Any,
    mode: UpdateConfigMode = UpdateConfigMode.REPLACE,
    add_first: bool = False,
):
    """Update a generic attribute in the config file.

    Arguments:

    - **config_yaml (dict)**:<br />
        Configuration dictionary.

    - **attribute_name (str)**:<br />
        Name of the attribute to update.

    - **attribute_value (Any)**:<br />
        Value of the attribute to update.

    - **mode (str, optional)**:<br />
        Mode to update the config file.
    """

    from af_pipeline.constants.af_constants import ConfigYaml

    af_input_jobs = config_yaml.get(ConfigYaml.AF_INPUT_JOBS, {})

    for job_cycle, job_sets in af_input_jobs.items():
        for idx, job_set in enumerate(job_sets):
            if attribute_name in job_set:
                af_input_jobs[job_cycle][idx][attribute_name] = (
                    attribute_value[job_cycle][idx]
                )
            elif add_first:
                af_input_jobs[job_cycle][idx] = {
                    attribute_name: attribute_value[job_cycle][idx],
                    **job_set,
                }
            else:
                af_input_jobs[job_cycle][idx] = {
                    **job_set,
                    attribute_name: attribute_value[job_cycle][idx],
                }

    updated_config = update_config(
        config_yaml=config_yaml,
        updates={ConfigYaml.AF_INPUT_JOBS: af_input_jobs},
        mode=mode
    )

    return updated_config

def update_config(
    config_yaml: dict,
    updates: dict = None,
    mode: UpdateConfigMode = UpdateConfigMode.REPLACE,
):
    """Update config file with a new field or update an existing field.

    Arguments:

    - **config_yaml (dict)**:<br />
        Configuration dictionary.

    - **updates (dict, optional)**:<br />
        Fields to update in the config file.

    - **mode (str, optional)**:<br />
        Mode to update the config file. ("append" or "replace").

    """

    update_fields = list(updates.keys()) if updates else []

    if len(update_fields) == 0:

        print("No fields to update in config")
        return None

    existing_fields = list(config_yaml.keys())

    for field in update_fields:

        add_field = False

        if field in existing_fields:
            if mode == UpdateConfigMode.REPLACE:
                config_yaml[field] = updates[field]
            elif mode == UpdateConfigMode.SOFT_REPLACE:
                #! only update if the field is not already set
                if config_yaml[field] is None or config_yaml[field] == "":
                    config_yaml[field] = updates[field]
            else:
                raise ValueError(
                    "Invalid mode. Use one of the following: "\
                    f": {list(UpdateConfigMode)}"
                )

        else:
            print(f"{field} not found in config")
            print("Adding field to config")
            add_field = True

        if add_field:
            config_yaml[field] = updates[field]
            add_field = False

    print(f"Config dict updated with {update_fields}")
    return config_yaml

def time_it(func):
    """ Decorator to measure the execution time of a function.

    Add `@time_it` above any function definition to measure its execution time.

    Args:
        func (function): Function to be decorated.

    Returns:
        wrapper (function): Decorated function with execution time measurement.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        print(f"Function {func.__name__} took: {end_time - start_time:.4f} seconds")
        return result
    return wrapper

def generate_seeds(num_seeds: int, set_seed: int=RANDOM_SEED) -> List[int]:
    """Generate `model_seeds`."""

    random.seed(set_seed)
    model_seeds = random.sample(range(1, SEED_MULTIPLIER * num_seeds), num_seeds)

    return model_seeds

def chain_id_gen():
    """ Generator to sequentially generate 52 alphabets to use as Chain IDs

    TODO: Extend to more than 52 chains if needed

    Yields:
        `i (str)`: Chain ID

    Example:

        >>> gen = chain_id_gen()
        >>> _chains= []
        >>> for _ in range(52):
        ...     _chains.append(next(gen))
        >>> print("".join(_chains))
        ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz
    """

    for i in (list(string.ascii_uppercase)):
        yield i
    for i in (list(string.ascii_lowercase)):
        yield i

def create_mask(
    partition_dict: Dict,
    hide_interactions: str | MaskedInteractionType = MaskedInteractionType.INTRA_PART,
    masked_value: bool | int | float | str | MaskedInteractionValue = MaskedInteractionValue.MASKED_V,
    unmasked_value: bool | int | float | str | MaskedInteractionValue = MaskedInteractionValue.UNMASKED_V,
):
    """ Create a binary 2D mask.

    Create a binary 2D mask for selecting only inter_part or intra_part interactions. \n
    The mask is created by setting the values of the intra_part or
    inter_part interactions to masked values. \n

    if `hide_interactions=="intra_part"`, intra_part interactions masked.

    if `hide_interactions=="inter_part"`, inter_part interactions masked.

    ## Arguments:

    - **partition_dict (Dict)**:<br />
        `{"chain_id":chain_length, ..., "total": total_length}`

    - **hide_interactions (str | MaskedInteractionType, optional):**:<br />
        Hide "intra_part" or "inter_part" interactions. Defaults to "intra_part".

    - **masked_value (bool | int | float | str | MaskedInteractionValue, optional):**:<br />
        Value to set for the masked interactions. Defaults to 1.

    - **unmasked_value (bool | int | float | str | MaskedInteractionValue, optional):**:<br />
        Value to set for the unmasked interactions. Defaults to 0.

    ## Returns:

    - **np.ndarray**:<br />
        Binary 2D mask for selecting only "inter_part" or "intra_part" interactions.

    ## Examples:

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

    assert hide_interactions in list(MaskedInteractionType), (
        f"hide_interactions should be one of the following: "\
        f"{list(MaskedInteractionType)}"
    )

    assert masked_value != unmasked_value, \
        "masked_value and unmasked_value should be different."

    assert all(
        isinstance(part_len, int) and part_len > 0
        for part_len in partition_dict.values()
    ), "All chain lengths in partition_dict should be positive integers."

    if MiscStrEnum.TOTAL not in partition_dict:
        warnings.warn(
            "The partition_dict does not contain 'total'. "
            "Assuming the total length is the sum of all chain lengths.",
            UserWarning,
        )

    sys_len = partition_dict.get(MiscStrEnum.TOTAL, sum(partition_dict.values()))

    mask_ = np.full((sys_len, sys_len), unmasked_value)

    prev = 0
    for chain in partition_dict:
        if chain == MiscStrEnum.TOTAL:
            continue
        l = partition_dict[chain]
        curr = prev + l
        mask_[prev:curr:, prev:curr] = masked_value
        prev += l

    if hide_interactions == MaskedInteractionType.INTRA_PART:
        return mask_

    elif hide_interactions == MaskedInteractionType.INTER_PART:
        new_mask_ = np.full((sys_len, sys_len), unmasked_value)
        new_mask_[mask_ == unmasked_value] = masked_value

        return new_mask_

def symmetrize_matrix(matrix: np.ndarray | None) -> np.ndarray:
    """ Symmetrize a matrix by averaging it with its transpose

    ## Arguments:

    - **matrix (np.ndarray | None)**:<br />
        Input matrix to be symmetrized.

    ## Returns:

    - **np.ndarray**:<br />
        Symmetrized matrix.

    ## Examples:

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
    """ Fill up the blanks in a list.

    ## Arguments:

    - **li (list)**:<br />
        List with missing numbers.

    ## Returns:

    - **list**:<br />
        List with all the missing numbers filled up between the minimum and
        maximum values.

    ## Examples:

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
    """ Returns a residue range string from a list of residue numbers.

    ## Arguments:

    - **res_range (list)**:<br />
        List of residue numbers, e.g., [1, 2, 3, 5, 6, 7]

    - **as_list (bool, optional):**:<br />
        If `True`, returns a list of residue range strings, e.g., ["1-3", "5-7"].
        If `False`, returns a single string with residue ranges separated by commas,
        e.g., "1-3,5-7". Defaults to `False`.

    ## Returns:

    - **str | list**:<br />
        Residue range string or list of residue range strings.

    ## Examples:

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
                f"{start}{RES_SEPARATOR}{prev}"
            ) if start != prev else ranges.append(str(start))

            start = prev = num

    if start == prev:
        ranges.append(str(start))

    else:
        ranges.append(f"{start}{RES_SEPARATOR}{prev}")

    if as_list:
        return ranges

    else:
        return RES_SPLITTER.join(ranges)

def get_res_range_from_key(
    res_range: str,
    return_type: ReturnType = ReturnType.LIST,
) -> list | set:
    """ Convert a residue range string to a list of residue numbers

    ## Arguments:

    - **res_range (str)**:<br />
        Residue range string, e.g., "1-3,5-7"

    - **return_type (ReturnType, optional):**:<br />
        Return type, either `list` or `set`. Defaults to `list`.

    ## Returns:

    - **list | set**:<br />
        List or set of residue numbers, e.g., [1, 2, 3, 5, 6, 7] or {1, 2, 3, 5, 6, 7}.

    ## Examples:

    >>> get_res_range_from_key("1-3,5-7")
    [1, 2, 3, 5, 6, 7]
    >>> get_res_range_from_key("1-3,5-7", return_type="set")
    {1, 2, 3, 5, 6, 7}
    """

    res_range_list = []

    for res_range in res_range.split(RES_SPLITTER):

        if f"{RES_SEPARATOR}" in res_range:
            start, end = map(int, res_range.split(RES_SEPARATOR))
            res_range_list.extend(list(range(start, end+1)))

        else:
            res_range_list.append(int(res_range))

    if return_type == ReturnType.SET:
        return set(res_range_list)

    return res_range_list

def convert_false_to_true(
    arr: np.ndarray | list,
    threshold:int=2
) -> np.ndarray:
    """ Convert False values in a binary array to True
    if the patch length is less than or equal to a threshold \n
    A patch is defined as a sequence of consecutive False values

    ## Arguments:

    - **arr (np.ndarray | list)**:<br />
        Binary array with False values.

    - **threshold (int, optional):**:<br />
        Threshold for patch length to convert False to True. Defaults to 2.

    ## Returns:

    - **np.ndarray**:<br />
        Binary array with False values converted to True if the patch length is
        less than or equal to threshold.
    """
    if isinstance(arr, list):
        arr = np.array(arr)

    where_false = list(np.argwhere(arr == False).flatten())

    false_patches = get_key_from_res_range(where_false, as_list=True)

    for patch in false_patches:

        patch = patch.split(RES_SEPARATOR)
        if len(patch) == 1:
            patch.append(patch[0])

        start = int(patch[0])
        end = int(patch[1])

        if end - start + 1 <= threshold:
            arr[start:end + 1] = True

    return arr

def get_duplicate_indices(
    my_li: list,
    return_type: str | ReturnType = ReturnType.LIST,
    keep_which: str | MiscStrEnum = MiscStrEnum.FIRST,
) -> list | dict:
    """ Get the indices of duplicate elements in a list.

    ## Arguments:

    - **my_li (list)**:<br />
        List of elements to check for duplicates.

    - **return_type (str | ReturnType, optional):**:<br />
        Output type, either "list" or "dict".
        Defaults to "list".

    - **keep_which (str | MiscStrEnum, optional):**:<br />
        - If "first", keeps the first occurrence of the duplicate element.
        - If "last", keeps the last occurrence of the duplicate element.
        - If None, keeps all occurrences of the duplicate element.
        Defaults to "first".

    ## Returns:

    - **list | dict**:<br />
        - If return_type is "list", returns a list of indices of duplicate
        elements.
        - If return_type is "dict", returns a dictionary with residue IDs
        as keys and duplicate indices as values.
        - first or last occurrence of the duplicate element is excluded
        from the output list or dict based on the keep_which parameter.

    ## Examples:
    >>> my_li = [3, 2, 3, 4, 2, 5, 1]
    >>> get_duplicate_indices(my_li)
    [2, 4]
    >>> get_duplicate_indices(my_li, return_type="dict")
    {3: [2], 2: [4]}
    >>> get_duplicate_indices(my_li, keep_which="last")
    [0, 1]
    >>> get_duplicate_indices(my_li, keep_which="last", return_type="dict")
    {3: [0], 2: [1]}
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
        if keep_which == MiscStrEnum.FIRST:
            indices.pop(0)
        elif keep_which == MiscStrEnum.LAST:
            indices.pop(-1)
        duplicate_indices.extend(indices)

    if return_type == ReturnType.LIST:
        return duplicate_indices

    elif return_type == ReturnType.DICT:

        duplicate_indices = {}

        for token_id in repeated_tokens:
            indices = [i for i, x in enumerate(my_li) if x == token_id]
            if keep_which == MiscStrEnum.FIRST:
                indices.pop(0)
            elif keep_which == MiscStrEnum.LAST:
                indices.pop(-1)
            duplicate_indices[token_id] = indices

        return duplicate_indices

def update_list(
    li: list,
    idxs_to_update: dict[Any, list[int]],
    replace_with_avg: bool = False,
    idxs_to_keep: dict = {},
) -> list:
    """ Update a list by replacing specified indices.

    ## Arguments:

    - **li (list)**:<br />
        User-defined list to be updated.

    - **idxs_to_update (dict[Any, list[int]])**:<br />
        Dictionary of indices to update.
        e.g. {token1: [idx1, idx2, ...], token2: [idx3, idx4, ...]}
            For each key, a list of indices is provided. In the output list, the
            values corresponding to these indices will be replaced with the
            average if `replace_with_avg` is True.

    - **replace_with_avg (bool, optional):**:<br />
        If `True`, replaces the specified indices with the average value.\n
        If `False`, replaces the specified indices with the values at the
        specified indices in `idxs_to_keep`.\n
        If `idxs_to_keep` is not provided, the first index in each
        list of `idxs_to_update` will be used to replace the indices.

    - **idxs_to_keep (dict, optional):**:<br />
        Dictionary of indices to keep.
        If provided, the specified indices will be replaced with the
        values at the specified indices in `idxs_to_keep`.

    ## Returns:

    - **list**:<br />
        Updated list with specified indices replaced.

    ## Examples:

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
) -> np.ndarray:
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

    ## Arguments:

    - **matrix (np.ndarray | None)**:<br />
        Square matrix to be updated.

    - **idxs_to_update (dict)**:<br />
        Dictionary of indices to update.

    - **replace_with_avg (bool, optional):**:<br />
        If `True`, replaces the rows and columns with the average value.\n
        If `False`, replaces the rows and columns with the values at the
        specified indices in `idxs_to_keep`.\n
        If `idxs_to_keep` is not provided, the first index in each list of
        `idxs_to_update` will be used to replace the columns and rows.

    - **idxs_to_keep (dict, optional):**:<br />
        Dictionary of indices to keep.

    ## Returns:

    - **np.ndarray**:<br />
        Updated square matrix with specified rows and columns replaced.

    ## Examples:

    >>> matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> idxs_to_update = {1: [0, 1]}
    >>> idxs_to_keep = {1: 0}
    >>> update_matrix_row_col(
    ... matrix, idxs_to_update, False, idxs_to_keep
    ... )
    array([[1., 3.],
           [7., 9.]])
    >>> matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> idxs_to_update = {1: [0, 1]}
    >>> update_matrix_row_col(
    ... matrix, idxs_to_update, True, {}
    ... )
    array([[3. , 4.5],
           [7.5, 9. ]])
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

    matrix = copy.deepcopy(matrix).astype(float)
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

def extract_protein_chain_mapping(
    protein_chain_mapping: dict,
) -> dict[str, str]:
    """ Extract the protein chain mapping from the provided dictionary.

    For e.g., if the user provides the following mapping:
    ```python
    {
        "ProteinA": ["A", "B"],
        "ProteinB": ["C"]
    }
    ```
    The function will return the following dictionary:
    ```python
    {
        "A": "ProteinA",
        "B": "ProteinA",
        "C": "ProteinB"
    }
    ```

    ## Arguments:

    - **protein_chain_mapping (dict)**:<br />
        Protein-to-chain map.

    ## Returns:

    - **protein_chain_map (dict)**:<br />
        Dictionary with chain IDs as keys and protein names as values.

    ## Examples:

    >>> protein_chain_mapping = {
    ... "ProteinA" : ["A", "B"],
    ... "ProteinB": ["C"]
    ... }
    >>> sorted(extract_protein_chain_mapping(protein_chain_mapping).items())
    [('A', 'ProteinA'), ('B', 'ProteinA'), ('C', 'ProteinB')]
    """

    chain_protein_map = {}

    if len(protein_chain_mapping) == 0:
        return chain_protein_map

    for protein_name, chain_ids in protein_chain_mapping.items():
        for chain_id in chain_ids:
            if chain_id not in chain_protein_map:
                chain_protein_map[chain_id] = protein_name

    return chain_protein_map

def generate_cmap(
    n: int,
    scheme: str | ColorMapScheme = ColorMapScheme.SOFT_WARM
) -> list[str]:
    """ Generate a list of n colors based on the specified color scheme.
    (modfied from chatgpt)

    ## Arguments:

    - **n (int)**:<br />
        Number of colors to generate.

    - **scheme (str | ColorMapScheme, optional):**:<br />
        Color scheme to use for generating colors. Choose from the following:
        - `standard`: Standard RGB colors.
        - `non_bright`: Muted colors with reduced brightness.
        - `earth_tone`: Colors inspired by natural earth tones.
        - `cool_tone`: Colors with a cool, calming palette.
        - `soft_warm`: Colors with a soft, warm palette.
        - `contrasting_non_bright`: Alternating between darker and lighter muted tones.
        Defaults to `soft_warm`.

    ## Returns:

    - **list**:<br />
        List of n colors in hexadecimal format.
    """

    import random
    import time

    colors = set()
    start = time.time()

    assert n > 0; "Number of colors must be greater than 0"
    if n != 2:
        assert scheme != ColorMapScheme.BINARY; "Binary scheme is only valid for 2 colors"

    if n == 2 and scheme == ColorMapScheme.BINARY:
        return [BinaryColorMap.ZERO_C, BinaryColorMap.ONE_C]

    while len(colors) < n:
        if scheme == ColorMapScheme.STANDARD:
            r = random.randint(0, 255)
            g = random.randint(0, 255)
            b = random.randint(0, 255)

        elif scheme == ColorMapScheme.NON_BRIGHT:
            r = random.randint(50, 180)
            g = random.randint(50, 180)
            b = random.randint(50, 180)

        elif scheme == ColorMapScheme.EARTH_TONE:
            r = random.randint(100, 180)
            g = random.randint(60, 140)
            b = random.randint(40, 120)

        elif scheme == ColorMapScheme.COOL_TONE:
            r = random.randint(50, 120)
            g = random.randint(100, 180)
            b = random.randint(120, 255)

        elif scheme == ColorMapScheme.SOFT_WARM:
            r = random.randint(180, 255)
            g = random.randint(130, 200)
            b = random.randint(90, 160)

        elif scheme == ColorMapScheme.CONTRASTING_NON_BRIGHT:
            if len(colors) % 2 == 0:  # Alternate between darker and lighter muted tones
                r = random.randint(40, 120)
                g = random.randint(40, 120)
                b = random.randint(40, 120)
            else:
                r = random.randint(140, 200)
                g = random.randint(140, 200)
                b = random.randint(140, 200)

        else:
            raise ValueError(
                "Invalid scheme. Choose from "\
                f"{list(ColorMapScheme)}"
            )

        color = "#{:02x}{:02x}{:02x}".format(r, g, b)
        colors.add(color)

        if time.time() - start > 10:
            break

    return list(colors)

if __name__ == "__main__":

    import doctest
    doctest.testmod()