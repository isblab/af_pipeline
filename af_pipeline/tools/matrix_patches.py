import numpy as np
import pandas as pd
from collections import defaultdict
from af_pipeline.utils.misc_utils import (
    get_key_from_res_range,
    get_res_range_from_key,
)
from af_pipeline.constants.af_constants import MiscStrEnum

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

        Arguments:

        - **matrix (np.ndarray)**:<br />
            Binary matrix where rows and columns represent different objects
            (e.g., chains in a protein complex).

        - **row_obj (str)**:<br />
            Identifier for the rows in the matrix

        - **col_obj (str)**:<br />
            Identifier for the columns in the matrix

        Returns:

        - **patches (dict)**:<br />
            Dictionary of interacting patches

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

        assert np.isin(self.matrix, [0, 1]).all() and np.any(self.matrix), (
            f"Matrix must be binary and non-empty, got {np.unique(self.matrix)}"
        )

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
        """Get the indices of 1s in a binary matrix rowwise or columnwise.

        Arguments:

        - **matrix (np.ndarray)**:<br />
            Binary matrix where rows and columns represent different objects
            (e.g., chains in a protein complex).

        - **axis (int, optional)**:<br />
            0 for rowwise, 1 for columnwise.

        Returns:

        - **one_sets (dict)**:<br />
            `{k:v}` where `v` is a set of indices of 1s for key `k`.

        Example:

            >>> matrix = np.array([
            ... [1, 0, 1],
            ... [0, 1, 0],
            ... [1, 1, 0]
            ... ])
            >>> MatrixPatches.get_one_sets_from_matrix(matrix=matrix, axis=0)
            {0: {np.int64(0), np.int64(2)}, 1: {np.int64(1)}, 2: {np.int64(0), np.int64(1)}}
            >>> MatrixPatches.get_one_sets_from_matrix(matrix=matrix, axis=1)
            {0: {np.int64(0), np.int64(2)}, 1: {np.int64(1), np.int64(2)}, 2: {np.int64(0)}}
        """

        assert np.isin(matrix, [0, 1]).all() and np.any(matrix), (
            f"Matrix must be binary and non-empty, got {np.unique(matrix)}"
        )

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
        """Add the subsets of the sets in list_of_sets to the one_sets.

        Arguments:

        - **one_sets (dict)**:<br />
            `{k:v}` where `v` is a set of indices of 1s for key `k`.

        Returns:

        - **new_one_sets (dict)**:<br />
            `{k:v}` where `v` is a list of sets of indices of 1s for key `k`
            each set is a subset of the original set and is present in
            the values of `one_sets`.

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
        each subset only contains consecutive indices.

        Arguments:

        - **one_sets (dict)**:<br />
            `{k:v}` where `v` is a set of indices of 1s for key `k`.

        Returns:

        - **new_one_sets (dict)**:<br />
            dictionary of lists of lists where each list contains the
            indices of 1s.

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
        each subset only contains consecutive indices.

        Arguments:

        - **one_set (set | list)**:<br />
            Set of indices of 1s.

        Returns:

        - **sub_sets (list)**:<br />
            List of lists where each list contains the indices of 1s.

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
        """Convert a dictionary to a pandas DataFrame.

        Arguments:

        - **one_sets (dict)**:<br />
            Dictionary to convert.

        - **columns (list)**:<br />
            Column names.

        Returns:

        - df (pd.DataFrame)**:<br />
            DataFrame with the dictionary keys as first column and values
            as second column in columns.

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
        """Group a DataFrame by a column and aggregate another column.

        Arguments:

        - **df (pd.DataFrame)**:<br />
            DataFrame with groupby_col and agg_col.

        - **groupby_col (str)**:<br />
            Column to group by (each value is a set).

        - **agg_col (str)**:<br />
            Column to aggregate (each value is a string).

        Returns:

        - **df_group (pd.DataFrame)**:<br />
            Grouped DataFrame with both columns as a set.

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
        df_group[agg_col] = df_group[agg_col].astype(object)
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
        into a new DataFrame with interacting residues ranges without duplicates.

        Arguments:

        - **df1 (pd.DataFrame)**:<br />
            DataFrame 1.

        - **df2 (pd.DataFrame)**:<br />
            DataFrame 2.

        - **colname_1 (str)**:<br />
            Column name 1.

        - **colname_2 (str)**:<br />
            Column name 2.

        Returns:

        - **new_df (pd.DataFrame)**:<br />
            Combined DataFrame of interacting residues ranges without
            duplicates.

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
        new_df.drop_duplicates(inplace=True, keep=MiscStrEnum.FIRST)
        new_df.reset_index(drop=True, inplace=True)

        return new_df

    @staticmethod
    def remove_subset_rows(
        df: pd.DataFrame,
        colname_1: str,
        colname_2: str
    ):
        """Remove rows that are subsets of other rows.
        (from chatgpt)

        Arguments:

        - **df (pd.DataFrame)**:<br />
            DataFrame with columns `colname_1` and `colname_2`.

        - **colname_1 (str)**:<br />
            column name 1.

        - **colname_2 (str)**:<br />
            column name 2.

        Returns:

        - **filtered_df (pd.DataFrame)**:<br />
            DataFrame with subset rows removed.

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
        """Check if row is a subset of other_row for two specified columns.

        Arguments:

        - **row (pd.Series)**:<br />
            Row to check if it is a subset of other_row.

        - **other_row (pd.Series)**:<br />
            Row to check against.

        - **colname_1 (str)**:<br />
            Column name 1.

        - **colname_2 (str)**:<br />
            Column name 2.

        Returns:

        - **(bool)**:<br />
            `True` if row is a subset of `other_row`, `False` otherwise.

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