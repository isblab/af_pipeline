"""
[file_utils](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/utils/file_utils)
===================

- Utility functions for handling files, including reading and writing json, fasta and pickle files.

"""

import yaml
import json
from typing import Any

def read_yaml(
    file_path: str
) -> Any:
    """Read a yaml file and return a dictionary.

    ## Arguments:

    - **file_path (str)**:<br />
        Path to yaml file.

    ## Returns:
    - **data (Any)**:<br />
        Data from yaml file.
    """

    with open(file_path, "r") as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    return data

def write_json(
    file_path: str,
    data: Any,
    indent: int | None = None,
):
    """Write data to a json file.

    ## Arguments:

    - **file_path (str)**:<br />
        Path to json file.

    - **data (Any)**:<br />
        Data to write.
    """

    with open(file_path, "w") as f:
        json.dump(data, f, indent=indent)


def read_json(file_path: str) -> Any:
    """Load a json file.

    ## Arguments:

    - **file_path (str)**:<br />
        Path to json file.

    ## Returns:

    - **data (Any)**:<br />
        Data from json file.
    """

    with open(file_path, "r") as f:
        data = json.load(f)

    return data

def read_pkl(file_path: str) -> Any:
    """Load a pickle file.

    ## Arguments:

    - **file_path (str)**:<br />
        Path to pickle file.

    ## Returns:

    - **data (Any)**:<br />
        Data from pickle file.
    """

    import pickle as pkl

    with open(file_path, "rb") as f:
        data = pkl.load(f)

    return data


def read_fasta(fasta_file: str) -> dict:
    """
    Read a fasta file and return a dictionary of sequences.

    ## Arguments:

    - **fasta_file (str)**:<br />
        Path to fasta file.

    ## Returns:

    - **all_sequences (dict)**:<br />
        `{sequence_header: sequence}`.
    """

    all_sequences = {}
    seq_id = ""

    with open(fasta_file, "r") as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith(">"):
            seq_id = line[1:].strip()
        else:
            seq = line.strip()
            all_sequences[seq_id] = (
                seq
                if seq_id not in all_sequences
                else all_sequences[seq_id] + seq
            )

    return all_sequences