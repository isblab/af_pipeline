import json

def write_json(
    file_path: str,
    data
):
    """Write data to a json file

    Args:

        file_path (str):
            Path to json file

        data (dict):
            Data to write
    """

    with open(file_path, "w") as f:
        json.dump(data, f)


def read_json(file_path: str):
    """Load a json file

    Args:

        file_path (str):
            Path to json file

    Returns:

        data:
            Data from json file
    """

    with open(file_path, "r") as f:
        data = json.load(f)

    return data


def read_fasta(fasta_file: str) -> dict:
    """
    Read a fasta file and return a dictionary of sequences

    Args:

        fasta_file (str):
            Path to fasta file

    Returns:

        all_sequences (dict):
            Dictionary in the format {sequence_header: sequence}
    """

    all_sequences = {}

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
