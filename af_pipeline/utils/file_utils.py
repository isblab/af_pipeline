import ruamel.yaml
import json
import yaml
from typing import Any

class NonAliasingRTRepresenter(ruamel.yaml.representer.RoundTripRepresenter):
    def ignore_aliases(self, data):
        return True

def write_json(
    file_path: str,
    data
):
    """Write data to a json file.

    Arguments:

    - **file_path (str)**:<br />
        Path to json file.

    - **data (dict)**:<br />
        Data to write.
    """

    with open(file_path, "w") as f:
        json.dump(data, f)


def read_json(file_path: str) -> dict | list:
    """Load a json file.

    Arguments:

    - **file_path (str)**:<br />
        Path to json file.

    Returns:

    - **data (dict|list)**:<br />
        Data from json file.
    """

    with open(file_path, "r") as f:
        data = json.load(f)

    return data

def read_pkl(file_path: str) -> Any:
    """Load a pickle file.

    Arguments:

    - **file_path (str)**:<br />
        Path to pickle file.

    Returns:

    - **data (any)**:<br />
        Data from pickle file.
    """

    import pickle as pkl

    with open(file_path, "rb") as f:
        data = pkl.load(f)

    return data


def read_fasta(fasta_file: str) -> dict:
    """
    Read a fasta file and return a dictionary of sequences.

    Arguments:

    - **fasta_file (str)**:<br />
        Path to fasta file.

    Returns:

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

def update_config(
    input_file: str,
    updates: dict = None,
    mode: str = "replace",
):
    """Update config file with a new field or update an existing field.

    Arguments:

    - **input_file (str)**:<br />
        Path to input config file.

    - **updates (dict, optional)**:<br />
        Fields to update in the config file.

    - **mode (str, optional)**:<br />
        Mode to update the config file. ("append" or "replace").

    """

    yaml = ruamel.yaml.YAML()
    yaml.Representer = NonAliasingRTRepresenter

    update_fields = list(updates.keys()) if updates else []

    if len(update_fields) == 0:

        print("No fields to update in config")
        return None

    yaml.preserve_quotes = True

    with open(input_file, "r") as f:
        config_yaml = yaml.load(f)

    existing_fields = list(config_yaml.keys())

    for field in update_fields:

        add_field = False

        if field in existing_fields:
            if mode == "replace":
                config_yaml[field] = updates[field]
            elif mode == "soft_replace":
                #! only update if the field is not already set
                if config_yaml[field] is None or config_yaml[field] == "":
                    config_yaml[field] = updates[field]
            elif mode == "append":
                #! need to change this, not working as expected
                config_yaml[field].update(updates[field])
            else:
                raise ValueError("Invalid mode. Use 'replace' or 'append")

        else:
            print(f"{field} not found in config")
            print("Adding field to config")
            add_field = True

        if add_field:
            config_yaml[field] = updates[field]
            add_field = False

    with open(input_file, "w") as f:
        yaml.dump(config_yaml, f)

    print(f"Config file updated with {update_fields}")

def update_job_names_in_config(
    input_file: str,
    job_set_names: dict,
    mode: str = "replace",
):
    """ Update the job names in the config file.

    Arguments:

    - **input_file (str)**:<br />
        Path to input config file.

    - **job_set_names (dict)**:<br />
        Dictionary of job cycle names.

    - **mode (str, optional)**:<br />
        Mode to update the config file.
    """

    config_yaml = yaml.load(open(input_file), Loader=yaml.FullLoader)

    af_input_jobs = config_yaml.get("af_input_jobs", {})

    for job_cycle, job_set_names in job_set_names.items():
        if job_cycle in af_input_jobs:
            for idx, job_set in enumerate(af_input_jobs[job_cycle]):
                if "name" in job_set:
                    af_input_jobs[job_cycle][idx]["name"] = job_set_names[idx]
                else:
                    af_input_jobs[job_cycle][idx] = {
                        "name": job_set_names[idx],
                        **job_set,
                    }
        else:
            print(f"{job_cycle} not found in af_input_jobs")

    update_config(
        input_file=input_file,
        updates={"af_input_jobs": af_input_jobs},
        mode=mode
    )

def update_af_offsets_in_config(
    input_file: str,
    af_offsets: dict,
    mode: str = "replace",
):
    """Update the `af_offsets` in the config file.

    Arguments:

    - **input_file (str)**:<br />
        Path to input config file.

    - **af_offsets (dict)**:<br />
        Dictionary of job cycle offsets.

    - **mode (str, optional)**:<br />
        Mode to update the config file.
    """

    config_yaml = yaml.load(open(input_file), Loader=yaml.FullLoader)

    af_input_jobs = config_yaml.get("af_input_jobs", {})

    for job_cycle, job_set_offsets in af_offsets.items():
        if job_cycle in af_input_jobs:
            for idx, job_set in enumerate(af_input_jobs[job_cycle]):
                if "af_offset" in job_set:
                    af_input_jobs[job_cycle][idx]["af_offset"] = job_set_offsets[idx]
                else:
                    af_input_jobs[job_cycle][idx] = {
                        **job_set,
                        "af_offset": job_set_offsets[idx],
                    }
        else:
            print(f"{job_cycle} not found in af_input_jobs")

    update_config(
        input_file=input_file,
        updates={"af_input_jobs": af_input_jobs},
        mode=mode
    )