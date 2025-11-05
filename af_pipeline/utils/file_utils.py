import ruamel.yaml
import json
import yaml

class NonAliasingRTRepresenter(ruamel.yaml.representer.RoundTripRepresenter):
    def ignore_aliases(self, data):
        return True

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
    """Update config file with a new field or update an existing field

    Args:
        input_file (str): Path to input config file
        updates (dict, optional): Fields to update in the config file.
            Defaults to None.
        mode (str, optional): Mode to update the config file.
            Defaults to "replace". ("append" or "replace")
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
    """Update the af_offsets in the config file

    Args:
        input_file (str): Path to input config file
        af_offsets (dict): Dictionary of job cycle offsets
        mode (str, optional): Mode to update the config file.
            Defaults to "replace". ("append" or "replace")
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