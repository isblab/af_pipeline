import json
import os
from af_pipeline.tools.misc_tools import get_key_from_res_range

def save_rigid_bodies_txt(
    output_dir: str,
    domains: list,
    protein_chain_map: dict,
    file_name: str = "rigid_bodies",
):
    """ Save rigid bodies to a text file.

    This function writes the rigid bodies information to a text file in a human-readable format.

    The output file will contain the rigid body index, chain ID, protein name (if available), and the residue range.

    Args:

        output_dir (str):
            Directory where the output file will be saved.

        domains (list):
            List of dictionaries, where each dictionary represents a rigid body.

        protein_chain_map (dict):
            A mapping of chain IDs to protein names.

        file_name (str, optional):
            Name of the output file without extension.
            Defaults to "rigid_bodies".
    """

    file_name += ".txt"
    output_path = os.path.join(output_dir, file_name)

    with open(output_path, "w") as f:

        for idx, rb_dict in enumerate(domains):
            f.write(f"Rigid Body {idx}\n")

            for chain_id, res_list in rb_dict.items():

                protein_name = protein_chain_map.get(chain_id, None)

                if len(res_list) > 0:
                    if protein_name:
                        f.write(
                            f"{protein_name}_{chain_id}: {get_key_from_res_range(res_range=res_list)}\n"
                        )
                    else:
                        f.write(
                            f"{chain_id}:{get_key_from_res_range(res_range=res_list)}\n"
                        )

            f.write("\n")

def save_rigid_bodies_json(
    output_dir: str,
    domains: list,
    protein_chain_map: dict,
    file_name: str = "rigid_bodies",
):
    """ Save rigid bodies to a JSON file.

    This function writes the rigid bodies information to a JSON file.

    For per-atom tokens JSON format is recommended over text format.

    Args:

        output_dir (str):
            Directory where the output file will be saved.

        domains (list):
            List of dictionaries, where each dictionary represents a rigid body.

        protein_chain_map (dict):
            A mapping of chain IDs to protein names.

        file_name (str, optional):
            Name of the output file without extension.
            Defaults to "rigid_bodies".
    """

    file_name += ".json"
    output_path = os.path.join(output_dir, file_name)

    rigid_bodies = []
    for idx, rb_dict in enumerate(domains):
        ch_dict = {}
        for chain_id, res_num_list in rb_dict.items():
            protein_name = protein_chain_map.get(chain_id, None)
            if not protein_name:
                protein_name = "Unknown"
            ch_dict[chain_id] = {
                "protein": protein_name,
                "residues": []
            }
            for atom_name, res_num in res_num_list:
                ch_dict[chain_id]["residues"].append(
                    (atom_name, res_num)
                )
        rigid_bodies.append(ch_dict)

    with open(output_path, "w") as f:
        json.dump(rigid_bodies, f, indent=4)