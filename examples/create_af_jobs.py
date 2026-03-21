import os
import yaml
from argparse import ArgumentParser
from af_pipeline.af_input.alphafold3 import AlphaFoldServer
from af_pipeline.af_input.alphafold2 import AlphaFold2
from af_pipeline.af_input.colabfold import ColabFold
from af_pipeline.utils.misc_utils import add_attribute
from af_pipeline.utils.file_utils import read_fasta, write_json
from pprint import pprint
from af_pipeline.constants.af_constants import (
    ConfigYaml,
    AFInputJobFields,
)

if __name__ == "__main__":

    args = ArgumentParser()

    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/af_input_jobs",
        help="output directory for alphafold input json files",
    )
    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/config.yaml",
        help="input yaml file containing the target proteins",
    )
    args.add_argument(
        "-p",
        "--protein_sequences",
        type=str,
        required=False,
        default="./input/protein_sequences.fasta",
        help="fasta file containing all protein sequences",
    )
    args.add_argument(
        "-n",
        "--nucleotide_sequences",
        type=str,
        required=False,
        default="./input/nucleic_acid_sequences.fasta",
        help="fasta file containing dna or rna sequences",
    )
    args = args.parse_args()

    config_path = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)

    config_yaml = yaml.load(open(config_path), Loader=yaml.FullLoader)
    protein_uniprot_map = config_yaml.get(ConfigYaml.PROTEIN_UNIPROT_MAP, None)
    input_dict = config_yaml.get(ConfigYaml.AF_INPUT_JOBS, None)

    protein_sequences = read_fasta(args.protein_sequences)
    nucleic_acid_sequences = (
        read_fasta(args.nucleotide_sequences)
        if args.nucleotide_sequences else None
    )

    # nucleic_acid_sequences or proteins is not a required argument to AFInput
    # Although, headers in protein sequences should match the entity names in
    # the input yaml file if the proteins are not provided

    af_input = AlphaFoldServer(
        input_dict=input_dict,
        protein_sequences=protein_sequences,
        nucleic_acid_sequences=nucleic_acid_sequences,
        entities_map=protein_uniprot_map,
    )

    job_cycles, job_set_names, af_offsets, cycle_seeds = af_input.create_af3_job_cycles()
    # pprint(job_cycles)
    AlphaFoldServer.write_job_files(
        job_cycles=job_cycles,
        output_dir=output_dir,
        num_jobs_per_file=20,
    )

    # This replaces/adds the job names in the config file
    updated_config = add_attribute(
        config_yaml=config_yaml,
        attribute_name=AFInputJobFields.JOB_SET_NAME,
        attribute_value=job_set_names,
        mode="replace",
        add_first=True,
    )

    # This replaces/adds the af_offsets in the config file
    updated_config = add_attribute(
        config_yaml=updated_config,
        attribute_name=AFInputJobFields.AF_OFFSET,
        attribute_value=af_offsets,
        mode="replace",
    )

    updated_config = add_attribute(
        config_yaml=updated_config,
        attribute_name=AFInputJobFields.MODEL_SEEDS,
        attribute_value=cycle_seeds,
        mode="replace",
    )

    write_json(
        file_path=os.path.join(os.path.dirname(output_dir), "af_input_jobs.json"),
        data=updated_config,
    )

    # For AlphaFold2
    af_input = AlphaFold2(
        input_dict=input_dict,
        protein_sequences=protein_sequences,
        entities_map=protein_uniprot_map,
    )

    job_cycles = af_input.create_af2_job_cycles()
    # pprint(job_cycles)
    AlphaFold2.write_job_files(
        job_cycles=job_cycles,
        output_dir=args.output,
    )

    # For ColabFold
    # af_input = ColabFold(
    #     protein_sequences=protein_sequences,
    #     input_dict=input_dict,
    #     entities_map=protein_uniprot_map,
    # )

    # job_cycles = af_input.create_colabfold_job_cycles()
    # AlphaFold2.write_job_files(
    #     job_cycles=job_cycles,
    #     output_dir=args.output,
    # )