import yaml
from argparse import ArgumentParser
from af_pipeline.af_input.alphafold3 import AlphaFoldServer
from af_pipeline.af_input.alphafold2 import AlphaFold2
from af_pipeline.af_input.colabfold import ColabFold
from af_pipeline.utils.file_utils import (
    read_fasta,
    update_config,
    update_job_names_in_config,
    update_af_offsets_in_config
)
from pprint import pprint

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
        "-t",
        "--pred_type",
        type=str,
        required=False,
        default="AF3",
        help="prediction type (AF2/AF3/ColabFold)",
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

    config_yaml = yaml.load(open(args.input), Loader=yaml.FullLoader)
    protein_uniprot_map = config_yaml.get("proteins", None)
    input_yml = config_yaml.get("af_input_jobs", None)

    protein_sequences = read_fasta(args.protein_sequences)
    nucleic_acid_sequences = (
        read_fasta(args.nucleotide_sequences)
        if args.nucleotide_sequences else None
    )

    # nucleic_acid_sequences or proteins is not a required argument to AFInput
    # Although, headers in protein sequences should match the entity names in
    # the input yaml file if the proteins are not provided

    af_input = AlphaFoldServer(
        input_dict=input_yml,
        protein_sequences=protein_sequences,
        nucleic_acid_sequences=nucleic_acid_sequences,
        entities_map=protein_uniprot_map,
    )

    job_cycles, job_set_names, af_offsets = af_input.create_af3_job_cycles()
    # pprint(job_cycles)
    af_input.write_job_files(
        job_cycles=job_cycles,
        output_dir=args.output,
        num_jobs_per_file=20,
    )

    # This replaces/adds the job names in the config file
    # update_job_names_in_config(
    #     input_file=args.input,
    #     job_set_names=job_set_names,
    #     mode="replace",
    # )

    # This replaces/adds the af_offsets in the config file
    # update_af_offsets_in_config(
    #     input_file=args.input,
    #     af_offsets=af_offsets,
    #     mode="replace",
    # )

    # For AlphaFold2
    af_input = AlphaFold2(
        input_yml=input_yml,
        protein_sequences=protein_sequences,
        entities_map=protein_uniprot_map,
    )

    job_cycles = af_input.create_af2_job_cycles()
    # pprint(job_cycles)
    af_input.write_job_files(
        job_cycles=job_cycles,
        output_dir=args.output,
    )

    # For ColabFold
    # af_input = ColabFold(
    #     protein_sequences=protein_sequences,
    #     input_yml=input_yml,
    #     entities_map=protein_uniprot_map,
    # )

    # job_cycles = af_input.create_colabfold_job_cycles()
    # af_input.write_job_files(
    #     job_cycles=job_cycles,
    #     output_dir=args.output,
    # )