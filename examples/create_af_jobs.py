import os
from argparse import ArgumentParser
from af_pipeline.af_input.alphafold3 import AlphaFoldServer
from af_pipeline.af_input.alphafold2 import AlphaFold2
from af_pipeline.af_input.colabfold import ColabFold
from af_pipeline.utils.file_utils import read_fasta, read_yaml
from af_pipeline.constants import af_constants

if __name__ == "__main__":

    af_constants.RES_RANGE_SEP = "to"

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
    args.add_argument(
        "-t",
        "--pred_type",
        type=str,
        required=False,
        default="AF3",
        choices=["AF3", "AF2", "ColabFold"],
        help="Type of predictions to create input jobs for. (Default: AF3)",
    )
    args = args.parse_args()

    config_path = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)

    config_dict = read_yaml(config_path)

    protein_sequences = read_fasta(args.protein_sequences)
    nucleic_acid_sequences = (
        read_fasta(args.nucleotide_sequences)
        if args.nucleotide_sequences else None
    )

    if args.pred_type == "AF3":

        af_input = AlphaFoldServer(
            config_dict=config_dict,
            protein_sequences=protein_sequences,
            nucleic_acid_sequences=nucleic_acid_sequences, # optional
        )

        af_input.write_job_files(output_dir=output_dir, num_jobs_per_file=20)

    elif args.pred_type == "AF2":

        af_input = AlphaFold2(
            config_dict=config_dict,
            protein_sequences=protein_sequences,
        )

        af_input.write_job_files(output_dir=output_dir)

    elif args.pred_type == "ColabFold":

        af_input = ColabFold(
            config_dict=config_dict,
            protein_sequences=protein_sequences,
        )

        af_input.write_job_files(output_dir=output_dir)