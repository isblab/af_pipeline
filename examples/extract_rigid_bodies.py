import os
import yaml
from argparse import ArgumentParser
from af_pipeline.rigid_bodies.rigid_bodies import RigidBodies

if __name__ == "__main__":

    args = ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/config.yaml",
        help="Path to input yaml file",
    )

    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/rigid_bodies",
        help="Output directory for rigid body extraction results",
    )

    args.add_argument(
        "--plddt_cutoff",
        type=float,
        required=False,
        default=70,
        help="PLDDT cutoff for confident residues",
    )

    args.add_argument(
        "--plddt_cutoff_idr",
        type=float,
        required=False,
        default=50,
        help="PLDDT cutoff for IDRs",
    )

    args.add_argument(
        "--pae_cutoff",
        type=float,
        required=False,
        default=12,
        help="PAE cutoff for confident residues",
    )

    args.add_argument(
        "--pae_power",
        type=float,
        required=False,
        default=1,
        help="Power to raise PAE to",
    )

    args.add_argument(
        "--resolution",
        type=float,
        required=False,
        default=0.5,
        help="Resolution of the PAE matrix",
    )

    args.add_argument(
        "--library",
        type=str,
        required=False,
        default="igraph",
        help="Library to use for domain prediction",
    )

    args.add_argument(
        "--num_res",
        type=int,
        required=False,
        default=1,
        help="Minimum number of residues in a domain",
    )

    args.add_argument(
        "--num_proteins",
        type=int,
        required=False,
        default=1,
        help="Minimum number of proteins in a domain",
    )

    args.add_argument(
        "--apply_plddt_filter",
        action="store_true",
        required=False,
        help="Apply PLDDT filter to the residues",
    )

    args.add_argument(
        "--idr_chains",
        type=str,
        required=False,
        default="",
        help="Comma separated list of chains to be considered as IDRs",
    )

    args.add_argument(
        "--protein_chain_map",
        nargs="+",
        required=False,
        default=[],
        help="List of protein:chain pairs to be considered for rigid body extraction. Format: protein1:chain1,protein2:chain2,...",
    )

    args = args.parse_args()

    assert args.plddt_cutoff >= 0 and args.plddt_cutoff <= 100, \
        f"Invalid calue for PLDDT cutoff {args.plddt_cutoff}"

    assert args.pae_cutoff >= 0, \
        f"Invalid value for PAE cutoff {args.pae_cutoff}"

    assert args.pae_power >= 0, \
        f"Invalid value for PAE power {args.pae_power}"

    assert args.resolution > 0, \
        f"Invalid value for resolution {args.resolution}"

    assert args.library in ["igraph", "networkx", "label_propagation"], \
        f"Invalid library {args.library}"

    assert args.num_res > 0, \
        f"Invalid value for num_res {args.num_res}"

    assert args.num_proteins > 0, \
        f"Invalid value for num_proteins {args.num_proteins}"

    config_yaml = yaml.load(open(args.input), Loader=yaml.FullLoader)
    input_dict = config_yaml.get("best_af3_predictions", None)
    idr_chains = args.idr_chains.split(",") if args.idr_chains else []

    for _, pred_to_analyse in enumerate(input_dict):

        af_offset = pred_to_analyse.get("af_offset", None)
        structure_path = pred_to_analyse.get("structure_path", None)
        data_path = pred_to_analyse.get("data_path", None)

        rigid_bodies_extractor = RigidBodies(
            data_file_path=data_path,
            structure_file_path=structure_path,
            af_offset=af_offset,
            idr_chains=idr_chains,
            rep_atom_dict={},
            average_token_pae=False,
            average_token_plddt=False,
            metric_level="per_token",
        )

        rigid_bodies_extractor.plddt_cutoff = args.plddt_cutoff
        rigid_bodies_extractor.plddt_cutoff_idr = args.plddt_cutoff_idr
        rigid_bodies_extractor.pae_cutoff = args.pae_cutoff
        rigid_bodies_extractor.pae_power = args.pae_power
        rigid_bodies_extractor.resolution = args.resolution
        rigid_bodies_extractor.library = args.library

        domains = rigid_bodies_extractor.extract_rigid_bodies(
            num_res=args.num_res,
            num_proteins=args.num_proteins,
            plddt_filter=args.apply_plddt_filter,
        )

        rigid_bodies_extractor.save_rigid_bodies(
            domains=domains,
            output_dir=args.output,
            output_format="txt",
            save_structure=True,
            structure_file_type="pdb",
            no_plddt_filter_for_structure=False,
            pae_plot=False,
            rb_assessment={
                "as_average": True,
                "symmetric_pae": True,
            },
            protein_chain_map=None
        )

        print(f"Rigid bodies for {os.path.basename(structure_path)} saved in {args.output}")