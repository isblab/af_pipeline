import os
import yaml
from argparse import ArgumentParser
from af_pipeline.initialize import Initialize
from af_pipeline.rigid_bodies.rigid_bodies import RigidBodies
from af_pipeline.constants.af_constants import ConfigYaml
from af_pipeline.utils.file_utils import read_json

if __name__ == "__main__":

    args = ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/config.json",
        help="Path to input json file",
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
        "--min_res",
        type=int,
        required=False,
        default=1,
        help="Minimum number of residues in a domain",
    )

    args.add_argument(
        "--min_proteins",
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
        help="List of protein:chain pairs to be considered for rigid body extraction. Format: protein1:chain1 protein2:chain2",
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

    assert args.min_res > 0, \
        f"Invalid value for min_res {args.min_res}"

    assert args.min_proteins > 0, \
        f"Invalid value for min_proteins {args.min_proteins}"

    # config_yaml = yaml.load(open(args.input), Loader=yaml.FullLoader)
    config_yaml = read_json(args.input)
    input_dict = config_yaml.get(ConfigYaml.best_pred, None)
    idr_chains = args.idr_chains.split(",") if args.idr_chains else []
    protein_chain_map = {
        protein: chain for pair in args.protein_chain_map for protein, chain in [pair.split(":")]
    }

    for pred_head, pred_to_analyse in input_dict.items():

        af_offset = pred_to_analyse.get("af_offset", None)
        structure_path = pred_to_analyse.get("structure_path", None)
        data_path = pred_to_analyse.get("data_path", None)

        initialize = Initialize(
            data_file_path=data_path,
            structure_file_path=structure_path,
            af_offset=af_offset,
            rep_atom_dict={},
            average_token_pae=False,
            average_token_plddt=False,
            metric_level="per_token",
            use_fast_cif_parser=False,
        )

        rigid_bodies_extractor = RigidBodies(
            library=args.library,
            pae_cutoff=args.pae_cutoff,
            pae_power=args.pae_power,
            resolution=args.resolution,
            plddt_cutoff=args.plddt_cutoff,
            plddt_cutoff_idr=args.plddt_cutoff_idr,
            idr_chains=idr_chains,
        )

        rigid_bodies_extractor.set_attributes_from(
            instance=initialize,
        )

        domains = rigid_bodies_extractor.extract_rigid_bodies(
            pae_matrix=rigid_bodies_extractor.pae,
            min_res=args.min_res,
            min_proteins=args.min_proteins,
            plddt_filter=args.apply_plddt_filter,
        )

        rigid_bodies_extractor.save_rigid_bodies(
            domains=domains,
            output_dir=args.output,
            rb_out_fmt="txt",
            save_structure=True,
            rb_struct_fmt="pdb",
            filter_struct_by_plddt=True,
            protein_chain_map=protein_chain_map
        )

        rigid_bodies_extractor.assess_rigid_bodies(
            domains=domains,
            output_dir=args.output,
            protein_chain_map=protein_chain_map,
            symmetric_pae=True,
            as_average=True,
        )

        rigid_bodies_extractor.show_rigid_bodies_on_pae_matrix(
            domains=domains,
            output_dir=args.output,
        )

        print(f"Rigid bodies for {os.path.basename(structure_path)} saved in {args.output}")