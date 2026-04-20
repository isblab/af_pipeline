import os
from argparse import ArgumentParser
from af_pipeline.interaction.interaction import Interaction
from af_pipeline.parser.initialize import Initialize
from af_pipeline.utils.file_utils import read_json
from af_pipeline.constants.af_constants import (
    BestPredictionFields,
    MetricLevel,
    PlotType
)

def get_chain_protein_map(
    protein_chain_map: dict[str, list[str]]
) -> dict[str, str]:

    chain_protein_map = {}
    for protein_name, chain_ids in protein_chain_map.items():
        for chain_id in chain_ids:
            chain_protein_map[chain_id] = protein_name
    return chain_protein_map

if __name__ == "__main__":

    args = ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./output/best_af_predictions.json",
        help="Path to input json file with best AF3 predictions",
    )

    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/interacting_patches",
        help="Path to output directory",
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
        "--idr_chains",
        type=str,
        required=False,
        default="",
        help="Comma separated list of chains to be considered as IDRs",
    )

    args.add_argument(
        "--contact_threshold",
        type=float,
        required=False,
        default=8.0,
        help="Distance cutoff in Angstroms for defining contacts (default: 8.0)",
    )

    args.add_argument(
        "--interaction_pae_cutoff",
        type=float,
        required=False,
        default=10.0,
        help="PAE cutoff for defining interactions (default: 10.0)",
    )

    args.add_argument(
        "--overwrite",
        action='store_true',
        default=False,
        help="Overwrite existing interacting patches output (default: False)",
    )

    args = args.parse_args()
    if not os.path.exists(args.input):
        raise FileNotFoundError(
            "Could not find best_af_predictions.json"
            "Please run `rank_af_predictions.py` to obtain it."
        )
    best_preds = read_json(args.input)
    idr_chains = args.idr_chains.split(",") if args.idr_chains else []
    os.makedirs(args.output, exist_ok=True)

    for pred_, pred_to_analyse in best_preds.items():
        if "monomer" in pred_.lower():
            continue
        structure_path = pred_to_analyse.get(BestPredictionFields.STRUCTURE_PATH)
        data_path = pred_to_analyse.get(BestPredictionFields.DATA_PATH)
        af_offset = pred_to_analyse.get(BestPredictionFields.AF_OFFSET, {})
        entity_chain_map = pred_to_analyse.get(BestPredictionFields.ENTITY_CHAIN_MAP, {})
        pred_dir_name = os.path.basename(
            os.path.dirname(os.path.dirname(structure_path))
        )
        output_dir = os.path.join(args.output, pred_dir_name)

        if os.path.exists(output_dir) and args.overwrite is False:
            print(f"Interacting patches for {pred_dir_name} already exist in {output_dir}. Skipping...")
            continue

        initialize = Initialize(
            data_file_path=data_path,
            structure_file_path=structure_path,
            af_offset=af_offset,
            rep_atom_dict={},
            average_token_pae=False,
            average_token_plddt=False,
            metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
            use_fast_cif_parser=False,
        )

        af_interaction = Interaction(
            contact_threshold=args.contact_threshold,
            plddt_cutoff=args.plddt_cutoff,
            pae_cutoff=args.interaction_pae_cutoff,
            plddt_cutoff_idr=args.plddt_cutoff_idr,
            idr_chains=idr_chains,
            save_plot=False,
            save_table=True,
            setup_instance=initialize,
        )

        regions_of_interest_ = af_interaction.create_regions_of_interest()

        for region_of_interest in regions_of_interest_:

            chain_ids = list(region_of_interest.keys())
            af_interaction.save_ppair_interaction(
                region_of_interest=region_of_interest,
                output_dir=output_dir,
                save_plot=True,
                plot_type=PlotType.STATIC,
                p1_name=entity_chain_map.get(chain_ids[0], chain_ids[0]),
                p2_name=entity_chain_map.get(chain_ids[1], chain_ids[1]),
                concat_residues=True,
                contact_probability=True,
            )
        print("-------------------------------")