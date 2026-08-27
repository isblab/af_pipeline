import os
import warnings
from pprint import pprint
from argparse import ArgumentParser
from af_pipeline.rank_predictions.rank_af import get_job_set_dirs
from af_pipeline.utils.file_utils import read_json, write_json
from af_pipeline.rank_predictions.rank_af import RankAF3, RankAF2
from af_pipeline.constants.af_constants import ConfigYaml
from af_pipeline.constants import af_constants

if __name__ == "__main__":

    af_constants.RES_RANGE_SEP = "to" # set to match the separator used in create_af_jobs.py for AF3 predictions

    args = ArgumentParser()

    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output",
        help="Output directory for ranked alphafold predictions",
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

    args.add_argument(
        "--pred_dirs",
        nargs="+",
        required=False,
        default=[
            "./input/AF_predictions/AF3",
            # "./input/AF_predictions/ColabFold/Act1_1_1to375_Cdc3_1_11to120"
            # "./input/AF_predictions/AF3/Collagen_triple_helix_targets",
            # "./input/AF_predictions/AF3/Actin_profilin_targets/Act1_1_1-375_Cdc3_1_1-127",
            # "./input/AF_predictions/AF3/Actin_profilin_targets/actin_profilin_complex",
            # "./input/AF_predictions/AF3/Actin_profilin_targets"
        ],
    )

    args = args.parse_args()

    # Get the best prediction for each job set directory
    best_predictions = {}

    job_set_dirs = set()

    # scan through all directories and subdirectories to find job set directories
    for pred_dir in args.pred_dirs:
        if not os.path.exists(pred_dir):
            warnings.warn(f"Prediction directory {pred_dir} does not exist. Skipping...")
            continue
        job_set_dirs.update(
            get_job_set_dirs(pred_dir=pred_dir, pred_type=args.pred_type)
        )

    for job_set_dir in list(set(job_set_dirs)):

        if args.pred_type == "AF3":
            ranker = RankAF3(
                job_set_dir=job_set_dir,
                af_input_jobs=None,
                soft_match=True,
            )
            best_pred_info = ranker.extract_af3_best_pred_data()

        if len(best_pred_info) > 0:
            best_predictions.update(best_pred_info)

    os.makedirs(args.output, exist_ok=True)
    write_json(
        file_path=os.path.join(args.output, "best_af_predictions.json"),
        data=best_predictions,
    )
    print(f"Best predictions for {len(best_predictions)} job sets written to {os.path.join(args.output, 'best_af_predictions.json')}")