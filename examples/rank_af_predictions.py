import os
import yaml
import warnings
from pprint import pprint
from argparse import ArgumentParser
from af_pipeline.rank_predictions.rank_af import get_job_set_dirs
from af_pipeline.utils.file_utils import read_json, write_json
from af_pipeline.rank_predictions.rank_af import RankAF3JobSet
from af_pipeline.constants.af_constants import ConfigYaml
from af_pipeline.constants import af_constants

if __name__ == "__main__":

    af_constants.RES_RANGE_SEP = "to"

    args = ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./output/af_input_jobs.json",
        help="Path to input yaml file containing the target proteins and their uniprot ids",
    )

    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output",
        help="Output directory for ranked alphafold predictions",
    )

    args.add_argument(
        "--pred_dirs",
        nargs="+",
        required=False,
        default=[
            "./input/AF_predictions/AF3",
            # "./input/AF_predictions/AF3/Collagen_triple_helix_targets",
            # "./input/AF_predictions/AF3/Actin_profilin_targets/Act1_1_1-375_Cdc3_1_1-127",
            # "./input/AF_predictions/AF3/Actin_profilin_targets/actin_profilin_complex",
            # "./input/AF_predictions/AF3/Actin_profilin_targets"
        ],
    )

    args = args.parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(
            "Could not find af_input_jobs.json"
            "Please run `create_af_jobs.py` to obtain it."
        )

    config_dict = read_json(args.input)

    job_set_dirs = set()

    # scan through all directories and subdirectories to find job set directories
    for pred_dir in args.pred_dirs:
        print(f"Scanning directory: {pred_dir}")
        if not os.path.exists(pred_dir):
            warnings.warn(f"Prediction directory {pred_dir} does not exist. Skipping...")
            continue
        job_set_dirs.update(get_job_set_dirs(pred_dir))

    pprint(job_set_dirs)

    # Get the best prediction for each job set directory
    best_predictions = {}

    for job_set_dir in list(set(job_set_dirs)):
        ranker = RankAF3JobSet(
            job_set_dir=job_set_dir,
            try_af_offset_from_path=False,
        )
        ranker.add_job_set_id(
            af_input_jobs=config_dict.get(ConfigYaml.AF_INPUT_JOBS, {}),
            soft_match=True,
        )
        best_pred_info = ranker.extract_af3_best_pred_data(
            af_input_jobs=config_dict.get(ConfigYaml.AF_INPUT_JOBS, {})
        )
        if len(best_pred_info) > 0:
            # print("Updating config yaml with ranked predictions")
            # update_config(
            #     input_file=args.input,
            #     updates={f"{ranker.cycle_name}_job_{str(ranker.job_set_name)}": best_pred_info},
            #     mode="replace",
            # )
            best_predictions.update(best_pred_info)

    os.makedirs(args.output, exist_ok=True)
    write_json(
        file_path=os.path.join(args.output, "best_af_predictions.json"),
        data=best_predictions,
    )

    # # Add the best predictions to the config file
    # if len(best_predictions) > 0:
    #     update_config(
    #         input_file=args.input,
    #         updates={ConfigYaml.best_pred: best_predictions},
    #         mode="replace",
    #     )