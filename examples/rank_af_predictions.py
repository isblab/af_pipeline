import warnings
from af_pipeline.rank_predictions.rank_af import (
    get_job_set_dirs,
    is_valid_af_master_dir,
    is_valid_cycle_dir,
    is_valid_job_set_dir,
)
from argparse import ArgumentParser
import yaml
from utils import update_config
import os
from pprint import pprint
from af_pipeline.rank_predictions.rank_af import RankAF3JobSet

if __name__ == "__main__":

    args = ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/config.yaml",
        help="Path to input yaml file containing the target proteins and their uniprot ids",
    )

    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/ranked_af3_predictions",
        help="Output directory for ranked predictions",
    )

    args.add_argument(
        "-k",
        "--pred_keys",
        nargs="+",
        default=["Lb2Cas12a_RNA_DNA_complex_8I54:1", "Actin_profilin:1,2"],
        help="Keys of AF3 predictions to rank. If not provided, all predictions will be ranked.",
    )

    args = args.parse_args()

    config_yaml = yaml.load(open(args.input), Loader=yaml.FullLoader)

    job_set_dirs = []

    af_job_set_dirs = config_yaml.get("af_job_set_dirs", None)

    if af_job_set_dirs is None:
        warnings.warn(
            "No AF job set directories found in the config yaml. "
            "Please provide 'af_job_set_dirs' in the config yaml."
        )

    if isinstance(af_job_set_dirs, list):
        for af_job_set_dir in af_job_set_dirs:
            print(is_valid_job_set_dir(af_job_set_dir))
            job_set_dirs.extend(get_job_set_dirs(af_job_set_dir))

    # af_master_dirs = config_yaml.get("af_master_dirs", None)

    # if af_master_dirs is None:
    #     warnings.warn(
    #         "No AF master directories found in the config yaml. "
    #         "Please provide 'af_master_dirs' in the config yaml."
    #     )

    # if isinstance(af_master_dirs, list):
    #     for af_master_dir in af_master_dirs:
    #         print(is_valid_af_master_dir(af_master_dir))
    #         job_set_dirs.extend(get_job_set_dirs(af_master_dir))

    pprint(job_set_dirs)

    best_predictions = []

    for job_set_dir in list(set(job_set_dirs)):
        ranker = RankAF3JobSet(
            job_set_dir=job_set_dir,
        )
        ranker.add_job_set_id(af_input_jobs=config_yaml.get("af_input_jobs", {}))
        best_pred_info = ranker.extract_af3_best_pred_data(
            af_input_jobs=config_yaml.get("af_input_jobs", {})
        )

        if len(best_pred_info) > 0:
            # print("Updating config yaml with ranked predictions")
            # update_config(
            #     input_file=args.input,
            #     updates={f"{ranker.cycle_name}_job_{str(ranker.job_set_name)}": best_pred_info},
            #     mode="replace",
            # )
            best_predictions.extend(best_pred_info)

    if len(best_predictions) > 0:
        update_config(
            input_file=args.input,
            updates={"best_af3_predictions": best_predictions},
            mode="replace",
        )