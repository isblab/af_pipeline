import os
import yaml
import warnings
from pprint import pprint
from argparse import ArgumentParser
from af_pipeline.rank_predictions.rank_af import (
    get_job_set_dirs,
    is_valid_af_master_dir,
    is_valid_cycle_dir,
    is_valid_job_set_dir,
)
from af_pipeline.utils.file_utils import update_config
from af_pipeline.rank_predictions.rank_af import RankAF3JobSet
from af_pipeline.constants.af_constants import ConfigYaml

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

    args = args.parse_args()

    config_yaml = yaml.load(open(args.input), Loader=yaml.FullLoader)

    job_set_dirs = []

    # There are three ways to scan through directories to get job set directories
    # One can choose the way that suits them best
    # If you need to find best predictions for all job sets across cycles, use 2.
    # If you need to find best predictions for a specific cycle, use 3.
    # If you need to find best predictions for a specific job set directory, use 1.

    # 1. Scan through af_job_set_dirs to get all job set directories
    af_job_set_dirs = config_yaml.get(ConfigYaml.job_set, None)

    if af_job_set_dirs is None:
        warnings.warn(
            "No AF job set directories found in the config yaml. "
            "Please provide 'af_job_set_dirs' in the config yaml."
        )

    if isinstance(af_job_set_dirs, list):
        for af_job_set_dir in af_job_set_dirs:
            print(is_valid_job_set_dir(af_job_set_dir))
            job_set_dirs.extend(get_job_set_dirs(af_job_set_dir))

    # 2. Scan through af_cycle_dirs to get all job set directories
    af_master_dirs = config_yaml.get(ConfigYaml.master, None)

    if af_master_dirs is None:
        warnings.warn(
            "No AF master directories found in the config yaml. "
            "Please provide 'af_master_dirs' in the config yaml."
        )

    if isinstance(af_master_dirs, list):
        for af_master_dir in af_master_dirs:
            print(is_valid_af_master_dir(af_master_dir))
            job_set_dirs.extend(get_job_set_dirs(af_master_dir))

    # 3. Scan through af_cycle_dirs to get all job set directories
    af_cycle_dirs = config_yaml.get(ConfigYaml.cycle, None)

    if af_cycle_dirs is None:
        warnings.warn(
            "No AF cycle directories found in the config yaml. "
            "Please provide 'af_cycle_dirs' in the config yaml."
        )

    if isinstance(af_cycle_dirs, list):
        for af_cycle_dir in af_cycle_dirs:
            print(is_valid_cycle_dir(af_cycle_dir))
            job_set_dirs.extend(get_job_set_dirs(af_cycle_dir))

    pprint(job_set_dirs)

    # Get the best prediction for each job set directory
    best_predictions = {}

    for job_set_dir in list(set(job_set_dirs)):
        ranker = RankAF3JobSet(
            job_set_dir=job_set_dir,
        )
        ranker.add_job_set_id(
            af_input_jobs=config_yaml.get(ConfigYaml.input, {}),
            soft_match=True,
        )
        best_pred_info = ranker.extract_af3_best_pred_data(
            af_input_jobs=config_yaml.get(ConfigYaml.input, {})
        )
        if len(best_pred_info) > 0:
            # print("Updating config yaml with ranked predictions")
            # update_config(
            #     input_file=args.input,
            #     updates={f"{ranker.cycle_name}_job_{str(ranker.job_set_name)}": best_pred_info},
            #     mode="replace",
            # )
            best_predictions.update(best_pred_info)

    # Add the best predictions to the config file
    if len(best_predictions) > 0:
        update_config(
            input_file=args.input,
            updates={ConfigYaml.best_pred: best_predictions},
            mode="replace",
        )