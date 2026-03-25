"""
Rank AlphaFold Predictions
===================================
- Rank AlphaFold predictions based on various metrics and select the best model for each job set.
"""

import os
import warnings
from collections import defaultdict
from af_pipeline.utils.file_utils import read_json
from af_pipeline.utils.misc_utils import chain_id_gen
from af_pipeline.constants import af_constants
from af_pipeline.constants.af_constants import (
    AF3Metrics,
    AF3SummaryConfidenceFields,
    AFInputJobFields,
    BestPredictionFields,
    FileFormat,
)

def get_directory_level(pred_dir:str) -> int | None:
    """ Get the level of the prediction directory

    **level 0** is the job directory -> this will have 5 models\n
    **level 1** is the job set directory -> this will have `n` directories given `n` model seeds\n
    **level 2** is the job cycle directory -> this will have `m` directories given `m` job sets\n
    **level 3** is the prediction directory -> this has all the job cycles

    ## Arguments:

    - **pred_dir (str)**:<br />
        Path to the prediction directory

    ## Returns:

    - **int | None**:<br />
        `int`: Level of the prediction directory
        `None`: If the directory does not exist or is not valid
    """
    if not os.path.isdir(pred_dir):
        print(f"Prediction directory {pred_dir} does not exist.")
        return None

    sub_items = os.listdir(pred_dir)
    sub_dirs = [item for item in sub_items if os.path.isdir(os.path.join(pred_dir, item))]
    sub_files = [item for item in sub_items if os.path.isfile(os.path.join(pred_dir, item))]

    for level in range(4):
        if (
            any(f".{FileFormat.CIF}" in file or f".{FileFormat.PDB}" in file
            for file in sub_files) or len(sub_dirs) == 0
        ):
            return level

        pred_dir = os.path.join(pred_dir, sub_dirs[0])
        sub_items = os.listdir(pred_dir)
        sub_dirs = [item for item in sub_items if os.path.isdir(os.path.join(pred_dir, item))]
        sub_files = [item for item in sub_items if os.path.isfile(os.path.join(pred_dir, item))]

    # return None if no valid level is found
    return None

def is_valid_job_dir(job_dir:str) -> tuple[bool, str]:
    """ Check if the job directory is valid

    A valid job directory is one that is at **level 0** and contains at least one
    prediction file (either .cif or .pdb).

    ## Arguments:

    - **job_dir (str)**:<br />
        Path to the job directory

    ## Returns:

    - **bool**:<br />
        `True` if the job directory is valid, `False` otherwise
    """

    dir_level = get_directory_level(job_dir)

    if dir_level is None or dir_level != 0:
        msg = f"{job_dir} is not at the expected level 0 for a job directory."
        return False, msg

    elif dir_level == 0:
        return True, f"{job_dir} is a valid job directory."

def is_valid_job_set_dir(job_set_dir:str) -> tuple[bool, str, list, list]:
    """ Check if the job set directory is valid

    A valid job set directory is one that is at **level 1** and all its subdirectories
    are valid job directories (i.e., at **level 0** and contains prediction files).

    ## Arguments:

    - **job_set_dir (str)**:<br />
        Path to the job set directory

    ## Returns:

    - **bool**:<br />
        `True` if the job set directory is valid, `False` otherwise
    """

    dir_level = get_directory_level(job_set_dir)

    return_val = True
    valid_job_dirs = []
    invalid_job_dirs = []

    if dir_level is None or dir_level != 1:
        msg = f"{job_set_dir} is not at the expected level 1 for a job set directory."
        return False, msg, valid_job_dirs, invalid_job_dirs

    job_dirs = [
        item for item in os.listdir(job_set_dir)
        if os.path.isdir(os.path.join(job_set_dir, item))
    ]

    if len(job_dirs) == 0:
        msg = f"No job directories found in {job_set_dir}."
        return False, msg, valid_job_dirs, invalid_job_dirs

    for job_dir in job_dirs:

        is_valid_job, _msg = is_valid_job_dir(os.path.join(job_set_dir, job_dir))

        if is_valid_job:
            valid_job_dirs.append(job_dir)
        else:
            invalid_job_dirs.append(job_dir)

        return_val = return_val and is_valid_job

    if return_val:
        msg = f"{job_set_dir} is a valid job set directory."
    else:
        msg = f"{job_set_dir} is not a valid job set directory. "

    return return_val, msg, valid_job_dirs, invalid_job_dirs

def is_valid_job_cycle_dir(job_cycle_dir:str) -> tuple[bool, str, list, list]:
    """ Check if the job cycle directory is valid

    A valid job cycle directory is one that is at **level 2** and all its subdirectories
    are valid job set directories (i.e., at **level 1**).

    ## Arguments:

    - **job_cycle_dir (str)**:<br />
        Path to the job cycle directory

    ## Returns:

    - **bool**:<br />
        `True` if the job cycle directory is valid, `False` otherwise
    """

    dir_level = get_directory_level(job_cycle_dir)

    return_val = True
    valid_job_set_dirs = []
    invalid_job_set_dirs = []

    if dir_level is None or dir_level != 2:
        msg = f"{job_cycle_dir} is not at the expected level 2 for a job cycle directory."
        return False, msg, valid_job_set_dirs, invalid_job_set_dirs

    job_set_dirs = [
        item for item in os.listdir(job_cycle_dir)
        if os.path.isdir(os.path.join(job_cycle_dir, item))
    ]

    if len(job_set_dirs) == 0:
        msg = f"No job set directories found in {job_cycle_dir}."
        return False, msg, valid_job_set_dirs, invalid_job_set_dirs

    for job_set_dir in job_set_dirs:
        is_valid_job_set, _msg, _, _ = is_valid_job_set_dir(
            os.path.join(job_cycle_dir, job_set_dir)
        )
        if is_valid_job_set:
            valid_job_set_dirs.append(job_set_dir)
        else:
            invalid_job_set_dirs.append(job_set_dir)
        return_val = return_val or is_valid_job_set

    if return_val:
        msg = f"{job_cycle_dir} is a valid job cycle directory."
    else:
        msg = f"{job_cycle_dir} is not a valid job cycle directory. "

    return return_val, msg, valid_job_set_dirs, invalid_job_set_dirs

def is_valid_af_master_dir(af_master_dir:str) -> tuple[bool, str, list, list]:
    """ Check if the AF master directory is valid

    A valid AF master directory is one that is at **level 3** and all its subdirectories
    are valid job cycle directories (i.e., at **level 2**).

    ## Arguments:

    - **af_master_dir (str)**:<br />
        Path to the AF master directory

    ## Returns:

    - **bool**:<br />
        `True` if the AF master directory is valid, `False` otherwise
    """

    dir_level = get_directory_level(af_master_dir)

    return_val = True
    valid_job_cycle_dirs = []
    invalid_job_cycle_dirs = []

    if dir_level is None or dir_level != 3:
        msg = f"{af_master_dir} is not at the expected level 3 for an AF master directory."
        return False, msg, valid_job_cycle_dirs, invalid_job_cycle_dirs

    job_cycle_dirs = [
        item for item in os.listdir(af_master_dir)
        if os.path.isdir(os.path.join(af_master_dir, item))
    ]

    if len(job_cycle_dirs) == 0:
        msg = f"No job cycle directories found in {af_master_dir}."
        return False, msg, valid_job_cycle_dirs, invalid_job_cycle_dirs

    for cycle_dir in job_cycle_dirs:
        is_valid_cycle, _msg, _, _ = is_valid_job_cycle_dir(os.path.join(af_master_dir, cycle_dir))
        if is_valid_cycle:
            valid_job_cycle_dirs.append(cycle_dir)
        else:
            invalid_job_cycle_dirs.append(cycle_dir)
        return_val = return_val or is_valid_cycle

    if return_val:
        msg = f"{af_master_dir} is a valid AF master directory."
    else:
        msg = f"{af_master_dir} is not a valid AF master directory."

    return return_val, msg, valid_job_cycle_dirs, invalid_job_cycle_dirs

def get_job_set_dirs(pred_dir:str) -> set:
    """ Get the job set directories from the input prediction directory.

    ## Arguments:

    - **pred_dir (str)**:<br />
        Path to the prediction directory

    ## Returns:

    - **list**:<br />
        List of job set directories found in the prediction directory
    """
    if not os.path.isdir(pred_dir):
        raise ValueError(f"Prediction directory {pred_dir} does not exist.")

    job_set_dirs = set()

    is_valid_af_master, msg1, valid_job_cycle_dirs, _ = is_valid_af_master_dir(pred_dir)

    if is_valid_af_master:

        for job_cycle_dir in valid_job_cycle_dirs:
            job_cycle_path = os.path.join(pred_dir, job_cycle_dir)
            is_valid_job_cycle, _, valid_job_set_dirs, _ = is_valid_job_cycle_dir(job_cycle_path)
            job_set_dirs.update([
                os.path.join(job_cycle_path, job_set_dir)
                for job_set_dir in valid_job_set_dirs
            ])

        return job_set_dirs

    is_valid_job_cycle, _, valid_job_set_dirs, _ = is_valid_job_cycle_dir(pred_dir)

    if is_valid_job_cycle:
        job_set_dirs.update([
            os.path.join(pred_dir, job_set_dir)
            for job_set_dir in valid_job_set_dirs
        ])

        return job_set_dirs

    is_valid_job_set, _, _, _ = is_valid_job_set_dir(pred_dir)

    if is_valid_job_set:
        job_set_dirs.add(pred_dir)
        return job_set_dirs

    else:
        raise ValueError(f"Prediction directory {pred_dir} is not a valid AF \
            master directory, job cycle directory, or job set directory.")

class RankAF3JobSet:
    """ Class to rank AF3 predictions for a given job set directory """

    job_set_dir: str
    """ Path to the job set directory"""

    job_set_name: str
    """ Name of the job set"""

    job_cycle_name: str
    """ Name of the job cycle"""

    job_set_id: int
    """ ID of the job set in the job cycle (1-indexed)"""

    try_af_offset_from_path: bool
    """ Whether to try to extract AF offset from the structure path"""

    def __init__(
        self,
        job_set_dir:str,
        try_af_offset_from_path: bool = False,
    ):
        self.job_set_dir = os.path.abspath(job_set_dir)
        self.job_cycle_name = os.path.basename(os.path.dirname(self.job_set_dir))
        self.job_set_name = os.path.basename(self.job_set_dir)
        self.try_af_offset_from_path = try_af_offset_from_path

    def extract_af3_best_pred_data(self, af_input_jobs: dict| None = None) -> list:
        """ Extract AF3 model metrics and paths to the best model and data for a
        given prediction directory

        ## Arguments:

        - **af_input_jobs (dict | None, optional):**:<br />
            Dictionary containing AF input data

        ## Returns:

        - **list**:<br />
            List of dictionaries containing paths to the best model and
            data for a given prediction directory along with offsets
        """

        best_pred_info = {}

        _best_seed, _ranking_score, structure_path, _best_model_idx = self.rank_seeds()

        if structure_path is None:
            warnings.warn(
                f"No valid AF3 predictions found for {self.job_cycle_name} job {self.job_set_name}. "
                "Skipping ranking of AF3 predictions."
            )
            return []

        data_path = self.get_data_path_from_structure_path(structure_path)

        af_offset = self.extract_af_offset_from_af_input_jobs(af_input_jobs)

        if self.try_af_offset_from_path:
            print("Trying to extract AF offset from structure path...")

            try:
                af_offset = self.extract_af_offset_from_path(structure_path)

            except ValueError as e:
                print(f"Error extracting AF offset from path {structure_path}: {e}")

        key = os.path.basename(os.path.dirname(os.path.dirname(structure_path)))

        best_pred_info[key] = {
            BestPredictionFields.STRUCTURE_PATH: structure_path,
            BestPredictionFields.DATA_PATH: data_path,
            BestPredictionFields.AF_OFFSET: af_offset
        }

        return best_pred_info

    def rank_seeds(self) -> tuple:
        """ Rank the AF3 predictions based on the ranking score

        ## Returns:

        - **tuple**:<br />
            (Model seed, ranking score, model path, model index) corresponding to the best model
        """

        af3_metrics_per_seed = self.get_af3_model_metrics_per_seed()

        if af3_metrics_per_seed is None:
            return None, None, None, None

        ranking = []
        for seed, metrics_list in af3_metrics_per_seed.items():
            for metrics in metrics_list:
                ranking.append(
                    (
                        seed,
                        metrics[AF3Metrics.RANKING_SCORE],
                        metrics[AF3Metrics.IPTM],
                        metrics[AF3Metrics.PTM],
                        metrics[AF3Metrics.FRACTION_DISORDERED],
                        metrics[AF3Metrics.MODEL_PATH],
                        metrics[AF3Metrics.MODEL_IDX],
                    )
                )

        # Sort by ranking_score
        ranking.sort(
            key=lambda item: (
                item[1],  # ranking_score
                item[6]*(-1), # model index
                item[2], # iptm
                item[3],  # ptm
                item[4],  # fraction_disordered
                item[0],  # seed
            ),
            reverse=True,
        )

        if len(ranking) == 0:
            warnings.warn(
                f"No valid AF3 predictions found in {self.job_set_dir}. "
                "Skipping ranking of AF3 predictions."
            )
            return None, None, None, None

        best_model = ranking[0]
        best_model_seed = best_model[0]
        best_ranking_score = best_model[1]
        best_model_path = best_model[5]
        best_model_idx = best_model[6]

        for idx, model in enumerate(ranking):
            print(
                f"Seed: {model[0]}, "
                f"Ranking Score: {model[1]:.2f}, "
                f"iptm: {model[2]}, "
                f"ptm: {model[3]}, "
                f"Fraction Disordered: {model[4]}, "
                f"Model Index: {model[6]}"
            )

        return best_model_seed, best_ranking_score, best_model_path, best_model_idx

    def get_af3_model_metrics_per_seed(self) -> dict | None:
        """ Get AF3 model metrics per seed

        Reads the AF3 job set directory and extracts the model metrics for
        each seed.
        The directory structure is expected to be as follows:
        ```
        job_set_dir/
        └── seed/
            ├── summary_confidences_0.json
            ├── summary_confidences_1.json
            ├── summary_confidences_2.json
            ├── summary_confidences_3.json
            ├── summary_confidences_4.json
            ├── job_request.json
            ├── model_0.cif
            ├── model_1.cif
            ├── model_2.cif
            ├── model_3.cif
            └── model_4.cif
        ```

        Each seed directory will contain 5 summary confidence files, one for each model.

        The function will return a dictionary with the seed as the key and a list of
        model metrics as the value. Each model metric will be a dictionary containing
        the following
        keys:

            - ranking_score
            - iptm
            - ptm
            - fraction_disordered
            - model_path
            - model_idx

        ## Returns:

        - **dict | None**:<br />
            Dictionary containing the best AF3 model metrics per seed
        """

        if not os.path.exists(self.job_set_dir):
            warnings.warn(
                f"Prediction directory {self.job_set_dir} does not exist. "
                "Skipping ranking of AF3 predictions."
            )
            return None

        seed_prediction_dirs = [
            os.path.join(self.job_set_dir, pred)
            for pred in os.listdir(self.job_set_dir)
            if os.path.isdir(os.path.join(self.job_set_dir, pred))
        ]

        if len(seed_prediction_dirs) == 0:
            warnings.warn(
                f"No seed prediction directories found in {self.job_set_dir}. "
                "Skipping ranking of AF3 predictions."
            )
            return None

        af3_metrics_per_seed = defaultdict(list)

        for af3_seed_prediction_dir in seed_prediction_dirs:

            model_metrics = {}

            summary_confidences_files = [
                os.path.join(af3_seed_prediction_dir, af3_file)
                for af3_file in os.listdir(af3_seed_prediction_dir)
                if "summary_confidences" in af3_file
            ]

            job_request_file = [
                os.path.join(af3_seed_prediction_dir, af3_file)
                for af3_file in os.listdir(af3_seed_prediction_dir)
                if "job_request" in af3_file
            ]

            # assert len(summary_confidences_files) == 5, \
            # f"There should be 5 summary_confidences files. \
            #     Found: {len(summary_confidences_files)}"

            if len(summary_confidences_files) != 5:
                warnings.warn(
                    f"There should be 5 summary_confidences files in {af3_seed_prediction_dir}. "
                    f"Found: {len(summary_confidences_files)}. "
                    "Skipping this."
                )
                continue

            # assert len(job_request_file) == 1
            # "There should be 1 job_request file"

            if len(job_request_file) != 1:
                warnings.warn(
                    f"There should be 1 job_request file in {af3_seed_prediction_dir}. "
                    f"Found: {len(job_request_file)}. "
                    "Skipping this."
                )
                continue

            for summary_confidence_path in summary_confidences_files:

                model_idx = summary_confidence_path.split("_")[-1].split(".")[0]
                af3_summary_confidences = self.parse_af3_summary_confidences(
                    af3_summary_confidence_file=summary_confidence_path
                )
                model_metrics[int(model_idx)] = af3_summary_confidences

            # choosing the best model based on ranking score
            model_seed = self.get_seed_from_job_request(job_request_file[0])

            for i, model_metric in model_metrics.items():

                af3_metrics_per_seed[model_seed].append(
                    {
                        AF3Metrics.RANKING_SCORE: model_metric[AF3Metrics.RANKING_SCORE],
                        AF3Metrics.IPTM: model_metric[AF3Metrics.IPTM],
                        AF3Metrics.PTM: model_metric[AF3Metrics.PTM],
                        AF3Metrics.FRACTION_DISORDERED: model_metric[AF3Metrics.FRACTION_DISORDERED],
                        AF3Metrics.MODEL_PATH: [
                            os.path.join(af3_seed_prediction_dir, af3_file)
                            for af3_file in os.listdir(af3_seed_prediction_dir)
                            if f"model_{i}" in af3_file and af3_file.endswith(FileFormat.CIF)
                        ][0],
                        AF3Metrics.MODEL_IDX: i,
                    }
                )

        return af3_metrics_per_seed

    @staticmethod
    def parse_af3_summary_confidences(af3_summary_confidence_file: str) -> dict:
        """ Get AF3 model metrics from summary confidence file\n
        Reads the summary confidence file and extracts the required metrics.\n
        The summary confidence file is expected to be in JSON format with the
        following keys:

            - fraction_disordered
            - has_clash
            - iptm
            - num_recycles
            - ptm
            - ranking_score

        ## Arguments:

        - **af3_summary_confidence_file (str)**:<br />
            Path to AF3 summary confidence file

        ## Returns:

        - **dict**:<br />
            Dictionary containing AF3 model metrics
        """

        metric_data = read_json(af3_summary_confidence_file)

        required_data = {
            AF3SummaryConfidenceFields.FRACTION_DISORDERED: metric_data.get(AF3SummaryConfidenceFields.FRACTION_DISORDERED),
            AF3SummaryConfidenceFields.HAS_CLASH: metric_data.get(AF3SummaryConfidenceFields.HAS_CLASH),
            AF3SummaryConfidenceFields.IPTM: metric_data.get(AF3SummaryConfidenceFields.IPTM),
            AF3SummaryConfidenceFields.NUM_RECYCLES: metric_data.get(AF3SummaryConfidenceFields.NUM_RECYCLES),
            AF3SummaryConfidenceFields.PTM: metric_data.get(AF3SummaryConfidenceFields.PTM),
            AF3SummaryConfidenceFields.RANKING_SCORE: metric_data.get(AF3SummaryConfidenceFields.RANKING_SCORE),
        }

        return required_data

    @staticmethod
    def get_seed_from_job_request(af3_job_request_file: str) -> int:
        """ Get seed from job request file

        ## Arguments:

        - **af3_job_request_file (str)**:<br />
            Path to job request file

        ## Returns:

        - **int**:<br />
            Seed of the model
        """

        job_request_data = read_json(af3_job_request_file)

        return job_request_data[0].get(AFInputJobFields.MODEL_SEEDS)[0]

    @staticmethod
    def get_data_path_from_structure_path(structure_path: str) -> str:
        """ Get the data path from the structure path

        Specific to AF3, where the structure file is in MMCIF format
        and the data file is in JSON format.

        example structure path -\n
            "/path/to/AF3_pred/p1_1_1to100_p2_2_101to200/model_0001.cif"

        The output will be a path to the data file in the format -\n
            "/path/to/AF3_pred/p1_1_1to100_p2_2_101to200/full_data_0001.json"

        ## Arguments:

        - **structure_path (str)**:<br />
            Path to the structure file

        ## Returns:

        - **str**:<br />
            Path to the data file
        """

        return structure_path.replace(
            f".{FileFormat.CIF}", f".{FileFormat.JSON}"
        ).replace("model_", "full_data_")

    def extract_af_offset_from_af_input_jobs(self, af_input_jobs: dict| None = None) -> dict:
        """ Extract the offsets for AF3 predictions from the AF jobs
        dictionary in the config yaml file.

        ## Arguments:

        - **af_input_jobs (dict | None, optional):**:<br />
            Dictionary containing AF input data

        ## Returns:

        - **dict**:<br />
            Dictionary containing AF3 offsets for each chain in the format -
            ```
            {
                "A": [start, end],
                "B": [start, end]
            }
            ```
        """
        af_offset = {}

        self.try_af_offset_from_path = False

        if af_input_jobs is None:
            warnings.warn("No AF jobs found. Skipping extraction of offsets.")
            return af_offset

        # atoz = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        # chain_count = 0

        job_sets = af_input_jobs.get(self.job_cycle_name, None)

        if job_sets is None:
            warnings.warn(
                f"No job sets found for job cycle {self.job_cycle_name} in AF jobs. "
                "Please check the input yaml file."
            )
            self.try_af_offset_from_path = True
            return af_offset

        assert isinstance(self.job_set_id, int), \
            "You haven't added job id yet. Please call add_job_set_id() first."

        if self.job_set_id < 1 or self.job_set_id > len(job_sets):
            self.try_af_offset_from_path = True
            return af_offset

        if AFInputJobFields.AF_OFFSET in job_sets[self.job_set_id-1]:
            # If af_offset is already present in the job set, use it
            af_offset = job_sets[self.job_set_id-1][AFInputJobFields.AF_OFFSET]
            if len(af_offset) == 0:
                warnings.warn(
                    f"No AF offset found for job {self.job_set_id} in AF jobs. "
                    "Please check the input yaml file."
                )
                self.try_af_offset_from_path = True
            else:
                self.try_af_offset_from_path = False
            return af_offset

        else:
            self.try_af_offset_from_path = True

        return af_offset

    @staticmethod
    def extract_af_offset_from_path(structure_path: str) -> dict:
        """ Extract the offset for AF3 prediction from the structure path

        example structure paths\n
            "/path/to/AF3_pred/p1_1_1to100_p2_2_101-200_1234/model_0001.cif"
            "/path/to/AF3_pred/p1_1_1to100_p2_2_101-200/model_0001.cif"

        The directory name is expected to be in the format\n
            p1_copy1_1-100_p2_copy2_101-200_seed

        where `p1`, `p2` are the protein names, `copy1`, `copy2` are the number of copies,
        and `1-100`, `101-200` are the residue ranges.

        The output will be a dictionary with the `chain_id` as key and
        `residue_range` as values. For example -
        ```
        {
            "A": [1, 100],
            "B": [101, 200],
            "C": [101, 200]
        }
        ```

        ## Arguments:

        - **structure_path (str)**:<br />
            Path to the structure file

        ## Returns:

        - **dict**:<br />
            Dictionary containing the offset for each chain in the structure
        """

        dirname = os.path.basename(os.path.dirname(structure_path))
        af_offset = {}

        assert len(dirname.split("_")) >= 3
        "Invalid directory name for AF3 prediction"

        assert len(dirname.split("_")) % 3 < 27
        "Invalid directory name for AF3 prediction, too many chains"

        if len(dirname.split("_")[0:-1]) % 3 == 0:
            seed_is_present = True

        elif len(dirname.split("_")[0:-1]) % 3 == 2:
            seed_is_present = False

        else:
            raise ValueError(
                "Invalid directory name for AF3 prediction, "
                "should be in the format p1_copy1_1to100_p2_copy2_101to200_seed "
                "or p1_copy1_1to100_p2_copy2_101to200"
            )

        if seed_is_present:
            p_c_r = dirname.split("_")[0:-1]  # Exclude the seed

        else:
            p_c_r = dirname.split("_") # Seed is not present in the path

        # atoz = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        chain_count = 0

        chainGen = chain_id_gen()
        for i in range(len(p_c_r) // 3):

            copy_number = p_c_r[i * 3 + 1]
            res_range = p_c_r[i * 3 + 2]

            for _ in range(int(copy_number)):

                chain_id = next(chainGen)
                try:
                    start, end = res_range.split(af_constants.RES_RANGE_SEP)
                except ValueError:
                    try:
                        start, end = res_range.split("-")
                    except ValueError:
                        raise ValueError(
                            f"Invalid residue range format in {structure_path}. "
                            "Expected format is 'STARTtEND' or 'START-END'."
                        )
                af_offset[chain_id] = [int(start), int(end)]
                chain_count += 1

        return af_offset

    def add_job_set_id(self, af_input_jobs: dict | None = None, soft_match = False) -> int:
        """ Get the `job_id` from the AF input jobs dictionary

        ## Arguments:

        - **af_input_jobs (dict | None, optional):**:<br />
            Dictionary containing AF input data

        - **soft_match (bool, optional):**:<br />
            If True, allow partial matching of job set names

        ## Returns:

        - **int**:<br />
            Job id for the current job set
        """
        if af_input_jobs is None:
            warnings.warn("No AF jobs found. Skipping extraction of job id.")
            return -1

        job_sets = af_input_jobs.get(self.job_cycle_name, None)

        if job_sets is None:
            return -1

        # assert isinstance(self.job_set_name, str), "You haven't added job set name yet."

        job_set_id = -1 # in case job set name is not found

        for idx, job in enumerate(job_sets):
            if job.get(AFInputJobFields.JOB_SET_NAME, "") == self.job_set_name:
                job_set_id = idx + 1  # Job IDs are 1-indexed

        if job_set_id == -1 and soft_match:
            for idx, job in enumerate(job_sets):
                if (
                    self.job_set_name in job.get(AFInputJobFields.JOB_SET_NAME, "")
                    or job.get(AFInputJobFields.JOB_SET_NAME, "") in self.job_set_name
                ):
                    job_set_id = idx + 1

        self.job_set_id = job_set_id
        return job_set_id