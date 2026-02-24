"""
Rank AlphaFold Predictions
===================================
- This module provides functionality to rank AF3 predictions based on various metrics
- Currently supports ranking of AlphaFold 3 (AF3) predictions from AlphaFold server.
- Note the following directory structure for AF3 predictions:
```
master_directory
    └── cycle/
        └── job_set/
            └── model_seed/
                ├── summary_confidences_0.json
                ├── summary_confidences_1.json
                ├── summary_confidences_2.json
                ├── summary_confidences_3.json
                ├── summary_confidences_4.json
                ├── job_request.json
                ├── full_data_0.json
                ├── full_data_1.json
                ├── full_data_2.json
                ├── full_data_3.json
                ├── full_data_4.json
                ├── model_0.cif
                ├── model_1.cif
                ├── model_2.cif
                ├── model_3.cif
                └── model_4.cif
```
- The best way to use this module is to provide the master directory containing all cycles.
- To rank a subset of predictions, one can also provide a specific cycle directory or job set directory.
"""

import os
import warnings
from collections import defaultdict
from af_pipeline.utils.file_utils import read_json
from af_pipeline.utils.misc_utils import chain_id_gen
from af_pipeline.constants.af_constants import RES_RANGE_SEP

def get_directory_level(pred_dir:str) -> int | None:
    """ Get the level of the prediction directory

    **level 0** is the job directory -> this will have 5 models\n
    **level 1** is the job set directory -> this will have `n` directories given `n` model seeds\n
    **level 2** is the cycle directory -> this will have `m` directories given `m` job sets\n
    **level 3** is the prediction directory -> this has all the cycles

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
        if any(".cif" in file or ".pdb" in file for file in sub_files) or len(sub_dirs) == 0:
            return level

        pred_dir = os.path.join(pred_dir, sub_dirs[0])
        sub_items = os.listdir(pred_dir)
        sub_dirs = [item for item in sub_items if os.path.isdir(os.path.join(pred_dir, item))]
        sub_files = [item for item in sub_items if os.path.isfile(os.path.join(pred_dir, item))]

    # return None if no valid level is found
    return None

def is_valid_job_dir(job_dir:str) -> bool:
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
        print(f"Job directory {job_dir} is not at the expected level 0.")
        return False

    elif dir_level == 0:
        return True

def is_valid_job_set_dir(job_set_dir:str) -> bool:
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

    if dir_level is None or dir_level != 1:
        print(f"Job set directory {job_set_dir} is not at the expected level 1.")
        return False

    job_dirs = [
        item for item in os.listdir(job_set_dir)
        if os.path.isdir(os.path.join(job_set_dir, item))
    ]

    if len(job_dirs) == 0:
        print(f"No job directories found in {job_set_dir}.")
        return False

    for job_dir in job_dirs:
        if not is_valid_job_dir(os.path.join(job_set_dir, job_dir)):
            return False

    return True

def is_valid_cycle_dir(cycle_dir:str) -> bool:
    """ Check if the cycle directory is valid

    A valid cycle directory is one that is at **level 2** and all its subdirectories
    are valid job set directories (i.e., at **level 1**).

    ## Arguments:

    - **cycle_dir (str)**:<br />
        Path to the cycle directory

    ## Returns:

    - **bool**:<br />
        `True` if the cycle directory is valid, `False` otherwise
    """

    dir_level = get_directory_level(cycle_dir)

    if dir_level is None or dir_level != 2:
        print(f"Cycle directory {cycle_dir} is not at the expected level 2.")
        return False

    job_set_dirs = [
        item for item in os.listdir(cycle_dir)
        if os.path.isdir(os.path.join(cycle_dir, item))
    ]

    if len(job_set_dirs) == 0:
        print(f"No job set directories found in {cycle_dir}.")
        return False

    for job_set_dir in job_set_dirs:
        if not is_valid_job_set_dir(os.path.join(cycle_dir, job_set_dir)):
            return False

    return True

def is_valid_af_master_dir(af_master_dir:str) -> bool:
    """ Check if the AF master directory is valid

    A valid AF master directory is one that is at **level 3** and all its subdirectories
    are valid cycle directories (i.e., at **level 2**).

    ## Arguments:

    - **af_master_dir (str)**:<br />
        Path to the AF master directory

    ## Returns:

    - **bool**:<br />
        `True` if the AF master directory is valid, `False` otherwise
    """

    dir_level = get_directory_level(af_master_dir)

    if dir_level is None or dir_level != 3:
        print(f"AF master directory {af_master_dir} is not at the expected level 3.")
        return False

    cycle_dirs = [
        item for item in os.listdir(af_master_dir)
        if os.path.isdir(os.path.join(af_master_dir, item))
    ]

    if len(cycle_dirs) == 0:
        print(f"No cycle directories found in {af_master_dir}.")
        return False

    for cycle_dir in cycle_dirs:

        if not is_valid_cycle_dir(os.path.join(af_master_dir, cycle_dir)):
            return False

    print(f"AF master directory {af_master_dir} is valid.")
    return True

def get_job_set_dirs(pred_dir:str) -> list:
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

    job_set_dirs = []

    if is_valid_af_master_dir(pred_dir):
        cycle_dirs = [
            item for item in os.listdir(pred_dir)
            if os.path.isdir(os.path.join(pred_dir, item))
        ]

        for cycle_dir in cycle_dirs:
            cycle_path = os.path.join(pred_dir, cycle_dir)
            job_set_dirs.extend([
                os.path.join(cycle_path, job_set_dir)
                for job_set_dir in os.listdir(cycle_path)
                if os.path.isdir(os.path.join(cycle_path, job_set_dir))
            ])

    elif is_valid_cycle_dir(pred_dir):
        job_set_dirs = [
            os.path.join(pred_dir, job_set_dir)
            for job_set_dir in os.listdir(pred_dir)
            if os.path.isdir(os.path.join(pred_dir, job_set_dir))
        ]

    elif is_valid_job_set_dir(pred_dir):
        job_set_dirs = [pred_dir]

    else:
        raise ValueError(f"Prediction directory {pred_dir} is not a valid AF \
            master directory, cycle directory, or job set directory.")

    return job_set_dirs


class RankAF3JobSet:
    """ Class to rank AF3 predictions for a given job set directory """

    job_set_dir: str
    """ Path to the job set directory"""

    job_set_name: str
    """ Name of the job set"""

    cycle_name: str
    """ Name of the cycle"""

    job_set_id: int
    """ ID of the job set in the cycle (1-indexed)"""

    try_af_offset_from_path: bool
    """ Whether to try to extract AF offset from the structure path"""

    def __init__(
        self,
        job_set_dir:str,
        try_af_offset_from_path: bool = False,
    ):
        self.job_set_dir = os.path.abspath(job_set_dir)
        self.cycle_name = os.path.basename(os.path.dirname(self.job_set_dir))
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
                f"No valid AF3 predictions found for {self.cycle_name} job {self.job_set_name}. "
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
            "structure_path": structure_path,
            "data_path": data_path,
            "af_offset": af_offset
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
                        metrics["ranking_score"],
                        metrics["iptm"],
                        metrics["ptm"],
                        metrics["fraction_disordered"],
                        metrics["model_path"],
                        metrics["model_idx"],
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
                        "ranking_score": model_metric["ranking_score"],
                        "iptm": model_metric["iptm"],
                        "ptm": model_metric["ptm"],
                        "fraction_disordered": model_metric["fraction_disordered"],
                        "model_path": [
                            os.path.join(af3_seed_prediction_dir, af3_file)
                            for af3_file in os.listdir(af3_seed_prediction_dir)
                            if f"model_{i}" in af3_file and af3_file.endswith("cif")
                        ][0],
                        "model_idx": i,
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
            "fraction_disordered": metric_data.get("fraction_disordered"),
            "has_clash": metric_data.get("has_clash"),
            "iptm": metric_data.get("iptm"),
            "num_recycles": metric_data.get("num_recycles"),
            "ptm": metric_data.get("ptm"),
            "ranking_score": metric_data.get("ranking_score"),
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

        return job_request_data[0].get("modelSeeds")[0]

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

        return structure_path.replace(".cif", ".json").replace("model_", "full_data_")

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

        job_sets = af_input_jobs.get(self.cycle_name, None)

        if job_sets is None:
            warnings.warn(
                f"No job sets found for cycle {self.cycle_name} in AF jobs. "
                "Please check the input yaml file."
            )
            self.try_af_offset_from_path = True
            return af_offset

        assert isinstance(self.job_set_id, int), \
            "You haven't added job id yet. Please call add_job_set_id() first."

        if self.job_set_id < 1 or self.job_set_id > len(job_sets):
            self.try_af_offset_from_path = True
            return af_offset

        if "af_offset" in job_sets[self.job_set_id-1]:
            # If af_offset is already present in the job set, use it
            af_offset = job_sets[self.job_set_id-1]["af_offset"]
            if len(af_offset) == 0:
                warnings.warn(
                    f"No AF offset found for job {self.job_set_id} in AF jobs. "
                    "Please check the input yaml file."
                )
                self.try_af_offset_from_path = True
            else:
                self.try_af_offset_from_path = False
            return af_offset

        # af_entities = job_sets[self.job_set_id-1].get("entities", None)

        # if af_entities is None:
        #     raise ValueError(
        #         f"No entities found for job {self.job_set_id} in AF jobs. "
        #         "Please check the input yaml file."
        #     )

        # for entity in af_entities:

        #     for _entity_count in range(entity.get("count", 1)):

        #         res_range = entity.get("range", None)

        #         if res_range:
        #             af_offset[atoz[chain_count]] = [
        #                 entity["range"][0],
        #                 entity["range"][1],
        #             ]

        #         else:
        #             warnings.warn(f"Entity {entity['name']} does not have a range.")
        #             self.try_af_offset_from_path = True

        #         chain_count += 1

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
                    start, end = res_range.split(RES_RANGE_SEP)
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

        job_sets = af_input_jobs.get(self.cycle_name, None)

        if job_sets is None:
            return -1

        # assert isinstance(self.job_set_name, str), "You haven't added job set name yet."

        job_set_id = -1 # in case job set name is not found

        for idx, job in enumerate(job_sets):
            if job.get("job_set_name", "") == self.job_set_name:
                job_set_id = idx + 1  # Job IDs are 1-indexed

        if job_set_id == -1 and soft_match:
            for idx, job in enumerate(job_sets):
                if (
                    self.job_set_name in job.get("job_set_name", "")
                    or job.get("job_set_name", "") in self.job_set_name
                ):
                    job_set_id = idx + 1

        self.job_set_id = job_set_id
        return job_set_id