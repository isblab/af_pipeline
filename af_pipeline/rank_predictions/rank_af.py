import os
import warnings
from pathlib import Path
from collections import defaultdict
from af_pipeline.utils.file_utils import read_json
from af_pipeline.utils.misc_utils import chain_id_gen
from af_pipeline.constants import af_constants
from af_pipeline.constants.af_constants import (
    AF3Metrics,
    AF3SummaryConfidenceFields,
    AFInputEntityFields,
    AFInputJobFields,
    BestPredictionFields,
    FileFormat,
)

def is_valid_job_dir(job_dir: str) -> bool:
    """ Check if a directory is a valid job directory.

    A valid job directory is defined as a directory that contains at least one
    structure file.

    ## Arguments:

    - **job_dir (str)**:<br />
        Path to the job directory

    ## Returns:

    - **bool**:<br />
        True if the directory is a valid job directory, False otherwise
    """

    files = [f for f in Path(job_dir).iterdir() if f.is_file()]

    if any(f.suffix in [f".{FileFormat.CIF}", f".{FileFormat.PDB}"] for f in files):
        return True
    else:
        return False

def get_job_set_dirs(pred_dir:str) -> set:
    """ Get the job set directories from the prediction directory.

    ## Arguments:

    - **pred_dir (str)**:<br />
        Path to the prediction directory

    ## Returns:

    - **set**:<br />
        Set of valid job set directories
    """

    if not os.path.isdir(pred_dir):
        raise ValueError(f"Prediction directory {pred_dir} does not exist.")

    job_set_dirs = set()

    directories = [Path(pred_dir)] + [d for d in Path(pred_dir).rglob('*') if d.is_dir()]

    for directory in directories:

        subdirectories = [d for d in directory.iterdir() if d.is_dir()]
        if len(subdirectories) == 0:
            continue

        if all([is_valid_job_dir(d) for d in subdirectories]):
            print(f"Found job set directory: {directory}")
            job_set_dirs.add(str(directory))

    if len(job_set_dirs) == 0:
        raise ValueError(f"No valid job set directories found in {pred_dir}.")

    return job_set_dirs

def extract_entity_chain_mapping(
    job_set_id: int = -1,
    af_input_jobs: list| None = None,
    structure_path: str | None = None,
    mapping_type: str = "chain_to_entity",
) -> dict:
    """ Extract the entity-chain mapping for a given job set ID from the AF
    input jobs in the config dictionary or from the structure path.

    ## Arguments:

    - **job_set_id (int, optional):**:<br />
        ID of the job set in af_input_jobs (1-indexed). This will be used to
        extract the entity-chain mapping from the AlphaFold input jobs in the config
        dictionary.

    - **af_input_jobs (list | None, optional):**:<br />
        List of AlphaFold input jobs from the config dictionary. This will be used to
        extract the entity-chain mapping for the given job_set_id.

    - **structure_path (str | None, optional):**:<br />
        Path to the structure file. This will be used to extract the
        entity-chain mapping if af_input_jobs is not provided.

    - **mapping_type (str, optional):**:<br />
        Type of mapping to extract. Should be either "chain_to_entity" or
        "entity_to_chain".

    ## Returns:

    - **dict**:<br />
        Dictionary containing the entity-chain mapping.
        - If mapping_type is "chain_to_entity", the keys will be chain IDs
            and the values will be entity names.
        - If mapping_type is "entity_to_chain", the keys will be entity names
            and the values will be lists of chain IDs.
    """

    assert not (af_input_jobs is None and structure_path is None), (
        "Either af_input_jobs or structure_path should be provided."
    )
    assert mapping_type in ["chain_to_entity", "entity_to_chain"], (
        "Invalid mapping type. Should be 'chain_to_entity' or 'entity_to_chain'."
    )

    mapping = extract_entity_chain_mapping_from_af_input_jobs(
        job_set_id=job_set_id,
        af_input_jobs=af_input_jobs,
        mapping_type=mapping_type,
    )

    if len(mapping) == 0:
        try:
            mapping = extract_entity_chain_mapping_from_path(
                structure_path=structure_path,
                mapping_type=mapping_type,
            )
        except ValueError as e:
            print(f"Error extracting entity-chain mapping from path {structure_path}: {e}")

    return mapping

def extract_entity_chain_mapping_from_af_input_jobs(
    job_set_id: int = -1,
    af_input_jobs: list| None = None,
    mapping_type: str = "chain_to_entity",
) -> dict:
    """ Extract the entity-chain mapping for a given job set ID from the AF
    input jobs in the config dictionary.

    ## Arguments:

    - **job_set_id (int, optional):**:<br />
        ID of the job set in af_input_jobs (1-indexed). This will be used to
        extract the entity-chain mapping from the AlphaFold input jobs in the config
        dictionary.

    - **af_input_jobs (list | None, optional):**:<br />
        List of AlphaFold input jobs from the config dictionary. This will be used to
        extract the entity-chain mapping for the given job_set_id.

    - **mapping_type (str, optional):**:<br />
        Type of mapping to extract. Should be either "chain_to_entity" or
        "entity_to_chain".

    ## Returns:

    - **dict**:<br />
        Dictionary containing the entity-chain mapping.
        - If mapping_type is "chain_to_entity", the keys will be chain IDs
            and the values will be entity names.
        - If mapping_type is "entity_to_chain", the keys will be entity names
            and the values will be lists of chain IDs.
    """

    assert mapping_type in ["chain_to_entity", "entity_to_chain"], (
        "Invalid mapping type. Should be 'chain_to_entity' or 'entity_to_chain'."
    )

    mapping = {}

    if af_input_jobs is None:
        warnings.warn("No AF jobs found. Skipping extraction of entity-chain mapping.")
        return mapping

    if job_set_id < 1 or job_set_id > len(af_input_jobs):
        return mapping

    job_set: dict = af_input_jobs[job_set_id-1]

    if mapping_type == "chain_to_entity":
        chainGen = chain_id_gen()
        entities = job_set.get(AFInputJobFields.ENTITIES, [])
        for entity in entities:
            entity_name = entity[AFInputEntityFields.NAME]
            entity_count = entity.get(AFInputEntityFields.COUNT, 1)
            for _ in range(entity_count):
                chain_id = next(chainGen)
                mapping[chain_id] = entity_name

    else:
        mapping = defaultdict(list)
        chainGen = chain_id_gen()
        entities = job_set.get(AFInputJobFields.ENTITIES, [])
        for entity in entities:
            entity_name = entity[AFInputEntityFields.NAME]
            entity_count = entity.get(AFInputEntityFields.COUNT, 1)
            for _ in range(entity_count):
                chain_id = next(chainGen)
                mapping[entity_name].append(chain_id)

    return mapping

def extract_entity_chain_mapping_from_path(
    structure_path: str,
    mapping_type: str = "chain_to_entity"
) -> dict:
    """ Extract the entity-chain mapping for a given structure path.

    ## Arguments:

    - **structure_path (str)**:<br />
        Path to the structure file. This will be used to extract the
        entity-chain mapping based on the directory name.

    - **mapping_type (str, optional):**:<br />
        Type of mapping to extract. Should be either "chain_to_entity" or
        "entity_to_chain".

    ## Returns:

    - **dict**:<br />
        Dictionary containing the entity-chain mapping.
        - If mapping_type is "chain_to_entity", the keys will be chain IDs
            and the values will be entity names.
        - If mapping_type is "entity_to_chain", the keys will be entity names
            and the values will be lists of chain IDs.
    """

    assert mapping_type in ["chain_to_entity", "entity_to_chain"], (
        "Invalid mapping type. Should be 'chain_to_entity' or 'entity_to_chain'."
    )

    dirname = os.path.basename(os.path.dirname(structure_path))
    mapping = {} if mapping_type == "chain_to_entity" else defaultdict(list)

    if len(dirname.split("_")[0:-1]) % 3 == 0:
        seed_is_present = True

    elif len(dirname.split("_")[0:-1]) % 3 == 2:
        seed_is_present = False

    else:
        raise ValueError(
            "Invalid directory name for AF3 prediction, "
            f"should be in the format p1_copy1_1{af_constants.RES_RANGE_SEP}100_p2_copy2_101{af_constants.RES_RANGE_SEP}200_seed "
            f"or p1_copy1_1{af_constants.RES_RANGE_SEP}100_p2_copy2_101{af_constants.RES_RANGE_SEP}200"
        )

    if seed_is_present:
        p_c_r = dirname.split("_")[0:-1]  # Exclude the seed

    else:
        p_c_r = dirname.split("_") # Seed is not present in the path

    chainGen = chain_id_gen()
    for i in range(len(p_c_r) // 3):
        entity_name = p_c_r[i * 3]
        entity_count = p_c_r[i * 3 + 1]

        for _ in range(int(entity_count)):

            chain_id = next(chainGen)
            if mapping_type == "chain_to_entity":
                mapping[chain_id] = entity_name
            else:
                mapping[entity_name].append(chain_id)

    return mapping

def extract_af_offset(
    job_set_id: int = -1,
    af_input_jobs: list| None = None,
    structure_path: str | None = None,
) -> dict:
    """ Extract the AF offset for a given job set ID from the AlphaFold input jobs
    in the config dictionary or from the structure prediction path.

    ## Arguments:

    - **job_set_id (int, optional):**:<br />
        ID of the job set in af_input_jobs (1-indexed). This will be used to
        extract the AF offset from the AlphaFold input jobs in the config dictionary.

    - **af_input_jobs (list | None, optional):**:<br />
        List of AlphaFold input jobs from the config dictionary. This will be used to
        extract the AF offset for the given job_set_id.

    - **structure_path (str | None, optional):**:<br />
        Path to the structure file. This will be used to extract the AF offset
        if af_input_jobs is not provided.

    ## Returns:

    - **dict**:<br />
        Dictionary containing the AF offset for each chain in the format -
        ```
        {
            "A": [start, end],
            "B": [start, end]
        }
        ```
    """

    assert not (af_input_jobs is None and structure_path is None), (
        "Either af_input_jobs or structure_path should be provided."
    )

    af_offset = extract_af_offset_from_af_input_jobs(
        job_set_id=job_set_id,
        af_input_jobs=af_input_jobs,
    )

    if len(af_offset) == 0:
        try:
            af_offset = extract_af_offset_from_path(
                structure_path=structure_path,
            )
        except ValueError as e:
            print(f"Error extracting AF offset from path {structure_path}: {e}")

    return af_offset

def extract_af_offset_from_af_input_jobs(
    job_set_id: int = -1,
    af_input_jobs: list| None = None
) -> dict:
    """ Extract the offsets for AF3 predictions from the AF jobs
    dictionary in the config yaml file.

    ## Arguments:

    - **job_set_id (int, optional):**:<br />
        ID of the job set in af_input_jobs (1-indexed). This will be used to
        extract the AF offset from the af_input_jobs.

    - **af_input_jobs (dict | None, optional):**:<br />
        Dictionary containing AlphaFold input data

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

    if af_input_jobs is None:
        warnings.warn("No AF jobs found. Skipping extraction of offsets.")
        return af_offset

    if job_set_id < 1 or job_set_id > len(af_input_jobs):
        return af_offset

    job_set: dict = af_input_jobs[job_set_id-1]
    af_offset = job_set.get(AFInputJobFields.AF_OFFSET, {})
    if len(af_offset) == 0:
        warnings.warn(
            f"No AF offset found for job {job_set_id} in AF jobs. "
            "Please check the input yaml file."
        )

    return af_offset

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
    print(dirname)
    af_offset = {}

    if len(dirname.split("_")[0:-1]) % 3 == 0:
        seed_is_present = True

    elif len(dirname.split("_")[0:-1]) % 3 == 2:
        seed_is_present = False

    else:
        raise ValueError(
            "Invalid directory name for AF3 prediction, "
            f"should be in the format p1_copy1_1{af_constants.RES_RANGE_SEP}100_p2_copy2_101{af_constants.RES_RANGE_SEP}200_seed "
            f"or p1_copy1_1{af_constants.RES_RANGE_SEP}100_p2_copy2_101{af_constants.RES_RANGE_SEP}200"
        )

    if seed_is_present:
        p_c_r = dirname.split("_")[0:-1]  # Exclude the seed

    else:
        p_c_r = dirname.split("_") # Seed is not present in the path

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
                        f"Expected format is 'START{af_constants.RES_RANGE_SEP}END' or 'START-END'."
                    )
            af_offset[chain_id] = [int(start), int(end)]

    return af_offset

def assign_job_set_id(
    job_set_name: str,
    af_input_jobs: list | None = None,
    soft_match = False
) -> int:
    """ Get the job set ID for the current job set based on the job set name.

    ## Arguments:

    - **job_set_name (str):**:<br />
        Name of the job set. This will be used to match the job set name in the
        AF input jobs in the config dictionary and extract the corresponding job
        set ID.

    - **af_input_jobs (dict | None, optional):**:<br />
        Dictionary containing AlphaFold input data

    - **soft_match (bool, optional):**:<br />
        If True, allow partial matching of job set names

    ## Returns:

    - **int**:<br />
        Job id for the current job set
    """
    job_set_id = -1

    if af_input_jobs is None:
        warnings.warn("No AF jobs found. Skipping extraction of job id.")
        return job_set_id

    idx = 0

    _condition = {
        True: lambda x, y: x in y or y in x,
        False: lambda x, y: x == y,
    }

    for job_set in af_input_jobs:

        lookup_name: str = job_set.get(AFInputJobFields.JOB_SET_NAME, "")

        if _condition[soft_match](job_set_name.lower(), lookup_name.lower()):
            job_set_id = idx + 1
            break

    return job_set_id

class RankAF2JobSet:
    """ Class to rank AF2 or ColabFold predictions for a given job set directory """

    job_set_dir: str
    """ Path to the job set directory"""

    job_set_name: str
    """ Name of the job set"""

    job_set_id: int
    """ ID of the job set in af_input_jobs (1-indexed)"""

    af_input_jobs: list | None
    """ List of AlphaFold input jobs from the config dictionary"""

    data_format: str = "json"
    """ Format of the AlphaFold2 output data file. Should be either "json" or "pkl".
    Default is "json"."""

    def __init__(
        self,
        job_set_dir:str,
        job_set_id: int = -1,
        af_input_jobs: list| None = None,
        soft_match: bool = False,
        **kwargs,
    ):
        self.job_set_dir = os.path.abspath(job_set_dir)
        self.job_set_name = os.path.basename(self.job_set_dir)
        self.job_set_id = job_set_id
        self.data_format = kwargs.get("data_format", "json")
        if job_set_id != -1:
            self.job_set_id = job_set_id
        else:
            self.job_set_id = assign_job_set_id(
                job_set_name=self.job_set_name,
                af_input_jobs=af_input_jobs,
                soft_match=soft_match,
            )
        self.af_input_jobs = af_input_jobs

    def extract_af2_best_pred_data(self) -> list:
        """ Extract the paths to the best model and data for a given AlphaFold2
        predictions directory and return them in a dictionary along with the
        residue offsets and entity-chain mapping.

        ## Returns:

        - **dict**:<br />
            Dictionary containing paths to the best model and data for a given
            prediction directory along with offsets and entity-chain mapping.
        """

        best_pred_info = {}

        ranking_debug_json = os.path.join(self.job_set_dir, "ranking_debug.json")
        ranking_debug = read_json(ranking_debug_json)
        order = ranking_debug.get("order", [])
        # scores = ranking_debug.get("iptm+ptm", [])

        structure_path = os.path.join(self.job_set_dir, f"relaxed_{order[0]}.pdb")
        if not os.path.exists(structure_path):
            structure_path = os.path.join(self.job_set_dir, f"unrelaxed_{order[0]}.pdb")

        data_path = os.path.join(self.job_set_dir, f"result_{order[0]}.{self.data_format}")

        af_offset = extract_af_offset(
            job_set_id=self.job_set_id,
            af_input_jobs=self.af_input_jobs,
            structure_path=structure_path,
        )

        mapping = extract_entity_chain_mapping(
            job_set_id=self.job_set_id,
            af_input_jobs=self.af_input_jobs,
            structure_path=structure_path,
            mapping_type="chain_to_entity",
        )

        key = os.path.basename(os.path.dirname(os.path.dirname(structure_path)))

        best_pred_info[key] = {
            BestPredictionFields.STRUCTURE_PATH: structure_path,
            BestPredictionFields.DATA_PATH: data_path,
            BestPredictionFields.AF_OFFSET: af_offset,
            BestPredictionFields.ENTITY_CHAIN_MAP: mapping,
        }

        return best_pred_info

    def extract_colabfold_best_pred_data(self) -> dict:
        """ Extract the paths to the best model and data for a given ColabFold
        predictions directory and return them in a dictionary along with the
        residue offsets and entity-chain mapping.

        ## Returns:

        - **dict**:<br />
            Dictionary containing paths to the best model and data for a given
            prediction directory along with offsets and entity-chain mapping.
        """

        best_pred_info = {}

        structure_path = self.get_best_colabfold_model_path()
        data_path = structure_path.replace("_unrelaxed_", "_scores_").replace(".pdb", ".json")

        key = os.path.basename(os.path.dirname(os.path.dirname(structure_path)))

        af_offset = extract_af_offset(
            job_set_id=self.job_set_id,
            af_input_jobs=self.af_input_jobs,
            structure_path=structure_path,
        )
        print(self.job_set_id)
        mapping = extract_entity_chain_mapping(
            job_set_id=self.job_set_id,
            af_input_jobs=self.af_input_jobs,
            structure_path=structure_path,
            mapping_type="chain_to_entity",
        )

        best_pred_info[key] = {
            BestPredictionFields.STRUCTURE_PATH: structure_path,
            BestPredictionFields.DATA_PATH: data_path,
            BestPredictionFields.AF_OFFSET: af_offset,
            BestPredictionFields.ENTITY_CHAIN_MAP: mapping,
        }

        return best_pred_info

    def get_best_colabfold_model_path(self) -> str:
        """ Obtain the path of the rank 1 model from the ColabFold predictions directory.

        ## Returns:

        - **str**:<br />
            Path to the best ColabFold model file.
        """

        model_paths = list(Path(self.job_set_dir).glob("*.pdb"))
        if len(model_paths) == 0:
            raise ValueError(f"No model files found in {self.job_set_dir}.")

        model_ranks = []
        for model_path in model_paths:
            model_name = model_path.stem
            rank_str = model_name.split("relaxed_rank_")[-1].split("_")[0]
            try:
                rank = int(rank_str)
                model_ranks.append((model_path, rank))
            except ValueError:
                print(f"Could not extract rank from model name {model_name}. Skipping this file.")

        model_ranks.sort(key=lambda x: x[1])
        structure_path = str(model_ranks[0][0])

        return structure_path

class RankAF3JobSet:
    """ Class to rank AF3 predictions for a given job set directory """

    job_set_dir: str
    """ Path to the job set directory"""

    job_set_name: str
    """ Name of the job set"""

    job_set_id: int
    """ ID of the job set in af_input_jobs (1-indexed)"""

    af_input_jobs: list | None
    """ List of AlphaFold input jobs from the config dictionary"""

    def __init__(
        self,
        job_set_dir:str,
        job_set_id: int = -1,
        af_input_jobs: list| None = None,
        soft_match: bool = False,
    ):
        self.job_set_dir = os.path.abspath(job_set_dir)
        self.job_set_name = os.path.basename(self.job_set_dir)
        if job_set_id != -1:
            self.job_set_id = job_set_id
        else:
            self.job_set_id = assign_job_set_id(
                job_set_name=self.job_set_name,
                af_input_jobs=af_input_jobs,
                soft_match=soft_match,
            )
        self.af_input_jobs = af_input_jobs

    def extract_af3_best_pred_data(self) -> dict:
        """ Extract AF3 model metrics and paths to the best model and data for a
        given prediction directory

        ## Returns:

        - **dict**:<br />
            Dictionary containing paths to the best model and data for a given
            prediction directory along with offsets and entity-chain mapping.
        """

        best_pred_info = {}

        _best_seed, _ranking_score, structure_path, _best_model_idx = self.rank_seeds()

        if structure_path is None:
            warnings.warn(
                f"No valid AF3 predictions found for {self.job_set_name}. "
                "Skipping ranking of AF3 predictions."
            )
            return []

        data_path = self.get_data_path_from_structure_path(structure_path)

        af_offset = extract_af_offset(
            job_set_id=self.job_set_id,
            af_input_jobs=self.af_input_jobs,
            structure_path=structure_path,
        )

        mapping = extract_entity_chain_mapping(
            job_set_id=self.job_set_id,
            af_input_jobs=self.af_input_jobs,
            structure_path=structure_path,
            mapping_type="chain_to_entity",
        )

        key = os.path.basename(os.path.dirname(os.path.dirname(structure_path)))

        best_pred_info[key] = {
            BestPredictionFields.STRUCTURE_PATH: structure_path,
            BestPredictionFields.DATA_PATH: data_path,
            BestPredictionFields.AF_OFFSET: af_offset,
            BestPredictionFields.ENTITY_CHAIN_MAP: mapping,
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

        print("\n" + self.job_set_name + f"\n{"-"*len(self.job_set_name)}")
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

            if len(summary_confidences_files) != 5:
                warnings.warn(
                    f"There should be 5 summary_confidences files in {af3_seed_prediction_dir}. "
                    f"Found: {len(summary_confidences_files)}. "
                    "Skipping this."
                )
                continue

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

        metric_data: dict = read_json(af3_summary_confidence_file)

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

        job_request_data: dict = read_json(af3_job_request_file)

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
