import os
import pytest
from af_pipeline.rank_predictions.rank_af import (
    RankAF2JobSet,
    RankAF3JobSet,
    extract_af_offset_from_path,
    get_job_set_dirs,
    assign_job_set_id,
    extract_af_offset_from_af_input_jobs,
    extract_entity_chain_mapping_from_af_input_jobs,
    extract_entity_chain_mapping_from_path,
)

dummy_job_set_dir = "path/to/job_set_1"
dummy_struct_path = "path/to/job_set_1/protA_1_1-10_protB_1_5-15_1/model_1.cif"

af3_pred_dir = "tests/test_data/af_predictions/af3"
af2_pred_dir = "tests/test_data/af_predictions/af2/job_set_1"
colabfold_pred_dir = "tests/test_data/af_predictions/colabfold/jobset_min"

af_input_jobs = [
    {
        "job_set_name": "job_set_1",
        "modelSeeds": [1, 2],
        "entities": [
            {
                "name": "protA",
                "type": "proteinChain",
                "range": [1, 10],
            },
            {
                "name": "protB",
                "type": "proteinChain",
                "range": [5, 15],
            },
        ],
        "af_offset": {
            "A": [1, 10],
            "B": [5, 15],
        }
    }
]

af_input_jobs2 = [
    {
        "job_set_name": "job_set_1",
        "modelSeeds": [1, 47],
        "entities": [
            {
                "name": "protA",
                "type": "proteinChain",
                "range": [1, 5],
            },
            {
                "name": "protB",
                "type": "proteinChain",
                "range": [11, 15],
            },
        ],
        "af_offset": {
            "A": [1, 5],
            "B": [11, 15],
        }
    }
]

af_input_jobs3 = [
    {
        "job_set_name": "jobset_min",
        "modelSeeds": [1, 0],
        "entities": [
            {
                "name": "protA",
                "type": "proteinChain",
                "range": [1, 19],
            },
            {
                "name": "protB",
                "type": "proteinChain",
                "range": [1, 28],
            },
        ],
        "af_offset": {
            "A": [1, 19],
            "B": [1, 28],
        }
    }
]

@pytest.fixture
def rank_af3_job_set():
    return RankAF3JobSet(
        job_set_dir=dummy_job_set_dir,
        af_input_jobs=af_input_jobs,
    )

@pytest.fixture
def rank_af3_job_set_real():
    return RankAF3JobSet(
        job_set_dir=af3_pred_dir+"/prota_1_1t5_protb_1_11t15",
        af_input_jobs=af_input_jobs2,
    )

@pytest.fixture
def rank_colabfold_job_set_real():
    return RankAF2JobSet(
        job_set_dir=colabfold_pred_dir,
        af_input_jobs=af_input_jobs3,
    )

def test_assign_job_set_id():

    job_set_id = assign_job_set_id(
        job_set_name="job_set_1",
        af_input_jobs=af_input_jobs,
        soft_match=False
    )
    assert job_set_id == 1

    af_input_jobs[0]["job_set_name"] = "job_set_1_subset"
    job_set_id = assign_job_set_id(
        job_set_name="job_set_1",
        af_input_jobs=af_input_jobs,
        soft_match=True
    )
    assert job_set_id == 1

    af_input_jobs[0]["job_set_name"] = "job_set_2"
    job_set_id = assign_job_set_id(
        job_set_name="job_set_1",
        af_input_jobs=af_input_jobs,
        soft_match=True
    )
    assert job_set_id == -1

def test_extract_af_offset_from_af_input_jobs():

    af_input_jobs[0]["job_set_name"] = "job_set_1"
    job_set_id = assign_job_set_id(
        job_set_name="job_set_1",
        af_input_jobs=af_input_jobs,
        soft_match=False
    )
    af_offset = extract_af_offset_from_af_input_jobs(
        job_set_id=job_set_id,
        af_input_jobs=af_input_jobs
    )
    expected_af_offset = {
        "A": [1, 10],
        "B": [5, 15],
    }
    assert af_offset == expected_af_offset; "AF offset does not match expected value"

def test_extract_af_offset_from_path1():

    af_offset = extract_af_offset_from_path(
        structure_path=dummy_struct_path
    )
    expected_af_offset = {
        "A": [1, 10],
        "B": [5, 15],
    }
    assert af_offset == expected_af_offset; "AF offset does not match expected value"

def test_get_job_set_dirs():

    job_set_dirs = get_job_set_dirs(pred_dir=af3_pred_dir, pred_type="AF3")
    expected_dirs = [af3_pred_dir + "/prota_1_1t5_protb_1_11t15"]
    assert set(job_set_dirs) == set(expected_dirs), \
        "Job set directories do not match expected value"

def test_extract_af_offset_from_path2():

    structure_path = "prota_1_1t5_protb_1_11t15_47/fold_prota_1_1t5_protb_1_11t15_47_model_0.cif"

    af_offset = extract_af_offset_from_path(structure_path)
    assert af_offset == {
        "A": [1, 5],
        "B": [11, 15],
    }

def test_extract_entity_chain_mapping_from_af_input_jobs():

    job_set_id = assign_job_set_id(
        job_set_name="job_set_1",
        af_input_jobs=af_input_jobs,
        soft_match=False
    )
    entity_chain_mapping = extract_entity_chain_mapping_from_af_input_jobs(
        job_set_id=job_set_id,
        af_input_jobs=af_input_jobs,
        mapping_type="chain_to_entity",
    )
    expected_mapping = {
        "A": "protA",
        "B": "protB",
    }
    assert entity_chain_mapping == expected_mapping, "Entity-chain mapping does not match expected value."
    entity_chain_mapping = extract_entity_chain_mapping_from_af_input_jobs(
        job_set_id=job_set_id,
        af_input_jobs=af_input_jobs,
        mapping_type="entity_to_chain",
    )
    expected_mapping = {
        "protA": ["A"],
        "protB": ["B"],
    }
    assert entity_chain_mapping == expected_mapping, "Chain-entity mapping does not match expected value."

def test_extract_entity_chain_mapping_from_path():

    structure_path = "prota_1_1t5_protb_1_11t15/prota_1_1t5_protb_1_11t15_47/fold_prota_1_1t5_protb_1_11t15_47_model_0.cif"

    entity_chain_mapping = extract_entity_chain_mapping_from_path(
        structure_path=structure_path,
        mapping_type="chain_to_entity",
    )
    expected_mapping = {
        "A": "prota",
        "B": "protb",
    }
    assert entity_chain_mapping == expected_mapping, "Entity-chain mapping does not match expected value."
    entity_chain_mapping = extract_entity_chain_mapping_from_path(
        structure_path=structure_path,
        mapping_type="entity_to_chain",
    )
    expected_mapping = {
        "prota": ["A"],
        "protb": ["B"],
    }
    assert entity_chain_mapping == expected_mapping, "Chain-entity mapping does not match expected value."

def test_extract_af3_best_pred_data(rank_af3_job_set_real: RankAF3JobSet):

    best_pred_info = rank_af3_job_set_real.extract_af3_best_pred_data()

    assert os.path.basename(best_pred_info["prota_1_1t5_protb_1_11t15"]["structure_path"]) == (
        "fold_prota_1_1t5_protb_1_11t15_1_model_0.cif"
    )

    best_pred_info = rank_af3_job_set_real.extract_af3_best_pred_data()

    assert best_pred_info["prota_1_1t5_protb_1_11t15"]["af_offset"] == {
        "A": [1, 5],
        "B": [11, 15],
    }

def test_extract_colabfold_best_pred_data(rank_colabfold_job_set_real: RankAF2JobSet):

    best_pred_info = rank_colabfold_job_set_real.extract_colabfold_best_pred_data()

    assert os.path.basename(best_pred_info["jobset_min"]["structure_path"]) == (
        "jobset_min_ed465_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_001.pdb"
    )

    best_pred_info = rank_colabfold_job_set_real.extract_colabfold_best_pred_data()

    assert best_pred_info["jobset_min"]["af_offset"] == {
        "A": [1, 19],
        "B": [1, 28],
    }