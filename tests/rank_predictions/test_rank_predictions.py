import os
import pytest
from af_pipeline.rank_predictions.rank_af import (
    RankAF3JobSet,
    get_job_set_dirs
)

dummy_job_set_dir = "path/to/cycle_1/job_set_1"
dummy_struct_path = "path/to/cycle_1/job_set_1/protA_1_1-10_protB_1_5-15_1/model_1.cif"

pred_dir = "tests/test_data/af_predictions/af3/dummy_cycle"

af_input_jobs = {
    "cycle_1": [
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
}

af_input_jobs2 = {
    "cycle_1": [
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
}



@pytest.fixture
def rank_af3_job_set():
    return RankAF3JobSet(
        job_set_dir=dummy_job_set_dir,
    )

@pytest.fixture
def rank_af3_job_set_real():
    return RankAF3JobSet(
        job_set_dir=pred_dir+"/prota_1_1t5_protb_1_11t15",
        try_af_offset_from_path=True,
    )

def test_add_job_set_id(rank_af3_job_set: RankAF3JobSet):

    job_set_id = rank_af3_job_set.add_job_set_id(
        af_input_jobs=af_input_jobs,
        soft_match=False
    )
    assert job_set_id == 1

    af_input_jobs["cycle_1"][0]["job_set_name"] = "job_set_1_subset"
    job_set_id = rank_af3_job_set.add_job_set_id(
        af_input_jobs=af_input_jobs,
        soft_match=True
    )
    assert job_set_id == 1

    af_input_jobs["cycle_1"][0]["job_set_name"] = "job_set_2"
    job_set_id = rank_af3_job_set.add_job_set_id(
        af_input_jobs=af_input_jobs,
        soft_match=True
    )
    assert job_set_id == -1

def test_extract_af_offset_from_af_input_jobs(rank_af3_job_set: RankAF3JobSet):

    af_input_jobs["cycle_1"][0]["job_set_name"] = "job_set_1"
    rank_af3_job_set.add_job_set_id(
        af_input_jobs=af_input_jobs,
        soft_match=False
    )
    af_offset = rank_af3_job_set.extract_af_offset_from_af_input_jobs(
        af_input_jobs=af_input_jobs
    )
    expected_af_offset = {
        "A": [1, 10],
        "B": [5, 15],
    }
    assert af_offset == expected_af_offset; "AF offset does not match expected value"

def test_extract_af_offset_from_path(rank_af3_job_set: RankAF3JobSet):

    af_offset = rank_af3_job_set.extract_af_offset_from_path(
        structure_path=dummy_struct_path
    )
    expected_af_offset = {
        "A": [1, 10],
        "B": [5, 15],
    }
    assert af_offset == expected_af_offset; "AF offset does not match expected value"

def test_get_job_set_dirs():

    job_set_dirs = get_job_set_dirs(pred_dir=pred_dir)
    expected_dirs = [
        "tests/test_data/af_predictions/af3/dummy_cycle/prota_1_1t5_protb_1_11t15",
    ]
    assert set(job_set_dirs) == set(expected_dirs), \
        "Job set directories do not match expected value"

def test_extract_af_offset_from_path(rank_af3_job_set_real: RankAF3JobSet):

    structure_path = "tests/test_data/af_predictions/af3/dummy_cycle/prota_1_1t5_protb_1_11t15/tests/test_data/af_predictions/af3/dummy_cycle/prota_1_1t5_protb_1_11t15/prota_1_1t5_protb_1_11t15_47/fold_prota_1_1t5_protb_1_11t15_47_model_0.cif"

    af_offset = rank_af3_job_set_real.extract_af_offset_from_path(structure_path)
    assert af_offset == {
        "A": [1, 5],
        "B": [11, 15],
    }

def test_extract_af3_best_pred_data(rank_af3_job_set_real: RankAF3JobSet):

    best_pred_info = rank_af3_job_set_real.extract_af3_best_pred_data()

    assert os.path.basename(best_pred_info["prota_1_1t5_protb_1_11t15"]["structure_path"]) == (
        "fold_prota_1_1t5_protb_1_11t15_1_model_0.cif"
    )

    best_pred_info = rank_af3_job_set_real.extract_af3_best_pred_data(
        af_input_jobs=af_input_jobs2
    )

    assert best_pred_info["prota_1_1t5_protb_1_11t15"]["af_offset"] == {
        "A": [1, 5],
        "B": [11, 15],
    }