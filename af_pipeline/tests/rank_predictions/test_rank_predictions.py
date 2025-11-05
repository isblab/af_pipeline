import pytest
from af_pipeline.rank_predictions.rank_af import RankAF3JobSet

dummy_job_set_dir = "path/to/cycle_1/job_set_1"
dummy_struct_path = "path/to/cycle_1/job_set_1/protA_1_1-10_protB_1_5-15_1/model_1.cif"

af_input_jobs = {
    "cycle_1": [
        {
            "name": "job_set_1",
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

@pytest.fixture
def rank_af3_job_set():
    return RankAF3JobSet(
        job_set_dir=dummy_job_set_dir,
    )

def test_add_job_set_id(rank_af3_job_set: RankAF3JobSet):

    job_set_id = rank_af3_job_set.add_job_set_id(
        af_input_jobs=af_input_jobs,
        soft_match=False
    )
    assert job_set_id == 1

    af_input_jobs["cycle_1"][0]["name"] = "job_set_1_subset"
    job_set_id = rank_af3_job_set.add_job_set_id(
        af_input_jobs=af_input_jobs,
        soft_match=True
    )
    assert job_set_id == 1

    af_input_jobs["cycle_1"][0]["name"] = "job_set_2"
    job_set_id = rank_af3_job_set.add_job_set_id(
        af_input_jobs=af_input_jobs,
        soft_match=True
    )
    assert job_set_id == -1

def test_extract_af_offset_from_af_input_jobs(rank_af3_job_set: RankAF3JobSet):

    af_input_jobs["cycle_1"][0]["name"] = "job_set_1"
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