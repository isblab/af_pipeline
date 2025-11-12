import os
import pytest
from af_pipeline.af_input.alphafold2 import AlphaFold2
from af_pipeline.constants.af_constants import RES_RANGE_SEP

test_out_dir = "tests/test_output"

protein_sequences_by_id = {
    "A12345": "MKTAYIAKQRQISFVKSHFSRQDILDLI",
    "B67890": "GAVLILLLVAVAVVAGVAA",
}
entities_map = {
    "protA": "A12345",
    "protB": "B67890",
}
input_dict = {
    "cycle1": [{
        "job_set_name": "jobset_min",
        "modelSeeds": [0, 1],
        "entities": [{
            "name": "protB",
            "type": "proteinChain",
        },{
            "name": "dnaA",
            "type": "dnaSequence",
        },{
            "name": "protA",
            "type": "proteinChain",
            "count": 2,
        },]
    }],
}

input_dict2 = {
    "cycle1": [{
        "modelSeeds": [0, 1],
        "entities": [{
            "name": "protB",
            "type": "proteinChain",
        },{
            "name": "protA",
            "type": "proteinChain",
            "count": 2,
            "range": [1, 10],
        },]
    }],
}

@pytest.fixture
def alphafold2():
    return AlphaFold2(
        input_dict=input_dict,
        protein_sequences=protein_sequences_by_id,
        entities_map=entities_map,
    )

@pytest.fixture
def alphafold2_no_job_names():
    return AlphaFold2(
        input_dict=input_dict2,
        protein_sequences=protein_sequences_by_id,
        entities_map=entities_map,
    )

def test_create_alphafold2_job_cycles(
    alphafold2: AlphaFold2,
    alphafold2_no_job_names: AlphaFold2,
):

    expected_seq_to_add = {
        f"protB_1_1{RES_RANGE_SEP}19": "GAVLILLLVAVAVVAGVAA",
        f"protA_1_1{RES_RANGE_SEP}28": "MKTAYIAKQRQISFVKSHFSRQDILDLI",
        f"protA_2_1{RES_RANGE_SEP}28": "MKTAYIAKQRQISFVKSHFSRQDILDLI",
    }

    job_cycles = alphafold2.create_af2_job_cycles()

    assert "cycle1" in job_cycles, "Job cycle 'cycle1' not found in output."

    job_list = job_cycles["cycle1"]
    assert len(job_list) == 1, "Expected one job in 'cycle1'."

    for seq_to_add, job_name in job_list:
        assert job_name == "jobset_min", "Job name mismatch."
        for key in expected_seq_to_add:
            assert seq_to_add[key] == expected_seq_to_add[key], \
                f"Sequence for {key} mismatch."

    job_cycles = alphafold2_no_job_names.create_af2_job_cycles()

    expected_seq_to_add = {
        f"protB_1_1{RES_RANGE_SEP}19": "GAVLILLLVAVAVVAGVAA",
        f"protA_1_1{RES_RANGE_SEP}10": "MKTAYIAKQR",
        f"protA_2_1{RES_RANGE_SEP}10": "MKTAYIAKQR",
    }

    for seq_to_add, job_name in job_cycles["cycle1"]:
        assert job_name == (
            f"protB_1_1{RES_RANGE_SEP}19_protA_2_1{RES_RANGE_SEP}10"
        ), "Job name mismatch."
        for key in expected_seq_to_add:
            assert seq_to_add[key] == expected_seq_to_add[key], \
                f"Sequence for {key} mismatch."

def test_write_job_files(alphafold2: AlphaFold2):

    job_cycles = alphafold2.create_af2_job_cycles()

    alphafold2.write_job_files(
        job_cycles=job_cycles,
        output_dir=os.path.join(test_out_dir, "af2_input_jobs"),
    )

    assert os.path.isfile(os.path.join(
        test_out_dir, "af2_input_jobs", "cycle1", "jobset_min.fasta"
    ))