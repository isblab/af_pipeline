import os
import pytest
from af_pipeline.af_input.colabfold import ColabFold

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
        },]
    }],
}

@pytest.fixture
def colabfold():
    return ColabFold(
        input_dict=input_dict,
        protein_sequences=protein_sequences_by_id,
        entities_map=entities_map,
    )

def test_create_colabfold_job_cycles(colabfold: ColabFold):

    expected_fasta_dict = {
        "jobset_min": "GAVLILLLVAVAVVAGVAA:\nMKTAYIAKQRQISFVKSHFSRQDILDLI",
    }

    job_cycles = colabfold.create_colabfold_job_cycles()

    assert "cycle1" in job_cycles, "Job cycle 'cycle1' not found in output."

    job_list = job_cycles["cycle1"]
    assert len(job_list) == 1, "Expected one job in 'cycle1'."

    for fasta_dict, job_name in job_list:
        assert job_name == "jobset_min", "Job name mismatch."
        assert fasta_dict[job_name] == expected_fasta_dict[job_name], \
            "FASTA content mismatch."

