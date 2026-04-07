import pytest
from af_pipeline.af_input.colabfold import ColabFold
from af_pipeline.constants.af_constants import (
    AFInputJobFields,
    EntityType,
    AFInputEntityFields,
    ConfigYaml,
)

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
        AFInputJobFields.JOB_SET_NAME: "jobset_min",
        AFInputJobFields.MODEL_SEEDS: [0, 1],
        AFInputJobFields.ENTITIES: [{
            AFInputEntityFields.NAME: "protB",
            AFInputEntityFields.TYPE: EntityType.PROTEIN_CHAIN,
        },{
            AFInputEntityFields.NAME: "dnaA",
            AFInputEntityFields.TYPE: EntityType.DNA_SEQUENCE,
        },{
            AFInputEntityFields.NAME: "protA",
            AFInputEntityFields.TYPE: EntityType.PROTEIN_CHAIN,
        },]
    }],
}
config_dict = {
    ConfigYaml.AF_INPUT_JOBS: input_dict,
    ConfigYaml.PROTEIN_UNIPROT_MAP: entities_map,
}

@pytest.fixture
def colabfold():
    return ColabFold(
        config_dict=config_dict,
        protein_sequences=protein_sequences_by_id,
    )

def test_create_colabfold_job_cycles(colabfold: ColabFold):

    expected_fasta_dict = {
        "jobset_min": "GAVLILLLVAVAVVAGVAA:\nMKTAYIAKQRQISFVKSHFSRQDILDLI",
    }

    colabfold.create_colabfold_job_cycles()

    assert "cycle1" in colabfold.job_cycles, "Job cycle 'cycle1' not found in output."

    job_list = colabfold.job_cycles["cycle1"]
    assert len(job_list) == 1, "Expected one job in 'cycle1'."

    for fasta_dict, job_name in job_list:
        assert job_name == "jobset_min", "Job name mismatch."
        assert fasta_dict[job_name] == expected_fasta_dict[job_name], \
            "FASTA content mismatch."

