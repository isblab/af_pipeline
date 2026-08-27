import os
import pytest
import tempfile
from af_pipeline.af_input.alphafold2 import AF2Config
from af_pipeline.constants.af_constants import (
    RES_RANGE_SEP,
    ConfigYaml,
    AFInputJobFields,
    AFInputEntityFields,
    EntityType,
)

test_out_dir = tempfile.mkdtemp()


protein_sequences_by_id = {
    "A12345": "MKTAYIAKQRQISFVKSHFSRQDILDLI",
    "B67890": "GAVLILLLVAVAVVAGVAA",
}
entities_map = {
    "protA": "A12345",
    "protB": "B67890",
}
input_job_sets = [{
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
        AFInputEntityFields.COUNT: 2,
    },]
}]
config_dict = {
    ConfigYaml.PROTEIN_UNIPROT_MAP: entities_map,
    ConfigYaml.AF_INPUT_JOBS: input_job_sets,
}

input_job_sets2 = [{
    AFInputJobFields.MODEL_SEEDS: [0, 1],
    AFInputJobFields.ENTITIES: [{
        AFInputEntityFields.NAME: "protB",
        AFInputEntityFields.TYPE: EntityType.PROTEIN_CHAIN,
    },{
        AFInputEntityFields.NAME: "protA",
        AFInputEntityFields.TYPE: EntityType.PROTEIN_CHAIN,
        AFInputEntityFields.COUNT: 2,
        AFInputEntityFields.RANGE: [1, 10],
    },]
}]
config_dict2 = {
    ConfigYaml.PROTEIN_UNIPROT_MAP: entities_map,
    ConfigYaml.AF_INPUT_JOBS: input_job_sets2,
}

@pytest.fixture
def alphafold2():
    return AF2Config(
        config_dict=config_dict,
        protein_sequences=protein_sequences_by_id,
    )

@pytest.fixture
def alphafold2_no_job_names():
    return AF2Config(
        config_dict=config_dict2,
        protein_sequences=protein_sequences_by_id,
    )

def test_write_job_files(alphafold2: AF2Config):

    alphafold2.write_job_files(
        output_dir=os.path.join(test_out_dir, "af2_input_jobs"),
    )

    assert os.path.isfile(os.path.join(
        test_out_dir, "af2_input_jobs", "jobset_min", "jobset_min.fasta"
    ))