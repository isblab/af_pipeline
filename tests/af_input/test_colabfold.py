import os
import pytest
import tempfile
from af_pipeline.af_input.colabfold import ColabFold
from af_pipeline.constants.af_constants import (
    AFInputJobFields,
    EntityType,
    AFInputEntityFields,
    ConfigYaml,
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
    },]
}]
config_dict = {
    ConfigYaml.AF_INPUT_JOBS: input_job_sets,
    ConfigYaml.PROTEIN_UNIPROT_MAP: entities_map,
}

@pytest.fixture
def colabfold():
    return ColabFold(
        config_dict=config_dict,
        protein_sequences=protein_sequences_by_id,
    )

def test_write_job_files(colabfold: ColabFold):

    colabfold.write_job_files(
        output_dir=os.path.join(test_out_dir, "cf_input_jobs"),
    )

    assert os.path.isfile(os.path.join(
        test_out_dir, "cf_input_jobs", "jobset_min", "jobset_min.fasta"
    ))