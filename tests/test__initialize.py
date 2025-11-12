from af_pipeline._initialize import _Initialize
import pytest
import Bio.PDB.Structure
from af_pipeline.parser.structure_parser import StructureParser
from af_pipeline.parser.data_parser import DataParser

struct_path1 = "tests/test_data/af_predictions/af3/fold_dummy_job_2/fold_dummy_job_2_model_0.cif"
data_path_1 = "tests/test_data/af_predictions/af3/fold_dummy_job_2/fold_dummy_job_2_full_data_0.json"

@pytest.fixture
def initialize1():
    return _Initialize(
        data_file_path=data_path_1,
        structure_file_path=struct_path1,
        af_offset={},
        rep_atom_dict={},
        average_token_pae=True,
        average_token_plddt=True,
        metric_level="per_token",
    )

@pytest.fixture
def initialize2():
    return _Initialize(
        data_file_path=data_path_1,
        structure_file_path=struct_path1,
        af_offset={},
        rep_atom_dict={},
        average_token_pae=True,
        average_token_plddt=True,
        metric_level="representative_token",
    )

def test_initialize_properties(
    initialize1: _Initialize,
    initialize2: _Initialize,
):
    """Test the _Initialize properties."""

    for initializer in [initialize1, initialize2]:

        assert isinstance(initializer.structure, Bio.PDB.Structure.Structure), \
            "Structure should be a Bio.PDB.Structure.Structure instance."

        assert isinstance(initializer.data_parser, DataParser), \
            "data_parser should be a DataParser instance."

        assert isinstance(initializer.structure_parser, StructureParser), \
            "structure_parser should be a StructureParser instance."

    with pytest.raises(Exception):
        _ = _Initialize(
            data_file_path=data_path_1,
            structure_file_path=struct_path1,
            af_offset={},
            rep_atom_dict={},
            average_token_pae=True,
            average_token_plddt=True,
            metric_level="invalid_level",
        )