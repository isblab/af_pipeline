import os
import pytest
from af_pipeline.parser.initialize import Initialize
from af_pipeline.rigid_bodies.rigid_bodies import RigidBodies

struct_path1 = "tests/test_data/af_predictions/af3/fold_dummy_job_1/fold_dummy_job_1_model_0.cif"
data_path_1 = "tests/test_data/af_predictions/af3/fold_dummy_job_1/fold_dummy_job_1_full_data_0.json"

intializer = Initialize(
    data_file_path=data_path_1,
    structure_file_path=struct_path1,
    af_offset={},
    rep_atom_dict={},
    average_token_pae=False,
    average_token_plddt=False,
    metric_level="representative_token",
)

@pytest.fixture
def rigid_bodies1():

    return RigidBodies(
        library="networkx",
        plddt_cutoff=70.0,
        pae_cutoff=2.0,
        setup_instance=intializer,
    )

@pytest.fixture
def rigid_bodies2():

    return RigidBodies(
        library="igraph",
        plddt_cutoff=70.0,
        pae_cutoff=2.0,
        setup_instance=intializer,
    )

def test_extract_rigid_bodies(
    rigid_bodies1: RigidBodies,
    rigid_bodies2: RigidBodies,
):
    """Test the extract_rigid_bodies method."""

    rigid_bodies_networkx = rigid_bodies1.extract_rigid_bodies(
        pae_matrix=rigid_bodies1.pae,
        min_res=1,
        min_proteins=1,
        plddt_filter=True,
    )

    assert isinstance(rigid_bodies_networkx, list), "Rigid bodies should be returned as a list."
    expected_rigid_bodies = [
        {'A': [('CB', 1), ('CB', 2), ('CB', 3), ('CB', 4), ('CB', 5)]},
        {'B': [('CB', 1), ('CB', 2), ('CB', 3), ('CB', 4), ('CA', 5)]}
    ]
    assert rigid_bodies_networkx == expected_rigid_bodies, "Extracted rigid bodies do not match expected output."

    rigid_bodies_igraphx = rigid_bodies2.extract_rigid_bodies(
        pae_matrix=rigid_bodies2.pae,
        min_res=1,
        min_proteins=1,
        plddt_filter=True,
    )
    print(rigid_bodies_igraphx)

    expected_rigid_bodies = [
        {'A': [('CB', 1), ('CB', 2), ('CB', 3), ('CB', 4), ('CB', 5)],
        'B': [('CB', 1), ('CB', 2), ('CB', 3), ('CB', 4), ('CA', 5)]}
    ]
    assert isinstance(rigid_bodies_igraphx, list), "Rigid bodies should be returned as a list."
    assert rigid_bodies_igraphx == expected_rigid_bodies, "Extracted rigid bodies do not match expected output."

def test_save_rigid_bodies(
    rigid_bodies1: RigidBodies,
):
    """Test the save_rigid_bodies method."""

    # Test saving with txt output
    rigid_bodies1.save_rigid_bodies(
        domains=[{'A': [('CB', 1), ('CB', 2), ('CB', 3), ('CB', 4), ('CB', 5)]}],
        output_dir="tests/test_output",
        rb_out_fmt="txt",
        save_structure=True,
        rb_struct_fmt="cif",
        filter_struct_by_plddt=True,
    )

    assert os.path.exists("tests/test_output/af3_rigid_bodies.txt"), "Rigid bodies txt file was not created."
    assert os.path.exists("tests/test_output/rigid_body_0.cif"), "Rigid body structure file was not created."

    rigid_bodies1.save_rigid_bodies(
        domains=[{'A': [('CB', 1), ('CB', 2), ('CB', 3), ('CB', 4), ('CB', 5)]}],
        output_dir="tests/test_output",
        rb_out_fmt="json",
        save_structure=True,
        rb_struct_fmt="pdb",
        filter_struct_by_plddt=False,
    )

    assert os.path.exists("tests/test_output/af3_rigid_bodies.json"), "Rigid bodies json file was not created."
    assert os.path.exists("tests/test_output/rigid_body_0.pdb"), "Rigid body structure file was not created."

def test_show_rigid_bodies_on_pae_matrix(
    rigid_bodies1: RigidBodies,
):
    """Test the show_rigid_bodies_on_pae_matrix method."""

    rigid_bodies1.show_rigid_bodies_on_pae_matrix(
        domains=[{'A': [('CB', 1), ('CB', 2), ('CB', 3), ('CB', 4), ('CB', 5)]}],
        output_dir="tests/test_output",
    )

    assert os.path.exists("tests/test_output/rigid_body_0.png"), "PAE matrix with rigid bodies was not created."

def test_assess_rigid_bodies(
    rigid_bodies1: RigidBodies,
):
    """Test the assess_rigid_bodies method."""

    rigid_bodies1.assess_rigid_bodies(
        domains=[{'A': [('CB', 1), ('CB', 2), ('CB', 3), ('CB', 4), ('CB', 5)]}],
        output_dir="tests/test_output",
        protein_chain_map={'A': 'protein1'},
        symmetric_pae=True,
        as_average=True,
    )

    assert os.path.exists("tests/test_output/rigid_body_0_assessment.xlsx"), "Rigid body assessment file was not created."