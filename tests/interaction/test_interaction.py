import os
import numpy as np
import pytest
from af_pipeline.parser.initialize import Initialize
from af_pipeline.interaction.interaction import Interaction

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
def interaction():
    return Interaction(
        contact_threshold=8.0,
        plddt_cutoff=70.0,
        pae_cutoff=2.0,
        setup_instance=intializer,
    )

def test_create_regions_of_interest(
    interaction: Interaction,
):
    """Test the create_regions_of_interest method."""

    regions_of_interest = interaction.create_regions_of_interest()

    assert isinstance(regions_of_interest, list), "Regions of interest should be returned as a list."
    expected_regions = [{'A': [1, 5], 'B': [1, 5]}]
    assert regions_of_interest == expected_regions, "Regions of interest do not match expected output."


def test_get_interaction_data(
    interaction: Interaction,
):
    """Test the Interaction class."""

    interaction_map, plddt1, plddt2, avg_pae = interaction.get_interaction_data(
        region_of_interest={"A": [1, 5], "B": [1, 5]},
    )

    assert isinstance(interaction_map, np.ndarray), "Interaction map should be a dictionary."
    assert interaction_map.shape == (5, 5), "Interaction map should have shape (5, 5)."
    expected_interaction_map = np.array([
        [1,0,0,0,0],
        [0,1,1,1,0],
        [1,1,1,1,0],
        [0,0,1,1,1],
        [0,0,1,1,1]
    ])
    assert np.array_equal(interaction_map, expected_interaction_map), "Interaction map does not match expected output."
    assert isinstance(plddt1, dict), "pLDDT values for chain 1 should be a dictionary."
    assert isinstance(plddt2, dict), "pLDDT values for chain 2 should be a dictionary."
    expected_plddt1 = {'A': np.array([95.01, 96.05, 96.94, 96.92, 97.15])}
    expected_plddt2 = {'B': np.array([95.86, 97.21, 97.5 , 97.19, 96.98])}

    for k, v in expected_plddt1.items():
        assert k in plddt1, f"Chain {k} not found in pLDDT values for chain 1."
        assert np.array_equal(plddt1[k], v), f"pLDDT values for chain {k} do not match expected output."

    for k, v in expected_plddt2.items():
        assert k in plddt2, f"Chain {k} not found in pLDDT values for chain 2."
        assert np.array_equal(plddt2[k], v), f"pLDDT values for chain {k} do not match expected output."

    assert isinstance(avg_pae, np.ndarray), "Average PAE should be a numpy array."
    assert avg_pae.shape == (5, 5), "Average PAE should have shape (5, 5)."
    expected_avg_pae = np.array([
        [2.5,2.0,2.4,2.6,3.5],
        [2.45,1.85,1.7,1.85,2.6],
        [2.95,1.8,1.5,1.6,2.45],
        [3.4,2.0,1.8,1.85,2.5],
        [3.85,2.4,2.3,2.0,2.45]
    ])
    assert np.isclose(avg_pae, expected_avg_pae).all(), "Average PAE values do not match expected output."

def test_get_confident_interaction_map(
    interaction: Interaction,
):
    """Test the get_confident_interaction_map method."""

    confident_interaction_map = interaction.get_confident_interaction_map(
        region_of_interest={"A": [1, 5], "B": [1, 5]},
    )
    assert isinstance(confident_interaction_map, np.ndarray), "Confident interaction map should be a numpy array."
    assert confident_interaction_map.shape == (5, 5), "Confident interaction map should have shape (5, 5)."
    expected_confident_interaction_map = np.array([
        [0,0,0,0,0],
        [0,1,1,1,0],
        [0,1,1,1,0],
        [0,0,1,1,0],
        [0,0,0,1,0]
    ])
    assert np.array_equal(confident_interaction_map, expected_confident_interaction_map), "Confident interaction map does not match expected output."

def test_get_interacting_patches(
    interaction: Interaction,
):
    """Test the get_interacting_patches method."""

    confident_interaction_map = interaction.get_confident_interaction_map(
        region_of_interest={"A": [1, 5], "B": [1, 5]},
    )

    interacting_patches = interaction.get_interacting_patches(
        contact_map=confident_interaction_map,
        region_of_interest={"A": [1, 5], "B": [1, 5]},
    )

    assert isinstance(interacting_patches, dict), "Interacting patches should be returned as a dictionary."
    expected_interacting_patches = {
        0: {"A": np.array([2, 3]), "B": np.array([2, 3, 4])},
        1: {"A": np.array([2, 3, 4]), "B": np.array([3, 4])},
        2: {"A": np.array([2, 3, 4, 5]), "B": np.array([4])},
    }
    for patch_id, chains in expected_interacting_patches.items():
        assert patch_id in interacting_patches, f"Patch ID {patch_id} not found in interacting patches."
        for chain_id, residues in chains.items():
            assert chain_id in interacting_patches[patch_id], f"Chain ID {chain_id} not found in patch {patch_id}."
            assert np.array_equal(interacting_patches[patch_id][chain_id], residues), f"Residues for chain {chain_id} in patch {patch_id} do not match expected output."