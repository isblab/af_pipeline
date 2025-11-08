import pytest
from af_pipeline.tools.structure_tools import RenumberResidues

offset = {"A": [10, 13], "B": [3, 5]}
region_of_interest = {"A": [10, 12]}
token_chain_ids = ["A", "A", "A", "A", "A", "B", "B", "B"]
token_res_ids = [1, 2, 3, 3, 4, 1, 2, 3]
token_atom_names = ["CB", "CB", "O", "CB", "CB", "O", "N", "C"]

@pytest.fixture
def renumber_residues1():
    return RenumberResidues()

@pytest.fixture
def renumber_residues2():
    return RenumberResidues(offset=offset)

def test_renumber_chain_res_num(
    renumber_residues1: RenumberResidues,
    renumber_residues2: RenumberResidues,
):
    chain_res_num = renumber_residues1.renumber_chain_res_num(
        chain_res_num=5,
        chain_id="A",
    )
    assert chain_res_num == 5; "Failed to renumber without offset"

    chain_res_num = renumber_residues2.renumber_chain_res_num(
        chain_res_num=5,
        chain_id="A",
    )
    assert chain_res_num == 14; "Failed to renumber with offset for chain A"

def test_renumber_region_of_interest(
    renumber_residues1: RenumberResidues,
    renumber_residues2: RenumberResidues,
):
    renumbered_roi1 = renumber_residues1.renumber_region_of_interest(
        region_of_interest=region_of_interest,
    )
    expected_roi = {"A": [10, 12]}
    assert renumbered_roi1 == expected_roi; "Failed to renumber ROI without offset"

    renumbered_roi2 = renumber_residues2.renumber_region_of_interest(
        region_of_interest=region_of_interest,
    )
    expected_roi = {"A": [1, 3]}
    assert renumbered_roi2 == expected_roi; "Failed to renumber ROI with offset for chain A"

def test_residue_map(
    renumber_residues1: RenumberResidues,
    renumber_residues2: RenumberResidues,
):
    idx_to_num, num_to_idx = renumber_residues1.residue_map(
        token_chain_ids=token_chain_ids,
        token_res_ids=token_res_ids,
        token_atom_names=token_atom_names,
    )
    expected_idx_to_num = {
        0: {"chain_id": "A", "token_num": 1, "atom_name": "CB"},
        1: {"chain_id": "A", "token_num": 2, "atom_name": "CB"},
        2: {"chain_id": "A", "token_num": 3, "atom_name": "O"},
        3: {"chain_id": "A", "token_num": 3, "atom_name": "CB"},
        4: {"chain_id": "A", "token_num": 4, "atom_name": "CB"},
        5: {"chain_id": "B", "token_num": 1, "atom_name": "O"},
        6: {"chain_id": "B", "token_num": 2, "atom_name": "N"},
        7: {"chain_id": "B", "token_num": 3, "atom_name": "C"},
    }
    expected_num_to_idx = {
        "A": {
            1: {"CB": 0},
            2: {"CB": 1},
            3: {"O": 2, "CB": 3},
            4: {"CB": 4},
        },
        "B": {
            1: {"O": 5},
            2: {"N": 6},
            3: {"C": 7}
        },
    }
    assert idx_to_num == expected_idx_to_num; "Failed to create idx_to_num mapping without offset"
    assert num_to_idx == expected_num_to_idx; "Failed to create num_to_idx mapping without offset"

    idx_to_num, num_to_idx = renumber_residues2.residue_map(
        token_chain_ids=token_chain_ids,
        token_res_ids=token_res_ids,
        token_atom_names=token_atom_names,
    )
    expected_idx_to_num = {
        0: {"chain_id": "A", "token_num": 10, "atom_name": "CB"},
        1: {"chain_id": "A", "token_num": 11, "atom_name": "CB"},
        2: {"chain_id": "A", "token_num": 12, "atom_name": "O"},
        3: {"chain_id": "A", "token_num": 12, "atom_name": "CB"},
        4: {"chain_id": "A", "token_num": 13, "atom_name": "CB"},
        5: {"chain_id": "B", "token_num": 3, "atom_name": "O"},
        6: {"chain_id": "B", "token_num": 4, "atom_name": "N"},
        7: {"chain_id": "B", "token_num": 5, "atom_name": "C"},
    }
    expected_num_to_idx = {
        "A": {
            10: {"CB": 0},
            11: {"CB": 1},
            12: {"O": 2, "CB": 3},
            13: {"CB": 4},
        },
        "B": {
            3: {"O": 5},
            4: {"N": 6},
            5: {"C": 7}
        },
    }
    assert idx_to_num == expected_idx_to_num; "Failed to create idx_to_num mapping with offset"
    assert num_to_idx == expected_num_to_idx; "Failed to create num_to_idx mapping with offset"