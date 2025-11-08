from af_pipeline.tools.structure_tools import RenumberResidues
from af_pipeline.parser1.structure_parser1 import StructureParser
import pytest

struct_file1 = "./tests/data/af_predictions/af3/dp1_2_1021to1950_4/fold_dp1_2_1021to1950_4_model_0.cif"

af_offset={
    "A": [1021, 1950],
    "B": [1021, 1950],
}

@pytest.fixture
def renumber():
    return RenumberResidues(
        af_offset=af_offset,
    )

def test_renumber_structure(renumber: RenumberResidues):
    struct_parser = StructureParser(
        struct_file_path=struct_file1,
    )
    my_struct = struct_parser.get_structure(
        parser=struct_parser.get_parser(),
    )
    renumbered_struct = renumber.renumber_structure(
        structure=my_struct,
    )

    for model in renumbered_struct:
        for chain in model:
            chain_id = chain.id[0]
            start, end = af_offset[chain_id]
            renumbered_res = list(range(start, end + 1))
            res_idx = 0
            for residue in chain:
                assert residue.id[1] == renumbered_res[res_idx], \
                    f"Residue {residue.id[1]} in chain {chain_id} is not renumbered correctly."
                res_idx += 1

def test_renumber_chain_res_num(renumber: RenumberResidues):
    chain_id = "A"
    chain_res_num = 1
    renumbered_res_num = renumber.renumber_chain_res_num(
        chain_res_num=chain_res_num,
        chain_id=chain_id,
    )
    assert renumbered_res_num == af_offset[chain_id][0], \
        f"Renumbered residue number {renumbered_res_num} does not match expected {af_offset[chain_id][0]} for chain {chain_id}."

def test_renumber_region_of_interest(renumber: RenumberResidues):
    region_of_interest = {
        "A": [1021, 1031],
    }
    renumbered_region = renumber.renumber_region_of_interest(
        region_of_interest=region_of_interest,
    )
    expected_region = {
        "A": [1, 11],
    }
    assert renumbered_region == expected_region, \
        f"Renumbered region {renumbered_region} does not match expected {expected_region}."

def test_residue_map():
    token_chain_ids = ["A", "A", "A", "B", "B", "C"]
    token_res_ids = [1, 1, 2, 1, 1, 1]
    af_offset = {
        "A": [1021, 1022],
    }
    renumber = RenumberResidues(
        af_offset=af_offset,
    )
    idx_to_num, num_to_idx = renumber.residue_map(
        token_chain_ids=token_chain_ids,
        token_res_ids=token_res_ids,
    )
    expected_idx_to_num = {
        0: {"chain_id": "A", "token_num": 1021},
        1: {"chain_id": "A", "token_num": 1021},
        2: {"chain_id": "A", "token_num": 1022},
        3: {"chain_id": "B", "token_num": 1},
        4: {"chain_id": "B", "token_num": 1},
        5: {"chain_id": "C", "token_num": 1},
    }
    expected_num_to_idx = {
        "A": {1021: [0, 1], 1022:[2]},
        "B": {1: [3, 4]},
        "C": {1: [5]},
    }
    print("Index to number mapping:", idx_to_num)
    print("Number to index mapping:", num_to_idx)
    assert idx_to_num == expected_idx_to_num, \
        f"Index to number mapping does not match expected."
    assert num_to_idx == expected_num_to_idx, \
        f"Number to index mapping does not match expected."