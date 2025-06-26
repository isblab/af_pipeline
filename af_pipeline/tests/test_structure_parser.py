import pytest
from Bio.PDB import PDBParser, MMCIFParser
import Bio.PDB.Structure
from af_pipeline.constants.af_constants import *
from af_pipeline.parser1.structure_parser import StructureParser
from af_pipeline.parser1.data_parser import DataParser
import numpy as np

struct_path1 = "./tests/data/af_predictions/af3/lb2cas12a_rna_dna_complex_1/fold_lb2cas12a_rna_dna_complex_1_model_0.cif"

data_path_1 = "./tests/data/af_predictions/af3/lb2cas12a_rna_dna_complex_1/fold_lb2cas12a_rna_dna_complex_1_full_data_0.json"

struct_path2 = "./tests/data/af_predictions/afdb/AF-Q49B96-F1-model_v4.pdb"

data1 = DataParser(data_file_path=data_path_1).get_data_dict()

struct_parser1 = StructureParser(
    struct_file_path=struct_path1
)

my_struct1 = struct_parser1.get_structure(
    parser = MMCIFParser()
)

struct_parser2 = StructureParser(
    struct_file_path=struct_path2,
)

my_struct2 = struct_parser2.get_structure(
    parser = PDBParser()
)

def test_get_structure():

    parser = MMCIFParser(QUIET=True)

    structure = struct_parser1.get_structure(
        parser=parser,
    )
    assert isinstance(structure, Bio.PDB.Structure.Structure), \
        "Expected a Bio.PDB.MMCIFParser object."

    parser = PDBParser(QUIET=True)

    structure = struct_parser2.get_structure(
        parser=parser,
    )
    assert isinstance(structure, Bio.PDB.Structure.Structure), \
        "Expected a Bio.PDB.PDBParser object."

def test_extract_perresidue_quantity1():

    struct_parser1.only_representative = False

    quantities = ["res_pos", "coords", "plddt"]

    for model in my_struct1:
        for chain in model:
            chain_id = chain.id[0]
            for residue in chain:

                result_respos, _ = struct_parser1.extract_perresidue_quantity(
                    residue=residue,
                    quantity=quantities[0],
                )
                result_coords, _ = struct_parser1.extract_perresidue_quantity(
                    residue=residue,
                    quantity=quantities[1],
                )
                result_plddt, _ = struct_parser1.extract_perresidue_quantity(
                    residue=residue,
                    quantity=quantities[2],
                )

                _condition = (
                    residue.xtra.get("is_modified", False) == True
                    or residue.xtra.get("entityType", None) is None
                    or residue.xtra.get("is_ligand", False) == True
                    or residue.xtra.get("is_ion", False) == True
                )

                if _condition:
                    for idx, atom in enumerate(residue):
                        # print(result_respos[idx], residue.id[1])
                        # print(result_coords[idx], atom.coord)
                        # print(result_plddt[idx], atom.bfactor)

                        assert result_respos[idx] == residue.id[1], \
                            f"Residue position mismatch: \
                            {result_respos[idx]} != {residue.id[1]} \
                            for residue {residue.get_resname()} \
                            in chain {chain_id}"

                        assert np.allclose(
                            result_coords[idx],
                            atom.coord,
                            rtol=1e-2,
                            atol=1e-2
                        ), f"Coordinates mismatch for atom {atom.name}: \
                            {result_coords[idx]} != {atom.coord}"

                        assert result_plddt[idx] == atom.bfactor, \
                            f"Plddt mismatch for atom {atom.name}: \
                            {result_plddt[idx]} != {atom.bfactor}"

                else:
                    # print(result_respos[0], residue.id[1])
                    # print(result_coords[0], residue[_[0]].coord)
                    # print(result_plddt[0], residue[_[0]].bfactor)

                    assert len(result_respos) == 1, \
                        f"Expected 1 residue position, got {len(result_respos)}"
                    assert len(result_coords) == 1, \
                        f"Expected 1 coordinate, got {len(result_coords)}"
                    assert len(result_plddt) == 1, \
                        f"Expected 1 plddt value, got {len(result_plddt)}"

                    assert result_respos[0] == residue.id[1], \
                        f"Residue position mismatch: \
                        {result_respos} != {residue.id[1]} \
                        for residue {residue.get_resname()} \
                        in chain {chain_id}"

                    assert np.allclose(
                        result_coords[0],
                        residue[_[0]].coord,
                        rtol=1e-2,
                        atol=1e-2
                    )

                    assert result_plddt[0] == residue[_[0]].bfactor, \
                        f"Plddt mismatch for atom {_[0]}: \
                        {result_plddt[0]} != {residue[_[0]].bfactor}"

def test_extract_perresidue_quantity2():

    struct_parser1.only_representative = True

    quantities = ["res_pos", "coords", "plddt"]

    for model in my_struct1:
        for chain in model:
            chain_id = chain.id[0]
            for residue in chain:

                result_respos, _ = struct_parser1.extract_perresidue_quantity(
                    residue=residue,
                    quantity=quantities[0],
                )
                result_coords, _ = struct_parser1.extract_perresidue_quantity(
                    residue=residue,
                    quantity=quantities[1],
                )
                result_plddt, _ = struct_parser1.extract_perresidue_quantity(
                    residue=residue,
                    quantity=quantities[2],
                )

                assert len(result_respos) == 1, \
                    f"Expected 1 residue position, got {len(result_respos)}"
                assert len(result_coords) == 1, \
                    f"Expected 1 coordinate, got {len(result_coords)}"
                assert len(result_plddt) == 1, \
                    f"Expected 1 plddt value, got {len(result_plddt)}"

                assert result_respos[0] == residue.id[1], \
                    f"Residue position mismatch: \
                    {result_respos} != {residue.id[1]} \
                    for residue {residue.get_resname()} \
                    in chain {chain_id}"

                assert np.allclose(
                    result_coords[0],
                    residue[_[0]].coord,
                    rtol=1e-2,
                    atol=1e-2
                )

                assert result_plddt[0] == residue[_[0]].bfactor, \
                    f"Plddt mismatch for atom {_[0]}: \
                    {result_plddt[0]} != {residue[_[0]].bfactor}"

def test_get_token_chain_res_ids():

    struct_parser1.only_representative = False
    token_chain_ids, token_res_ids = struct_parser1.get_token_chain_res_ids(
        structure=my_struct1,
    )

    expected_token_chain_ids = data1["token_chain_ids"]
    expected_token_res_ids = data1["token_res_ids"]

    assert token_chain_ids == expected_token_chain_ids, \
        f"Token chain IDs mismatch: \
        {token_chain_ids} != {expected_token_chain_ids}"

    assert token_res_ids == expected_token_res_ids, \
        f"Token residue IDs mismatch: \
        {token_res_ids} != {expected_token_res_ids}"

def test_get_rep_coordinates():

    coords_list = struct_parser1.get_rep_coordinates(
        structure=my_struct1,
    )
    assert isinstance(coords_list, list) and len(coords_list) > 0, \
        "Expected a list of coordinates."

def get_rep_plddt():

    plddt_list = struct_parser1.get_rep_plddt(
        structure=my_struct1,
    )
    assert isinstance(plddt_list, list) and len(plddt_list) > 0, \
        "Expected a list of plddt values."