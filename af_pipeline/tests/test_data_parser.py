import pytest
from af_pipeline.parser1.data_parser import DataParser
import numpy as np

data_path_1 = "./tests/data/af_predictions/dummy_files/dummy_data_afdb.json"

data_path_2 = "./tests/data/af_predictions/dummy_files/dummy_data_af3_1.json"

data1 = DataParser(data_file_path=data_path_1).get_data_dict()
data2 = DataParser(data_file_path=data_path_2).get_data_dict()

@pytest.fixture
def data_parser1():
    return DataParser(
        data_file_path=data_path_1,
    )

@pytest.fixture
def data_parser2():
    return DataParser(
        data_file_path=data_path_2,
    )

def test_get_data_dict_afdb(data_parser1: DataParser):
    """Test the DataParser.get_data_dict method for AFDB data."""

    data = data_parser1.get_data_dict()
    assert isinstance(data, dict), "Data should be a Dictionary."
    assert "predicted_aligned_error" in data, \
        "Data should contain 'predicted_aligned_error' key."

def test_get_data_dict_af3(data_parser2: DataParser):
    """Test the DataParser.get_data_dict method for AF3 data."""

    data = data_parser2.get_data_dict()
    assert isinstance(data, dict), "Data should be a dictionary."
    assert "pae" in data, \
        "Data should contain 'pae' key."
    assert "token_res_ids" in data, \
        "Data should contain 'token_res_ids' key."
    assert "token_chain_ids" in data, \
        "Data should contain 'token_chain_ids' key."
    assert "contact_probs" in data, \
        "Data should contain 'contact_probs' key."

def test_get_token_chain_ids():
    """Test the DataParser.get_token_chain_ids method."""

    token_chain_ids = DataParser.get_token_chain_ids(data1)

    assert token_chain_ids is None, \
        "Token chain IDs should be None for AFDB .json files."

    token_chain_ids = DataParser.get_token_chain_ids(data2)
    assert isinstance(token_chain_ids, list), \
        "Token chain IDs should be a list for AF3 .json files."

def test_get_token_res_ids():
    """Test the DataParser.get_token_res_ids method."""

    token_res_ids = DataParser.get_token_res_ids(data1)

    assert token_res_ids is None, \
        "Token residue IDs should be None for AFDB .json files."

    token_res_ids = DataParser.get_token_res_ids(data2)
    assert isinstance(token_res_ids, list), \
        "Token residue IDs should be a list for AF3 .json files."

def test_get_pae():
    """Test the DataParser.get_pae method."""

    pae = DataParser.get_pae(data1)

    assert isinstance(pae, np.ndarray), \
        "PAE should be a numpy array for AFDB .json files."

    pae = DataParser.get_pae(data2)
    assert isinstance(pae, np.ndarray), \
        "PAE should be a numpy array for AF3 .json files."

def test_get_contact_probs_mat():
    """Test the DataParser.get_contact_probs_mat method."""

    contact_probs = DataParser.get_contact_probs_mat(data1)

    assert contact_probs is None, \
        "Contact probabilities should be None for AFDB .json files."

    contact_probs = DataParser.get_contact_probs_mat(data2)
    assert isinstance(contact_probs, np.ndarray), \
        "Contact probabilities should be a numpy array for AF3 .json files."

def test_get_atom_chain_ids():
    """Test the DataParser.get_atom_chain_ids method."""

    atom_chain_ids = DataParser.get_atom_chain_ids(data1)

    assert atom_chain_ids is None, \
        "Atom chain IDs should be None for AFDB .json files."

    atom_chain_ids = DataParser.get_atom_chain_ids(data2)
    assert isinstance(atom_chain_ids, list), \
        "Atom chain IDs should be a list for AF3 .json files."

def test_get_atom_plddts():
    """Test the DataParser.get_atom_plddts method."""

    atom_plddts = DataParser.get_atom_plddts(data1)

    assert atom_plddts is None, \
        "Atom plddts should be None for AFDB .json files."

    atom_plddts = DataParser.get_atom_plddts(data2)
    assert isinstance(atom_plddts, list), \
        "Atom plddts should be a list for AF3 .json files."
