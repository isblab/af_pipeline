import pytest
import Bio.PDB.Structure
from af_pipeline.constants.af_constants import MetricLevel
from af_pipeline.parser.structure_parser import StructureParser
from af_pipeline.parser.data_parser import DataParser
import numpy as np


struct_path1 = "tests/test_data/af_predictions/af3/fold_dummy_job_1/fold_dummy_job_1_model_0.cif"
data_path_1 = "tests/test_data/af_predictions/af3/fold_dummy_job_1/fold_dummy_job_1_full_data_0.json"

struct_path2 = "tests/test_data/af_predictions/afdb/A0A075B6S0/AF-A0A075B6S0-F1-model_v6.pdb"
data_path2 = "tests/test_data/af_predictions/afdb/A0A075B6S0/AF-A0A075B6S0-F1-predicted_aligned_error_v6.json"

struct_path3 = "dummy_path"

struct_path4 = "tests/test_data/af_predictions/af3/fold_dummy_job_ptm_1/fold_dummy_job_ptm_1_model_0.cif"
data_path_4 = "tests/test_data/af_predictions/af3/fold_dummy_job_ptm_1/fold_dummy_job_ptm_1_full_data_0.json"

struct_path5 = "tests/test_data/af_predictions/af3/fold_dummy_job_2/fold_dummy_job_2_model_0.cif"
data_path_5 = "tests/test_data/af_predictions/af3/fold_dummy_job_2/fold_dummy_job_2_full_data_0.json"

@pytest.fixture
def struct_parser1():
    return StructureParser(
        structure_file_path=struct_path1,
        preserve_header_footer=True,
    )

@pytest.fixture
def struct_parser2():
    return StructureParser(
        structure_file_path=struct_path2,
        preserve_header_footer=False,
    )

@pytest.fixture
def struct_parser3():
    return StructureParser(
        structure_file_path=struct_path3,
        preserve_header_footer=False,
    )

@pytest.fixture
def struct_parser4():
    return StructureParser(
        structure_file_path=struct_path4,
        preserve_header_footer=False,
    )

@pytest.fixture
def struct_parser5():
    return StructureParser(
        structure_file_path=struct_path5,
        preserve_header_footer=True,
    )

@pytest.fixture
def my_struct1(struct_parser1: StructureParser):
    return struct_parser1.get_structure_obj()

@pytest.fixture
def my_struct2(struct_parser2: StructureParser):
    return struct_parser2.get_structure_obj()

@pytest.fixture
def my_struct4(struct_parser4: StructureParser):
    return struct_parser4.get_structure_obj()

@pytest.fixture
def my_struct5(struct_parser5: StructureParser):
    return struct_parser5.get_structure_obj()

def test_structure_type(
    struct_parser1: StructureParser,
    struct_parser2: StructureParser,
    struct_parser3: StructureParser,
):

    assert struct_parser1.structure_type == "cif", \
        f"Structure type should be 'cif', got {struct_parser1.structure_type}."

    assert struct_parser2.structure_type == "pdb", \
        f"Structure type should be 'pdb', got {struct_parser2.structure_type}."

    with pytest.raises(Exception):
        _ = struct_parser3.structure_type

def test_get_structure(
    struct_parser1: StructureParser,
    struct_parser2: StructureParser,
):

    structure = struct_parser1.get_structure_obj()
    assert isinstance(structure, Bio.PDB.Structure.Structure), \
        f"Expected a Bio.PDB.Structure.Structure object. \
          Got {structure.__class__} instead."

    structure = struct_parser2.get_structure_obj()
    assert isinstance(structure, Bio.PDB.Structure.Structure), \
        f"Expected a Bio.PDB.Structure.Structure object. \
          Got {structure.__class__} instead."

def test_get_residues(
    my_struct1: Bio.PDB.Structure.Structure,
    my_struct2: Bio.PDB.Structure.Structure,
):

    expected_res_chain_ids1 = [
        (1, 'A'), (2, 'A'), (3, 'A'), (4, 'A'), (5, 'A'),
        (1, 'B'), (2, 'B'), (3, 'B'), (4, 'B'), (5, 'B'),
    ]
    for idx, (residue, chain_id) in enumerate(StructureParser.get_residues(my_struct1)):
        (exp_res_pos, exp_chain_id) = expected_res_chain_ids1[idx]
        assert (residue.id[1], chain_id) == (exp_res_pos, exp_chain_id), \
            f"Expected residue position and chain ID {(exp_res_pos, exp_chain_id)}, \
              got {(residue.id[1], chain_id)} instead."

    expected_res_chain_ids2 = [
        (1, 'A'), (2, 'A'), (3, 'A'), (4, 'A'), (5, 'A'),
        (6, 'A'), (7, 'A'), (8, 'A'), (9, 'A'), (10, 'A'),
        (11, 'A'), (12, 'A'), (13, 'A'), (14, 'A'), (15, 'A'),
        (16, 'A'),
    ]

    for idx, (residue, chain_id) in enumerate(StructureParser.get_residues(my_struct2)):
        (exp_res_pos, exp_chain_id) = expected_res_chain_ids2[idx]
        assert (residue.id[1], chain_id) == (exp_res_pos, exp_chain_id), \
            f"Expected residue position and chain ID {(exp_res_pos, exp_chain_id)}, \
              got {(residue.id[1], chain_id)} instead."

def test_extract_perresidue_quantities(
    my_struct1: Bio.PDB.Structure.Structure,
):

    quantities = [
        "res_pos",
        "res_name",
        "coord",
        "plddt",
        "chain_id",
        "entity_type",
        "rep_atom",
        "rep_atom_local_idx",
    ]

    residue_objs = []

    for model in my_struct1:
        for chain in model:
            for residue in chain:
                residue_objs.append(residue)

    expected_plddts = [
        95.01, 96.05, 96.94, 96.92, 97.15,
        95.86, 97.21, 97.50, 97.19, 96.98
    ]

    expected_coords = np.array([
        np.array([6.040, 5.294, 1.426,]),
        np.array([3.826, 0.732, 4.299,]),
        np.array([-0.053, 3.002, 0.605,]),
        np.array([-3.025, 0.727, 5.079,]),
        np.array([-6.065, 0.717, 0.266,]),
        np.array([4.656, 2.930,  -3.371,]),
        np.array([4.598, -2.373, -1.090,]),
        np.array([-0.779, -0.994, -1.186,]),
        np.array([-1.076, -3.751, 3.668,]),
        np.array([-5.201, -4.716, 2.310,]),
    ])

    expected_rep_atoms = ["CB"]*9 + ["CA"]

    res_idx = 0
    for residue in residue_objs:

        quantities_dict = StructureParser.extract_perresidue_quantities(
            residue=residue,
            quantities=quantities,
            rep_atom=None,
        )

        assert set(quantities_dict.keys()) == set(quantities), \
            f"Expected quantities keys {quantities}, \
              but got {list(quantities_dict.keys())} instead."

        assert quantities_dict["res_pos"] == residue.id[1], \
            f"Residue position does not match for {residue.id}."

        assert quantities_dict["res_name"] == residue.resname, \
            f"Residue name does not match for {residue.full_id}."

        assert np.allclose(
            quantities_dict["coord"],
            expected_coords[res_idx]
        ), f"Coordinates do not match for {residue.full_id}."

        assert quantities_dict["plddt"] == expected_plddts[res_idx], \
            f"Plddt does not match for {residue.full_id}."

        assert quantities_dict["chain_id"] == residue.get_parent().id[0], \
            f"Chain ID does not match for {residue.full_id}."

        assert quantities_dict["entity_type"] == "proteinChain", \
            f"Entity type does not match for {residue.full_id}."

        assert quantities_dict["rep_atom"] == expected_rep_atoms[res_idx], \
            f"Representative atom does not match for {residue.full_id}."

        res_idx += 1

def test_get_token_chain_ids(
    my_struct1: Bio.PDB.Structure.Structure,
    my_struct2: Bio.PDB.Structure.Structure,
    my_struct4: Bio.PDB.Structure.Structure,
):
    """Test the StructureParser.get_token_chain_ids method."""

    token_chain_ids1 = StructureParser.get_token_chain_ids(
        structure=my_struct1,
        rep_atom_dict={},
        metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
    )
    expected_chain_ids1 = ['A'] * 5 + ['B'] * 5

    assert token_chain_ids1 == expected_chain_ids1, \
        "Token chain IDs do not match."

    token_chain_ids2 = StructureParser.get_token_chain_ids(
        structure=my_struct2,
        rep_atom_dict={},
        metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
    )
    expected_chain_ids2 = ['A'] * 16

    assert token_chain_ids2 == expected_chain_ids2, \
        "Token chain IDs do not match."

    token_chain_ids4 = StructureParser.get_token_chain_ids(
        structure=my_struct4,
        rep_atom_dict={},
        metric_level=MetricLevel.PER_TOKEN,
    )
    expected_chain_ids1 = ['A'] * 5 + ['B'] * 14

    assert token_chain_ids4 == expected_chain_ids1, \
        "Token chain IDs do not match with `only_representative` set to `False`."

def test_get_token_res_ids(
    my_struct1: Bio.PDB.Structure.Structure,
    my_struct2: Bio.PDB.Structure.Structure,
    my_struct4: Bio.PDB.Structure.Structure,
):
    """Test the StructureParser.get_token_res_ids method."""

    token_res_ids1 = StructureParser.get_token_res_ids(
        structure=my_struct1,
        rep_atom_dict={},
        metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
    )
    expected_res_ids1 = [
        1, 2, 3, 4, 5,
        1, 2, 3, 4, 5,
    ]

    assert token_res_ids1 == expected_res_ids1, \
        "Token residue IDs do not match."

    token_res_ids2 = StructureParser.get_token_res_ids(
        structure=my_struct2,
        rep_atom_dict={},
        metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
    )
    expected_res_ids2 = list(range(1, 17))

    assert token_res_ids2 == expected_res_ids2, \
        "Token residue IDs do not match."

    token_res_ids4 = StructureParser.get_token_res_ids(
        structure=my_struct4,
        rep_atom_dict={},
        metric_level=MetricLevel.PER_TOKEN,
    )
    expected_res_ids1 = [
        1, 2, 3, 4, 5,
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 5,
    ]

    assert token_res_ids4 == expected_res_ids1, \
        "Token res IDs do not match with `only_representative` set to `False`."

def test_get_plddt(
    my_struct1: Bio.PDB.Structure.Structure,
    my_struct2: Bio.PDB.Structure.Structure,
    my_struct4: Bio.PDB.Structure.Structure,
    my_struct5: Bio.PDB.Structure.Structure,
):
    """Test the StructureParser.get_plddt method."""

    plddts1 = StructureParser.get_plddt(
        structure=my_struct1,
        rep_atom_dict={},
        average_token_plddt=False,
        metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
    )
    expected_plddts = [
        95.01, 96.05, 96.94, 96.92, 97.15, # chain A
        95.86, 97.21, 97.50, 97.19, 96.98, # chain B
    ]

    assert plddts1 == expected_plddts, "pLDDT values do not match."

    plddts2 = StructureParser.get_plddt(
        structure=my_struct2,
        rep_atom_dict={},
        average_token_plddt=False,
        metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
    )
    expected_plddts = [
        56.72, 63.31, 67.69, 68.88, 73.88, 74.94, 76.44, 74.00, 73.62, 76.00,
        84.44, 84.38, 87.25, 89.06, 87.69, 88.06
    ]

    assert plddts2 == expected_plddts, "pLDDT values do not match."

    plddts4 = StructureParser.get_plddt(
        structure=my_struct4,
        rep_atom_dict={},
        average_token_plddt=False,
        metric_level=MetricLevel.PER_TOKEN,
    )
    expected_plddts = [
        87.26, 90.17, 94.25, 93.06, 93.43,
        94.85, 97.08, 96.15, 96.90, 96.47, 97.58, 95.97, 96.06, 92.66, 92.04, 92.49, 96.45, 94.35, 93.57,
    ]

    assert plddts4 == expected_plddts, \
        "pLDDT values do not match with `only_representative` set to `False`."

    plddts5 = StructureParser.get_plddt(
        structure=my_struct5,
        rep_atom_dict={},
        average_token_plddt=False,
        metric_level=MetricLevel.PER_ATOM,
    )
    expected_plddts = [
        95.39, 96.60, 96.72, 95.02, 95.09, 91.42, 89.83, 80.08, 93.22, 94.31,
        94.15, 92.33, 93.96, 88.29, 85.01, 79.72, 79.79, 94.47, 94.25, 94.29,
        92.40, 93.20, 87.98, 89.30, 95.52, 94.53, 92.01, 89.52, 92.43, 91.68,
        91.81, 84.87, 88.50,
    ]
    assert plddts5 == expected_plddts, \
        "pLDDT values do not match with `per_atom` set to `True`."

    plddts5 = StructureParser.get_plddt(
        structure=my_struct5,
        rep_atom_dict={},
        average_token_plddt=True,
        metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
    )
    expected_plddts = [
        np.mean([95.39, 96.60, 96.72, 95.02, 95.09, 91.42, 89.83, 80.08,]),
        np.mean([93.22, 94.31, 94.15, 92.33, 93.96, 88.29, 85.01, 79.72, 79.79,]),
        np.mean([94.47, 94.25, 94.29, 92.40, 93.20, 87.98, 89.30,]),
        np.mean([95.52, 94.53, 92.01, 89.52, 92.43, 91.68, 91.81, 84.87, 88.50,]),
    ]

    assert plddts5 == expected_plddts, \
        "pLDDT values do not match with `only_representative` and \
         `average_token_plddt` set to `True`."

    plddts4 = StructureParser.get_plddt(
        structure=my_struct4,
        rep_atom_dict={"SEP": "P"},
        average_token_plddt=False,
        metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
    )
    expected_plddts = [
        87.26, 90.17, 94.25, 93.06, 93.43,
        94.85, 96.06, 96.45, 94.35, 93.57,
    ]
    assert plddts4 == expected_plddts, \
        "pLDDT vals do not match with `only_representative` set to `True` and \
         `rep_atom_dict` set to include 'P' as representative atom in 'SEP'."

    with pytest.raises(Exception):

        _ = StructureParser.get_plddt(
            structure=my_struct5,
            per_atom=False,
            rep_atom_dict={},
            average_token_plddt=True,
            only_representative=None,
        )

    with pytest.raises(Exception):

        _ = StructureParser.get_plddt(
            structure=my_struct4,
            per_atom=False,
            rep_atom_dict={"SEP": "A"},
            average_token_plddt=False,
            only_representative=True,
        )

    with pytest.raises(Exception):

        _ = StructureParser.get_plddt(
            structure=my_struct4,
            per_atom=False,
            rep_atom_dict={"SEP": 123},
            average_token_plddt=False,
            only_representative=True,
        )

def test_get_coordinates(
    my_struct1: Bio.PDB.Structure.Structure,
    my_struct4: Bio.PDB.Structure.Structure,
    my_struct5: Bio.PDB.Structure.Structure,
):
    """Test the StructureParser.get_coordinates method."""

    coords1 = StructureParser.get_coordinates(
        structure=my_struct1,
        rep_atom_dict={},
        metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
    )
    expected_coords = np.array([
        np.array([6.040, 5.294, 1.426,]),
        np.array([3.826, 0.732, 4.299,]),
        np.array([-0.053, 3.002, 0.605,]),
        np.array([-3.025, 0.727, 5.079,]),
        np.array([-6.065, 0.717, 0.266,]),
        np.array([4.656, 2.930,  -3.371,]),
        np.array([4.598, -2.373, -1.090,]),
        np.array([-0.779, -0.994, -1.186,]),
        np.array([-1.076, -3.751, 3.668,]),
        np.array([-5.201, -4.716, 2.310,]),
    ])
    assert np.allclose(coords1, expected_coords), "Coordinates do not match."

    coords4 = StructureParser.get_coordinates(
        structure=my_struct4,
        rep_atom_dict={},
        metric_level=MetricLevel.PER_TOKEN,
    )
    expected_coords = np.array([
        np.array([-0.016, 7.704, 3.096 ]),
        np.array([3.279, 3.212, 1.945 ]),
        np.array([-1.516, 4.197, -1.290]),
        np.array([1.778, -0.932, -1.844]),
        np.array([-1.768, -3.591, 0.090 ]),
        np.array([4.856, 0.195, 1.727 ]),
        np.array([2.680, -0.572, 3.945 ]),
        np.array([1.491, -0.429, 4.654 ]),
        np.array([1.912, -0.261, 6.066 ]),
        np.array([2.288, -1.518, 6.583 ]),
        np.array([0.697, 0.773, 4.157 ]),
        np.array([1.304, 1.847, 3.979 ]),
        np.array([3.869, -1.856, 6.647 ]),
        np.array([4.247, -2.518, 5.391 ]),
        np.array([3.983, -2.806, 7.837 ]),
        np.array([4.640, -0.580, 6.889 ]),
        np.array([-2.466, 1.206, 2.522 ]),
        np.array([-1.804, 5.402, 6.287 ]),
        np.array([-5.142, 6.511, 3.441 ]),
    ])
    assert np.allclose(coords4, expected_coords), \
        "Coordinates do not match for struct_parser4 with \
         `only_representative` set to `False`."

    coords5 = StructureParser.get_coordinates(
        structure=my_struct5,
        rep_atom_dict={},
        metric_level=MetricLevel.PER_ATOM,
    )
    expected_coords = np.array([
        np.array([-2.745, -1.474, 7.791]),
        np.array([-1.840, -2.199, 6.874]),
        np.array([-1.342, -1.251, 5.796]),
        np.array([-2.095, -0.402, 5.309]),
        np.array([-2.563, -3.379, 6.221]),
        np.array([-1.677, -4.182, 5.282]),
        np.array([-2.525, -5.587, 4.554]),
        np.array([-2.814, -6.606, 6.005]),
        np.array([-0.074, -1.363, 5.447]),
        np.array([0.511, -0.557, 4.376]),
        np.array([1.170, -1.466, 3.350]),
        np.array([2.041, -2.266, 3.692]),
        np.array([1.529, 0.414 , 4.961]),
        np.array([0.891, 1.427 , 5.888]),
        np.array([1.900, 2.399 , 6.450]),
        np.array([3.115, 2.234 , 6.214]),
        np.array([1.476, 3.329 , 7.156]),
        np.array([0.736, -1.348, 2.140]),
        np.array([1.282, -2.162, 1.053]),
        np.array([1.772, -1.233, -0.052]),
        np.array([1.006, -0.414, -0.550]),
        np.array([0.229, -3.118, 0.492]),
        np.array([0.821, -3.964, -0.624]),
        np.array([-0.316, -4.007, 1.611]),
        np.array([3.019, -1.397, -0.413]),
        np.array([3.611, -0.541, -1.430]),
        np.array([3.986, -1.350, -2.675]),
        np.array([4.207, -2.566, -2.575]),
        np.array([4.867, 0.146, -0.891]),
        np.array([4.612, 1.013, 0.328]),
        np.array([3.736, 2.530, -0.107]),
        np.array([2.271, 2.354, 0.872]),
        np.array([4.146, -0.754, -3.760]),
    ])
    assert np.allclose(coords5, expected_coords), \
        "Coordinates do not match with `per_atom` set to `True`."

def test_get_token_atom_names(
    my_struct1: Bio.PDB.Structure.Structure,
    my_struct4: Bio.PDB.Structure.Structure,
):
    """Test the StructureParser.get_token_atom_names method."""

    atom_names1 = StructureParser.get_token_atom_names(
        structure=my_struct1,
        rep_atom_dict={},
        metric_level=MetricLevel.REPRESENTATIVE_TOKEN,
    )
    expected_atom_names1 = ["CB"]*9 + ["CA"]
    assert atom_names1 == expected_atom_names1, "Token atom names do not match."

    atom_names4 = StructureParser.get_token_atom_names(
        structure=my_struct4,
        rep_atom_dict={},
        metric_level=MetricLevel.PER_TOKEN,
    )
    expected_atom_names4 = ["CB"]*6 + [
        "N", "CA", "CB", "OG", "C", "O", "P", "O1P", "O2P", "O3P",
    ] + ["CB"]*2 + ["CA"]
    assert atom_names4 == expected_atom_names4, \
        "Token atom names do not match with `only_representative` set to `False`."
