import os
import pytest
import tempfile
from af_pipeline.af_input.alphafold3 import (
    AFServerConfig, AFJobSet, AFSequence, Entity
)
from af_pipeline.constants.af_constants import (
    MAX_TEMPLATE_DATE,
    RES_RANGE_SEP,
    ConfigYaml,
    AFInputJobFields,
    AFInputEntityFields,
    EntityType,
)

test_out_dir = tempfile.mkdtemp()

###############################################################################
protein_sequences_by_name = {
    "protA": "MKTAYIAKQRQISFVKSHFSRQDILDLI",
    "protB": "GAVLILLLVAVAVVAGVAA",
}

entities_map = {
    "protA": "A12345",
    "protB": "B67890",
}

protein_sequences_by_id = {
    "A12345": "MKTAYIAKQRQISFVKSHFSRQDILDLI",
    "B67890": "GAVLILLLVAVAVVAGVAA",
}

nucleic_acid_sequences = {
    "dnaA": "ATGCGTACGTAGCTAGCTAGCTAGCTA",
    "rnaA": "AUGCGUACGUAGCUAGCUAGCUAGCUA",
}
###############################################################################
entity_info_protein_min = {
    "name": "protB",
    "type": EntityType.PROTEIN_CHAIN,
}

entity_info_protein_max = {
    AFInputEntityFields.NAME: "protA",
    AFInputEntityFields.TYPE: EntityType.PROTEIN_CHAIN,
    AFInputEntityFields.COUNT: 2,
    AFInputEntityFields.RANGE: [2, 25],
    AFInputEntityFields.USE_STRUCTURE_TEMPLATE: True,
    AFInputEntityFields.MAX_TEMPLATE_DATE: "2023-01-01",
    AFInputEntityFields.GLYCANS: [["BMA", 5]],
    AFInputEntityFields.MODIFICATIONS: [["CCD_HY3", 11]],
}

entity_info_dna_min = {
    AFInputEntityFields.NAME: "dnaA",
    AFInputEntityFields.TYPE: EntityType.DNA_SEQUENCE,
}

entity_info_dna_max = {
    AFInputEntityFields.NAME: "dnaA",
    AFInputEntityFields.TYPE: EntityType.DNA_SEQUENCE,
    AFInputEntityFields.COUNT: 1,
    AFInputEntityFields.RANGE: [3, 20],
    AFInputEntityFields.MODIFICATIONS: [["CCD_6OG", 3], ["CCD_6MA", 5]],
}

entity_info_rna_min = {
    AFInputEntityFields.NAME: "rnaA",
    AFInputEntityFields.TYPE: EntityType.RNA_SEQUENCE,
}

entity_info_rna_max = {
    AFInputEntityFields.NAME: "rnaA",
    AFInputEntityFields.TYPE: EntityType.RNA_SEQUENCE,
    AFInputEntityFields.COUNT: 1,
    AFInputEntityFields.RANGE: [3, 20],
    AFInputEntityFields.MODIFICATIONS: [["CCD_5MC", 4], ["CCD_5MU", 6]],
}

entity_info_ligand = {
    AFInputEntityFields.NAME: "CCD_ATP",
    AFInputEntityFields.TYPE: EntityType.LIGAND,
    AFInputEntityFields.COUNT: 2,
}

entity_info_ion = {
    AFInputEntityFields.NAME: "MG",
    AFInputEntityFields.TYPE: EntityType.ION,
    AFInputEntityFields.COUNT: 2,
}
###############################################################################
job_set_info_min = {
    AFInputJobFields.JOB_SET_NAME: "jobset_min",
    AFInputJobFields.MODEL_SEEDS: [0, 1],
    AFInputJobFields.ENTITIES: [
        entity_info_protein_min,
        entity_info_dna_min,
        entity_info_rna_min,
        entity_info_ligand,
        entity_info_ion,
    ]
}

job_set_info_min_no_name = {
    AFInputJobFields.MODEL_SEEDS: [0, 1],
    AFInputJobFields.ENTITIES: [
        entity_info_protein_min,
        entity_info_dna_min,
        entity_info_rna_min,
        entity_info_ligand,
        entity_info_ion,
    ]
}

job_set_info_max = {
    AFInputJobFields.JOB_SET_NAME: "jobset_max",
    AFInputJobFields.MODEL_SEEDS: [0, 1],
    AFInputJobFields.ENTITIES: [
        entity_info_protein_max,
        entity_info_dna_max,
        entity_info_rna_max,
        entity_info_ligand,
        entity_info_ion,
    ]
}

job_set_info_max_no_name = {
    AFInputJobFields.MODEL_SEEDS: [0, 1],
    AFInputJobFields.ENTITIES: [
        entity_info_protein_max,
        entity_info_dna_max,
        entity_info_rna_max,
        entity_info_ligand,
        entity_info_ion,
    ]
}
###############################################################################
job_sets_list = [
    job_set_info_min,
    job_set_info_max,
]
###############################################################################
input_job_sets = [job_sets_list[0]]
config_dict = {
    ConfigYaml.AF_INPUT_JOBS: input_job_sets,
    ConfigYaml.PROTEIN_UNIPROT_MAP: entities_map,
}

###############################################################################
# Entity class Fixtures and Tests
###############################################################################
@pytest.fixture
def pentity1():
    return Entity(
        entity_info=entity_info_protein_min,
        protein_sequences=protein_sequences_by_name,
    )

@pytest.fixture
def pentity2():
    return Entity(
        entity_info=entity_info_protein_max,
        protein_sequences=protein_sequences_by_name,
    )

@pytest.fixture
def pentity3():
    return Entity(
        entity_info=entity_info_protein_max,
        protein_sequences=protein_sequences_by_id,
        entities_map=entities_map,
    )

@pytest.fixture
def dentity1():
    return Entity(
        entity_info=entity_info_dna_min,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

@pytest.fixture
def dentity2():
    return Entity(
        entity_info=entity_info_dna_max,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

@pytest.fixture
def rentity1():
    return Entity(
        entity_info=entity_info_rna_min,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

@pytest.fixture
def rentity2():
    return Entity(
        entity_info=entity_info_rna_max,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

@pytest.fixture
def lentity():
    return Entity(
        entity_info=entity_info_ligand,
    )

@pytest.fixture
def ientity():
    return Entity(
        entity_info=entity_info_ion,
    )

###############################################################################

def test_get_template_settings(
    pentity1: Entity, pentity2: Entity, pentity3: Entity,
    dentity1: Entity, dentity2: Entity,
    rentity1: Entity, rentity2: Entity,
    lentity: Entity, ientity: Entity,
):
    """Test Entity.get_template_settings method for proteinChain entity."""

    assert_msg = "Template settings do not match expected values."

    # Test for proteinChain entity
    template_settings = pentity1.get_template_settings()
    expected_settings = {
        "useStructureTemplate": True,
        "maxTemplateDate": MAX_TEMPLATE_DATE,
        "templates": [],
    }
    assert template_settings == expected_settings, assert_msg

    template_settings = pentity2.get_template_settings()
    expected_settings = {
        "useStructureTemplate": True,
        "maxTemplateDate": "2023-01-01",
        "templates": [],
    }
    assert template_settings == expected_settings, assert_msg

    template_settings = pentity3.get_template_settings()
    assert template_settings == expected_settings, assert_msg

    # Test for dnaSequence entity
    template_settings = dentity1.get_template_settings()
    expected_settings = {}
    assert template_settings == expected_settings, assert_msg

    template_settings = dentity2.get_template_settings()
    assert template_settings == expected_settings, assert_msg

    # Test for rnaSequence entity
    template_settings = rentity1.get_template_settings()
    assert template_settings == expected_settings, assert_msg

    template_settings = rentity2.get_template_settings()
    assert template_settings == expected_settings, assert_msg

    # Test for ligand entity
    template_settings = lentity.get_template_settings()
    assert template_settings == expected_settings, assert_msg

    # Test for ion entity
    template_settings = ientity.get_template_settings()
    assert template_settings == expected_settings, assert_msg

def test_get_real_sequence(
    pentity1: Entity, pentity2: Entity, pentity3: Entity,
    dentity1: Entity, dentity2: Entity,
    rentity1: Entity, rentity2: Entity,
    lentity: Entity, ientity: Entity,
):
    """Test Entity.get_real_sequence method for proteinChain entity."""

    assert_msg = "Real sequence does not match expected sequence."

    # Test for proteinChain entity
    real_sequence = pentity1.get_real_sequence()
    expected_sequence = "GAVLILLLVAVAVVAGVAA"
    assert real_sequence == expected_sequence, assert_msg

    real_sequence = pentity2.get_real_sequence()
    expected_sequence = "MKTAYIAKQRQISFVKSHFSRQDILDLI"
    assert real_sequence == expected_sequence, assert_msg

    real_sequence = pentity3.get_real_sequence()
    assert real_sequence == expected_sequence, assert_msg

    # Test for dnaSequence entity
    real_sequence = dentity1.get_real_sequence()
    expected_sequence = "ATGCGTACGTAGCTAGCTAGCTAGCTA"
    assert real_sequence == expected_sequence, assert_msg

    real_sequence = dentity2.get_real_sequence()
    assert real_sequence == expected_sequence, assert_msg

    # Test for rnaSequence entity
    real_sequence = rentity1.get_real_sequence()
    expected_sequence = "AUGCGUACGUAGCUAGCUAGCUAGCUA"
    assert real_sequence == expected_sequence, assert_msg

    real_sequence = rentity2.get_real_sequence()
    assert real_sequence == expected_sequence, assert_msg

    # Test for ligand entity
    real_sequence = lentity.get_real_sequence()
    expected_sequence = ""
    assert real_sequence == expected_sequence, assert_msg

    # Test for ion entity
    real_sequence = ientity.get_real_sequence()
    assert real_sequence == expected_sequence, assert_msg

def test_get_glycans(
    pentity1: Entity, pentity2: Entity,
    dentity1: Entity,
    rentity1: Entity,
    lentity: Entity, ientity: Entity,
):
    """Test Entity.get_glycans method for proteinChain entity."""

    assert_msg = "Glycans do not match expected glycans."

    # Test for proteinChain entity
    glycans = pentity1.get_glycans()
    expected_glycans = []
    assert glycans == expected_glycans, assert_msg

    glycans = pentity2.get_glycans()
    expected_glycans = [{
        "residues": "BMA",
        "position": 5 - 2 + 1,
    }]
    assert glycans == expected_glycans, assert_msg

    # Test for dnaSequence entity
    glycans = dentity1.get_glycans()
    expected_glycans = []
    assert glycans == expected_glycans, assert_msg

    # Test for rnaSequence entity
    glycans = rentity1.get_glycans()
    assert glycans == expected_glycans, assert_msg

    # Test for ligand entity
    glycans = lentity.get_glycans()
    assert glycans == expected_glycans, assert_msg

    # Test for ion entity
    glycans = ientity.get_glycans()
    assert glycans == expected_glycans, assert_msg

def test_get_modifications_protein(
    pentity1: Entity, pentity2: Entity,
    dentity1: Entity, dentity2: Entity,
    rentity1: Entity, rentity2: Entity,
    lentity: Entity, ientity: Entity,
):
    """Test Entity.get_modifications method for proteinChain entity."""

    assert_msg = "Modifications do not match expected modifications."

    # Test for proteinChain entity
    modifications = pentity1.get_modifications()
    expected_modifications = []
    assert modifications == expected_modifications, assert_msg

    modifications = pentity2.get_modifications()
    expected_modifications = [{
        "ptmType": "CCD_HY3",
        "ptmPosition": 11 - 2 + 1,
    }]
    assert modifications == expected_modifications, assert_msg

    # Test for dnaSequence entity
    modifications = dentity1.get_modifications()
    expected_modifications = []
    assert modifications == expected_modifications, assert_msg

    modifications = dentity2.get_modifications()
    expected_modifications = [
        {
            "modificationType": "CCD_6OG",
            "basePosition": 3 - 3 + 1,
        },
        {
            "modificationType": "CCD_6MA",
            "basePosition": 5 - 3 + 1,
        },
    ]
    assert modifications == expected_modifications, assert_msg

    # Test for rnaSequence entity
    modifications = rentity1.get_modifications()
    expected_modifications = []
    assert modifications == expected_modifications, assert_msg

    modifications = rentity2.get_modifications()
    expected_modifications = [
        {
            "modificationType": "CCD_5MC",
            "basePosition": 4 - 3 + 1,
        },
        {
            "modificationType": "CCD_5MU",
            "basePosition": 6 - 3 + 1,
        },
    ]
    assert modifications == expected_modifications, assert_msg

    # Test for ligand entity
    modifications = lentity.get_modifications()
    expected_modifications = []
    assert modifications == expected_modifications, assert_msg

    # Test for ion entity
    modifications = ientity.get_modifications()
    assert modifications == expected_modifications, assert_msg

def test_get_entity_range(
    pentity1: Entity, pentity2: Entity,
    dentity1: Entity, dentity2: Entity,
    rentity1: Entity, rentity2: Entity,
    lentity: Entity, ientity: Entity,
):
    """Test Entity.get_entity_range method for proteinChain entity."""

    assert_msg = "Entity range does not match expected range."

    # Test for proteinChain entity
    entity_range = pentity1.get_entity_range()
    expected_range = (1, 19)
    assert entity_range == expected_range, assert_msg

    entity_range = pentity2.get_entity_range()
    expected_range = (2, 25)
    assert entity_range == expected_range, assert_msg

    # Test for dnaSequence entity
    entity_range = dentity1.get_entity_range()
    expected_range = (1, 27)
    assert entity_range == expected_range, assert_msg

    entity_range = dentity2.get_entity_range()
    expected_range = (3, 20)
    assert entity_range == expected_range, assert_msg

    # Test for rnaSequence entity
    entity_range = rentity1.get_entity_range()
    expected_range = (1, 27)
    assert entity_range == expected_range, assert_msg

    entity_range = rentity2.get_entity_range()
    expected_range = (3, 20)
    assert entity_range == expected_range, assert_msg

    # Test for ligand entity
    entity_range = lentity.get_entity_range()
    expected_range = (1, 1)
    assert entity_range == expected_range, assert_msg

    # Test for ion entity
    entity_range = ientity.get_entity_range()
    assert entity_range == expected_range, assert_msg

###############################################################################
# AFSequence class Fixtures and Tests
###############################################################################
@pytest.fixture
def paf_sequence1():
    return AFSequence(
        entity_info=entity_info_protein_min,
        protein_sequences=protein_sequences_by_name,
    )

@pytest.fixture
def paf_sequence2():
    return AFSequence(
        entity_info=entity_info_protein_max,
        protein_sequences=protein_sequences_by_name,
    )

@pytest.fixture
def paf_sequence3():
    return AFSequence(
        entity_info=entity_info_protein_max,
        protein_sequences=protein_sequences_by_id,
        entities_map=entities_map,
    )

@pytest.fixture
def daf_sequence1():
    return AFSequence(
        entity_info=entity_info_dna_min,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

@pytest.fixture
def daf_sequence2():
    return AFSequence(
        entity_info=entity_info_dna_max,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

@pytest.fixture
def raf_sequence1():
    return AFSequence(
        entity_info=entity_info_rna_min,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

@pytest.fixture
def raf_sequence2():
    return AFSequence(
        entity_info=entity_info_rna_max,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

@pytest.fixture
def laf_sequence():
    return AFSequence(
        entity_info=entity_info_ligand,
    )

@pytest.fixture
def iaf_sequence():
    return AFSequence(
        entity_info=entity_info_ion,
    )

###############################################################################
def test_create_af_sequence(
    paf_sequence1: AFSequence, paf_sequence2: AFSequence,
    paf_sequence3: AFSequence,
    daf_sequence1: AFSequence, daf_sequence2: AFSequence,
    raf_sequence1: AFSequence, raf_sequence2: AFSequence,
    laf_sequence: AFSequence, iaf_sequence: AFSequence,
):
    """Test AFSequence creation for different entity types."""

    assert_msg = "AFSequence does not match expected outcome."

    af_seq = paf_sequence1.create_af_sequence()
    expected_af_seq = {
        "proteinChain": {
            "sequence": "GAVLILLLVAVAVVAGVAA",
            "glycans": [],
            "modifications": [],
            "count": 1,
            "maxTemplateDate": MAX_TEMPLATE_DATE,
            "useStructureTemplate": True,
            "templates": [],
        }
    }
    assert af_seq == expected_af_seq, assert_msg

    af_seq = paf_sequence2.create_af_sequence()
    expected_af_seq = {
        "proteinChain": {
            "sequence": "KTAYIAKQRQISFVKSHFSRQDIL",
            "glycans": [{
                "residues": "BMA",
                "position": 5 - 2 + 1,
            }],
            "modifications": [{
                "ptmType": "CCD_HY3",
                "ptmPosition": 11 - 2 + 1,
            }],
            "count": 2,
            "maxTemplateDate": "2023-01-01",
            "useStructureTemplate": True,
            "templates": [],
        }
    }
    assert af_seq == expected_af_seq, assert_msg

    af_seq = paf_sequence3.create_af_sequence()
    assert af_seq == expected_af_seq, assert_msg

    # Test for dnaSequence entity
    af_seq = daf_sequence1.create_af_sequence()
    expected_af_seq = {
        "dnaSequence": {
            "sequence": "ATGCGTACGTAGCTAGCTAGCTAGCTA",
            "modifications": [],
            "count": 1,
        }
    }
    assert af_seq == expected_af_seq, assert_msg

    af_seq = daf_sequence2.create_af_sequence()
    expected_af_seq = {
        "dnaSequence": {
            "sequence": "GCGTACGTAGCTAGCTAG",
            "modifications": [
                {
                    "modificationType": "CCD_6OG",
                    "basePosition": 3 - 3 + 1,
                },
                {
                    "modificationType": "CCD_6MA",
                    "basePosition": 5 - 3 + 1,
                },
            ],
            "count": 1,
        }
    }
    assert af_seq == expected_af_seq, assert_msg

    # Test for rnaSequence entity
    af_seq = raf_sequence1.create_af_sequence()
    expected_af_seq = {
        "rnaSequence": {
            "sequence": "AUGCGUACGUAGCUAGCUAGCUAGCUA",
            "modifications": [],
            "count": 1,
        }
    }
    assert af_seq == expected_af_seq, assert_msg

    af_seq = raf_sequence2.create_af_sequence()
    expected_af_seq = {
        "rnaSequence": {
            "sequence": "GCGUACGUAGCUAGCUAG",
            "modifications": [
                {
                    "modificationType": "CCD_5MC",
                    "basePosition": 4 - 3 + 1,
                },
                {
                    "modificationType": "CCD_5MU",
                    "basePosition": 6 - 3 + 1,
                },
            ],
            "count": 1,
        }
    }
    assert af_seq == expected_af_seq, assert_msg

    # Test for ligand entity
    af_seq = laf_sequence.create_af_sequence()
    expected_af_seq = {
        "ligand": {
            "count": 2,
            "ligand": "CCD_ATP",
        }
    }
    assert af_seq == expected_af_seq, assert_msg

    # Test for ion entity
    af_seq = iaf_sequence.create_af_sequence()
    expected_af_seq = {
        "ion": {
            "count": 2,
            "ion": "MG",
        }
    }
    assert af_seq == expected_af_seq, assert_msg

###############################################################################
# AFJobSet class Fixtures and Tests
###############################################################################
@pytest.fixture
def af_jobset1():
    return AFJobSet(
        job_set_info=job_set_info_min,
        protein_sequences=protein_sequences_by_name,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

@pytest.fixture
def af_jobset2():
    return AFJobSet(
        job_set_info=job_set_info_max,
        protein_sequences=protein_sequences_by_name,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

@pytest.fixture
def af_jobset3():
    return AFJobSet(
        job_set_info=job_set_info_min_no_name,
        protein_sequences=protein_sequences_by_id,
        nucleic_acid_sequences=nucleic_acid_sequences,
        entities_map=entities_map,
    )

@pytest.fixture
def af_jobset4():
    return AFJobSet(
        job_set_info=job_set_info_max_no_name,
        protein_sequences=protein_sequences_by_id,
        nucleic_acid_sequences=nucleic_acid_sequences,
        entities_map=entities_map,
    )

###############################################################################
def test_create_job_set(
    af_jobset1: AFJobSet, af_jobset2: AFJobSet, af_jobset3: AFJobSet, af_jobset4: AFJobSet,
):

    assert_msg = "AFJobSet does not match expected outcome."
    assert_msg2 = "AFoffsets for chains do not match expected outcome."

    af_jobset = af_jobset1.create_job_set()
    expected_af_jobset = {
        "name": "jobset_min",
        "modelSeeds": [0, 1],
        "sequences": [
            {
                "proteinChain": {
                    "sequence": "GAVLILLLVAVAVVAGVAA",
                    "glycans": [],
                    "modifications": [],
                    "count": 1,
                    "maxTemplateDate": MAX_TEMPLATE_DATE,
                    "useStructureTemplate": True,
                    "templates": [],
                }
            },
            {
                "dnaSequence": {
                    "sequence": "ATGCGTACGTAGCTAGCTAGCTAGCTA",
                    "modifications": [],
                    "count": 1,
                }
            },
            {
                "rnaSequence": {
                    "sequence": "AUGCGUACGUAGCUAGCUAGCUAGCUA",
                    "modifications": [],
                    "count": 1,
                }
            },
            {
                "ligand": {
                    "count": 2,
                    "ligand": "CCD_ATP",
                }
            },
            {
                "ion": {
                    "count": 2,
                    "ion": "MG",
                }
            },
        ],
    }
    expected_af_offsets = {
        "A": [1, 19],
        "B": [1, 27],
        "C": [1, 27],
        "D": [1, 1],
        "E": [1, 1],
        "F": [1, 1],
        "G": [1, 1],
    }
    assert af_jobset == expected_af_jobset, assert_msg
    assert af_jobset1.job_set_af_offset == expected_af_offsets, assert_msg2

    af_jobset = af_jobset3.create_job_set()
    expected_af_jobset["name"] = f"protB_1_1{RES_RANGE_SEP}19_dnaA_1_1{RES_RANGE_SEP}27_rnaA_1_1{RES_RANGE_SEP}27_CCD_ATP_2_1{RES_RANGE_SEP}1_MG_2_1{RES_RANGE_SEP}1"
    assert af_jobset == expected_af_jobset, assert_msg

    af_jobset = af_jobset2.create_job_set()
    expected_af_jobset = {
        "name": "jobset_max",
        "modelSeeds": [0, 1],
        "sequences": [
            {
                "proteinChain": {
                    "sequence": "KTAYIAKQRQISFVKSHFSRQDIL",
                    "glycans": [{
                        "residues": "BMA",
                        "position": 5 - 2 + 1,
                    }],
                    "modifications": [{
                        "ptmType": "CCD_HY3",
                        "ptmPosition": 11 - 2 + 1,
                    }],
                    "count": 2,
                    "maxTemplateDate": "2023-01-01",
                    "useStructureTemplate": True,
                    "templates": [],
                }
            },
            {
                "dnaSequence": {
                    "sequence": "GCGTACGTAGCTAGCTAG",
                    "modifications": [
                        {
                            "modificationType": "CCD_6OG",
                            "basePosition": 3 - 3 + 1,
                        },
                        {
                            "modificationType": "CCD_6MA",
                            "basePosition": 5 - 3 + 1,
                        },
                    ],
                    "count": 1,
                }
            },
            {
                "rnaSequence": {
                    "sequence": "GCGUACGUAGCUAGCUAG",
                    "modifications": [
                        {
                            "modificationType": "CCD_5MC",
                            "basePosition": 4 - 3 + 1,
                        },
                        {
                            "modificationType": "CCD_5MU",
                            "basePosition": 6 - 3 + 1,
                        },
                    ],
                    "count": 1,
                }
            },
            {
                "ligand": {
                    "count": 2,
                    "ligand": "CCD_ATP",
                }
            },
            {
                "ion": {
                    "count": 2,
                    "ion": "MG",
                }
            },
        ],
    }
    expected_af_offsets = {
        "A": [2, 25],
        "B": [2, 25],
        "C": [3, 20],
        "D": [3, 20],
        "E": [1, 1],
        "F": [1, 1],
        "G": [1, 1],
        "H": [1, 1],
    }
    assert af_jobset == expected_af_jobset, assert_msg
    assert af_jobset2.job_set_af_offset == expected_af_offsets, assert_msg2

    af_jobset = af_jobset4.create_job_set()
    expected_af_jobset["name"] = f"protA_2_2{RES_RANGE_SEP}25_dnaA_1_3{RES_RANGE_SEP}20_rnaA_1_3{RES_RANGE_SEP}20_CCD_ATP_2_1{RES_RANGE_SEP}1_MG_2_1{RES_RANGE_SEP}1"
    assert af_jobset == expected_af_jobset, assert_msg

###############################################################################
# AlphaFoldServer class Fixtures and Tests
###############################################################################
@pytest.fixture
def alphafoldserver():
    return AFServerConfig(
        config_dict=config_dict,
        protein_sequences=protein_sequences_by_name,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

def test_write_job_files(alphafoldserver: AFServerConfig):

    alphafoldserver.write_job_files(
        output_dir=os.path.join(test_out_dir, "af3_input_jobs"),
        num_jobs_per_file=20,
    )

    assert os.path.isfile(os.path.join(
        test_out_dir, "af3_input_jobs", "jobset_min", "jobset_min_set_0.json"
    ))