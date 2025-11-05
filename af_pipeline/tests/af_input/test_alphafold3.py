import pytest
from af_pipeline.af_input.alphafold3 import (
    AlphaFoldServer, AFCycle, AFJobSet, AFSequence, Entity
)
from af_pipeline.constants.af_constants import MAX_TEMPLATE_DATE

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
    "type": "proteinChain",
}

entity_info_protein_max = {
    "name": "protA",
    "type": "proteinChain",
    "count": 2,
    "range": [2, 25],
    "useStructureTemplate": True,
    "maxTemplateDate": "2023-01-01",
    "glycans": [["BMA", 5]],
    "modifications": [["CCD_HY3", 11]],
}

entity_info_dna_min = {
    "name": "dnaA",
    "type": "dnaSequence",
}

entity_info_dna_max = {
    "name": "dnaA",
    "type": "dnaSequence",
    "count": 1,
    "range": [3, 20],
    "modifications": [["CCD_6OG", 3], ["CCD_6MA", 5]],
}

entity_info_rna_min = {
    "name": "rnaA",
    "type": "rnaSequence",
}

entity_info_rna_max = {
    "name": "rnaA",
    "type": "rnaSequence",
    "count": 1,
    "range": [3, 20],
    "modifications": [["CCD_5MC", 4], ["CCD_5MU", 6]],
}

entity_info_ligand = {
    "name": "CCD_ATP",
    "type": "ligand",
    "count": 2,
}

entity_info_ion = {
    "name": "MG",
    "type": "ion",
    "count": 2,
}
###############################################################################
job_set_info_min = {
    "name": "jobset_min",
    "modelSeeds": [0, 1],
    "entities": [
        entity_info_protein_min,
        entity_info_dna_min,
        entity_info_rna_min,
        entity_info_ligand,
        entity_info_ion,
    ]
}

job_set_info_min_no_name = {
    "modelSeeds": [0, 1],
    "entities": [
        entity_info_protein_min,
        entity_info_dna_min,
        entity_info_rna_min,
        entity_info_ligand,
        entity_info_ion,
    ]
}

job_set_info_max = {
    "name": "jobset_max",
    "modelSeeds": [0, 1],
    "entities": [
        entity_info_protein_max,
        entity_info_dna_max,
        entity_info_rna_max,
        entity_info_ligand,
        entity_info_ion,
    ]
}

job_set_info_max_no_name = {
    "modelSeeds": [0, 1],
    "entities": [
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
input_dict = {
    "cycle1": [job_sets_list[0]],
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
    }
    assert template_settings == expected_settings, assert_msg

    template_settings = pentity2.get_template_settings()
    expected_settings = {
        "useStructureTemplate": True,
        "maxTemplateDate": "2023-01-01",
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
    expected_af_jobset["name"] = "protB_1_1-19_dnaA_1_1-27_rnaA_1_1-27_CCD_ATP_2_1-1_MG_2_1-1"
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
    expected_af_jobset["name"] = "protA_2_2-25_dnaA_1_3-20_rnaA_1_3-20_CCD_ATP_2_1-1_MG_2_1-1"
    assert af_jobset == expected_af_jobset, assert_msg

###############################################################################
# AFCycle class Fixtures and Tests
###############################################################################
@pytest.fixture
def af_cycle():
    return AFCycle(
        job_sets_list=job_sets_list,
        protein_sequences=protein_sequences_by_name,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

###############################################################################
def test_update_cycle(af_cycle: AFCycle):
    """Test AFCycle.update_cycle method."""

    assert_msg = "Updated cycle does not match expected outcome."

    af_cycle.update_cycle()
    expected_job_list = [
        {
            "name": "jobset_min_0",
            "modelSeeds": [0],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "GAVLILLLVAVAVVAGVAA",
                        "glycans": [],
                        "modifications": [],
                        "count": 1,
                        "maxTemplateDate": MAX_TEMPLATE_DATE,
                        "useStructureTemplate": True,
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
        },
        {
            "name": "jobset_min_1",
            "modelSeeds": [1],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "GAVLILLLVAVAVVAGVAA",
                        "glycans": [],
                        "modifications": [],
                        "count": 1,
                        "maxTemplateDate": MAX_TEMPLATE_DATE,
                        "useStructureTemplate": True,
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
        },
        {
            "name": "jobset_max_0",
            "modelSeeds": [0],
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
        },
        {
            "name": "jobset_max_1",
            "modelSeeds": [1],
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
        },
    ]
    expected_job_set_names = [
        "jobset_min",
        "jobset_max",
    ]
    expected_job_set_af_offsets = [
        {"A": [1, 19],
         "B": [1, 27],
         "C": [1, 27],
         "D": [1, 1],
         "E": [1, 1],
         "F": [1, 1],
         "G": [1, 1]},
        {'A': [2, 25],
         'B': [2, 25],
         'C': [3, 20],
         'D': [3, 20],
         'E': [1, 1],
         'F': [1, 1],
         'G': [1, 1],
         'H': [1, 1]},
    ]

    assert af_cycle.job_list == expected_job_list, assert_msg
    assert af_cycle.job_set_names == expected_job_set_names, assert_msg
    assert af_cycle.job_set_af_offsets == expected_job_set_af_offsets, assert_msg

###############################################################################
# AlphaFold3 class Fixtures and Tests
###############################################################################
@pytest.fixture
def alphafoldserver():
    return AlphaFoldServer(
        input_dict=input_dict,
        protein_sequences=protein_sequences_by_name,
        nucleic_acid_sequences=nucleic_acid_sequences,
    )

def test_create_af3_job_cycles(alphafoldserver: AlphaFoldServer):
    """Test AlphaFoldServer.create_af3_job_cycles method."""

    assert_msg = "AF3 job cycles do not match expected outcome."

    job_cycles, job_set_names, af_offsets = alphafoldserver.create_af3_job_cycles()
    expected_af3_job_cycles = {
        "cycle1": [
            {
                "name": "jobset_min_0",
                "modelSeeds": [0],
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": "GAVLILLLVAVAVVAGVAA",
                            "glycans": [],
                            "modifications": [],
                            "count": 1,
                            "maxTemplateDate": MAX_TEMPLATE_DATE,
                            "useStructureTemplate": True,
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
            },
            {
                "name": "jobset_min_1",
                "modelSeeds": [1],
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": "GAVLILLLVAVAVVAGVAA",
                            "glycans": [],
                            "modifications": [],
                            "count": 1,
                            "maxTemplateDate": MAX_TEMPLATE_DATE,
                            "useStructureTemplate": True,
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
            },
        ],
    }
    expected_job_set_names = {"cycle1": ["jobset_min"]}
    expected_af_offsets = {
        "cycle1": [
            {
                "A": [1, 19],
                "B": [1, 27],
                "C": [1, 27],
                "D": [1, 1],
                "E": [1, 1],
                "F": [1, 1],
                "G": [1, 1],
            }
        ]
    }

    assert job_cycles == expected_af3_job_cycles, assert_msg
    assert job_set_names == expected_job_set_names, assert_msg
    assert af_offsets == expected_af_offsets, assert_msg