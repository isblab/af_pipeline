# some constants

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from af_pipeline.utils.file_utils import read_json, read_pkl
from typing import Final
from dataclasses import dataclass

RES_RANGE_SEP: Final[str] = "t"

SEED_MULTIPLIER: Final[int] = 100
""" Multiplier to generate model seeds."""

MAX_TEMPLATE_DATE = "2021-09-30"

JOB_LIMIT_PER_JSON = 100

AVAILABLE_PARSERS = {
    "pdb": PDBParser,
    "cif": MMCIFParser,
}

@dataclass(frozen=True)
class ConfigYaml:
    input: str = "af_input_jobs"
    master: str = "af_master_dirs"
    cycle: str = "af_cycle_dirs"
    job_set: str = "af_job_set_dirs"
    best_pred: str = "best_af3_predictions"

ALLOWED_STRUCTURE_FORMATS = list(AVAILABLE_PARSERS.keys())

AVAILABLE_DATA_READERS = {
    "pkl": read_pkl,
    "json": read_json,
}

ALLOWED_DATA_FORMATS = list(AVAILABLE_DATA_READERS.keys())

AVAILABLE_ATOM_QUANTITIES = [
    "coord",
    "plddt",
    "atom_name",
    "res_pos",
    "res_name",
    "chain_id",
    "entity_type",
    "atom_local_idx",
]

AVAILABLE_RESIDUE_QUANTITIES =[
    "res_pos",
    "res_name",
    "coord",
    "plddt",
    "chain_id",
    "entity_type",
    "atoms",
    "atom_local_idxs",
    "rep_atom",
    "rep_atom_local_idx",
]

PTM = [
    "CCD_SEP",
    "CCD_TPO",
    "CCD_PTR",
    "CCD_NEP",
    "CCD_HIP",
    "CCD_ALY",
    "CCD_MLY",
    "CCD_M3L",
    "CCD_MLZ",
    "CCD_2MR",
    "CCD_AGM",
    "CCD_MCS",
    "CCD_HYP",
    "CCD_HY3", #
    "CCD_LYZ",
    "CCD_AHB",
    "CCD_P1L",
    "CCD_SNN", #
    "CCD_SNC",
    "CCD_TRF",
    "CCD_KCR",
    "CCD_CIR", #
    "CCD_YHA",
]

DNA_MOD_PURINES = [
    "CCD_6OG",
    "CCD_6MA",
    "CCD_8OG",
]
DNA_MOD_PYRIMIDINES = [
    "CCD_5CM",
    "CCD_C34",
    "CCD_5HC",
    "CCD_1CC",
    "CCD_5FC",
]
DNA_MOD_OTHER = [
    "CCD_3DR",
]

DNA_MOD = DNA_MOD_PURINES + DNA_MOD_PYRIMIDINES + DNA_MOD_OTHER

RNA_MOD_PURINES = [
    "CCD_A2M",
    "CCD_MA6",
    "CCD_6MZ",
    "CCD_2MG",
    "CCD_OMG",
    "CCD_7MG",
]
RNA_MOD_PYRIMIDINES = [
    "CCD_PSU",
    "CCD_5MC",
    "CCD_OMC",
    "CCD_4OC",
    "CCD_5MU",
    "CCD_OMU",
    "CCD_UR3",
    "CCD_RSQ",
]

RNA_MOD = RNA_MOD_PURINES + RNA_MOD_PYRIMIDINES

LIGAND = [
    "CCD_ADP",
    "CCD_ATP",
    "CCD_AMP",
    "CCD_GTP",
    "CCD_GDP",
    "CCD_FAD",
    "CCD_NAD",
    "CCD_NAP",
    "CCD_NDP",
    "CCD_HEM",
    "CCD_HEC",
    "CCD_PLM",
    "CCD_OLA",
    "CCD_MYR",
    "CCD_CIT",
    "CCD_CLA",
    "CCD_CHL",
    "CCD_BCL",
    "CCD_BCB",
]

ION = [
    "MG",
    "ZN",
    "CL",
    "CA",
    "NA",
    "MN",
    "K",
    "FE",
    "CU",
    "CO",
]

ENTITY_TYPES = [
    "proteinChain",
    "dnaSequence",
    "rnaSequence",
    "ligand",
    "ion",
]

STD_RESIDUES = [
    "GLY",
    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "THR",
    "SER",
    "MET",
    "CYS",
    "PRO",
    "PHE",
    "TYR",
    "TRP",
    "HIS",
    "LYS",
    "ARG",
    "ASP",
    "GLU",
    "ASN",
    "GLN",
]

DNA_PURINES = [
    "DA",
    "DG"
]
RNA_PURINES = [
    "A",
    "G"
]
RNA_PYRIMIDINES = [
    "U",
    "C"
]
DNA_PYRIMIDINES = [
    "DC",
    "DT"
]

REP_ATOMS = {
    "proteinChain": "CB",
    "is_ca_only": "CA",
    "is_purine": "C4",
    "is_pyrimidine": "C2",
}

ALLOWED_LIGANDS = [ligand.split("_")[1] for ligand in LIGAND]
ALLOWED_PTMS = [ptm_ccd.split("_")[1] for ptm_ccd in PTM]
ALLOWED_DNA_MODS = [mod_dna.split("_")[1] for mod_dna in DNA_MOD]
ALLOWED_RNA_MODS = [mod_rna.split("_")[1] for mod_rna in RNA_MOD]

ALLOWED_PURINE_MODS = [mod_purine.split("_") for mod_purine in DNA_MOD_PURINES + RNA_MOD_PURINES]
ALLOWED_PYRIMIDINE_MODS = [mod_pyrimidine.split("_") for mod_pyrimidine in DNA_MOD_PYRIMIDINES + RNA_MOD_PYRIMIDINES]

PROTEIN_ENTITIES = ALLOWED_PTMS + STD_RESIDUES
DNA_ENTITIES = ALLOWED_DNA_MODS + DNA_PURINES + DNA_PYRIMIDINES
RNA_ENTITIES = ALLOWED_RNA_MODS + RNA_PURINES + RNA_PYRIMIDINES

PURINES_STD = DNA_PURINES + RNA_PURINES
PYRIMIDINES_STD = DNA_PYRIMIDINES + RNA_PYRIMIDINES


PURINES = PURINES_STD + ALLOWED_PURINE_MODS
PYRIMIDINES = PYRIMIDINES_STD + ALLOWED_PYRIMIDINE_MODS

ONLY_CA_RESIDUES = [
    "GLY",
    "HY3",
    "SNN",
    "CIR",
]

CHAINWISE_ASSESSMENT_COLUMNS = {
    True: [
        "Chain ID",
        "Protein Name",
        "Average pLDDT",
        "Average ipLDDT",
        "Interface Residues",
        "Chain Type",
    ],
    False: [
        "Chain ID",
        "Protein Name",
        "Residue Number",
        "Average pLDDT",
        "Average ipLDDT",
        "Chain Type",
    ],
}

CHAIN_PAIRWISE_ASSESSMENT_COLUMNS = {
    (True, True): [
        "Protein Name",
        "Chain Type",
        "Chain ID",
        "Interface Residues",
        "Number of contacts",
        "Average ipLDDT chain1",
        "Average ipLDDT chain2",
        "Average PAE",
        "Average iPAE",
        # "Minimum PAE",
    ],
    (True, False): [
        "Protein Name",
        "Chain Type",
        "Chain ID",
        "Interface Residues",
        "Number of contacts",
        "Average ipLDDT chain1",
        "Average ipLDDT chain2",
        "Average PAE ij",
        "Average PAE ji",
        # "Minimum PAE ij",
        # "Minimum PAE ji",
        "Average iPAE ij",
        "Average iPAE ji",
    ],
    (False, True): [
        "Protein Name 1",
        "Chain Type 1",
        "Chain ID 1",
        "Protein Name 2",
        "Chain Type 2",
        "Chain ID 2",
        "Residue 1",
        "Residue 2",
        "pLDDT 1",
        "pLDDT 2",
        "PAE",
    ],
    (False, False): [
        "Protein Name 1",
        "Chain Type 1",
        "Chain ID 1",
        "Protein Name 2",
        "Chain Type 2",
        "Chain ID 2",
        "pLDDT 1",
        "pLDDT 2",
        "Residue 1",
        "Residue 2",
        "PAE ij",
        "PAE ji",
    ]
}


OVERALL_ASSESSMENT_COLUMNS = {
    True: {
        "Number of Chains": "num_chains",
        "Number of Interacting Chain Pairs": "num_interacting_chain_pairs",
        "Interface Residues": "num_interface_residues",
        "Number of Contacts": "num_contacts",
        "Average ipLDDT": "avg_iplddt",
        "Average IDR ipLDDT": "avg_idr_iplddt",
        "Average iPAE": "avg_ipae",
    },
    False: {
        "Number of Chains": "num_chains",
        "Number of Interacting Chain Pairs": "num_interacting_chain_pairs",
        "Interface Residues": "num_interface_residues",
        "Number of Contacts": "num_contacts",
        "Average ipLDDT": "avg_iplddt",
        "Average IDR ipLDDT": "avg_idr_iplddt",
        "Average iPAE ij": "avg_ipae_ij",
        "Average iPAE ji": "avg_ipae_ji",
    }
}
