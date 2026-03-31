# some constants

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.PDBParser import PDBParser
from af_pipeline.utils.file_utils import read_json, read_pkl
from typing import Final
from dataclasses import dataclass
from enum import StrEnum, IntEnum, auto
from string import Template

RANDOM_SEED = 47

class ColorMapScheme(StrEnum):
    SOFT_WARM = auto()
    STANDARD = auto()
    NON_BRIGHT = auto()
    EARTH_TONE = auto()
    COOL_TONE = auto()
    CONTRASTING_NON_BRIGHT = auto()
    BINARY = auto()

class BinaryColorMap(StrEnum):
    ZERO_C = "black"
    ONE_C = "green"

class MiscStrEnum(StrEnum):
    TOTAL = auto()
    FIRST = auto()
    LAST = auto()
    UNKNOWN = auto()
    REVERSE_COMPLEMENT = "RevComp"

RES_SEPARATOR = "-"
RES_SPLITTER = ","

class UpdateConfigMode(StrEnum):
    REPLACE = auto()
    SOFT_REPLACE = auto()

class MaskedInteractionValue(IntEnum):
    MASKED_V = 1
    UNMASKED_V = 0

class MaskedInteractionType(StrEnum):
    INTRA_PART = auto()
    INTER_PART = auto()

class ReturnType(StrEnum):
    DICT = auto()
    DATAFRAME = auto()
    LIST = auto()
    ARRAY = auto()
    SET = auto()

class MetricLevel(StrEnum):
    PER_TOKEN = auto()
    REPRESENTATIVE_TOKEN = auto()
    PER_ATOM = auto()

VALID_AF3_METRIC_LEVELS = [
    MetricLevel.PER_TOKEN,
    MetricLevel.REPRESENTATIVE_TOKEN,
]

VALID_METRIC_LEVELS = list(MetricLevel)

class DataFileKeys(StrEnum):
    TOKEN_CHAIN_IDS = auto()
    TOKEN_RES_IDS = auto()
    ATOM_CHAIN_IDS = auto()
    ATOM_PLDDTS = auto()
    PAE = auto()
    PREDICTED_ALIGNED_ERROR = auto()
    CONTACT_PROBS = auto()

@dataclass
class InitializeConstants:
    average_token_pae = True
    average_token_plddt = True
    metric_level = MetricLevel.PER_TOKEN
    use_fast_cif_parser = False

class CommunityDetectionLibrary(StrEnum):
    IGRAPH = auto()
    NETWORKX = auto()
    LABEL_PROPAGATION = auto()

class FileFormat(StrEnum):
    PDB = auto()
    CIF = auto()
    TXT = auto()
    JSON = auto()
    PKL = auto()
    XLSX = auto()
    PNG = auto()
    HTML = auto()
    FASTA = auto()

@dataclass
class RigidBodiesConstants:
    valid_libraries = list(CommunityDetectionLibrary)
    library = CommunityDetectionLibrary.NETWORKX
    pae_cutoff = 12.0
    pae_power = 1
    plddt_cutoff = 70.0
    plddt_cutoff_idr = 50.0
    resolution = 0.5
    min_res = 1
    min_proteins = 1
    plddt_filter = True
    random_seed = RANDOM_SEED
    valid_rb_out_fmts = [
        FileFormat.TXT,
        FileFormat.JSON,
    ]
    rb_out_fmt = FileFormat.TXT
    valid_rb_struct_fmts = [
        FileFormat.PDB,
        FileFormat.CIF,
    ]
    rb_struct_fmt = FileFormat.CIF
    save_structure = True
    filter_struct_by_plddt = True
    rb_name_template = Template("rigid_body_${rb_idx}")
    rb_assessment_name_template = Template("rigid_body_${rb_idx}_assessment")

class KeywordArg(StrEnum):
    PLDDT_CUTOFF_IDR = auto()
    IDR_CHAINS = auto()
    SAVE_PLOT = auto()
    SAVE_TABLE = auto()
    RANDOM_SEED = auto()
    SETUP_INSTANCE = auto()
    PROTEIN_CHAIN_MAP = auto()

class PlotType(StrEnum):
    STATIC = auto()
    INTERACTIVE = auto()
    BOTH = auto()

@dataclass
class InteractionConstants:
    contact_threshold = 8.0 # Distance threshold in (Angstorm) to define a contact between residue pairs.
    plddt_cutoff = 70.0 # pLDDT cutoff to consider a confident prediction.
    pae_cutoff = 5.0 # PAE cutoff to consider a confident prediction.
    plddt_cutoff_idr = 50.0 # pLDDT cutoff for IDR chains.
    save_plot = False
    save_table = False
    plot_type = PlotType.STATIC
    valid_plot_types = [
        PlotType.STATIC,
        PlotType.INTERACTIVE,
        PlotType.BOTH
    ]

class ResidueMapKeys(StrEnum):
    CHAIN_ID = auto()
    TOKEN_NUM = auto()
    ATOM_NAME = auto()

class ResidueMapDepth(StrEnum):
    ATOM = auto()
    RESIDUE = auto()

class TokenLevel(StrEnum):
    ATOM = auto()
    RESIDUE = auto()

class ResidueDecoration(StrEnum):
    ENTITY_TYPE = "entityType"
    TOKEN_LEVEL = auto()
    IS_MODIFIED = auto()
    IS_CA_ONLY = auto()
    IS_PURINE = auto()
    IS_PYRIMIDINE = auto()

class AtomDecoration(StrEnum):
    IS_REPRESENTATIVE = auto()

class EntityType(StrEnum):
    PROTEIN_CHAIN = "proteinChain"
    DNA_SEQUENCE = "dnaSequence"
    RNA_SEQUENCE = "rnaSequence"
    LIGAND = auto()
    ION = auto()

class InteractionMapType(StrEnum):
    DISTANCE = auto()
    CONTACT = auto()

VALID_INTERACTION_MAP_TYPES = [
    InteractionMapType.DISTANCE,
    InteractionMapType.CONTACT,
]

# @dataclass
# class RenumberResiduesConstants:
#     valid_depths = [
#         ResidueMapDepth.ATOM,
#         ResidueMapDepth.RESIDUE,
#     ]

@dataclass
class StructureParserConstants:
    preserve_header_footer = False
    use_fast_cif_parser = False

RES_RANGE_SEP: Final[str] = "t"

SEED_MULTIPLIER: Final[int] = 100
""" Multiplier to generate model seeds."""

MAX_TEMPLATE_DATE = "2021-09-30"

JOB_LIMIT_PER_JSON = 100

AVAILABLE_PARSERS = {
    FileFormat.PDB: PDBParser,
    FileFormat.CIF: MMCIFParser,
    "fast_cif": FastMMCIFParser,
}

AF_JOB_FILE = Template("${fname}_set_${set_idx}")

class AFServerSequenceFields(StrEnum):
    SEQUENCE = "sequence"
    GLYCANS = auto()
    MODIFICATIONS = auto()
    COUNT = auto()

class AFInputJobFields(StrEnum):
    NAME = auto()
    JOB_NAME = auto()
    JOB_SET_NAME = auto()
    MODEL_SEEDS = "modelSeeds"
    ENTITIES = auto()
    AF_OFFSET = auto()
    SEQUENCES = auto()

class AFInputEntityFields(StrEnum):
    NAME = auto()
    IDENTIFIER = auto()
    SEQUENCE = auto()
    TYPE = auto()
    COUNT = auto()
    RANGE = auto()
    STRAND = auto()
    USE_STRUCTURE_TEMPLATE = "useStructureTemplate"
    MAX_TEMPLATE_DATE = "maxTemplateDate"
    GLYCANS = auto()
    MODIFICATIONS = auto()

class NucleicAcidStrand(StrEnum):
    DOUBLE = auto()
    SINGLE = auto()

class GlycanModificationFields(StrEnum):
    RESIDUES = auto()
    POSITION = auto()

class ProteinModificationFields(StrEnum):
    PTM_TYPE = "ptmType"
    PTM_POSITION = "ptmPosition"

class NucleicAcidModificationFields(StrEnum):
    MODIFICATION_TYPE = "modificationType"
    BASE_POSITION = "basePosition"

# @dataclass(frozen=True)
# class ConfigYaml:
#     input: str = "af_input_jobs"
#     master: str = "af_master_dirs"
#     cycle: str = "af_cycle_dirs"
#     job_set: str = "af_job_set_dirs"
#     best_pred: str = "best_af3_predictions"

class ConfigYaml(StrEnum):
    AF_INPUT_JOBS = auto()
    AF_MASTER_DIRS = auto()
    AF_CYCLE_DIRS = auto()
    AF_JOB_SET_DIRS = auto()
    BEST_AF3_PREDICTIONS = auto()
    PROTEIN_UNIPROT_MAP = auto()

class BestPredictionFields(StrEnum):
    STRUCTURE_PATH = auto()
    DATA_PATH = auto()
    AF_OFFSET = auto()

class AF3Metrics(StrEnum):
    RANKING_SCORE = auto()
    IPTM = auto()
    PTM = auto()
    FRACTION_DISORDERED = auto()
    MODEL_PATH = auto()
    MODEL_IDX = auto()

class AF3SummaryConfidenceFields(StrEnum):
    FRACTION_DISORDERED = auto()
    HAS_CLASH = auto()
    IPTM = auto()
    NUM_RECYCLES = auto()
    PTM = auto()
    RANKING_SCORE = auto()

ALLOWED_STRUCTURE_FORMATS = list(AVAILABLE_PARSERS.keys())

AVAILABLE_DATA_READERS = {
    FileFormat.PKL: read_pkl,
    FileFormat.JSON: read_json,
}

ALLOWED_DATA_FORMATS = list(AVAILABLE_DATA_READERS.keys())

class AtomQuantity(StrEnum):
    COORD = auto()
    PLDDT = auto()
    ATOM_NAME = auto()
    RES_POS = auto()
    RES_NAME = auto()
    CHAIN_ID = auto()
    ENTITY_TYPE = auto()
    ATOM_LOCAL_IDX = auto()

AVAILABLE_ATOM_QUANTITIES = [
    AtomQuantity.COORD,
    AtomQuantity.PLDDT,
    AtomQuantity.ATOM_NAME,
    AtomQuantity.RES_POS,
    AtomQuantity.RES_NAME,
    AtomQuantity.CHAIN_ID,
    AtomQuantity.ENTITY_TYPE,
    AtomQuantity.ATOM_LOCAL_IDX,
]

class ResidueQuantity(StrEnum):
    RES_POS = auto()
    RES_NAME = auto()
    COORD = auto()
    PLDDT = auto()
    CHAIN_ID = auto()
    ENTITY_TYPE = auto()
    ATOMS = auto()
    ATOM_LOCAL_IDXS = auto()
    REP_ATOM = auto()
    REP_ATOM_LOCAL_IDX = auto()

class QuantityLevel(StrEnum):
    REPRESENTATIVE_ATOM = auto()
    PER_ATOM = auto()
    AVERAGE_ATOM = auto()

AVAILABLE_RESIDUE_QUANTITIES =[
    ResidueQuantity.RES_POS,
    ResidueQuantity.RES_NAME,
    ResidueQuantity.COORD,
    ResidueQuantity.PLDDT,
    ResidueQuantity.CHAIN_ID,
    ResidueQuantity.ENTITY_TYPE,
    ResidueQuantity.ATOMS,
    ResidueQuantity.ATOM_LOCAL_IDXS,
    ResidueQuantity.REP_ATOM,
    ResidueQuantity.REP_ATOM_LOCAL_IDX,
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
        "Total Residues": "num_total_residues",
        "Sequence Coverage": "rb_coverage",
        "Number of Contacts": "num_contacts",
        "Average pLDDT": "avg_plddt",
        "Average ipLDDT": "avg_iplddt",
        "Average IDR ipLDDT": "avg_idr_iplddt",
        "Average iPAE": "avg_ipae",
    },
    False: {
        "Number of Chains": "num_chains",
        "Number of Interacting Chain Pairs": "num_interacting_chain_pairs",
        "Interface Residues": "num_interface_residues",
        "Total Residues": "num_total_residues",
        "Sequence Coverage": "rb_coverage",
        "Number of Contacts": "num_contacts",
        "Average pLDDT": "avg_plddt",
        "Average ipLDDT": "avg_iplddt",
        "Average IDR ipLDDT": "avg_idr_iplddt",
        "Average iPAE ij": "avg_ipae_ij",
        "Average iPAE ji": "avg_ipae_ji",
    }
}
