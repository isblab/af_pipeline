"""
[af_constants](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/constants/af_constants.py)
===================================

"""

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.PDBParser import PDBParser
from af_pipeline.utils.file_utils import read_json, read_pkl
from typing import Final
from dataclasses import dataclass
from enum import StrEnum, auto
from string import Template

RANDOM_SEED = 47

class ColorMapScheme(StrEnum):
    """ Color map schemes for plotting."""
    SOFT_WARM = auto()
    STANDARD = auto()
    NON_BRIGHT = auto()
    EARTH_TONE = auto()
    COOL_TONE = auto()
    CONTRASTING_NON_BRIGHT = auto()
    BINARY = auto()

class BinaryColorMap(StrEnum):
    """ Color map for binary contact maps."""
    ZERO_C = "black"
    ONE_C = "green"

class MiscStrEnum(StrEnum):
    """ Miscellaneous string enums."""
    TOTAL = auto()
    FIRST = auto()
    LAST = auto()
    UNKNOWN = auto()
    REVERSE_COMPLEMENT = "RevComp"

RES_SEPARATOR = "-"
""" Separator for residue ranges. For example, a residue range can be represented as `1-100`."""

RES_SPLITTER = ","
""" Separator for multiple residue ranges. For example, `1-100,150-200`."""

class MaskedInteractionType(StrEnum):
    """ Which interactions to mask in the 2D interaction map."""
    INTRA_PART = auto()
    INTER_PART = auto()

class ReturnType(StrEnum):
    """ Return types for various functions."""
    DICT = auto()
    DATAFRAME = auto()
    LIST = auto()
    ARRAY = auto()
    SET = auto()

class MetricLevel(StrEnum):
    """ Levels at which the confidence metrics are calculated."""
    PER_TOKEN = auto()
    REPRESENTATIVE_TOKEN = auto()
    PER_ATOM = auto()

VALID_AF3_METRIC_LEVELS = [
    MetricLevel.PER_TOKEN,
    MetricLevel.REPRESENTATIVE_TOKEN,
]
""" Valid metric levels for AlphaFold3 predictions."""

VALID_METRIC_LEVELS = list(MetricLevel)
""" Valid metric levels for confidence metrics."""

class DataFileKeys(StrEnum):
    """ Keys for the data stored in JSON or PKL files accompanying the predicted structures."""
    TOKEN_CHAIN_IDS = auto()
    TOKEN_RES_IDS = auto()
    ATOM_CHAIN_IDS = auto()
    ATOM_PLDDTS = auto()
    PAE = auto()
    PREDICTED_ALIGNED_ERROR = auto()
    CONTACT_PROBS = auto()

@dataclass
class InitializeConstants:
    """ Default values for initializing the Initialize instance."""
    average_token_pae = False
    average_token_plddt = False
    metric_level = MetricLevel.PER_TOKEN
    use_fast_cif_parser = False

class CommunityDetectionLibrary(StrEnum):
    """ Libraries for community detection used in af_pipeline.rigid_bodies.rigid_bodies.RigidBodies class."""
    IGRAPH = auto()
    NETWORKX = auto()
    LABEL_PROPAGATION = auto()

class FileFormat(StrEnum):
    """ File formats for input and output files."""
    PDB = auto()
    CIF = auto()
    TXT = auto()
    JSON = auto()
    PKL = auto()
    XLSX = auto()
    PNG = auto()
    HTML = auto()
    FASTA = auto()
    CSV = auto()

@dataclass
class RigidBodiesConstants:
    """ Default values for af_pipeline.rigid_bodies.rigid_bodies.RigidBodies class."""
    valid_libraries = list(CommunityDetectionLibrary)
    """ Libraries for community detection used in af_pipeline.rigid_bodies.rigid_bodies.RigidBodies class."""
    library = CommunityDetectionLibrary.NETWORKX
    """ Library for community detection. Default is NETWORKX. Options are:
    - IGRAPH: Use the igraph library for community detection.
    - NETWORKX: Use the networkx library for community detection.
    - LABEL_PROPAGATION: Use the label propagation algorithm for community detection."""
    pae_cutoff = 12.0
    """ PAE cutoff in Angstroms to define edges in the graph for community detection."""
    pae_power = 1
    """ Power to which the PAE values are raised when defining edge weights for community detection."""
    plddt_cutoff = 70.0
    """ pLDDT cutoff to consider a confident prediction when extracting rigid bodies."""
    plddt_cutoff_idr = 50.0
    """ pLDDT cutoff for IDR chains."""
    resolution = 0.5
    """ Resolution parameter for community detection algorithms that require it."""
    min_res = 1
    """ Minimum number of residues in a rigid body."""
    min_proteins = 1
    """ Minimum number of unique protein chains in a rigid body."""
    plddt_filter = True
    """ Whether to filter out residues with low pLDDT scores when extracting rigid bodies."""
    random_seed = RANDOM_SEED
    """ Random seed for reproducibility in community detection using fast label propagation."""
    valid_rb_out_fmts = [
        FileFormat.TXT,
        FileFormat.JSON,
    ]
    """ File formats for saving the extracted rigid bodies. Options are:
    - TXT: Save rigid bodies in a text file format.
    - JSON: Save rigid bodies in a JSON file format."""
    rb_out_fmt = FileFormat.TXT
    """ Default file format for saving the extracted rigid bodies. (Default is TXT)"""
    valid_rb_struct_fmts = [
        FileFormat.PDB,
        FileFormat.CIF,
    ]
    """ File formats for saving the extracted rigid body structures. Options are:
    - PDB: Save rigid body structures in PDB format.
    - CIF: Save rigid body structures in CIF format."""
    rb_struct_fmt = FileFormat.CIF
    """ Default file format for saving the extracted rigid body structures. (Default is CIF)"""
    save_structure = True
    """ Whether to save the extracted rigid body structures."""
    filter_struct_by_plddt = True
    """ Whether to filter the extracted rigid body structures by pLDDT scores."""
    rb_name_template = Template("rigid_body_${rb_idx}")
    """ Template for naming the extracted rigid bodies. The template should contain
    `${rb_idx}` which will be replaced by the rigid body index."""
    rb_assessment_name_template = Template("rigid_body_${rb_idx}_assessment")
    """ Template for naming the rigid body assessment files. The template should contain
    `${rb_idx}` which will be replaced by the rigid body index."""

class KeywordArg(StrEnum):
    """ Keyword arguments for various functions."""
    PLDDT_CUTOFF_IDR = auto()
    IDR_CHAINS = auto()
    SAVE_PLOT = auto()
    SAVE_TABLE = auto()
    RANDOM_SEED = auto()
    SETUP_INSTANCE = auto()
    PROTEIN_CHAIN_MAP = auto()
    PAE_POWER = auto()
    RESOLUTION = auto()

class PlotType(StrEnum):
    """ Types of plots. \n
    - STATIC: Static plots using matplotlib.
    - INTERACTIVE: Interactive plots using plotly.
    - BOTH: Both static and interactive plots."""
    STATIC = auto()
    INTERACTIVE = auto()
    BOTH = auto()

@dataclass
class InteractionConstants:
    """ Default values for af_pipeline.interaction.interaction.Interaction class."""
    contact_threshold = 8.0
    """ Distance threshold in (Angstorm) to define a contact between residue pairs."""
    plddt_cutoff = 70.0
    """ pLDDT cutoff to consider a confident prediction."""
    pae_cutoff = 5.0
    """ PAE cutoff to consider a confident prediction."""
    plddt_cutoff_idr = 50.0
    """ pLDDT cutoff for IDR chains."""
    save_plot = False
    """ Whether to save the interaction map plot."""
    save_table = False
    """ Whether to save the interaction map table."""
    valid_plot_types = [
        PlotType.STATIC,
        PlotType.INTERACTIVE,
        PlotType.BOTH
    ]
    """ Type of plot to generate.
    - STATIC: Static plots using matplotlib.
    - INTERACTIVE: Interactive plots using plotly.
    - BOTH: Both static and interactive plots."""
    plot_type = PlotType.STATIC
    """ Default plot type. (Default is STATIC)"""

class ResidueMapKeys(StrEnum):
    """ Keys for residue to index and index to residue mapping."""
    CHAIN_ID = auto()
    TOKEN_NUM = auto()
    ATOM_NAME = auto()

class ResidueMapDepth(StrEnum):
    """ Depth of the residue mapping."""
    ATOM = auto()
    RESIDUE = auto()

class TokenLevel(StrEnum):
    """ Token level."""
    ATOM = auto()
    RESIDUE = auto()

class ResidueDecoration(StrEnum):
    """ Decoration for residues in the Biopython structure object."""
    ENTITY_TYPE = "entityType"
    TOKEN_LEVEL = auto()
    IS_MODIFIED = auto()
    IS_CA_ONLY = auto()
    IS_PURINE = auto()
    IS_PYRIMIDINE = auto()

class AtomDecoration(StrEnum):
    """ Decoration for atoms in the Biopython structure object."""
    IS_REPRESENTATIVE = auto()

class EntityType(StrEnum):
    """ Type of entity (AlphaFold3 definition)."""
    PROTEIN_CHAIN = "proteinChain"
    DNA_SEQUENCE = "dnaSequence"
    RNA_SEQUENCE = "rnaSequence"
    LIGAND = auto()
    ION = auto()

ENTITY_TYPES = list(EntityType)

class InteractionMapType(StrEnum):
    """ Type of interaction map to generate."""
    DISTANCE = auto()
    CONTACT = auto()

VALID_INTERACTION_MAP_TYPES = [
    InteractionMapType.DISTANCE,
    InteractionMapType.CONTACT,
]
""" Valid interaction map types for generating interaction maps. Options are:
- DISTANCE: Generate distance maps.
- CONTACT: Generate contact maps."""

@dataclass
class StructureParserConstants:
    """ Default parameter values for af_pipeline.parser.structure_parser.StructureParser."""
    preserve_header_footer = False
    use_fast_cif_parser = False

RES_RANGE_SEP: Final[str] = "t"
""" Separator for residue range in AlphaFold job name."""

SEED_MULTIPLIER: Final[int] = 100
""" Multiplier to generate model seeds."""

MAX_TEMPLATE_DATE = "2021-09-30"
""" Maximum template date for AlphaFold 3 predictions."""

JOB_LIMIT_PER_JSON = 100
""" Maximum number of jobs to store in a single JSON file for AlphaFold Server."""

AVAILABLE_PARSERS = {
    FileFormat.PDB: PDBParser,
    FileFormat.CIF: MMCIFParser,
    "fast_cif": FastMMCIFParser,
}
""" Available parsers for parsing structure files."""

AF_JOB_FILE = Template("${fname}_set_${set_idx}")
""" Template for naming the AlphaFold job files. The template should contain
`${fname}` which will be replaced by the base name of the input file and `${set_idx}`
which will be replaced by the job set index."""

class AFServerSequenceFields(StrEnum):
    """ Fields related to 'sequences' in the AlphaFold Server input."""
    SEQUENCE = "sequence"
    GLYCANS = auto()
    MODIFICATIONS = auto()
    COUNT = auto()
    UNPAIREDMSA = "unpairedMsa"
    PAIREDMSA = "pairedMsa"
    TEMPLATES = "templates"

class AFInputJobFields(StrEnum):
    NAME = auto()
    JOB_NAME = auto()
    JOB_SET_NAME = auto()
    MODEL_SEEDS = "modelSeeds"
    ENTITIES = auto()
    AF_OFFSET = auto()
    SEQUENCES = auto()
    MSA = auto()

class AFInputEntityFields(StrEnum):
    """ Fields for entities in the config dictionary used by modules in af_pipeline.af_input"""
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
    PAIRED_MSA = "pairedMsa"
    UNPAIRED_MSA = "unpairedMsa"

class MSAFields(StrEnum):
    """ Fields for multiple sequence alignments (MSA) in the config dictionary used by modules in af_pipeline.af_input"""
    TYPE = auto()
    PATH = auto()

class MSAType(StrEnum):
    """ Types of multiple sequence alignments (MSA) for AlphaFold Server input."""
    PAIRED = "pairedMsa"
    UNPAIRED = "unpairedMsa"

class NucleicAcidStrand(StrEnum):
    """ Strand type for nucleic acids."""
    DOUBLE = auto()
    SINGLE = auto()

class GlycanModificationFields(StrEnum):
    """ Fields for glycan modifications in the config dictionary used by modules in af_pipeline.af_input"""
    RESIDUES = auto()
    POSITION = auto()

class ProteinModificationFields(StrEnum):
    """ Fields for protein modifications in the config dictionary used by modules in af_pipeline.af_input"""
    PTM_TYPE = "ptmType"
    PTM_POSITION = "ptmPosition"

class NucleicAcidModificationFields(StrEnum):
    """ Fields for nucleic acid modifications in the config dictionary used by modules in af_pipeline.af_input"""
    MODIFICATION_TYPE = "modificationType"
    BASE_POSITION = "basePosition"

class ConfigYaml(StrEnum):
    """ Keys for the config YAML file used by modules in af_pipeline.af_input."""
    AF_INPUT_JOBS = auto()
    PROTEIN_UNIPROT_MAP = auto()

class BestPredictionFields(StrEnum):
    """ Keys for the best AlphaFold prediction dictionary used by modules in af_pipeline.rank_predictions."""
    STRUCTURE_PATH = auto()
    DATA_PATH = auto()
    AF_OFFSET = auto()
    ENTITY_CHAIN_MAP = auto()

class AF3Metrics(StrEnum):
    """ Metrics for assessing the confidence of AlphaFold 3 predictions."""
    RANKING_SCORE = auto()
    IPTM = auto()
    PTM = auto()
    FRACTION_DISORDERED = auto()
    MODEL_PATH = auto()
    MODEL_IDX = auto()

class AF3SummaryConfidenceFields(StrEnum):
    """ Fields in the summary of confidence metrics for AlphaFold 3 predictions."""
    FRACTION_DISORDERED = auto()
    HAS_CLASH = auto()
    IPTM = auto()
    NUM_RECYCLES = auto()
    PTM = auto()
    RANKING_SCORE = auto()

ALLOWED_MSA_TYPES = list(MSAType)

ALLOWED_STRUCTURE_FORMATS = list(AVAILABLE_PARSERS.keys())

AVAILABLE_DATA_READERS = {
    FileFormat.PKL: read_pkl,
    FileFormat.JSON: read_json,
}

ALLOWED_DATA_FORMATS = list(AVAILABLE_DATA_READERS.keys())

class AtomQuantity(StrEnum):
    """ Quantities related to atoms in the Biopython structure object."""
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
    """ Quantities related to residues in the Biopython structure object."""
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
""" Allowed post-translational modifications (PTMs) by AlphaFold server. """

DNA_MOD_PURINES = [
    "CCD_6OG",
    "CCD_6MA",
    "CCD_8OG",
]
""" Allowed DNA modifications on purines by AlphaFold server. """

DNA_MOD_PYRIMIDINES = [
    "CCD_5CM",
    "CCD_C34",
    "CCD_5HC",
    "CCD_1CC",
    "CCD_5FC",
]
""" Allowed DNA modifications on pyrimidines by AlphaFold server. """
DNA_MOD_OTHER = [
    "CCD_3DR",
]
""" Allowed DNA modifications on other nucleotides by AlphaFold server. """

DNA_MOD = DNA_MOD_PURINES + DNA_MOD_PYRIMIDINES + DNA_MOD_OTHER

RNA_MOD_PURINES = [
    "CCD_A2M",
    "CCD_MA6",
    "CCD_6MZ",
    "CCD_2MG",
    "CCD_OMG",
    "CCD_7MG",
]
""" Allowed RNA modifications on purines by AlphaFold server. """

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
""" Allowed RNA modifications on pyrimidines by AlphaFold server. """

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
""" Allowed ligands by AlphaFold server. """

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
""" Allowed ions by AlphaFold server. """

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
""" Standard amino acid residues. """

DNA_PURINES = [
    "DA",
    "DG"
]
""" Standard DNA purines. """
RNA_PURINES = [
    "A",
    "G"
]
""" Standard RNA purines. """
RNA_PYRIMIDINES = [
    "U",
    "C"
]
""" Standard RNA pyrimidines. """
DNA_PYRIMIDINES = [
    "DC",
    "DT"
]
""" Standard DNA pyrimidines. """

REP_ATOMS = {
    "proteinChain": "CB",
    "is_ca_only": "CA",
    "is_purine": "C4",
    "is_pyrimidine": "C2",
}
""" Representative atoms for different entity types and residue types.
- For protein chains, the representative atom is CB (or CA for glycine).
- For nucleic acids, the representative atom is C4 for purines and C2 for pyrimidines. """

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

class ChainType(StrEnum):
    IDR = "IDR"
    R = "R"

class ChainAssessment(StrEnum):
    CHAIN_ID = "Chain ID"
    PROTEIN_NAME = "Protein Name"
    AVERAGE_PLDDT = "Average pLDDT"
    AVERAGE_IPLDDT = "Average ipLDDT"
    INTERFACE_RESIDUES = "Interface Residues"
    IS_INTERFACE_RESIDUE = "Is Interface Residue"
    CHAIN_TYPE = "Chain Type"
    RESIDUE_NUMBER = "Residue Number"

# as_average
CHAINWISE_ASSESSMENT_COLUMNS = {
    True: [
        ChainAssessment.CHAIN_ID,
        ChainAssessment.PROTEIN_NAME,
        ChainAssessment.AVERAGE_PLDDT,
        ChainAssessment.AVERAGE_IPLDDT,
        ChainAssessment.INTERFACE_RESIDUES,
        ChainAssessment.CHAIN_TYPE,
    ],
    False: [
        ChainAssessment.CHAIN_ID,
        ChainAssessment.PROTEIN_NAME,
        ChainAssessment.RESIDUE_NUMBER,
        ChainAssessment.AVERAGE_PLDDT,
        # ChainAssessment.AVERAGE_IPLDDT,
        ChainAssessment.IS_INTERFACE_RESIDUE,
        ChainAssessment.CHAIN_TYPE,
    ],
}

class ChainPairAssessment(StrEnum):
    CHAIN_ID = "Chain ID"
    CHAIN_ID_1 = "Chain ID 1"
    CHAIN_ID_2 = "Chain ID 2"
    PROTEIN_NAME = "Protein Name"
    PROTEIN_NAME_1 = "Protein Name 1"
    PROTEIN_NAME_2 = "Protein Name 2"
    CHAIN_TYPE = "Chain Type"
    CHAIN_TYPE_1 = "Chain Type 1"
    CHAIN_TYPE_2 = "Chain Type 2"
    INTERFACE_RESIDUES = "Interface Residues"
    NUMBER_OF_CONTACTS = "Number of contacts"
    AVERAGE_IPLDDT_CHAIN1 = "Average ipLDDT chain1"
    AVERAGE_IPLDDT_CHAIN2 = "Average ipLDDT chain2"
    AVERAGE_PAE = "Average PAE"
    AVERAGE_PAE_IJ = "Average PAE ij"
    AVERAGE_PAE_JI = "Average PAE ji"
    AVERAGE_IPAE = "Average iPAE"
    AVERAGE_IPAE_IJ = "Average iPAE ij"
    AVERAGE_IPAE_JI = "Average iPAE ji"
    RESIDUE_1 = "Residue 1"
    RESIDUE_2 = "Residue 2"
    PAE_IJ = "PAE ij"
    PAE_JI = "PAE ji"
    PAE = "PAE"
    PLDDT_1 = "pLDDT 1"
    PLDDT_2 = "pLDDT 2"

CHAIN_PAIRWISE_ASSESSMENT_COLUMNS = {
    (True, True): [
        ChainPairAssessment.PROTEIN_NAME,
        ChainPairAssessment.CHAIN_ID,
        ChainPairAssessment.CHAIN_TYPE,
        ChainPairAssessment.INTERFACE_RESIDUES,
        ChainPairAssessment.NUMBER_OF_CONTACTS,
        ChainPairAssessment.AVERAGE_IPLDDT_CHAIN1,
        ChainPairAssessment.AVERAGE_IPLDDT_CHAIN2,
        ChainPairAssessment.AVERAGE_PAE,
        ChainPairAssessment.AVERAGE_IPAE,
    ],
    (True, False): [
        ChainPairAssessment.PROTEIN_NAME,
        ChainPairAssessment.CHAIN_ID,
        ChainPairAssessment.CHAIN_TYPE,
        ChainPairAssessment.INTERFACE_RESIDUES,
        ChainPairAssessment.NUMBER_OF_CONTACTS,
        ChainPairAssessment.AVERAGE_IPLDDT_CHAIN1,
        ChainPairAssessment.AVERAGE_IPLDDT_CHAIN2,
        ChainPairAssessment.AVERAGE_PAE_IJ,
        ChainPairAssessment.AVERAGE_PAE_JI,
        ChainPairAssessment.AVERAGE_IPAE_IJ,
        ChainPairAssessment.AVERAGE_IPAE_JI,
    ],
    (False, True): [
        ChainPairAssessment.PROTEIN_NAME_1,
        ChainPairAssessment.CHAIN_ID_1,
        ChainPairAssessment.CHAIN_TYPE_1,
        ChainPairAssessment.PROTEIN_NAME_2,
        ChainPairAssessment.CHAIN_TYPE_2,
        ChainPairAssessment.CHAIN_ID_2,
        ChainPairAssessment.RESIDUE_1,
        ChainPairAssessment.RESIDUE_2,
        ChainPairAssessment.PLDDT_1,
        ChainPairAssessment.PLDDT_2,
        ChainPairAssessment.PAE,
    ],
    (False, False): [
        ChainPairAssessment.PROTEIN_NAME_1,
        ChainPairAssessment.CHAIN_ID_1,
        ChainPairAssessment.CHAIN_TYPE_1,
        ChainPairAssessment.PROTEIN_NAME_2,
        ChainPairAssessment.CHAIN_TYPE_2,
        ChainPairAssessment.CHAIN_ID_2,
        ChainPairAssessment.RESIDUE_1,
        ChainPairAssessment.RESIDUE_2,
        ChainPairAssessment.PLDDT_1,
        ChainPairAssessment.PLDDT_2,
        ChainPairAssessment.PAE_IJ,
        ChainPairAssessment.PAE_JI,
    ]
}

class OverallAssessment(StrEnum):
    NUMBER_OF_CHAINS = "Number of Chains"
    NUMBER_OF_INTERACTING_CHAIN_PAIRS = "Number of Interacting Chain Pairs"
    INTERFACE_RESIDUES = "Interface Residues"
    TOTAL_RESIDUES = "Total Residues"
    SEQUENCE_COVERAGE = "Sequence Coverage"
    NUMBER_OF_CONTACTS = "Number of Contacts"
    AVERAGE_PLDDT = "Average pLDDT"
    AVERAGE_IPLDDT = "Average ipLDDT"
    AVERAGE_IDR_IPLDDT = "Average IDR ipLDDT"
    AVERAGE_IPAE = "Average iPAE"
    AVERAGE_IPAE_IJ = "Average iPAE ij"
    AVERAGE_IPAE_JI = "Average iPAE ji"


OVERALL_ASSESSMENT_COLUMNS = {
    True: [
        OverallAssessment.NUMBER_OF_CHAINS, #: "num_chains",
        OverallAssessment.NUMBER_OF_INTERACTING_CHAIN_PAIRS, #: "num_interacting_chain_pairs",
        OverallAssessment.INTERFACE_RESIDUES, #: "num_interface_residues",
        OverallAssessment.TOTAL_RESIDUES, #: "num_total_residues",
        OverallAssessment.SEQUENCE_COVERAGE, #: "rb_coverage",
        OverallAssessment.NUMBER_OF_CONTACTS, #: "num_contacts",
        OverallAssessment.AVERAGE_PLDDT, #: "avg_plddt",
        OverallAssessment.AVERAGE_IPLDDT, #: "avg_iplddt",
        OverallAssessment.AVERAGE_IDR_IPLDDT, #: "avg_idr_iplddt",
        OverallAssessment.AVERAGE_IPAE, #: "avg_ipae",
    ],
    False: [
        OverallAssessment.NUMBER_OF_CHAINS, #: "num_chains",
        OverallAssessment.NUMBER_OF_INTERACTING_CHAIN_PAIRS, #: "num_interacting_chain_pairs",
        OverallAssessment.INTERFACE_RESIDUES, #: "num_interface_residues",
        OverallAssessment.TOTAL_RESIDUES, #: "num_total_residues",
        OverallAssessment.SEQUENCE_COVERAGE, #: "rb_coverage",
        OverallAssessment.NUMBER_OF_CONTACTS, #: "num_contacts",
        OverallAssessment.AVERAGE_PLDDT, #: "avg_plddt",
        OverallAssessment.AVERAGE_IPLDDT, #: "avg_iplddt",
        OverallAssessment.AVERAGE_IDR_IPLDDT, #: "avg_idr_iplddt",
        OverallAssessment.AVERAGE_IPAE_IJ, #: "avg_ipae_ij",
        OverallAssessment.AVERAGE_IPAE_JI, #: "avg_ipae_ji",
    ]
}


class MMSeqs2API:

    BASE_URL = "https://a3m.mmseqs.com"
    """ Base URL for the MMseqs2 API."""

    @staticmethod
    def get_ticket_url(use_pairing=False):
        """ Get the ticket URL for MMseqs2 API.

        Args:
            use_pairing (bool): Whether to use pairing mode. Default is False.

        Returns:
            str: The ticket URL.
        """

        if use_pairing:
            return f"{MMSeqs2API.BASE_URL}/ticket/pair"
        else:
            return f"{MMSeqs2API.BASE_URL}/ticket/msa"

    @staticmethod
    def get_status_url(ID):
        """ Get the status URL for MMseqs2 API.

        Args:
            ID (str): The ticket ID.

        Returns:
            str: The status URL.
        """

        return f"{MMSeqs2API.BASE_URL}/ticket/{ID}"

    @staticmethod
    def get_download_url(ID):
        """ Get the download URL for MMseqs2 API.

        Args:
            ID (str): The ticket ID.

        Returns:
            str: The download URL.
        """

        return f"{MMSeqs2API.BASE_URL}/result/download/{ID}"