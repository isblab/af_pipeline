# some constants

MAX_TEMPLATE_DATE = "2021-09-30"

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

