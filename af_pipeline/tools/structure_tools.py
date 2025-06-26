from collections import defaultdict
from typing import Dict
import warnings
import Bio
import Bio.PDB
import Bio.PDB.Structure
import Bio.PDB.Residue
from af_pipeline.constants.af_constants import *

def add_header_footer(
    structure: Bio.PDB.Structure.Structure,
    struct_file_path: str,
) -> Bio.PDB.Structure.Structure:
    """Add the header and footer information to the structure object.

    Args:

        structure (Bio.PDB.Structure.Structure):
            Biopython Structure object.

        struct_file_path (str):
            path to the structure file.

    Returns:

        structure (Bio.PDB.Structure.Structure):
            Biopython Structure object with 'header_footer'.
    """

    with open(struct_file_path, "r") as f:
        lines = f.readlines()

    header_info = []
    header_section = ""

    for line in lines:
        header_section += line

        if line.startswith("#"):
            header_info.append(header_section)
            header_section = ""

        if line.startswith("_atom_site"):
            break

    footer_info = []
    footer_section = ""

    for line in lines[::-1]:
        footer_section = line + footer_section

        if line.startswith("#"):
            footer_info.append(footer_section)
            footer_section = ""

        if line.startswith("ATOM") or line.startswith("HETATM"):
            break

    structure.header_footer = {
        "header": header_info,
        "footer": footer_info,
    }

    return structure

def decorate_residue(
    residue: Bio.PDB.Residue.Residue
):
    """Decorate the residue with entity type based on its symbol.

    #! decorate residue if modified

    Adds `entityType` to the residue's `xtra` attribute. \n

    Args:

        residue (Bio.PDB.Residue.Residue):
            Biopython residue object.
    """

    symbol = residue.get_resname()

    if symbol in PROTEIN_ENTITIES:
        residue.xtra["entityType"] = "proteinChain"

        if symbol in ALLOWED_PTMS:
            residue.xtra["is_modified"] = True

    elif symbol in DNA_ENTITIES:
        residue.xtra["entityType"] = "dnaSequence"

        if symbol in ALLOWED_DNA_MODS:
            residue.xtra["is_modified"] = True

    elif symbol in RNA_ENTITIES:
        residue.xtra["entityType"] = "rnaSequence"

        if symbol in ALLOWED_RNA_MODS:
            residue.xtra["is_modified"] = True

    elif symbol in ALLOWED_LIGANDS:
        residue.xtra["entitiyType"] = "ligand"
        residue.xtra["is_ligand"] = True


    elif symbol in ION:
        residue.xtra["entityType"] = "ion"
        residue.xtra["is_ion"] = True

    else:
        warnings.warn(
            f"""
            The residue {symbol} does not belong to any known entity types.
            Setting 'entityType' to None.
            """
        )
        residue.xtra["entityType"] = None


class RenumberResidues:
    """Class to renumber the residues based on the offset.

    Attributes:

        af_offset (Dict | None):
            Offset describing start and end residue number for each chain in
            the predicted structure.
    """

    def __init__(self, af_offset: Dict | None = None):
        """Initialize the RenumberResidues class.

        Args:

            af_offset (Dict | None, optional):
                Offset describing start and end residue number for each chain in the predicted structure.
                (Default: None) \n
                example: `{'A': [1, 100], 'B': [101, 200]}`. \n
                If None, no renumbering is done.
        """

        self.af_offset = af_offset

    def renumber_structure(
        self,
        structure: Bio.PDB.Structure.Structure,
    ):
        """Renumber the residues in the structure based on the offset.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

        Returns:
            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object with renumbered residues.
        """

        for model in structure:
            for chain in model:
                chain_id = chain.id
                for residue in chain:
                    h, num, ins = residue.id

                    num = self.renumber_chain_res_num(
                        chain_res_num=num,
                        chain_id=chain_id,
                    )

                    residue.id = (h, num, ins)

        return structure

    def renumber_chain_res_num(
        self,
        chain_res_num: int,
        chain_id: str,
    ):
        """Renumber the residue number based on the offset.

        Args:

            chain_res_num (int):
                Residue number in the predicted structure.

            af_offset (dict):
                Offset describing start and end residue number for each chain
                in the predicted structure.
                e.g. `{'A': [1, 100], 'B': [101, 200]}`. \n
                If `None`, no renumbering is done.

        Returns:

            chain_res_num (int):
                renumbered residue number
        """

        if self.af_offset is not None and chain_id in self.af_offset:
            chain_res_num += self.af_offset[chain_id][0] - 1

        return chain_res_num

    def renumber_region_of_interest(
        self,
        region_of_interest: Dict,
    ):
        """Offset the interacting region to the AF2/3 numbering.

        Interacting region defined by user is as per the original numbering
        (UniProt in case of proteins). \n
        However, if the prediction is done on a fragment of the protein, the
        numbering will be different. \n
        This function offsets the interacting region to the numbering of the
        predicted structure. \n
        By default, the offset is assumed to be 0.

        Args:

            region_of_interest (Dict):
                Dictionary containing the region of interest for each chain.

        Returns:

            renumbered_region_of_interest (Dict):
                Dictionary containing the renumbered region of interest for
                each chain.
        """

        renumbered_region_of_interest = {}

        for chain_id in region_of_interest:

            start, end = region_of_interest[chain_id]

            if self.af_offset and chain_id in self.af_offset:

                start = start - (self.af_offset[chain_id][0] - 1)
                end = end - (self.af_offset[chain_id][0] - 1)

            renumbered_region_of_interest[chain_id] = [start, end]

        return renumbered_region_of_interest

    def residue_map(
        self,
        token_chain_ids: list,
        token_res_ids: list
    ):
        """Create a map of residue indices to residue numbers and vice-versa.

        res_idx is essentially token index. \n
        res_num is the residue number. \n
        res_num = res_idx + 1 if af_offset is not provided. \n
        res_num = res_idx + af_offset if af_offset is provided. \n
        af_offset informs what is the starting residue number for each chain.

        Args:

            token_chain_ids (list):
                Tokenized chain IDs.

            token_res_ids (list):
                Tokenized residue IDs.

        Returns:

            tuple (idx_to_num, num_to_idx):
                idx_to_num (Dict):
                    Dictionary mapping residue indices to residue numbers.

                num_to_idx (Dict):
                    Dictionary mapping residue numbers to residue indices.
        """

        idx_to_num = {}
        num_to_idx = defaultdict(dict)

        for res_idx, (chain_id, res_num) in enumerate(
            zip(token_chain_ids, token_res_ids)
        ):

            res_num = self.renumber_chain_res_num(
                chain_res_num=res_num,
                chain_id=chain_id,
            )
            print(res_idx, res_num, chain_id)

            idx_to_num[res_idx] = {
                "chain_id": chain_id,
                "res_num": res_num,
            }

            if res_num not in num_to_idx[chain_id]:
                num_to_idx[chain_id][res_num] = [res_idx]
            else:
                num_to_idx[chain_id][res_num].append(res_idx)

        return idx_to_num, num_to_idx
