import warnings
import Bio
import Bio.PDB
import Bio.PDB.Structure
import Bio.PDB.Residue
import Bio.PDB.Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select, PDBIO
from collections import defaultdict
from typing import Dict, Any, overload
from af_pipeline.constants.af_constants import *


def has_per_atom_token(residue: Bio.PDB.Residue.Residue) -> bool:
    """ Check if the residue has per-atom token.

    Assuming the residue is decorated with the appropriate attributes.
    - `is_modified`: True if the residue is modified.
    - `is_ligand`: True if the residue is a ligand.
    - `entityType`: None if the residue is not a protein, DNA, RNA, or ligand.

    Args:

        residue (Bio.PDB.Residue.Residue):
            Biopython residue object.

    Returns:

        condition (bool):
            True if the residue has per-atom token, False otherwise.
    """

    condition = (
        residue.xtra.get("is_modified")
        or residue.xtra.get("is_ligand")
        or residue.xtra.get("entityType") is None
    )

    return condition

def save_structure_obj(
    structure: Bio.PDB.Structure.Structure,
    out_file: str,
    res_select_obj: Select = Select(),
    save_type: str = "cif",
    preserve_header_footer = False,
):
    """Save the selection in Biopython structure object as a PDB or CIF file.

    Args:

        structure (Bio.PDB.Structure.Structure):
            Biopython Structure object to save.

        out_file (str):
            Path to the output file where the structure will be saved.

        res_select_obj (Bio.PDB.Select):
            Biopython Select object to filter the residues to save.
            Defaults to `Select()` which saves all residues.

        save_type (str):
            Type of file to save the structure as.
            Can be either "pdb" or "cif". Defaults to "cif".

        preserve_header_footer (bool):
            If True, the header and footer information from the structure
            will be preserved in the saved file. Defaults to False.
    """

    if save_type == "pdb":

        io = PDBIO()
        io.set_structure(structure)
        io.save(out_file, res_select_obj)

        if preserve_header_footer:
            warnings.warn(
                """

                PDB files do not support headers and footers.
                Saving without headers and footers.
                """
            )

    elif save_type == "cif":

        io = MMCIFIO()
        io.set_structure(structure)
        io.save(out_file, res_select_obj)

        if preserve_header_footer:

            header_footer = structure.xtra.get("header_footer", {})

            if {"header", "footer"}.issubset(set(header_footer.keys())):

                with open(out_file, "r") as f:
                    lines = f.readlines()

                lines.insert(0, "\n")
                lines.append("\n")

                for header_line in header_footer["header"][::-1]:
                    lines.insert(0, f"{header_line}")

                for footer_line in header_footer["footer"][::-1]:
                    lines.append(f"{footer_line}")

                with open(out_file, "w") as f:
                    for line in lines:
                        f.write(line)

            else:
                warnings.warn(
                    """

                    No header or footer information found in the structure.
                    Saving without headers and footers.
                    """
                )

def add_header_footer(
    structure: Bio.PDB.Structure.Structure,
    structure_file_path: str,
) -> Bio.PDB.Structure.Structure:
    """Add the header and footer information to the structure object.

    While parsing the cif file using Biopython MMCIFParser,
    the header and footer information information is not retained.
    This function extracts the header and footer information from the
    structure file and adds it to the structure object. \n

    The information is stored in the `xtra` attribute of the structure
    object under the key `header_footer`. \n

    Args:

        structure (Bio.PDB.Structure.Structure):
            Biopython Structure object.

        structure_file_path (str):
            path to the structure file.

    Returns:

        structure (Bio.PDB.Structure.Structure):
            Biopython Structure object with 'header_footer'.
    """

    with open(structure_file_path, "r") as f:
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

    structure.xtra["header_footer"] = {
        "header": header_info,
        "footer": footer_info,
    }

    return structure

def decorate_residue(
    residue: Bio.PDB.Residue.Residue,
    xtra_field: str | None = None,
    xtra_value: Any = None,
):
    """Decorate the residue.

    This function adds `entityType` attribute to the residue's `xtra`
    attribute depending on the type of entity. \n

    Allowed entity types are:
    - proteinChain
    - dnaSequence
    - rnaSequence
    - ligand
    - ion
    - None (if the residue does not belong to any of the above types)

    It also adds `is_modified`, `is_ca_only`, `is_ligand`, `is_ion`,
    `is_purine`, and `is_pyrimidine` attributes to the residue's `xtra`
    attribute wherever applicable. \n

    Args:

        residue (Bio.PDB.Residue.Residue):
            Biopython residue object.
    """

    symbol = residue.get_resname()

    if symbol in PROTEIN_ENTITIES:
        residue.xtra["entityType"] = "proteinChain"

        if symbol in ALLOWED_PTMS:
            residue.xtra["is_modified"] = True

        if symbol in ONLY_CA_RESIDUES:
            residue.xtra["is_ca_only"] = True

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

            The residue "{symbol}" does not belong to any known entity types.
            Setting 'entityType' to None.
            """
        )
        residue.xtra["entityType"] = None

    if symbol in PURINES:
        residue.xtra["is_purine"] = True

    elif symbol in PYRIMIDINES:
        residue.xtra["is_pyrimidine"] = True

    if xtra_field is not None and xtra_value is not None:
        if xtra_field in residue.xtra:
            warnings.warn(
                f"""

                The field '{xtra_field}' already exists in the residue's xtra.
                Overwriting the value.
                """
            )
        residue.xtra[xtra_field] = xtra_value

def decorate_atom(
    atom: Bio.PDB.Atom.Atom,
    xtra_field: str | None = None,
    xtra_value: Any = None,
):
    symbol = atom.get_name()
    residue = atom.get_parent()

    if not isinstance(residue, Bio.PDB.Residue.Residue):
        raise TypeError(
            f"""

            Expected a Bio.PDB.Residue.Residue object, got {type(residue)}
            """
        )

    residue = residue.get_resname()

    if residue in PROTEIN_ENTITIES and residue not in ONLY_CA_RESIDUES:
        if symbol == "CB":
            atom.xtra["is_representative"] = True

    elif residue in ONLY_CA_RESIDUES:
        if symbol == "CA":
            atom.xtra["is_representative"] = True

    elif residue in PURINES:
        if symbol == "C4":
            atom.xtra["is_representative"] = True

    elif residue in PURINES:
        if symbol == "C2":
            atom.xtra["is_representative"] = True

    else:
        atom.xtra["is_representative"] = False

    if xtra_field is not None and xtra_value is not None:
        if xtra_field in atom.xtra:
            warnings.warn(
                f"""

                The field '{xtra_field}' already exists in the atom's xtra.
                Overwriting the value.
                """
            )
        atom.xtra[xtra_field] = xtra_value


class RenumberResidues:
    """Class to renumber the residues based on the offset.

    Attributes:

        af_offset (dict):
            Offset describing start and end residue number for each chain in
            the predicted structure.
    """

    def __init__(self, af_offset: dict = {}):
        """Initialize the RenumberResidues class.

        Args:

            af_offset (dict, optional):
                Offset describing start and end residue number for each chain
                in the predicted structure.
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

        if chain_id in self.af_offset:
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

            if chain_id in self.af_offset:

                start = start - (self.af_offset[chain_id][0] - 1)
                end = end - (self.af_offset[chain_id][0] - 1)

            renumbered_region_of_interest[chain_id] = [start, end]

        return renumbered_region_of_interest

    def residue_map(
        self,
        token_chain_ids: list,
        token_res_ids: list,
        token_atom_names: list,
    ):
        """Create a map of residue indices to residue numbers and vice-versa.

        res_idx is essentially token index. \n
        res_num is the residue number. \n
        res_num = res_idx + 1 if af_offset is not provided. \n
        res_num = res_idx + af_offset if af_offset is provided. \n
        af_offset informs what is the starting residue number for each chain.

        Args:

            token_chain_ids (list):
                Token chain IDs.

            token_res_ids (list):
                Token residue IDs.

            token_atom_names (list):
                Token atom names.

        Returns:

            tuple (idx_to_num, num_to_idx):
                idx_to_num (Dict):
                    Dictionary mapping residue indices to residue numbers.

                num_to_idx (Dict):
                    Dictionary mapping residue numbers to residue indices.
        """

        idx_to_num = {}
        num_to_idx = defaultdict(dict)

        for token_idx, (chain_id, res_num, atom_name) in enumerate(
            zip(token_chain_ids, token_res_ids, token_atom_names)
        ):

            res_num = self.renumber_chain_res_num(
                chain_res_num=res_num,
                chain_id=chain_id,
            )

            idx_to_num[token_idx] = {
                "chain_id": chain_id,
                "res_num": res_num,
                "atom_name": atom_name,
            }

            if res_num not in num_to_idx[chain_id]:
                num_to_idx[chain_id][res_num] = {
                    atom_name: token_idx
                }
            else:
                num_to_idx[chain_id][res_num][atom_name] = token_idx

        return idx_to_num, num_to_idx


class ResidueSelect(Select):
    """Class to select residues in the structure based on the input dictionary.

    Attributes:
        confident_residues (Dict):
            Dictionary containing the chain ID as key and a list of residue
            numbers as value. \n
            e.g. {"A": [1, 2, 3], "B": [4, 5, 6]}
    """

    def __init__(
        self,
        confident_residues: Dict
    ):
        """Initialize the ResidueSelect class \n

        Args:

            confident_residues (Dict):
                Dictionary containing the chain ID as key and a list of residue numbers as value.
        """
        self.confident_residues = confident_residues

    @overload
    def accept_residue(
        self,
        residue: Bio.PDB.Residue.Residue
    ) -> bool:
        """Accept the residue if it's in `self.confident_residues`.

        Args:

            residue (Bio.PDB.Residue.Residue):
                Biopython residue object.

        Returns:

            bool:
                True if the residue is in the `self.confident_residues`.
        """

        chain = residue.parent
        if not isinstance(chain, Chain):
            raise TypeError(
                f"Expected a Bio.PDB.Chain.Chain object, got {type(chain)}"
            )
        chain_id = chain.id

        return residue.id[1] in self.confident_residues[chain_id]

    @overload
    def accept_residue(self, residue):
        """Overload this to reject residues for output."""
        return 1

    def accept_residue(
        self,
        residue: Bio.PDB.Residue.Residue
    ) -> bool:
        """Accept the residue if it's in `self.confident_residues`.

        Args:

            residue (Bio.PDB.Residue.Residue):
                Biopython residue object.

        Returns:

            bool:
                True if the residue is in the `self.confident_residues`.
        """

        chain = residue.parent
        if not isinstance(chain, Chain):
            raise TypeError(
                f"Expected a Bio.PDB.Chain.Chain object, got {type(chain)}"
            )
        chain_id = chain.id

        return residue.id[1] in self.confident_residues[chain_id]