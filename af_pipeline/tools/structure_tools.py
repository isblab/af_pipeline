"""
Tools to work with structure files
==================================
- This module provides utility functions and classes to work with structure files (PDB or CIF)
- Uses `Biopython` for structure manipulation and parsing.
- The structure object refers to `Bio.PDB.Structure.Structure` object from Biopython.
- The residue object refers to `Bio.PDB.Residue.Residue` object from Biopython.
Hence, the "residue" term is used for both amino acids and nucleotides.
"""
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
from af_pipeline.constants.af_constants import (
    PROTEIN_ENTITIES,
    DNA_ENTITIES,
    RNA_ENTITIES,
    ION,
    PURINES,
    PYRIMIDINES,
    ALLOWED_PTMS,
    ALLOWED_DNA_MODS,
    ALLOWED_RNA_MODS,
    ALLOWED_LIGANDS,
    ONLY_CA_RESIDUES,
)


def has_per_atom_token(residue: Bio.PDB.Residue.Residue) -> bool:
    """ Check if the residue has per-atom token.

    TODO: add reference to tokenization in AF3.

    Assuming the `residue` is decorated with the appropriate attributes.
    - `is_modified`: `True` if the `residue` is modified.
    - `is_ligand`: `True` if the `residue` is a ligand.
    - `entityType`: `None` if the `residue` is not a protein, DNA, RNA, or ligand.

    Args:
        residue (Bio.PDB.Residue.Residue):
            Biopython residue object.

    Returns:
        `condition (bool)`:
            `True` if the residue has per-atom token, `False` otherwise.
    """

    condition = (
        residue.xtra.get("is_modified")
        or residue.xtra.get("is_ligand")
        or residue.xtra.get("entityType") is None
    )

    if isinstance(condition, bool) is False:
        raise TypeError(
            f"""

            Expected a boolean value for condition, got {type(condition)}
            """
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
            Biopython structure object to save.

        out_file (str):
            Path to the output file where the structure will be saved.

        res_select_obj (Bio.PDB.Select):
            Biopython Select object to filter the residues to save.
            Defaults to `Select()` which saves all residues.

        save_type (str):
            Type of file to save the structure as.
            Can be either "pdb" or "cif". Defaults to "cif".

        preserve_header_footer (bool):
            If `True`, the header and footer information from the structure
            will be preserved in the saved file. Defaults to `False`.
            Only applicable for CIF files.
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

    While parsing the cif file using Biopython `MMCIFParser`,
    the header and footer information information is not retained.
    This function extracts the header and footer information from the
    structure file and adds it to the structure object. \n

    The information is stored in the `xtra` attribute of the structure
    object under the key `header_footer`. \n

    Args:

        structure (Bio.PDB.Structure.Structure):
            Biopython structure object.

        structure_file_path (str):
            Path to the structure file.

    Returns:

        `structure (Bio.PDB.Structure.Structure)`:
            Biopython Structure object with `header_footer` added in `xtra`.
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
    attribute depending on the type of entity.

    Allowed entity types are:
    ```python
    - "proteinChain"
    - "dnaSequence"
    - "rnaSequence"
    - "ligand"
    - "ion"
    - None #(if the residue does not belong to any of the above types)
    ```

    It also adds the following flags to the residue's `xtra` wherever applicable -
    ```python
    - "is_modified" # True if the residue is a modified residue.
    - "is_ca_only" # True if the residue has only CA atom and no CB atom.
    - "is_ligand" # True if the residue is a ligand.
    - "is_ion" # True if the residue is an ion.
    - "is_purine" # True if the nucleotide is a purine.
    - "is_pyrimidine" # True if the nucleotide is a pyrimidine.
    ```

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
    """ Decorate the atom.

    This function adds `is_representative` attribute to the atom's `xtra`

    Args:
        atom (Bio.PDB.Atom.Atom): Biopython atom object.
        xtra_field (str | None, optional): field to add to the atom's `xtra`.
        xtra_value (Any, optional): value of the field to add to the atom's `xtra`.
    """

    symbol = atom.get_name()
    residue = atom.get_parent()

    if not isinstance(residue, Bio.PDB.Residue.Residue):
        raise TypeError(
            f"""

            Expected a Bio.PDB.Residue.Residue object, got {type(residue)}
            """
        )

    residue = residue.get_resname()

    if residue in PROTEIN_ENTITIES and residue not in ONLY_CA_RESIDUES and symbol == "CB":
        atom.xtra["is_representative"] = True

    elif residue in ONLY_CA_RESIDUES and symbol == "CA":
        atom.xtra["is_representative"] = True

    elif residue in PURINES and symbol == "C4":
        atom.xtra["is_representative"] = True

    elif residue in PURINES and symbol == "C2":
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

    If the input sequence to the AlphFold does not start from the first residue,
    the residue numbering in the predicted structure will be different from
    the original (UniProt) numbering.

    This class renumbers the residues in the structure based on the provided `offset`.

    The `offset` is a dictionary with chain IDs as keys and a list of two integers
    defining the `start` and `end` residue numbers for that chain in the predicted structure.

    NOTE: Although written with the purpose for AlphFold predicted structures,
    this class can be used to renumber residues in any structure file.
    """

    offset: dict
    """ Offset describing start and end residue number for each chain in
    the predicted structure.\n
    example: `{'A': [1, 100], 'B': [101, 200]}`."""

    def __init__(self, offset: dict = {}):

        self.offset = offset

    def renumber_structure(
        self,
        structure: Bio.PDB.Structure.Structure,
    ):
        """Renumber the residues in the structure based on the offset.

        Args:
            structure (Bio.PDB.Structure.Structure):
                Biopython structure object.

        Returns:
            `structure (Bio.PDB.Structure.Structure)`:
                Biopython structure object with renumbered residues.
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

        Given a residue index (1-indexed) in a chain, this function renumbers it
        based on the offset provided for that chain.

        Args:

            chain_res_num (int):
                Residue index (1-indexed) within the chain in the predicted structure.

            chain_id (str):
                Chain ID of the residue.

        Returns:
            `chain_res_num (int)`: Renumbered residue number
        """

        if chain_id in self.offset:
            chain_res_num += self.offset[chain_id][0] - 1

        return chain_res_num

    def renumber_region_of_interest(
        self,
        region_of_interest: Dict[str, list],
    ):
        """Offset the region of interest to the AF2/3 numbering.

        Region of interest is defined by the user is as per the original numbering
        (UniProt in case of proteins). \n
        However, if the prediction is done on a fragment of the protein, the
        residue numbering in the predicted structure will be different compared to
        the one provided by the user. \n
        This function renumbers the region of interest to the numbering of the
        predicted structure based on the `offset`. \n

        Args:
            region_of_interest (dict):
                Dictionary containing the region of interest for each chain.

        Returns:
            `renumbered_region_of_interest (dict)`:
                Dictionary containing the renumbered region of interest for
                each chain.
        """

        renumbered_region_of_interest = {}

        for chain_id in region_of_interest:

            start, end = region_of_interest[chain_id]

            if chain_id in self.offset:

                start = start - (self.offset[chain_id][0] - 1)
                end = end - (self.offset[chain_id][0] - 1)

            renumbered_region_of_interest[chain_id] = [start, end]

        return renumbered_region_of_interest

    def residue_map(
        self,
        token_chain_ids: list,
        token_res_ids: list,
        token_atom_names: list,
    ):
        """Create a map of residue indices to residue numbers and vice-versa.

        `res_idx` is essentially the token index. \n
        `res_num` is the residue number in the chain which will be same as UniProt
        numbering if the `offset` is provided correctly or full-length sequences
        were used in the prediction. \n
        `offset` informs what is the starting residue number for each chain.

        For example, if the input is as follows: \n
        ```python
        token_chain_ids = ['A', 'A', 'A', 'B', 'B']
        token_res_ids = [1, 1, 2, 1, 2]
        token_atom_names = ['N', 'CA', 'N', 'N', 'CA']
        ```
        The output will be: \n
        ```python
        idx_to_num = {
            0: {'chain_id': 'A', 'res_num': 1, 'atom_name': 'N'},
            1: {'chain_id': 'A', 'res_num': 1, 'atom_name': 'CA'},
            2: {'chain_id': 'A', 'res_num': 2, 'atom_name': 'N'},
            3: {'chain_id': 'B', 'res_num': 1, 'atom_name': 'N'},
            4: {'chain_id': 'B', 'res_num': 2, 'atom_name': 'CA'},
        }
        num_to_idx = {
            'A': {
                1: {'N': 0, 'CA': 1},
                2: {'N': 2},
            },
            'B': {
                1: {'N': 3},
                2: {'CA': 4},
            },
        }
        ```

        If the offset is provided as: \n
        ```python
        offset = {
            'A': [10, 11],
            'B': [30, 31],
        }
        ```
        The output will be: \n
        ```python
        idx_to_num = {
            0: {'chain_id': 'A', 'res_num': 10, 'atom_name': 'N'},
            1: {'chain_id': 'A', 'res_num': 10, 'atom_name': 'CA'},
            2: {'chain_id': 'A', 'res_num': 11, 'atom_name': 'N'},
            3: {'chain_id': 'B', 'res_num': 30, 'atom_name': 'N'},
            4: {'chain_id': 'B', 'res_num': 31, 'atom_name': 'CA'},
        }
        num_to_idx = {
            'A': {
                10: {'N': 0, 'CA': 1},
                11: {'N': 2},
            },
            'B': {
                30: {'N': 3},
                31: {'CA': 4},
            },
        }
        ```

        Args:
            token_chain_ids (list):
                Token chain IDs. Same as provided in the AF3 JSON output.

            token_res_ids (list):
                Token residue IDs. Same as provided in the AF3 JSON output.

            token_atom_names (list):
                Token atom names.

        Returns:
            `tuple (idx_to_num, num_to_idx)`:\n
                `idx_to_num (Dict)`:
                    Dictionary mapping residue indices to residue numbers.

                `num_to_idx (Dict)`:
                    Dictionary mapping residue numbers to residue indices.
        """

        idx_to_num = {}
        num_to_idx = defaultdict(dict)

        for token_idx, (chain_id, token_num, atom_name) in enumerate(
            zip(token_chain_ids, token_res_ids, token_atom_names)
        ):

            token_num = self.renumber_chain_res_num(
                chain_res_num=token_num,
                chain_id=chain_id,
            )

            idx_to_num[token_idx] = {
                "chain_id": chain_id,
                "token_num": token_num,
                "atom_name": atom_name,
            }

            if token_num not in num_to_idx[chain_id]:
                num_to_idx[chain_id][token_num] = {
                    atom_name: token_idx
                }
            else:
                num_to_idx[chain_id][token_num][atom_name] = token_idx

        return idx_to_num, num_to_idx


class ResidueSelect(Select):
    """Class to select residues in the structure based on the input dictionary.

    This is achieved by overriding the `Bio.PDB.Select.accept_residue` method
    from Biopython.
    """

    confident_residues: Dict
    """ Dictionary containing the chain ID as key and a list of residue
    numbers as value.\n
    e.g. `{"A": [1, 2, 3], "B": [4, 5, 6]}`
    """

    def __init__(
        self,
        confident_residues: Dict
    ):
        self.confident_residues = confident_residues

    @overload
    def accept_residue(
        self,
        residue: Bio.PDB.Residue.Residue
    ) -> bool:
        """Accept the residue if it's in `confident_residues`.

        Args:
            residue (Bio.PDB.Residue.Residue):
                Biopython residue object.

        Returns:
            `bool`: `True` if the residue is in the `confident_residues`.
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
        """Accept the residue if it's in `confident_residues`.

        Args:
            residue (Bio.PDB.Residue.Residue):
                Biopython residue object.

        Returns:
            `bool`: `True` if the residue is in the `confident_residues`.
        """

        chain = residue.parent
        if not isinstance(chain, Chain):
            raise TypeError(
                f"Expected a Bio.PDB.Chain.Chain object, got {type(chain)}"
            )
        chain_id = chain.id

        return residue.id[1] in self.confident_residues[chain_id]