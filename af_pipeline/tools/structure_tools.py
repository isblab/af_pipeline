"""
[structure_tools](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/tools/structure_tools.py)
==================================
- This module provides utility functions and classes to work with structure files (PDB or CIF).
- Uses `Biopython` for structure manipulation and parsing.
- The structure object refers to `Bio.PDB.Structure.Structure`.
- The residue object refers to `Bio.PDB.Residue.Residue`.<br />
  Hence, the "residue" term is used for amino acids, nucleotides, ions, ligands.
"""
from textwrap import dedent
import warnings
import Bio
import Bio.PDB
import Bio.PDB.Structure
import Bio.PDB.Residue
import Bio.PDB.Atom
import numpy as np
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
    AtomDecoration,
    EntityType,
    FileFormat,
    InteractionMapType,
    ResidueDecoration,
    ResidueMapDepth,
    VALID_INTERACTION_MAP_TYPES,
    ResidueMapKeys,
    TokenLevel,
)

def get_interaction_map(
    coords1: np.ndarray,
    coords2: np.ndarray,
    contact_threshold: float,
    map_type: InteractionMapType,
) -> np.ndarray:
    """ Create an interaction map, given the input coordinates.

    Returns a distance map or a contact map, based on the map_type specified.

    ## Arguments:

    - **coords1 (np.ndarray)**:<br />
        Coordinates of shape (N, 3) for the first set of atoms/residues.

    - **coords2 (np.ndarray)**:<br />
        Coordinates of shape (M, 3) for the second set of atoms/residues.

    - **contact_threshold (float)**:<br />
        Distance threshold to define a contact.
        Only used if `map_type` is `InteractionMapType.CONTACT`.

    - **map_type (InteractionMapType)**:<br />
        Type of interaction map to create.
        Can be either `InteractionMapType.DISTANCE` or `InteractionMapType.CONTACT`.

         - `InteractionMapType.DISTANCE`:
            Returns a distance map of shape (N, M)
            where each element (i, j) is the distance between the i-th atom/residue
            in `coords1` and the j-th atom/residue in `coords2`.

         - `InteractionMapType.CONTACT`:
            Returns a binary contact map of shape (N, M)
            where each element (i, j) is 1 if the distance between the i-th atom/residue
            in `coords1` and the j-th atom/residue in `coords2` is less than or equal to
            `contact_threshold`, and 0 otherwise.

    ## Returns:

    - **np.ndarray**:<br />
        The interaction map as a numpy array.
    """

    if map_type not in VALID_INTERACTION_MAP_TYPES:
        raise ValueError(
            f"Invalid map_type specified."\
            f"Expected one of {VALID_INTERACTION_MAP_TYPES}, got {map_type}"
        )

    distance_map = get_distance_map(coords1, coords2)

    if map_type == InteractionMapType.DISTANCE:
        return distance_map

    elif map_type == InteractionMapType.CONTACT:
        contact_map = get_contact_map(distance_map, contact_threshold)
        return contact_map

def get_distance_map(coords1: np.ndarray, coords2: np.ndarray):
    """ Create an all-v-all distance map.

    Returns a matrix of distances between all pairs of atoms/residues in the
    two sets of coordinates.

    ## Arguments:

    - **coords1 (np.ndarray)**:<br />
        Coordinates of shape (N, 3) for the first set of atoms/residues.

    - **coords2 (np.ndarray)**:<br />
        Coordinates of shape (M, 3) for the second set of atoms/residues.

    ## Returns:

    - **np.ndarray**:<br />
        A matrix of shape (N, M) where each element (i, j) is the Euclidean distance
        between the i-th atom/residue in `coords1` and the j-th atom/residue in `coords2`.
    """

    from scipy.spatial import distance_matrix

    distance_map = distance_matrix(coords1, coords2)

    return distance_map

def get_contact_map(distance_map: np.ndarray, contact_threshold: float):
    """ Given the distance map, create a binary contact map by thresholding distances.

    Returns a binary matrix, where 1 indicates a contact and 0 indicates no contact.

    ## Arguments:

    - **distance_map (np.ndarray)**:<br />
        A matrix of shape (N, M) where each element (i, j) is the Euclidean distance
        between the i-th atom/residue in `coords1` and the j-th atom/residue in `coords2`.

    - **contact_threshold (float)**:<br />
        Distance threshold to define a contact. If the distance between two atoms/residues
        is less than or equal to this threshold, they are considered to be in contact.

    ## Returns:

    - **np.ndarray**:<br />
        A binary matrix of shape (N, M) where each element (i, j) is 1 if the distance
        between the i-th atom/residue in `coords1` and the j-th atom/residue in `coords2`
        is less than or equal to `contact_threshold`, and 0 otherwise.
    """

    contact_map = np.where(distance_map <= contact_threshold, 1, 0)

    return contact_map

def get_token_level(residue: Bio.PDB.Residue.Residue) -> TokenLevel:
    """ Get the token level of the residue.

    Assuming the `residue` object is decorated with the appropriate attributes.
    - `tokenLevel`: `TokenLevel.ATOM` if the residue has per-atom token, `TokenLevel.RESIDUE` otherwise.

    Arguments:

    - **residue (Bio.PDB.Residue.Residue)**:<br />
        Biopython residue object.

    Returns:

    - **token_level (TokenLevel)**:<br />
        Token level of the residue.
    """

    token_level = residue.xtra.get(
        ResidueDecoration.TOKEN_LEVEL, TokenLevel.RESIDUE
    )

    return token_level

def has_modifications(structure: Bio.PDB.Structure.Structure) -> bool:
    """ Check if the structure has any modified residues.

    Assuming the `residue` objects are decorated with the appropriate attributes.
    - `is_modified`: `True` if the `residue` is modified.

    Arguments:

    - **structure (Bio.PDB.Structure.Structure)**:<br />
        Biopython structure object.

    Returns:

    - **condition (bool)**:<br />
        `True` if the structure has any modified residues, `False` otherwise.
    """

    for residue in structure.get_residues():
        if residue.xtra.get(ResidueDecoration.IS_MODIFIED, False):
            return True

    return False

def save_structure_obj(
    structure: Bio.PDB.Structure.Structure,
    out_file: str,
    res_select_obj: Select = None,
    save_type: str = FileFormat.CIF,
    preserve_header_footer = False,
):
    """Save the selection in Biopython structure object as a PDB or CIF file.

    Arguments:

    - **structure (Bio.PDB.Structure.Structure)**:<br />
        Biopython structure object to save.

    - **out_file (str)**:<br />
        Path to the output file where the structure will be saved.

    - **res_select_obj (Bio.PDB.Select, optional)**:<br />
        Biopython Select object to filter the residues to save.
        Defaults to `Select()` which saves all residues.

    - **save_type (str, optional)**:<br />
        Type of file to save the structure as.
        Can be either "pdb" or "cif".

    - **preserve_header_footer (bool, optional)**:<br />
        If `True`, the header and footer information from the structure
        will be preserved in the saved file.
        > [!NOTE]
        > The header and footer information can only be preserved for CIF files.
    """

    res_select_obj = res_select_obj or Select()

    if save_type == FileFormat.PDB:

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

    elif save_type == FileFormat.CIF:

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

    While parsing the `CIF` file using Biopython `MMCIFParser`,
    the header and footer information information is not retained.<br />
    This function extracts the header and footer information from the
    structure file and adds it to the structure object. \n

    The information is stored in the `xtra` attribute of the structure
    object under the key `header_footer`. \n

    > [!NOTE]
    > - `PDB` format does not have header and footer information.
    > - Using this function leads to parsing of the structure file twice. Once
    >   to get the structure object and once to extract the header and footer information.

    Arguments:

    - **structure (Bio.PDB.Structure.Structure)**:<br />
        Biopython structure object.

    - **structure_file_path (str)**:<br />
        Path to the structure file.

    Returns:

    - **structure (Bio.PDB.Structure.Structure)**:<br />
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
    - "is_modified" # boolean indicating if the residue is modified
    - "is_ca_only" # boolean indicating if the residue has only CA atom
    - "is_purine" # boolean indicating if the nucleotide is purine
    - "is_pyrimidine" # boolean indicating if the nucleotide is pyrimidine
    - "is_ion" # boolean indicating if the residue is an ion.
    - "is_ligand" # boolean indicating if the residue is a ligand.
    ```

    Arguments:

    - **residue (Bio.PDB.Residue.Residue)**:<br />
        Biopython residue object.
    """

    symbol = residue.get_resname()

    if symbol in PROTEIN_ENTITIES:
        residue.xtra[ResidueDecoration.ENTITY_TYPE] = EntityType.PROTEIN_CHAIN

        if symbol in ALLOWED_PTMS:
            residue.xtra[ResidueDecoration.IS_MODIFIED] = True

        if symbol in ONLY_CA_RESIDUES:
            residue.xtra[ResidueDecoration.IS_CA_ONLY] = True

    elif symbol in DNA_ENTITIES:
        residue.xtra[ResidueDecoration.ENTITY_TYPE] = EntityType.DNA_SEQUENCE

        if symbol in ALLOWED_DNA_MODS:
            residue.xtra[ResidueDecoration.IS_MODIFIED] = True
    elif symbol in RNA_ENTITIES:
        residue.xtra[ResidueDecoration.ENTITY_TYPE] = EntityType.RNA_SEQUENCE

        if symbol in ALLOWED_RNA_MODS:
            residue.xtra[ResidueDecoration.IS_MODIFIED] = True

    elif symbol in ALLOWED_LIGANDS:
        residue.xtra[ResidueDecoration.ENTITY_TYPE] = EntityType.LIGAND

    elif symbol in ION:
        residue.xtra[ResidueDecoration.ENTITY_TYPE] = EntityType.ION

    else:
        warnings.warn(dedent(f"""
            The residue "{symbol}" does not belong to any known entity types.
            Setting 'entityType' to None.""")
        )
        residue.xtra[ResidueDecoration.ENTITY_TYPE] = None

    if symbol in PURINES:
        residue.xtra[ResidueDecoration.IS_PURINE] = True

    elif symbol in PYRIMIDINES:
        residue.xtra[ResidueDecoration.IS_PYRIMIDINE] = True

    if xtra_field is not None and xtra_value is not None:
        if xtra_field in residue.xtra:
            warnings.warn(dedent(f"""
                The field '{xtra_field}' already exists in the residue's xtra.
                Overwriting the value.""")
            )
        residue.xtra[xtra_field] = xtra_value

    condition = (
        residue.xtra.get(ResidueDecoration.IS_MODIFIED)
        or residue.xtra.get(ResidueDecoration.ENTITY_TYPE) in [EntityType.LIGAND, None]
    )
    _get_token_level = {True: TokenLevel.ATOM, False: TokenLevel.RESIDUE}
    residue.xtra[ResidueDecoration.TOKEN_LEVEL] = _get_token_level[condition]

def decorate_atom(
    atom: Bio.PDB.Atom.Atom,
    xtra_field: str | None = None,
    xtra_value: Any = None,
):
    """ Decorate the atom.

    This function adds `is_representative` attribute to the atom's `xtra`

    Arguments:

    - **atom (Bio.PDB.Atom.Atom)**:<br />
        Biopython atom object.

    - **xtra_field (str | None, optional)**:<br />
        field to add to the atom's `xtra`.

    - **xtra_value (Any, optional)**:<br />
        value of the field to add to the atom's `xtra`.
    """

    symbol = atom.get_name()
    residue = atom.get_parent()

    if not isinstance(residue, Bio.PDB.Residue.Residue):
        raise TypeError(
            f"Expected a `Bio.PDB.Residue.Residue` object, got {type(residue)}"
        )

    residue = residue.get_resname()

    residue_attrs = (
        residue in PROTEIN_ENTITIES,
        residue in ONLY_CA_RESIDUES,
        residue in PURINES,
        residue in PYRIMIDINES,
        symbol
    )

    atom_xtra_dict = {
        (True, False, False, False, "CB"): True,
        (True, True, False, False, "CA"): True,
        (False, False, True, False, "C4"): True,
        (False, False, False, True, "C2"): True,
    }

    atom.xtra[AtomDecoration.IS_REPRESENTATIVE] = atom_xtra_dict.get(
        residue_attrs, False
    )

    if (
        xtra_field is not None
        and xtra_value is not None
        and xtra_field in atom.xtra
    ):
        warnings.warn(dedent(f"""
            The field '{xtra_field}' already exists in the atom's xtra.
            Overwriting the value.""")
        )
        atom.xtra[xtra_field] = xtra_value


class RenumberResidues:
    """Class to renumber the residues based on the offset.

    If the input sequence to the AlphFold does not start from the first residue,
    the residue numbering in the predicted structure will be different from
    the original (UniProt in case of proteins) numbering.

    This class renumbers the residues in the structure based on the provided `offset`.

    > [!NOTE]
    > Although written with the purpose for AlphFold predicted structures,
    > this class can be used to renumber residues in any structure file.
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

        Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython structure object.

        Returns:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython structure object with renumbered residues.
        """

        for model in structure:
            for residue in model.get_residues():
                chain_id = residue.parent.id
                h, num, ins = residue.id

                num = self.renumber_chain_res_num(
                    chain_res_num=num,
                    chain_id=chain_id,
                )

                residue.id = (h, num, ins)

        return structure

    def original_chain_res_num(
        self,
        chain_res_num: int,
        chain_id: str,
    ):
        """Get the original residue number based on the offset.

        Inverse of `renumber_chain_res_num`.

        Arguments:

        - **chain_res_num (int)**:<br />
            Residue index (1-indexed) within the chain in the predicted structure.

        - **chain_id (str)**:<br />
            Chain ID of the residue.

        Returns:

        - **chain_res_num (int)**:<br />
            Original residue number
        """

        if chain_id in self.offset:
            chain_res_num -= (self.offset[chain_id][0] - 1)

        return chain_res_num

    def renumber_chain_res_num(
        self,
        chain_res_num: int,
        chain_id: str,
    ):
        """Renumber the residue number based on the offset.

        Given a residue index (1-indexed) in a chain, this function renumbers it
        based on the offset provided for that chain.

        Arguments:

        - **chain_res_num (int)**:<br />
            Residue index (1-indexed) within the chain in the predicted structure.

        - **chain_id (str)**:<br />
            Chain ID of the residue.

        Returns:

        - **chain_res_num (int)**:<br />
            Renumbered residue number
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

        Arguments:

        - **region_of_interest (dict)**:
            Dictionary containing the region of interest for each chain.

        Returns:

        - **renumbered_region_of_interest (dict)**:
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
        depth: ResidueMapDepth = ResidueMapDepth.ATOM,
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

        Arguments:

        - **token_chain_ids (list)**:<br />
            Token chain IDs. Same as provided in the AF3 JSON output.

        - **token_res_ids (list)**:<br />
            Token residue IDs. Same as provided in the AF3 JSON output.

        - **token_atom_names (list)**:<br />
            Token atom names.

        Returns:

        - **(tuple)**:<br />

            - `idx_to_num (Dict)`:
                Dictionary mapping token indices to token numbers.

            - `num_to_idx (Dict)`:
                Dictionary mapping token numbers to token indices.
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
                ResidueMapKeys.CHAIN_ID: chain_id,
                ResidueMapKeys.TOKEN_NUM: token_num,
            }

            if depth == ResidueMapDepth.ATOM:
                idx_to_num[token_idx][ResidueMapKeys.ATOM_NAME] = atom_name

            if token_num not in num_to_idx[chain_id]:
                if depth == ResidueMapDepth.ATOM:
                    num_to_idx[chain_id][token_num] = {
                        atom_name: token_idx
                    }
                elif depth == ResidueMapDepth.RESIDUE:
                    num_to_idx[chain_id][token_num] = token_idx
            else:
                if depth == ResidueMapDepth.ATOM:
                    num_to_idx[chain_id][token_num][atom_name] = token_idx
                elif depth == ResidueMapDepth.RESIDUE:
                    warnings.warn(
                        f"""

                        Multiple tokens found for residue number {token_num}
                        in chain {chain_id}. Overwriting the previous token index.
                        """
                    )
                    num_to_idx[chain_id][token_num] = token_idx

        return idx_to_num, num_to_idx


class ResidueSelect(Select):
    """Class to select residues in the structure based on the input dictionary.

    This is achieved by overriding the `Bio.PDB.Select.accept_residue` method
    from Biopython.
    """

    confident_residues: Dict
    """ Dictionary containing the chain ID as key and a list of residue
    numbers as value.<br />
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
        ...

    @overload
    def accept_residue(self, residue):
        """Overload this to reject residues for output."""
        return 1

    def accept_residue(
        self,
        residue: Bio.PDB.Residue.Residue
    ) -> bool:
        """Accept the residue if it's in `confident_residues`.

        Arguments:

        - **residue (Bio.PDB.Residue.Residue)**:<br />
            Biopython residue object.

        Returns:

        - **(bool)**:<br />
            `True` if the residue is in the `confident_residues`.
        """

        chain = residue.parent
        if not isinstance(chain, Chain):
            raise TypeError(
                f"Expected a Bio.PDB.Chain.Chain object, got {type(chain)}"
            )
        chain_id = chain.id

        return residue.id[1] in self.confident_residues.get(chain_id, [])