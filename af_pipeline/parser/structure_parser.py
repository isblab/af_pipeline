"""
Structure Parser Module
=======================
- This module provides the StructureParser class to parse structure files[^struct_file]
(.pdb or .cif) using Biopython.


[^struct_file]: *The "structure file" in this context refers to the file that contains\
    the 3D coordinates of the predicted structure.*

<details>
    <summary>
        <b>Note to the maintainer:</b>
    </summary>
- The parser heavily relies on Biopython's PDB and MMCIF parsers.
- `xtra` attribute of Biopython Residue objects is used to store additional
  information such as.
```python
    - "entityType" # entity type
    - "is_modified" # boolean indicating if the residue is modified
    - "is_ca_only" # boolean indicating if the residue has only CA atom
    - "is_purine" # boolean indicating if the nucleotide is purine
    - "is_pyrimidine" # boolean indicating if the nucleotide is pyrimidine
    - "is_ion" # True if the residue is an ion.
    - "is_ligand" # True if the residue is a ligand.
```
- Following methods have equivalent implementations in `af_pipeline.parser.data_parser.DataParser`:
    - `get_token_chain_ids`
    - `get_token_res_ids`<br />
  This is because, the output of these methods are used to define the mapping between
  the tokens indices and token numbers or "residue"[^residue] numbers.
- The other methods that are expected to be used by the end-users are:
    - `get_plddt`
    - `get_coordinates`<br />
- All the above four methods can be used in multiple modes defined by the input
  arguments. However, only a subset of the modes are used in the current pipeline.

[^residue]: *The term "residue" refers to `Bio.PDB.Residue.Residue` object from Biopython. \
    Hence, the "residue" term is used for both amino acids and nucleotides.*

</details>
"""

import os
import warnings
import numpy as np
import Bio
import Bio.PDB
import Bio.PDB.Structure
import Bio.PDB.Residue
import Bio.PDB.Atom
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from typing import Any, Generator
from af_pipeline.constants.af_constants import (
    ALLOWED_STRUCTURE_FORMATS,
    AVAILABLE_PARSERS,
    AVAILABLE_ATOM_QUANTITIES,
    AVAILABLE_RESIDUE_QUANTITIES, REP_ATOMS
)
from af_pipeline.tools.structure_tools import (
    add_header_footer,
    decorate_residue,
    has_per_atom_token,
)

class StructureParser:
    """ Class to parse structure files (.pdb or .cif) using Biopython."""

    structure_file_path: str
    """ Path to the AF2/3 structure file (.pdb or .cif). """

    preserve_header_footer: bool
    """ If `True`, the header and footer information is preserved in the
    structure object.
    > [!NOTE]
    > `preserve_header_footer` is only applicable for .cif files.
    """

    def __init__(
        self,
        structure_file_path: str,
        preserve_header_footer: bool = False,
    ):

        self.structure_file_path = structure_file_path
        self.preserve_header_footer = preserve_header_footer

    @property
    def structure_type(self) -> str:
        """Structure file type based on the file extension.

        See af_pipeline.constants.af_constants.ALLOWED_STRUCTURE_FORMATS for
        supported formats.

        Returns:

        - **ext (str)**:<br />
            Structure file type based on the file extension.
        """

        ext = os.path.splitext(self.structure_file_path)[1].replace(".", "")

        if ext not in ALLOWED_STRUCTURE_FORMATS:
            raise Exception(
                f"""

                Unsupported file format: {ext}.
                Supported formats are {ALLOWED_STRUCTURE_FORMATS}.
                """
            )

        return ext

    @property
    def parser(self) -> PDBParser | MMCIFParser:
        """Parser (PDB/CIF) for the input structure file.

        Returns:

        - **parser (Bio.PDB.PDBParser | Bio.PDB.MMCIFParser)**:<br />
            Appropriate biopython parser for the structure file.
        """

        parser = AVAILABLE_PARSERS[self.structure_type]

        if callable(parser):
            parser = parser()
        else:
            raise TypeError(
                f"""
                Parser is expected to be in {AVAILABLE_PARSERS.values()}.
                Got {type(parser)} instead which is not callable.
                """
            )

        return parser

    def get_structure_obj(self) -> Structure:
        """Return the Biopython Structure object for the structure file.

        Returns:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython Structure object.
        """

        basename = os.path.basename(self.structure_file_path)

        structure = self.parser.get_structure(basename, self.structure_file_path)

        if not isinstance(structure, Structure):
            raise TypeError(
                f"""

                Expected a Bio.PDB.Structure.Structure object.
                Got {type(structure)} instead.
                """
            )

        if self.preserve_header_footer:
            structure = add_header_footer(
                structure=structure,
                structure_file_path=self.structure_file_path
            )

        # decorate residues with entity types and modifications
        for model in structure:
            for chain in model:
                for residue in chain:
                    decorate_residue(residue=residue)

        return structure

    @staticmethod
    def get_residues(
        structure: Structure
    ) -> Generator[tuple[Residue, str], Any, None]:
        """Get residues in the structure.

        Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython Structure object.

        Yields:

        - **(tuple)**:<br />
            Tuple containing:

            - residue (`Bio.PDB.Residue.Residue`)
            - chain ID (`str`).
        """

        for model in structure:
            for chain in model:
                chain_id = chain.id[0]
                for residue in chain:

                    yield residue, chain_id

    @staticmethod
    def get_atoms(
        structure: Structure
    ) -> Generator[tuple[Atom, Residue, str], Any, None]:
        """Get atoms in the structure.

        Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython Structure object.

        Yields:

        - **tuple**:<br />
            Tuple containing -

            - atom (`Bio.PDB.Atom.Atom`)
            - parent residue (`Bio.PDB.Residue.Residue`)
            - chain ID (`str`).
        """

        for model in structure:
            for chain in model:
                chain_id = chain.id[0]
                for residue in chain:
                    for atom in residue:

                        yield atom, residue, chain_id

    @staticmethod
    def sanity_check_quantities(
        quantities: list,
        requested_on: Atom | Residue
    ):

        if isinstance(requested_on, Atom):
            available_quantities = AVAILABLE_ATOM_QUANTITIES
        elif isinstance(requested_on, Residue):
            available_quantities = AVAILABLE_RESIDUE_QUANTITIES
        else:
            raise TypeError(
                f"""

                Expected requested_on to be Bio.PDB.Atom.Atom or
                Bio.PDB.Residue.Residue object.
                Got {type(requested_on)} instead.
                """
            )

        if set(quantities).issubset(set(available_quantities)) is False:
            raise ValueError(
                f"""

                Could not recognize the following quantities:
                {set(quantities).difference(set(available_quantities))}

                Allowed quantities are: {available_quantities}
                """
            )

    def sanity_check_not_orphan(requested_on: Atom | Residue):

        if isinstance(requested_on, Atom):

            assert (
                isinstance(requested_on.parent, Residue)
                and isinstance(requested_on.parent.parent, Chain)
            ), (
                f"""
                Expected a Bio.PDB.Residue.Residue as parent object and
                Bio.PDB.Chain.Chain as grand-parent object. Got
                {type(requested_on.parent)} & {type(requested_on.parent.parent)}
                instead.
                """
            )

        elif isinstance(requested_on, Residue):
            assert isinstance(requested_on.parent, Chain), (
                f"""

                Expected a Bio.PDB.Chain.Chain as parent object.
                Got {type(requested_on.parent)} instead.
                """
            )

    @staticmethod
    def extract_peratom_quantities(
        atom: Atom,
        quantities: list = ["coord"]
    ) -> dict[str, Any]:
        """ Extract per-atom quantities from a `Bio.PDB.Atom.Atom` object.

        Allowed quantities are:\n
        ```python
        - "coord" # XYZ atom coordinates (np.ndarray).
        - "plddt" # pLDDT value (float).
        - "atom_name" # Atom name (str). e.g. "CA","CB",etc.
        - "res_pos" # Residue position (int).
        - "res_name" # Residue name (str). e.g. "ALA", "GLY", etc.
        - "chain_id" # Chain ID (str).
        - "entity_type" # Entity type (str). e.g. "proteinChain"
        - "atom_local_idx" # Local index of the atom in the residue (int).
        ```

        Arguments:

        - **atom (Bio.PDB.Atom.Atom)**:<br />
            Biopython Atom object.

        - **quantities (list, optional)**:<br />
            List of quantities to extract.

        Returns:

        - **peratom_quantities (dict)**:<br />
            Dictionary with the following `key`:`value` pair.<br />
            - `quantity_name`:`quantity_value`
        """

        StructureParser.sanity_check_quantities(quantities, atom)
        StructureParser.sanity_check_not_orphan(atom)

        peratom_quantities = {}


        quantity_funcs = {
            "coord": lambda a: a.coord,
            "plddt": lambda a: a.bfactor,
            "atom_name": lambda a: a.name,
            "res_pos": lambda a: a.parent.id[1],
            "res_name": lambda a: a.parent.resname,
            "chain_id": lambda a: a.parent.parent.id[0],
            "entity_type": lambda a: a.parent.xtra.get("entityType", None),
            "atom_local_idx": lambda a: next(
                idx for idx, at in enumerate(a.parent) if at.name == a.name
            ),
        }

        for quantity in quantities:
            peratom_quantities[quantity] = quantity_funcs[quantity](atom)

        return peratom_quantities

    @staticmethod
    def extract_perresidue_quantities(
        residue: Residue,
        quantities: list = ["coord"],
        rep_atom: str | Bio.PDB.Atom.Atom | None = None,
    ) -> dict[str, Any]:
        """ Extract per-residue quantities from a residue object.

        Following quantities are allowed: \n
        ```python
        - "res_pos" # Residue position (int).
        - "res_name" # Residue name (str). e.g. "ALA", "GLY", etc.
        - "coord" # XYZ coordinates of the representative atom (numpy array).
        - "plddt" # pLDDT value of the representative atom (float).
        - "chain_id" # Chain ID (str).
        - "entity_type" # Entity type (str). e.g. "proteinChain"
        - "atoms" # Atom names in the residue (list).
        - "atom_local_idxs" # Local indices of atoms in the residue. (list)
        - "rep_atom" # Name of the representative atom (str).
        - "rep_atom_local_idx" # Local index of the representative atom (int)
        ```

        By default, the representative atom is determined by `get_rep_atom`
        method based on the entity type of the residue.<br />
        In case a specific representative atom is to be used, it can be passed
        using the `rep_atom` argument as a string (atom name e.g. "N") or as a
        `Bio.PDB.Atom.Atom` object.

        Arguments:

        - **residue (Bio.PDB.Residue.Residue)**:<br />
            Biopython Residue object.

        - **quantities (list, optional)**:<br />
            List of quantities to extract.

        - **rep_atom (str, optional)**:<br />
            Representative atom to use for quantities that depend on it.
            If `None`, the representative atom is determined based on the
            entity type of the residue.

        Returns:

        - **perresidue_quantities (dict)**:<br />
            Dictionary with the following `key`:`value` pair.
            - `quantity_name`:`quantity_value`
        """

        StructureParser.sanity_check_quantities(quantities, residue)
        StructureParser.sanity_check_not_orphan(residue)

        if rep_atom is None:
            rep_atom = StructureParser.get_rep_atom(residue=residue)

        elif isinstance(rep_atom, str):
            try:
                rep_atom = residue[rep_atom]

            except KeyError:
                raise KeyError(
                    f"""

                    Representative atom {rep_atom} not found in
                    residue {residue.resname} of chain {residue.parent.id[0]}.
                    """
                )

        elif isinstance(rep_atom, Bio.PDB.Atom.Atom):
            assert rep_atom.parent == residue, (
                f"""

                Representative atom {rep_atom.name} does not belong to
                residue {residue.resname} of chain {residue.parent.id[0]}.
                """
            )

        else:
            raise TypeError(
                f"""

                `rep_atom` should be a string or biopython atom object.
                Got {type(rep_atom)} instead.
                """
            )

        perresidue_quantities = {}

        quantity_funcs = {
            "res_pos": lambda ra: ra.parent.id[1],
            "res_name": lambda ra: ra.parent.resname,
            "chain_id": lambda ra: ra.parent.parent.id[0],
            "entity_type": lambda ra: ra.parent.xtra.get("entityType", None),
            "atoms": lambda ra: [atom.name for atom in ra.parent],
            "atom_local_idxs": lambda ra: list(range(len(ra.parent))),
            "coord": lambda ra: ra.coord,
            "plddt": lambda ra: ra.bfactor,
            "rep_atom": lambda ra: ra.name,
            "rep_atom_local_idx": lambda ra: next(
                idx for idx, atom in enumerate(ra.parent) if atom.name == ra.name
            ),
        }

        for quantity in quantities:

            perresidue_quantities[quantity] = quantity_funcs[quantity](rep_atom)

        return perresidue_quantities

    @staticmethod
    def get_rep_atom(
        residue: Residue
    ) -> Bio.PDB.Atom.Atom:
        """ Get the representative atom for a residue based on its entity type.

        The representative atoms for the most common entity types are:
        ```python
        - proteinChain # "CA" (if only CA is present) or "CB" (if present).
        - dnaSequence # "C4" for purines, "C2" for pyrimidines.
        - rnaSequence # "C4" for purines, "C2" for pyrimidines.
        - ion # Atom name is the residue name (e.g. "NA", "CL").
        - ligan # First atom in the residue.
        - unknown # First atom in the residue.
        ```

        > [!NOTE]
        > Choosing representative atom depends on the decoration of the residue
        > using `af_pipeline.tools.structure_tools.decorate_residue` method.
        > Make sure to decorate the residues before using this method.

        Arguments:

        - **residue (Bio.PDB.Residue.Residue)**:<br />
            Biopython Residue object.

        Returns:

        - **rep_atom (Bio.PDB.Atom.Atom)**:<br />
            Representative atom for the residue.
        """

        symbol = residue.get_resname()
        rep_atom = ""

        residue_attrs = (
            residue.xtra.get("entityType"),
            residue.xtra.get("is_ca_only", False),
            residue.xtra.get("is_purine", False),
            residue.xtra.get("is_pyrimidine", False),
        )

        representative_atom_dict = {
            ("proteinChain", True, False, False): residue[REP_ATOMS["is_ca_only"]],
            ("proteinChain", False, False, False): residue[REP_ATOMS["proteinChain"]],
            ("dnaSequence", False, True, False): residue[REP_ATOMS["is_purine"]],
            ("dnaSequence", False, False, True): residue[REP_ATOMS["is_pyrimidine"]],
            ("rnaSequence", False, True, False): residue[REP_ATOMS["is_purine"]],
            ("rnaSequence", False, False, True): residue[REP_ATOMS["is_pyrimidine"]],
            ("ion", False, False, False): residue[symbol],
            ("dnaSequence", False, False, False): residue.child_list[0],
            ("rnaSequence", False, False, False): residue.child_list[0],
            ("ligand", False, False, False): residue.child_list[0],
            (None, False, False, False): residue.child_list[0],
        }

        rep_atom = representative_atom_dict.get(residue_attrs, None)

        if rep_atom is None:
            raise Exception(
                f"""

                Representative atom could not be determined for
                residue {residue.resname} with attributes {residue_attrs}
                (entitType, is_ca_only, is_purin, is_pyrimidine).
                """
            )

        return rep_atom

    @staticmethod
    def get_token_atom_names(
        structure: Structure,
        rep_atom_dict: dict = {},
        only_representative: bool = False,
    ) -> list[str]:
        """Get token atom IDs for the structure.

        Token atom IDs is a list of atom names for each token. \n
        If the residue has per-atom tokens, it token atom IDs are the atom
        names for each atom in the residue. \n
        Otherwise, it is the representative atom name for the residue.

        > [!NOTE]
        > The default parameter settings replicate `token_atom_names` from
        > AlphaFold3 output JSON.

        Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython Structure object.

        - **rep_atom_dict (dict, optional)**:<br />
            Dictionary with residue names as keys and representative
            atoms as values.<br />
            If `only_representative` is `True`, this dictionary is used to get
            the representative atom for the specified residue.

        - **only_representative (bool, optional)**:<br />
            If `True`, returns only representative atoms for all residues
            irrespective of per-atom tokens.

        Returns:

        - **(list)**:<br />
            Token atom IDs.
        """

        token_atom_ids = []

        # (only_representative, has_per_atom_token)
        handle_dict = {
            "rep_atom": [
                (True, False),
                (True, True),
                (False, False),
            ],
            "per_atom": [
                (False, True),
            ],
        }

        for residue, _chain_id in StructureParser.get_residues(structure):

            handle = (only_representative, has_per_atom_token(residue))

            if handle in handle_dict["per_atom"]:
                for atom in residue:
                    quants = StructureParser.extract_peratom_quantities(
                        atom=atom,
                        quantities=["atom_name"]
                    )
                    token_atom_ids.append(quants["atom_name"])

            elif handle in handle_dict["rep_atom"]:
                rep_atom = rep_atom_dict.get(
                    residue.get_resname(),
                    StructureParser.get_rep_atom(residue=residue)
                )
                quants = StructureParser.extract_perresidue_quantities(
                    residue=residue,
                    quantities=["rep_atom"],
                    rep_atom=rep_atom
                )
                token_atom_ids.append(quants["rep_atom"])

        return token_atom_ids

    @staticmethod
    def get_token_chain_ids(
        structure: Structure,
        rep_atom_dict: dict = {},
        only_representative: bool = False,
    ) -> list[str]:
        """ Get token chain IDs for the structure.

        > [!NOTE]
        > The default parameter settings replicate `token_chain_ids` from
        > AlphaFold3 output JSON.

        Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython Structure object.

        - **rep_atom_dict (dict, optional)**:<br />
            Dictionary with residue names as keys and representative
            atoms as values.<br />
            If `only_representative` is `True`, this dictionary is used to get
            the representative atom for the specified residue.

        - **only_representative (bool, optional)**:<br />
            If `True`, returns only representative chain IDs for all residues
            irrespective of per-atom tokens.

        Returns:

        - **token_chain_ids (list)**:<br />
            Token chain IDs for the structure.
        """

        token_chain_ids = []

        # (only_representative, has_per_atom_token)
        handle_dict = {
            "rep_atom": [
                (True, False),
                (True, True),
                (False, False),
            ],
            "per_atom": [
                (False, True),
            ],
        }

        for residue, _chain_id in StructureParser.get_residues(structure):

            handle = (only_representative, has_per_atom_token(residue))

            if handle in handle_dict["per_atom"]:
                for atom in residue:
                    quants = StructureParser.extract_peratom_quantities(
                        atom=atom,
                        quantities=["chain_id"]
                    )
                    token_chain_ids.append(quants["chain_id"])

            elif handle in handle_dict["rep_atom"]:
                rep_atom = rep_atom_dict.get(
                    residue.get_resname(),
                    StructureParser.get_rep_atom(residue=residue)
                )
                quants = StructureParser.extract_perresidue_quantities(
                    residue=residue,
                    quantities=["chain_id"],
                    rep_atom=rep_atom
                )
                token_chain_ids.append(quants["chain_id"])

        return token_chain_ids

    @staticmethod
    def get_token_res_ids(
        structure: Structure,
        rep_atom_dict: dict = {},
        only_representative: bool = False,
    ) -> list[str]:
        """ Get token residue IDs for the structure.

        > [!NOTE]
        > The default parameter settings replicate `token_res_ids` from
        > AlphaFold3 output JSON.

        Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython Structure object.

        - **rep_atom_dict (dict, optional)**:<br />
            Dictionary with residue names as keys and representative
            atoms as values.<br />
            If `only_representative` is `True`, this dictionary is used to get
            the representative atom for the specified residue.

        - **only_representative (bool, optional)**:<br />
            If `True`, returns only representative residue IDs for all
            residues irrespective of per-atom tokens.

        Returns:

        - **token_res_ids (list)**:<br />
            Token residue IDs for the structure.
        """

        token_res_ids = []

        # (only_representative, has_per_atom_token)
        handle_dict = {
            "rep_atom": [
                (True, False),
                (True, True),
                (False, False),
            ],
            "per_atom": [
                (False, True),
            ],
        }

        for residue, _chain_id in StructureParser.get_residues(structure):

            handle = (only_representative, has_per_atom_token(residue))

            if handle in handle_dict["per_atom"]:
                for atom in residue:
                    quants = StructureParser.extract_peratom_quantities(
                        atom=atom,
                        quantities=["res_pos"]
                    )
                    token_res_ids.append(quants["res_pos"])

            elif handle in handle_dict["rep_atom"]:
                rep_atom = rep_atom_dict.get(
                    residue.get_resname(),
                    StructureParser.get_rep_atom(residue=residue)
                )
                quants = StructureParser.extract_perresidue_quantities(
                    residue=residue,
                    quantities=["res_pos"],
                    rep_atom=rep_atom
                )
                token_res_ids.append(quants["res_pos"])

        return token_res_ids

    @staticmethod
    def get_plddt(
        structure: Structure,
        per_atom: bool = False,
        rep_atom_dict: dict = {},
        average_token_plddt: bool = False,
        only_representative: bool = False,
    ) -> list[float]:
        """Get pLDDT values from the structure.

        > [!NOTE]
        > To replicate `atom_plddts` from AF3 JSON file, set `per_atom` to `True`.

        Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython Structure object.

        - **per_atom (bool)**:<br />
            If `True`, returns pLDDT values for each atom.<br />
            If `False`, returns pLDDT values for each residue/token.<br />
            > [!NOTE]
            > `per_atom` option **supersedes** all other options.

        - **rep_atom_dict (dict)**:<br />
            Dictionary with residue names as keys and representative
            atoms as values.<br />
            If `only_representative` is `True`, this dictionary is used to get
            the representative atom for the specified residue.

        - **average_token_plddt (bool)**:<br />
            If `True`, averages pLDDT values for all atoms in the residue
            with per-atom tokens and returns a single value per residue.<br />
            > [!NOTE]
            > `average_token_plddt` option has effect only if `only_representative`
            > is `True`.

        - **only_representative (bool)**:<br />
            If `True`, returns only representative pLDDT values for all
            residues.

        Returns:

        - **plddt_values (list)**:<br />
            List of pLDDT values.<br />
            If `per_atom` is `True`, contains pLDDT values for each atom.<br />
            If `per_atom` is `False`, contains pLDDT values for each residue or token.
        """

        plddt_values = []
        rep_atom = None

        if per_atom:
            for atom, _res, _ch_id in StructureParser.get_atoms(structure):
                quants = StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=["plddt"]
                )
                plddt_values.append(quants["plddt"])

            return plddt_values

        # (only_representative, average_token_plddt, has_per_atom_token)
        handle_dict = {
            "rep_atom_plddt": [
                (True, False, True),
                (True, False, False),
                (False, True, False),
                (False, False, False),
            ],
            "avg_atom_plddt": [
                (True, True, True),
                (True, True, False),
            ],
            "all_atom_plddt": [
                (False, True, True),
                (False, False, True),
            ]
        }

        for residue, _ch_id in StructureParser.get_residues(structure):

            handle = (
                only_representative,
                average_token_plddt,
                has_per_atom_token(residue)
            )

            if handle in handle_dict["rep_atom_plddt"]:

                rep_atom = rep_atom_dict.get(
                    residue.get_resname(),
                    StructureParser.get_rep_atom(residue=residue)
                )
                quants = StructureParser.extract_perresidue_quantities(
                    residue=residue,
                    quantities=["plddt"],
                    rep_atom=rep_atom,
                )
                plddt_values.append(quants["plddt"])

            elif handle in handle_dict["avg_atom_plddt"]:
                # Average pLDDT for all atoms in the residue
                atom_plddt_values = [
                    StructureParser.extract_peratom_quantities(
                        atom=atom,
                        quantities=["plddt"]
                    )["plddt"] for atom in residue
                ]
                plddt_values.append(np.mean(atom_plddt_values))

            elif handle in handle_dict["all_atom_plddt"]:

                for atom in residue:
                    quants = StructureParser.extract_peratom_quantities(
                        atom=atom,
                        quantities=["plddt"]
                    )
                    plddt_values.append(quants["plddt"])

            else:
                raise Exception(
                    f"""

                    Unexpected combination of flags:
                    only_representative={only_representative},
                    average_token_plddt={average_token_plddt},
                    has_per_atom_token={has_per_atom_token(residue)}.

                    Expected booleans for all three flags.
                    """
                )

        return plddt_values

    @staticmethod
    def get_coordinates(
        structure: Bio.PDB.Structure.Structure,
        per_atom: bool = False,
        rep_atom_dict: dict = {},
        only_representative: bool = False,
    ) -> list[np.ndarray]:
        """Get coordinates from the structure.

        Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython Structure object.

        - **per_atom (bool)**:<br />
            If `True`, returns coordinates for each atom.
            If `False`, returns coordinates for each residue or token.
            > [!NOTE]
            > `per_atom` option **supersedes** all other options.

        - **rep_atom_dict (dict)**:<br />
            Dictionary with residue names as keys and representative
            atoms as values.<br />
            If `only_representative` is `True`, this dictionary is used to get
            the representative atom for the specified residue.

        - **only_representative (bool)**:<br />
            If `True`, returns only representative coordinates for all
            residues.

        Returns:

        - **coords (list)**:<br />
            List of coordinates. <br />
            If `per_atom` is `True`, contains coordinates for each atom.<br />
            If `per_atom` is `False`, contains coordinates for each residue or token.
        """

        coords = []
        rep_atom = None

        if per_atom:

            for atom, _res, _ch_id in StructureParser.get_atoms(structure):
                quants = StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=["coord"]
                )
                coords.append(quants["coord"])

            return coords

        # (only_representative, has_per_atom_token)
        handle_dict = {
            "rep_atom_coord": [
                (True, True),
                (True, False),
                (False, False),
            ],
            "all_atom_coord": [
                (False, True),
            ]
        }

        for residue, _ch_id in StructureParser.get_residues(structure):

            handle = (only_representative, has_per_atom_token(residue))

            if handle in handle_dict["all_atom_coord"]:

                for atom in residue:
                    quants = StructureParser.extract_peratom_quantities(
                        atom=atom,
                        quantities=["coord"]
                    )
                    coords.append(quants["coord"])

            elif handle in handle_dict["rep_atom_coord"]:

                rep_atom = rep_atom_dict.get(
                    residue.get_resname(),
                    StructureParser.get_rep_atom(residue=residue)
                )
                quants = StructureParser.extract_perresidue_quantities(
                    residue=residue,
                    quantities=["coord"],
                    rep_atom=rep_atom,
                )
                coords.append(quants["coord"])

        return coords