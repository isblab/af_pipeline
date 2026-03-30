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
    - "token_level" # token level (per-atom or per-residue)
    - "is_modified" # boolean indicating if the residue is modified
    - "is_ca_only" # boolean indicating if the residue has only CA atom
    - "is_purine" # boolean indicating if the nucleotide is purine
    - "is_pyrimidine" # boolean indicating if the nucleotide is pyrimidine
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
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.PDBParser import PDBParser
from textwrap import dedent
from typing import Any, Generator
from af_pipeline.constants.af_constants import (
    ALLOWED_STRUCTURE_FORMATS,
    AVAILABLE_PARSERS,
    AVAILABLE_ATOM_QUANTITIES,
    AVAILABLE_RESIDUE_QUANTITIES,
    REP_ATOMS,
    VALID_AF3_METRIC_LEVELS,
    VALID_METRIC_LEVELS,
    EntityType,
    AtomQuantity,
    FileFormat,
    MetricLevel,
    ResidueDecoration,
    ResidueQuantity,
    TokenLevel,
    StructureParserConstants as SPCons,
)
from af_pipeline.tools.structure_tools import (
    add_header_footer,
    decorate_residue,
    get_token_level,
)

class StructureParser:
    """ Class to parse structure files (.pdb or .cif) using Biopython."""

    structure_file_path: str
    """ Path to the structure file (PDB or CIF). """

    preserve_header_footer: bool
    """ If `True`, the header and footer information is preserved in the
    structure object.
    > [!NOTE]
    > `preserve_header_footer` is only applicable for .cif files.
    """

    use_fast_cif_parser: bool
    """ If `True`, uses the FastMMCIFParser for parsing .cif files."""

    def __init__(
        self,
        structure_file_path: str,
        preserve_header_footer: bool = SPCons.preserve_header_footer,
        use_fast_cif_parser: bool = SPCons.use_fast_cif_parser,
    ):

        self.structure_file_path = structure_file_path
        self.preserve_header_footer = preserve_header_footer
        self.use_fast_cif_parser = use_fast_cif_parser

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

        if ext == FileFormat.CIF and self.use_fast_cif_parser:
            ext = "fast_cif"

        if ext not in ALLOWED_STRUCTURE_FORMATS:
            raise Exception(
                f"""

                Unsupported file format: {ext}.
                Supported formats are {ALLOWED_STRUCTURE_FORMATS}.
                """
            )

        return ext

    @property
    def parser(self) -> PDBParser | MMCIFParser | FastMMCIFParser:
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
    def sanity_check_metric_level(metric_level: MetricLevel):
        assert metric_level in VALID_METRIC_LEVELS, (
            f"""
            Metric level should be from {VALID_METRIC_LEVELS}.
            Got '{metric_level}' instead."""
        )

    @staticmethod
    def sanity_check_af3_metric_level(metric_level: MetricLevel):
        assert metric_level in VALID_AF3_METRIC_LEVELS, (
            f"""
            Metric level should be from {VALID_AF3_METRIC_LEVELS}.
            Got '{metric_level}' instead."""
        )

    @staticmethod
    def sanity_check_quantities(
        quantities: list[ResidueQuantity | AtomQuantity],
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

    @staticmethod
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
        quantities: list[AtomQuantity] = [AtomQuantity.COORD]
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
            AtomQuantity.COORD: lambda a: a.coord,
            AtomQuantity.PLDDT: lambda a: a.bfactor,
            AtomQuantity.ATOM_NAME: lambda a: a.name,
            AtomQuantity.RES_POS: lambda a: a.parent.id[1],
            AtomQuantity.RES_NAME: lambda a: a.parent.resname,
            AtomQuantity.CHAIN_ID: lambda a: a.parent.parent.id[0],
            AtomQuantity.ENTITY_TYPE: lambda a: a.parent.xtra.get(
                ResidueDecoration.ENTITY_TYPE, None
            ),
            AtomQuantity.ATOM_LOCAL_IDX: lambda a: next(
                idx for idx, at in enumerate(a.parent) if at.name == a.name
            ),
        }

        for quantity in quantities:
            peratom_quantities[quantity] = quantity_funcs[quantity](atom)

        return peratom_quantities

    @staticmethod
    def extract_perresidue_quantities(
        residue: Residue,
        quantities: list[ResidueQuantity] = [ResidueQuantity.COORD],
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
            ResidueQuantity.RES_POS: lambda ra: ra.parent.id[1],
            ResidueQuantity.RES_NAME: lambda ra: ra.parent.resname,
            ResidueQuantity.CHAIN_ID: lambda ra: ra.parent.parent.id[0],
            ResidueQuantity.ENTITY_TYPE: lambda ra: ra.parent.xtra.get(
                ResidueDecoration.ENTITY_TYPE, None
            ),
            ResidueQuantity.ATOMS: lambda ra: [atom.name for atom in ra.parent],
            ResidueQuantity.ATOM_LOCAL_IDXS: lambda ra: list(range(len(ra.parent))),
            ResidueQuantity.COORD: lambda ra: ra.coord,
            ResidueQuantity.PLDDT: lambda ra: ra.bfactor,
            ResidueQuantity.REP_ATOM: lambda ra: ra.name,
            ResidueQuantity.REP_ATOM_LOCAL_IDX: lambda ra: next(
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
            residue.xtra.get(ResidueDecoration.ENTITY_TYPE),
            residue.xtra.get(ResidueDecoration.IS_CA_ONLY, False),
            residue.xtra.get(ResidueDecoration.IS_PURINE, False),
            residue.xtra.get(ResidueDecoration.IS_PYRIMIDINE, False),
        )

        representative_atom_dict = {
            (EntityType.PROTEIN_CHAIN, True, False, False): residue.child_dict.get(REP_ATOMS[ResidueDecoration.IS_CA_ONLY]),
            (EntityType.PROTEIN_CHAIN, False, False, False): residue.child_dict.get(REP_ATOMS[EntityType.PROTEIN_CHAIN]),
            (EntityType.DNA_SEQUENCE, False, True, False): residue.child_dict.get(REP_ATOMS[ResidueDecoration.IS_PURINE]),
            (EntityType.DNA_SEQUENCE, False, False, True): residue.child_dict.get(REP_ATOMS[ResidueDecoration.IS_PYRIMIDINE]),
            (EntityType.RNA_SEQUENCE, False, True, False): residue.child_dict.get(REP_ATOMS[ResidueDecoration.IS_PURINE]),
            (EntityType.RNA_SEQUENCE, False, False, True): residue.child_dict.get(REP_ATOMS[ResidueDecoration.IS_PYRIMIDINE]),
            (EntityType.ION, False, False, False): residue.child_dict.get(symbol),
            (EntityType.DNA_SEQUENCE, False, False, False): residue.child_list[0],
            (EntityType.RNA_SEQUENCE, False, False, False): residue.child_list[0],
            (EntityType.LIGAND, False, False, False): residue.child_list[0],
            (None, False, False, False): residue.child_list[0],
        }

        rep_atom = representative_atom_dict.get(residue_attrs, None)

        if rep_atom is None:
            raise Exception(dedent(f"""
                Representative atom could not be determined for
                residue {residue.resname} with attributes {residue_attrs}
                (entityType, is_ca_only, is_purine, is_pyrimidine).
                """)
            )

        return rep_atom

    @staticmethod
    def get_token_atom_names(
        structure: Structure,
        rep_atom_dict: dict = {},
        metric_level: MetricLevel = MetricLevel.REPRESENTATIVE_TOKEN,
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
            If `metric_level` is `representative_token`, this dictionary is used
            to get the representative atom for the specified residue.

        - **metric_level (MetricLevel, optional)**:<br />
            Metric level for the parser, either "per_token" or "representative_token".

        Returns:

        - **(list)**:<br />
            Token atom IDs.
        """

        StructureParser.sanity_check_af3_metric_level(metric_level)
        token_atom_ids = []

        def _get_rep_atom_name(residue: Bio.PDB.Residue.Residue) -> list[str]:
            rep_atom = rep_atom_dict.get(
                residue.get_resname(),
                StructureParser.get_rep_atom(residue=residue)
            )
            quants = StructureParser.extract_perresidue_quantities(
                residue=residue,
                quantities=[ResidueQuantity.REP_ATOM],
                rep_atom=rep_atom
            )
            return [quants[ResidueQuantity.REP_ATOM]]

        def _get_per_atom_names(residue: Bio.PDB.Residue.Residue) -> list[str]:
            return [
                StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=[AtomQuantity.ATOM_NAME]
                )[AtomQuantity.ATOM_NAME] for atom in residue
            ]

        # (metric_level, token_level)
        handle_dict = {
            (MetricLevel.REPRESENTATIVE_TOKEN, TokenLevel.ATOM): _get_rep_atom_name,
            (MetricLevel.REPRESENTATIVE_TOKEN, TokenLevel.RESIDUE): _get_rep_atom_name,
            (MetricLevel.PER_TOKEN, TokenLevel.ATOM): _get_per_atom_names,
            (MetricLevel.PER_TOKEN, TokenLevel.RESIDUE): _get_rep_atom_name,
        }

        for residue, _chain_id in StructureParser.get_residues(structure):

            token_level = get_token_level(residue)
            worker_func = handle_dict.get((metric_level, token_level))

            if callable(worker_func):
                token_atom_ids.extend(worker_func(residue))

            else:
                raise Exception(dedent(f"""
                    Unexpected combination of flags:
                    {metric_level=} and {token_level=}.""")
                )

        return token_atom_ids

    @staticmethod
    def get_token_chain_ids(
        structure: Structure,
        rep_atom_dict: dict = {},
        metric_level: MetricLevel = MetricLevel.REPRESENTATIVE_TOKEN,
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
            If `metric_level` is `representative_token`, this dictionary is used
            to get the representative atom for the specified residue.

        - **metric_level (MetricLevel, optional)**:<br />
            Metric level for the parser, either "per_token" or "representative_token".

        Returns:

        - **token_chain_ids (list)**:<br />
            Token chain IDs for the structure.
        """

        StructureParser.sanity_check_af3_metric_level(metric_level)
        token_chain_ids = []

        def _get_rep_chain_id(residue: Bio.PDB.Residue.Residue) -> list[str]:
            rep_atom = rep_atom_dict.get(
                residue.get_resname(),
                StructureParser.get_rep_atom(residue=residue)
            )
            quants = StructureParser.extract_perresidue_quantities(
                residue=residue,
                quantities=[ResidueQuantity.CHAIN_ID],
                rep_atom=rep_atom
            )
            return [quants[ResidueQuantity.CHAIN_ID]]

        def _get_per_atom_chain_ids(residue: Bio.PDB.Residue.Residue) -> list[str]:
            return [
                StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=[AtomQuantity.CHAIN_ID]
                )[AtomQuantity.CHAIN_ID] for atom in residue
            ]

        # (metric_level, token_level)
        handle_dict = {
            (MetricLevel.REPRESENTATIVE_TOKEN, TokenLevel.ATOM): _get_rep_chain_id,
            (MetricLevel.REPRESENTATIVE_TOKEN, TokenLevel.RESIDUE): _get_rep_chain_id,
            (MetricLevel.PER_TOKEN, TokenLevel.ATOM): _get_per_atom_chain_ids,
            (MetricLevel.PER_TOKEN, TokenLevel.RESIDUE): _get_rep_chain_id,
        }

        for residue, _chain_id in StructureParser.get_residues(structure):

            token_level = get_token_level(residue)
            worker_func = handle_dict.get((metric_level, token_level))

            if callable(worker_func):
                token_chain_ids.extend(worker_func(residue))

            else:
                raise Exception(dedent(f"""
                    Unexpected combination of flags:
                    {metric_level=} and {token_level=}.""")
                )

        return token_chain_ids

    @staticmethod
    def get_token_res_ids(
        structure: Structure,
        rep_atom_dict: dict = {},
        metric_level: MetricLevel = MetricLevel.REPRESENTATIVE_TOKEN,
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
            If `metric_level` is `representative_token`, this dictionary is used
            to get the representative atom for the specified residue.

        - **metric_level (MetricLevel, optional)**:<br />
            Metric level for the parser, either "per_token" or "representative_token".

        Returns:

        - **token_res_ids (list)**:<br />
            Token residue IDs for the structure.
        """

        StructureParser.sanity_check_af3_metric_level(metric_level)
        token_res_ids = []

        def _get_rep_res_id(residue: Bio.PDB.Residue.Residue) -> list[str]:
            rep_atom = rep_atom_dict.get(
                residue.get_resname(),
                StructureParser.get_rep_atom(residue=residue)
            )
            quants = StructureParser.extract_perresidue_quantities(
                residue=residue,
                quantities=[ResidueQuantity.RES_POS],
                rep_atom=rep_atom
            )
            return [quants[ResidueQuantity.RES_POS]]

        def _get_per_atom_res_ids(residue: Bio.PDB.Residue.Residue) -> list[str]:
            return [
                StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=[AtomQuantity.RES_POS]
                )[AtomQuantity.RES_POS] for atom in residue
            ]

        # (metric_level, token_level)
        handle_dict = {
            (MetricLevel.REPRESENTATIVE_TOKEN, TokenLevel.ATOM): _get_rep_res_id,
            (MetricLevel.REPRESENTATIVE_TOKEN, TokenLevel.RESIDUE): _get_rep_res_id,
            (MetricLevel.PER_TOKEN, TokenLevel.ATOM): _get_per_atom_res_ids,
            (MetricLevel.PER_TOKEN, TokenLevel.RESIDUE): _get_rep_res_id,
        }

        for residue, _chain_id in StructureParser.get_residues(structure):

            token_level = get_token_level(residue)
            worker_func = handle_dict.get((metric_level, token_level))

            if callable(worker_func):
                token_res_ids.extend(worker_func(residue))

            else:
                raise Exception(dedent(f"""
                    Unexpected combination of flags:
                    {metric_level=} and {token_level=}.""")
                )

        return token_res_ids

    @staticmethod
    def get_plddt(
        structure: Structure,
        rep_atom_dict: dict = {},
        average_token_plddt: bool = False,
        metric_level: MetricLevel = MetricLevel.REPRESENTATIVE_TOKEN,
    ) -> list[float]:
        """Get pLDDT values from the structure.

        > [!NOTE]
        > To replicate `atom_plddts` from AF3 JSON file, set `metric_level` to "per_atom".

        Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython Structure object.

        - **rep_atom_dict (dict)**:<br />
            Dictionary with residue names as keys and representative
            atoms as values.<br />
            If `metric_level` is `representative_token`, this dictionary is used
            to get the representative atom for the specified residue.

        - **average_token_plddt (bool)**:<br />
            If `True`, averages pLDDT values for all atoms in the residue
            with per-atom tokens and returns a single value per residue.<br />
            > [!NOTE]
            > `average_token_plddt` option has effect only if `only_representative`
            > is `True`.

        - **metric_level (MetricLevel)**:<br />
            Metric level for the parser.
            "per_token" or "representative_token" or "per_atom".

        Returns:

        - **plddt_values (list)**:<br />
            List of pLDDT values.
        """

        StructureParser.sanity_check_metric_level(metric_level)
        plddt_values = []

        if metric_level == MetricLevel.PER_ATOM:

            plddt_values = [
                StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=[AtomQuantity.PLDDT]
                )[AtomQuantity.PLDDT]
                for atom, _res, _ch_id in StructureParser.get_atoms(structure)
            ]

            return plddt_values

        def _get_rep_plddt_val(residue: Bio.PDB.Residue.Residue) -> list[float]:
            rep_atom = rep_atom_dict.get(
                residue.get_resname(),
                StructureParser.get_rep_atom(residue=residue)
            )
            quants = StructureParser.extract_perresidue_quantities(
                residue=residue,
                quantities=[ResidueQuantity.PLDDT],
                rep_atom=rep_atom,
            )
            return [quants[ResidueQuantity.PLDDT]]

        def _get_avg_plddt_val(residue: Bio.PDB.Residue.Residue) -> list[float]:
            atom_plddt_values = [
                StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=[AtomQuantity.PLDDT]
                )[AtomQuantity.PLDDT] for atom in residue
            ]
            return [np.mean(atom_plddt_values)]

        def _get_per_atom_plddt_vals(residue: Bio.PDB.Residue.Residue) -> list[float]:
            return [
                StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=[AtomQuantity.PLDDT]
                )[AtomQuantity.PLDDT] for atom in residue
            ]

        # (metric_level, average_token_plddt, token_level)
        handle_dict = {
            # single plddt value per residue; token_level flag has no effect in this case

            # from a representative atom
            (MetricLevel.REPRESENTATIVE_TOKEN, False, TokenLevel.ATOM): _get_rep_plddt_val,
            (MetricLevel.REPRESENTATIVE_TOKEN, False, TokenLevel.RESIDUE): _get_rep_plddt_val,

            # average of all atoms in the residue
            (MetricLevel.REPRESENTATIVE_TOKEN, True, TokenLevel.ATOM): _get_avg_plddt_val,
            (MetricLevel.REPRESENTATIVE_TOKEN, True, TokenLevel.RESIDUE): _get_avg_plddt_val,

            # multiple plddt values per residue (one for each atom in the residue)

            # if token is residue-level
            (MetricLevel.PER_TOKEN, True, TokenLevel.ATOM): _get_avg_plddt_val,
            (MetricLevel.PER_TOKEN, False, TokenLevel.RESIDUE): _get_rep_plddt_val,

            # if token is atom-level
            (MetricLevel.PER_TOKEN, False, TokenLevel.ATOM): _get_per_atom_plddt_vals,
            (MetricLevel.PER_TOKEN, True, TokenLevel.ATOM): _get_per_atom_plddt_vals,
        }

        for residue, _ch_id in StructureParser.get_residues(structure):

            token_level = get_token_level(residue)
            handle = (metric_level, average_token_plddt, token_level)
            worker_func = handle_dict.get(handle)

            if callable(worker_func):
                plddt_values.extend(worker_func(residue))

            else:
                raise Exception(dedent(f"""
                    Unexpected combination of flags:
                    {metric_level=},
                    {average_token_plddt=},
                    {token_level=}.""")
                )

        return plddt_values

    @staticmethod
    def get_coordinates(
        structure: Bio.PDB.Structure.Structure,
        rep_atom_dict: dict = {},
        metric_level: MetricLevel = MetricLevel.REPRESENTATIVE_TOKEN,
    ) -> list[np.ndarray]:
        """Get coordinates from the structure.

        Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Biopython Structure object.

        - **rep_atom_dict (dict)**:<br />
            Dictionary with residue names as keys and representative
            atoms as values.<br />
            If `metric_level` is `representative_token`, this dictionary is used
            to get the representative atom for the specified residue.

        - **metric_level (MetricLevel)**:<br />
            Metric level for the parser.
            "per_token" or "representative_token" or "per_atom".

        Returns:

        - **coords (list)**:<br />
            List of coordinates.
        """

        StructureParser.sanity_check_metric_level(metric_level)
        coords = []

        if metric_level == MetricLevel.PER_ATOM:

            coords = [
                StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=[AtomQuantity.COORD]
                )[AtomQuantity.COORD]
                for atom, _res, _ch_id in StructureParser.get_atoms(structure)
            ]

            return coords

        def _get_rep_coords(residue: Bio.PDB.Residue.Residue) -> list[np.ndarray]:
            rep_atom = rep_atom_dict.get(
                residue.get_resname(),
                StructureParser.get_rep_atom(residue=residue)
            )
            quants = StructureParser.extract_perresidue_quantities(
                residue=residue,
                quantities=[ResidueQuantity.COORD],
                rep_atom=rep_atom,
            )
            return [quants[ResidueQuantity.COORD]]

        def _get_per_atom_coords(residue: Bio.PDB.Residue.Residue) -> list[np.ndarray]:
            return [
                StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=[AtomQuantity.COORD]
                )[AtomQuantity.COORD] for atom in residue
            ]

        # (metric_level, token_level)
        handle_dict = {
            (MetricLevel.REPRESENTATIVE_TOKEN, TokenLevel.ATOM): _get_rep_coords,
            (MetricLevel.REPRESENTATIVE_TOKEN, TokenLevel.RESIDUE): _get_rep_coords,
            (MetricLevel.PER_TOKEN, TokenLevel.ATOM): _get_per_atom_coords,
            (MetricLevel.PER_TOKEN, TokenLevel.RESIDUE): _get_rep_coords,
        }

        for residue, _ch_id in StructureParser.get_residues(structure):

            token_level = get_token_level(residue)
            worker_func = handle_dict.get((metric_level, token_level))

            if callable(worker_func):
                coords.extend(worker_func(residue))

            else:
                raise Exception(dedent(f"""
                    Unexpected combination of flags:
                    {metric_level=}, {token_level=}.""")
                )

        return coords