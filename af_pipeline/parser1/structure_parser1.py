import os
import warnings
import numpy as np
import Bio
import Bio.PDB
import Bio.PDB.Structure
import Bio.PDB.Residue
import Bio.PDB.Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from af_pipeline.constants.af_constants import *
from af_pipeline.tools.structure_tools import (
    add_header_footer,
    decorate_residue,
    decorate_atom,
    has_per_atom_token,
)

#TODO
# Remove which_parser argument/attribute

class StructureParser:
    """ Class to parse structure files (.pdb or .cif) using Biopython.

    Attributes:

        structure_file_path (str):
            Path to the AF2/3 structure file (.pdb or .cif).

        preserve_header_footer (bool):
            If True, the header and footer information is preserved in the
            structure object.
            This is only applicable for .cif files. \n
            (Default: False)

        which_parser (str):
            Which parser to use for the structure file. \n
            "biopython" for Biopython's MMCIFParser, \n
            "pdbe" for PDBe's CifFileReader (not implemented yet). \n
            (Default: "biopython")
    """

    def __init__(
        self,
        structure_file_path: str,
        preserve_header_footer: bool = False,
        which_parser: str = "biopython",
    ):
        """ Initialize the StructureParser.

        Args:

            structure_file_path (str):
                Path to the AF2/3 structure file (.pdb or .cif).

            preserve_header_footer (bool):
                If True, the header and footer information is preserved in the
                structure object.
                This is only applicable for .cif files. \n
                (Default: False)

            which_parser (str):
                Which parser to use for the structure file. \n
                "biopython" for Biopython's MMCIFParser, \n
                "pdbe" for PDBe's CifFileReader (not implemented yet). \n
                (Default: "biopython")
    """

        self.structure_file_path = structure_file_path
        self.preserve_header_footer = preserve_header_footer
        self.which_parser = which_parser


    def get_parser(self):
        """Get the required parser (PDB/CIF) for the input structure file.

        Args:

            structure_file_path (str):
                Path to the structure file (.pdb or .cif).

        Returns:

            parser (Bio.PDB.PDBParser | Bio.PDB.MMCIFParser):
                Parser object.
        """

        parser = None
        ext = os.path.splitext(self.structure_file_path)[1]

        if "pdb" in ext:

            parser = PDBParser()

            if self.preserve_header_footer:
                warnings.warn(
                    """

                    Header can only be preserved for CIF files.
                    Output will not contain header/footer information.
                    """
                )

        elif "cif" in ext:

            if self.which_parser == "biopython":
                parser = MMCIFParser()

            elif self.which_parser == "pdbe":
                raise NotImplementedError(
                    "PDBe parser is not implemented yet. "
                    "Please use Biopython parser."
                )

        else:
            raise Exception("Incorrect file format.. Suported .pdb/.cif only.")

        if not isinstance(parser, (PDBParser, MMCIFParser)):
            raise TypeError(
                f"Parser should be either PDBParser or MMCIFParser. "
                f"Got {type(parser)} instead."
            )

        return parser


    def get_structure(
        self,
        parser: PDBParser | MMCIFParser,
    ):
        """Return the Biopython Structure object for the structure file.

        Args:

            parser (Bio.PDB.PDBParser | Bio.PDB.MMCIFParser):
                Parser object.

        Returns:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.
        """

        basename = os.path.basename(self.structure_file_path)

        if (
            isinstance(parser, PDBParser)
            or isinstance(parser, MMCIFParser)
        ):

            structure = parser.get_structure(basename, self.structure_file_path)

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
                        for atom in residue:
                            decorate_atom(atom=atom)

        else:
            raise Exception(
                "Parser should be either PDBParser or MMCIFParser."
            )

        return structure

    @staticmethod
    def get_residues(structure: Bio.PDB.Structure.Structure):
        """Get residues in the structure.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

        Yields:

            tuple: (residue, chain_id)
                Tuple containing the residue (Bio.PDB.Residue.Residue) and its
                chain ID (str).
        """

        for model in structure:
            for chain in model:
                chain_id = chain.id[0]
                for residue in chain:

                    yield residue, chain_id

    @staticmethod
    def get_atoms(structure: Bio.PDB.Structure.Structure):
        """Get atoms in the structure.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

        Yields:

            tuple: (atom, residue, chain_id)
                Tuple containing the atom (Bio.PDB.Atom.Atom), its residue
                (Bio.PDB.Residue.Residue), and the chain ID (str).
        """

        for model in structure:
            for chain in model:
                chain_id = chain.id[0]
                for residue in chain:
                    for atom in residue:
                        yield atom, residue, chain_id

    @staticmethod
    def extract_peratom_quantities(
        atom: Bio.PDB.Atom.Atom,
        quantities: list = ["coord"]
    ):
        """ Extract per-atom quantities from a Bio.PDB.Atom.Atom object.

        Allowed quantities are:
        - "coord": XYZ atom coordinates (numpy array).
        - "plddt": pLDDT value (float).
        - "atom_name": Atom name (str). e.g. "CA","CB",etc.
        - "res_pos": Residue position (int).
        - "res_name": Residue name (str). e.g. "ALA", "GLY", etc.
        - "chain_id": Chain ID (str).
        - "entity_type": Entity type (str). e.g. "proteinChain"
        - "atom_local_idx": Local index of the atom in the residue (int).

        Args:

            atom (Bio.PDB.Atom.Atom):
                Biopython Atom object.

            quantities (list, optional):
                List of quantities to extract. Defaults to ["coord"].

        Returns:

            peratom_quantities (dict):
                Dictionary with extracted quantities. \n
                Keys are the quantity names, values are the corresponding data.
        """

        allowed_quantities = [
            "coord",
            "plddt",
            "atom_name",
            "res_pos",
            "res_name",
            "chain_id",
            "entity_type",
            "atom_local_idx",
        ]

        assert set(quantities).issubset(
            set(allowed_quantities)
        ), f"Allowed quantities are: {allowed_quantities}"

        peratom_quantities = {}

        residue = atom.get_parent()
        if not isinstance(residue, Residue):
            raise TypeError(
                f"""

                Expected a Bio.PDB.Residue.Residue object.
                Got {type(residue)} instead.
                """
            )

        chain = residue.get_parent()
        if not isinstance(chain, Chain):
            raise TypeError(
                f"""

                Expected a Bio.PDB.Chain.Chain object.
                Got {type(chain)} instead.
                """
            )

        for quantity in quantities:
            if quantity == "coord":
                peratom_quantities["coord"] = atom.coord

            elif quantity == "plddt":
                peratom_quantities["plddt"] = atom.bfactor

            elif quantity == "atom_name":
                peratom_quantities["atom_name"] = atom.get_name()

            elif quantity == "res_pos":
                peratom_quantities["res_pos"] = residue.id[1]

            elif quantity == "res_name":
                peratom_quantities["res_name"] = residue.get_resname()

            elif quantity == "chain_id":
                peratom_quantities["chain_id"] = chain.id[0]

            elif quantity == "entity_type":
                peratom_quantities["entity_type"] = residue.xtra.get("entityType", None)

            elif quantity == "atom_local_idx":
                for idx, a in enumerate(residue):
                    if a.get_name() == atom.get_name():
                        peratom_quantities["atom_local_idx"] = idx
                        break

        return peratom_quantities

    @staticmethod
    def extract_perresidue_quantities(
        residue: Bio.PDB.Residue.Residue,
        quantities: list = ["coord"],
        rep_atom: str | Bio.PDB.Atom.Atom | None = None,
    ):
        """ Extract per-residue quantities from a residue object.

        Following quantities are allowed:
        - "res_pos": Residue position (int).
        - "res_name": Residue name (str). e.g. "ALA", "GLY", etc.
        - "coord": XYZ coordinates of the representative atom (numpy array).
        - "plddt": pLDDT value of the representative atom (float).
        - "chain_id": Chain ID (str).
        - "entity_type": Entity type (str). e.g. "proteinChain"
        - "atoms": Atom names in the residue (list).
        - "atom_local_idxs": Local indices of atoms in the residue. (list)
        - "rep_atom": Name of the representative atom (str).
        - "rep_atom_local_idx": Local index of the representative atom (int)

        Args:

            residue (Bio.PDB.Residue.Residue):
                Biopython Residue object.

            quantities (list, optional):
                List of quantities to extract. Defaults to ["coord"].

            rep_atom (str, optional):
                Representative atom to use for quantities that depend on it.
                If None, the representative atom is determined based on the
                entity type of the residue. Defaults to None.

        Returns:
            perresidue_quantities (dict):
                Dictionary with extracted quantities. \n
                Keys are the quantity names, values are the corresponding data.
        """

        allowed_quantities = [
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

        assert set(quantities).issubset(
            set(allowed_quantities)
        ), f"Allowed quantities are: {allowed_quantities}"

        chain = residue.get_parent()
        if not isinstance(chain, Chain):
            raise TypeError(
                f"""

                Expected a Bio.PDB.Chain.Chain object.
                Got {type(chain)} instead.
                """
            )

        perresidue_quantities = {}

        for quantity in quantities:

            if quantity == "res_pos":
                perresidue_quantities["res_pos"] = residue.id[1]

            elif quantity == "res_name":
                perresidue_quantities["res_name"] = residue.get_resname()

            elif quantity == "chain_id":
                perresidue_quantities["chain_id"] = chain.id[0]

            elif quantity == "entity_type":
                perresidue_quantities["entity_type"] = residue.xtra.get(
                    "entityType", None
                )

            elif quantity == "atoms":
                perresidue_quantities["atoms"] = [
                    atom.get_name() for atom in residue
                ]

            elif quantity == "atom_local_idxs":
                perresidue_quantities["atom_local_idxs"] = list(
                    range(len(residue))
                )

            else:
                if rep_atom is None:
                    rep_atom = StructureParser.get_rep_atom(residue=residue)
                else:
                    if isinstance(rep_atom, str):
                        try:
                            rep_atom = residue[rep_atom]

                        except KeyError:
                            raise KeyError(
                                f"""

                                Representative atom {rep_atom} not found in
                                residue {residue.get_resname()}.
                                """
                            )

                    elif isinstance(rep_atom, Bio.PDB.Atom.Atom):
                        rep_atom = rep_atom

                    else:
                        raise TypeError(
                            f"""

                            rep_atom should be a string.
                            Got {type(rep_atom)} instead.
                            """
                        )

                if quantity == "coord":
                    perresidue_quantities["coord"] = rep_atom.coord

                elif quantity == "plddt":
                    perresidue_quantities["plddt"] = rep_atom.bfactor

                elif quantity == "rep_atom":
                    perresidue_quantities["rep_atom"] = rep_atom.get_name()

                elif quantity == "rep_atom_local_idx":
                    for idx, atom in enumerate(residue):
                        if atom.get_name() == rep_atom.get_name():
                            perresidue_quantities["rep_atom_local_idx"] = idx
                            break

        return perresidue_quantities


    @staticmethod
    def get_rep_atom(residue: Bio.PDB.Residue.Residue):
        """ Get the representative atom for a residue based on its entity type.

        The representative atoms for the most common entity types are:
        - Protein chain: "CA" (if only CA is present) or "CB" (if present).
        - DNA/RNA sequence: "C4" for purines, "C2" for pyrimidines.
        - Ions: Atom name is the residue name (e.g. "NA", "CL").
        - Ligands: First atom in the residue.
        - Unknown entity types: First atom in the residue.

        Args:

            residue (Bio.PDB.Residue.Residue):
                Biopython Residue object.

        Returns:
            rep_atom (Bio.PDB.Atom.Atom):
                Representative atom for the residue. \n
        """

        symbol = residue.get_resname()
        rep_atom = ""

        if residue.xtra.get("entityType") == "proteinChain":

            if residue.xtra.get("is_ca_only"):
                rep_atom = residue["CA"]

            else:
                rep_atom = residue["CB"]

        elif residue.xtra.get("entityType") in ["dnaSequence", "rnaSequence"]:

            if residue.xtra.get("is_purine"):
                rep_atom = residue["C4"]

            elif residue.xtra.get("is_pyrimidine"):
                rep_atom = residue["C2"]

            else:
                rep_atom = residue.child_list[0]  # Fallback to first atom

        elif residue.xtra.get("entityType") == "ion":
            rep_atom = residue[symbol]

        elif residue.xtra.get("entityType") == "ligand":
            rep_atom = residue.child_list[0]  # Fallback to first atom

        elif residue.xtra.get("entityType") is None:
            warnings.warn(
                f"""

                Unknown entity type for residue {symbol}.
                Using first atom or user-specified atom as representative atom.
                """
            )
            rep_atom = residue.child_list[0]  # Fallback to first atom

        else:
            raise NotImplementedError(
                f"""

                Entity type "{residue.xtra.get('entityType')}" not implemented.
                """
            )

        return rep_atom

    @staticmethod
    def get_token_atom_names(
        structure: Bio.PDB.Structure.Structure,
        rep_atom_dict: dict = {},
        only_representative: bool = False,
    ):
        """Get token atom IDs for the structure.

        Token atom IDs is a list of atom names for each token. \n
        If the residue has per-atom tokens, it token atom IDs are the atom
        names for each atom in the residue. \n
        Otherwise, it is the representative atom name for the residue.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

            rep_atom_dict (dict):
                Dictionary with residue names as keys and representative
                atoms as values. \n
                If only_representative is True, this dictionary is used to get
                the representative atom for the specified residue.

            only_representative (bool):
                If True, returns only representative atoms for all residues
                irrespective of per-atom tokens.

        Returns:

            list: Token atom IDs.
        """

        token_atom_ids = []

        for residue, _chain_id in StructureParser.get_residues(structure):

            if has_per_atom_token(residue) and only_representative is False:
                for atom in residue:
                    quants = StructureParser.extract_peratom_quantities(
                        atom=atom,
                        quantities=["atom_name"]
                    )
                    token_atom_ids.append(quants["atom_name"])
            else:
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
        structure: Bio.PDB.Structure.Structure,
        rep_atom_dict: dict = {},
        only_representative: bool = False,
    ):
        """ Get token chain IDs for the structure.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

            rep_atom_dict (dict, optional):
                Dictionary with residue names as keys and representative
                atoms as values. \n
                If only_representative is True, this dictionary is used to get
                the representative atom for the specified residue. \n
                Defaults to {}.

            only_representative (bool, optional):
                If True, returns only representative chain IDs for all residues
                irrespective of per-atom tokens. \n
                Defaults to False.

        Returns:

            token_chain_ids (list):
                Token chain IDs for the structure. \n
        """

        token_chain_ids = []

        for residue, _chain_id in StructureParser.get_residues(structure):

            if has_per_atom_token(residue) and only_representative is False:
                for atom in residue:
                    quants = StructureParser.extract_peratom_quantities(
                            atom=atom,
                            quantities=["chain_id"]
                        )
                    token_chain_ids.append(quants["chain_id"])
            else:
                rep_atom = rep_atom_dict.get(
                    residue.get_resname(),
                    StructureParser.get_rep_atom(residue=residue).id
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
        structure: Bio.PDB.Structure.Structure,
        rep_atom_dict: dict = {},
        only_representative: bool = False,
    ):
        """ Get token residue IDs for the structure.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

            rep_atom_dict (dict, optional):
                Dictionary with residue names as keys and representative
                atoms as values. \n
                If only_representative is True, this dictionary is used to get
                the representative atom for the specified residue. \n
                Defaults to {}.

            only_representative (bool, optional):
                If True, returns only representative residue IDs for all
                residues irrespective of per-atom tokens. \n
                Defaults to False.

        Returns:
            token_res_ids (list):
                Token residue IDs for the structure. \n
        """

        token_res_ids = []

        for residue, _chain_id in StructureParser.get_residues(structure):

            if has_per_atom_token(residue) and only_representative is False:
                for atom in residue:
                    quants = StructureParser.extract_peratom_quantities(
                        atom=atom,
                        quantities=["res_pos"]
                    )
                    token_res_ids.append(quants["res_pos"])

            else:
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
        structure: Bio.PDB.Structure.Structure,
        per_atom: bool = False,
        rep_atom_dict: dict = {},
        average_token_plddt: bool = False,
        only_representative: bool = False,
    ):
        """Get pLDDT values from the structure.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

            per_atom (bool):
                If True, returns pLDDT values for each atom.
                If False, returns pLDDT values for each residue/token.

            rep_atom_dict (dict):
                Dictionary with residue names as keys and representative
                atoms as values. \n
                If only_representative is True, this dictionary is used to get
                the representative atom for the specified residue. \n
                defaults to {}.

            average_token_plddt (bool):
                If True, averages pLDDT values for all atoms in the residue
                with per-atom tokens and returns a single value per residue.

            only_representative (bool):
                If True, returns only representative pLDDT values for all
                residues irrespective of per-atom tokens. \n
                Defaults to False.

        Returns:

            plddt_values (list):
                List of pLDDT values. \n
                If per_atom is True, it contains pLDDT values for each atom.
                If per_atom is False, it contains pLDDT values for each residue
                or token.
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

        else:
            for residue, _ch_id in StructureParser.get_residues(structure):

                if has_per_atom_token(residue):

                    if only_representative is False:
                        for atom in residue:
                            quants = StructureParser.extract_peratom_quantities(
                                atom=atom,
                                quantities=["plddt"]
                            )
                            plddt_values.append(quants["plddt"])

                    elif average_token_plddt is True:
                        # Average pLDDT for all atoms in the residue
                        atom_plddt_values = [
                            StructureParser.extract_peratom_quantities(
                                atom=atom,
                                quantities=["plddt"]
                            )["plddt"] for atom in residue
                        ]
                        plddt_values.append(np.mean(atom_plddt_values))

                    else:
                        rep_atom = rep_atom_dict.get(
                            residue.get_resname(),
                            StructureParser.get_rep_atom(residue=residue)
                        )

                else:
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

        return plddt_values

    @staticmethod
    def get_coordinates(
        structure: Bio.PDB.Structure.Structure,
        per_atom: bool = False,
        rep_atom_dict: dict = {},
        only_representative: bool = False,
    ):
        """Get coordinates from the structure.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

            per_atom (bool):
                If True, returns coordinates for each atom.
                If False, returns coordinates for each residue or token.

            rep_atom_dict (dict):
                Dictionary with residue names as keys and representative
                atoms as values. \n
                If only_representative is True, this dictionary is used to get
                the representative atom for the specified residue. \n
                Defaults to {}.

            only_representative (bool):
                If True, returns only representative coordinates for all
                residues irrespective of per-atom tokens. \n
                Defaults to False.

        Returns:

            coords (list):
                List of coordinates. \n
                If per_atom is True, it contains coordinates for each atom.
                If per_atom is False, it contains coordinates for each residue
                or token.
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

        else:
            for residue, _ch_id in StructureParser.get_residues(structure):

                if has_per_atom_token(residue):

                    if only_representative is False:

                        for atom in residue:
                            quants = StructureParser.extract_peratom_quantities(
                                atom=atom,
                                quantities=["coord"]
                            )
                            coords.append(quants["coord"])
                    else:
                        rep_atom = rep_atom_dict.get(
                            residue.get_resname(),
                            StructureParser.get_rep_atom(residue=residue)
                        )

                else:
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
