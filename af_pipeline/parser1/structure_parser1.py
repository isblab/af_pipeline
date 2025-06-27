import warnings
import Bio.PDB
import Bio.PDB.PDBParser
import Bio.PDB.Residue
import Bio.PDB.Structure
import os
from af_pipeline.constants.af_constants import *
import Bio
from Bio.PDB import PDBParser, MMCIFParser
from af_pipeline.tools.structure_tools import (
    add_header_footer,
    decorate_residue,
    decorate_atom,
)


class StructureParser:

    def __init__(
        self,
        struct_file_path: str,
        preserve_header_footer: bool = False,
        which_parser: str = "biopython",
    ):

        self.struct_file_path = struct_file_path
        self.preserve_header_footer = preserve_header_footer
        self.which_parser = which_parser.lower()


    def get_parser(self):
        """Get the required parser (PDB/CIF) for the input file.

        Args:

            struct_file_path (str):
                Path to the AF2/3 structure file (.pdb or .cif).

        Returns:

            parser (Bio.PDB.PDBParser | Bio.PDB.MMCIFParser):
                Parser object.
        """

        ext = os.path.splitext(self.struct_file_path)[1]

        if "pdb" in ext:
            parser = PDBParser()

            if self.preserve_header_footer:
                raise Exception("Header can only be preserved for CIF files.")

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

        return parser


    def get_structure(
        self,
        parser: Bio.PDB.PDBParser | Bio.PDB.MMCIFParser,
    ) -> Bio.PDB.Structure.Structure:
        """Return the Biopython Structure object for the structure file.

        Args:

            parser (Bio.PDB.PDBParser | Bio.PDB.MMCIFParser):
                Parser object.

        Returns:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.
        """

        basename = os.path.basename(self.struct_file_path)

        if (
            isinstance(parser, Bio.PDB.PDBParser)
            or isinstance(parser, Bio.PDB.MMCIFParser)
        ):

            structure = parser.get_structure(basename, self.struct_file_path)

            if self.preserve_header_footer:
                structure = add_header_footer(
                    structure=structure,
                    struct_file_path=self.struct_file_path
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

        peratom_quantities = {
            k: [] for k in quantities
        }

        for quantity in quantities:
            if quantity == "coord":
                peratom_quantities["coord"] = atom.coord

            elif quantity == "plddt":
                peratom_quantities["plddt"] = atom.bfactor

            elif quantity == "atom_name":
                peratom_quantities["atom_name"] = atom.get_name()

            elif quantity == "res_pos":
                peratom_quantities["res_pos"] = atom.get_parent().id[1]

            elif quantity == "res_name":
                peratom_quantities["res_name"] = atom.get_parent().get_resname()

            elif quantity == "chain_id":
                peratom_quantities["chain_id"] = atom.get_parent().get_parent().id[0]

            elif quantity == "entity_type":
                peratom_quantities["entity_type"] = atom.get_parent().xtra.get("entityType", None)

            elif quantity == "atom_local_idx":
                for idx, a in enumerate(atom.get_parent()):
                    if a.get_name() == atom.get_name():
                        peratom_quantities["atom_local_idx"] = idx
                        break

        return peratom_quantities

    @staticmethod
    def extract_perresidue_quantities(
        residue: Bio.PDB.Residue.Residue,
        quantities: list = ["coord"],
        rep_atom: str = None,
    ):

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

        perresidue_quantities = {
            k: None for k in quantities
        }

        for quantity in quantities:

            if quantity == "res_pos":
                perresidue_quantities["res_pos"] = residue.id[1]

            elif quantity == "res_name":
                perresidue_quantities["res_name"] = residue.get_resname()

            elif quantity == "chain_id":
                perresidue_quantities["chain_id"] = residue.get_parent().id[0]

            elif quantity == "entity_type":
                perresidue_quantities["entity_type"] = residue.xtra.get(
                    "entityType", None
                )

            elif quantity == "atoms":
                perresidue_quantities["atoms"] = [
                    atom.get_name() for atom in residue
                ]

            elif quantity == "atom_local_idxs":
                perresidue_quantities["atom_local_idxs"] = list(range(len(residue)))

            else:
                if rep_atom is None:
                    rep_atom = StructureParser.get_rep_atom(residue=residue)
                else:
                    try:
                        rep_atom = residue[rep_atom]
                    except KeyError:
                        raise KeyError(
                            f"""
                            Representative atom {rep_atom} not found in \
                            residue {residue.get_resname()}.
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
    def get_rep_atom(
        residue: Bio.PDB.Residue.Residue,
    ):

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

        else:
            warnings.warn(
                f"Unknown entity type for residue {symbol}. "
                "Skipping per-residue quantities extraction."
            )
            rep_atom = residue.child_list[0]  # Fallback to first atom

        return rep_atom

    @staticmethod
    def per_atom_token(
        residue: Bio.PDB.Residue.Residue,
    ):
        condition = (
            residue.xtra.get("is_modified")
            or residue.xtra.get("is_ligand")
            or residue.xtra.get("entityType") is None
        )

        return condition

    @staticmethod
    def get_token_chain_ids(
        structure: Bio.PDB.Structure.Structure,
        rep_atom_dict: dict | None = None,
        only_representative: bool = False,
    ):

        token_chain_ids = []

        for residue, _chain_id in StructureParser.get_residues(structure):

            if (
                StructureParser.per_atom_token(residue)
                and only_representative is False
            ):
                for atom in residue:
                    quants = StructureParser.extract_peratom_quantities(
                            atom=atom,
                            quantities=["chain_id"]
                        )
                    token_chain_ids.append(quants["chain_id"])
            else:
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

    @staticmethod
    def get_token_res_ids(
        structure: Bio.PDB.Structure.Structure,
        rep_atom_dict: dict = {},
        only_representative: bool = False,
    ):

        token_res_ids = []

        for residue, _chain_id in StructureParser.get_residues(structure):
            if (
                StructureParser.per_atom_token(residue)
                and only_representative is False
            ):
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
        rep_atom_dict: dict | None = None,
        only_representative: bool = False,
    ):
        """Get pLDDT values from the structure.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

            per_atom (bool):
                If True, returns pLDDT values for each atom.
                If False, returns pLDDT values for each residue.

        Returns:

            list: pLDDT values.
        """

        plddt_values = []

        if per_atom:
            for atom, _res, _ch_id in StructureParser.get_atoms(structure):
                quants = StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=["plddt"]
                )
                plddt_values.append(quants["plddt"])

        else:
            for residue, _ch_id in StructureParser.get_residues(structure):
                if (
                    StructureParser.per_atom_token(residue)
                    and only_representative is False
                ):
                    for atom in residue:
                        quants = StructureParser.extract_peratom_quantities(
                            atom=atom,
                            quantities=["plddt"]
                        )
                        plddt_values.append(quants["plddt"])
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
        rep_atom_dict: dict | None = None,
        only_representative: bool = False,
    ):
        """Get coordinates from the structure.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

            per_atom (bool):
                If True, returns coordinates for each atom.
                If False, returns coordinates for each residue.

        Returns:

            list: Coordinates.
        """

        coords = []

        if per_atom:
            for atom, _res, _ch_id in StructureParser.get_atoms(structure):
                quants = StructureParser.extract_peratom_quantities(
                    atom=atom,
                    quantities=["coord"]
                )
                coords.append(quants["coord"])

        else:
            for residue, _ch_id in StructureParser.get_residues(structure):
                if (
                    StructureParser.per_atom_token(residue)
                    and only_representative is False
                ):
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
                    quants = StructureParser.extract_perresidue_quantities(
                        residue=residue,
                        quantities=["coord"],
                        rep_atom=rep_atom,
                    )
                    coords.append(quants["coord"])

        return coords
