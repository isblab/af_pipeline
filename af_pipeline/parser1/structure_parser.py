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
)

class StructureParser:
    """Class to parse the AF2/3 predicted structures

    Attributes:

        struct_file_path (str):
            Path to the AF2/3 structure file (.pdb or .cif).

        preserve_header_footer (bool):
            If True, the header and footer information is preserved in the
            structure object.
            This is only applicable for CIF files. \n
            (Default: False)

        which_parser (str):
            Which parser to use for the CIF file. \n
            "biopython" for Biopython's MMCIFParser or PDBParser, \n
            "pdbe" for PDBe's CifFileReader (not implemented yet). \n
            (Default: "biopython")

        only_representative (bool):
            If True, only the representative atoms are considered for
            modified protein/dna/rna chains and ligands. \n
            (Default: False)
    """

    def __init__(
        self,
        struct_file_path: str,
        **kwargs,
    ):
        """ Initialize the StructureParser class.

        Args:

            struct_file_path (str):
                Path to the AF2/3 structure file (.pdb or .cif).

            **preserve_header_footer (bool, optional):
                If True, the header and footer information is preserved in
                the structure object.
                This is only applicable for CIF files. \n
                (Default: False)

            **which_parser (str, optional):
                Which parser to use for the CIF file. \n
                "biopython" for Biopython's MMCIFParser or PDBParser, \n
                "pdbe" for PDBe's CifFileReader (not implemented yet). \n
                (Default: "biopython")

            **only_representative (bool, optional):
                If True, only the representative atoms are considered for
                modified protein/dna/rna chains and ligands. \n
                (Default: False)
        """

        self.struct_file_path = struct_file_path
        self.preserve_header_footer = kwargs.get(
            "preserve_header_footer", False
        )
        self.which_parser = kwargs.get("which_parser", "biopython")
        self.only_representative = kwargs.get("only_representative", False)

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

    def extract_perresidue_quantity(
        self,
        residue: Bio.PDB.Residue.Residue,
        quantity: str,
    ):
        """Extract per-residue quantities from the Biopython residue object.

        This is exploiting the `xtra` tag in the Biopython residue object
        which is decorated with entity types and modifications by the
        `decorate_residue` function. \n

        Given the Biopython residue object, return the specified quantity: \n
            1. residue or nucleotide or ion position \n
            2. Co-ordinates of the representative atom or all atoms \n
            3. pLDDT value of the representative atom or all atoms. \n

        Depending on the `only_representative` flag, the function will
        either return the requested quantity for the representative
        atom of the residue or all atoms in the residue. \n

        The representative atom is determined based on the entity type
        of the residue and its symbol. \n
        - for `proteinChain`, the representative atom is "CB" or "CA" in case
        "CB" is not present.
        - for `dnaSequence` and `rnaSequence`, the representative atom is "C4"
        for purines and "C2" for pyrimidines.
        - for `ion`, the representative atom is the symbol of the ion.
        - for `ligand`, the first atom in the ligand is selected as the
        representative atom.

        For modified protein/dna/rna chains and ligands, per atom quantities
        can be extracted if `only_representative` is set to False. \n

        Args:

            residue (Bio.PDB.Residue.Residue):
                Biopython residue object.

            quantity (str):
                Quantity to extract. ("res_pos", "coords", "plddt")

            only_representative (bool):
                If True, only the representative atoms are considered for
                modified protein/dna/rna chains and ligands. \n
                (Default: False)

        Returns:

            extracted_quantity (list):
                The extracted quantity based on the specified `quantity`.
                - If `quantity=="res_pos"`, returns residue position.
                - If `quantity=="coords"`, returns coordinates of the
                    representative atom or all atoms.
                - If `quantity=="plddt"`, returns pLDDT value of the
                    representative atom or all atoms.
        """

        # Using representative atoms as specified by AF3.
        # https://github.com/google-deepmind/alphafold3/blob/main/src/alphafold3/model/features.py

        symbol = residue.get_resname()
        rep_atom = residue.child_list[0].get_name()

        # 1) Determine the representative atom based on the entity type
        if residue.xtra.get("entityType") == "proteinChain":

            # if the residue is standard residue
            # or user wants only representative atoms
            if self.only_representative or symbol in STD_RESIDUES:

                if (
                    "CB" in residue.child_dict
                    and symbol in PROTEIN_ENTITIES
                ):  # this includes modifications
                    rep_atom = "CB"

                elif (
                    "CB" not in residue.child_dict
                    and symbol in ONLY_CA_RESIDUES
                ):  # this includes modifications
                    rep_atom = "CA"

                else:
                    raise Exception(
                        f"""
                        Are you sure this is a protein chain?
                        residue {symbol} in chain {residue.parent.id}
                        does not have a Cb-atom or a Ca-atom.
                        """
                    )

            # if the residue is a post-translationally modified
            # return all atoms in the PTM residue
            elif symbol in ALLOWED_PTMS:
                rep_atom = residue.child_list

            else:
                raise Exception(
                    f"Unknown protein entity type for residue {symbol} "
                    f"in chain {residue.parent.id}."
                )

        elif residue.xtra.get("entityType") in ["dnaSequence", "rnaSequence"]:

            # if the dna/rna nucleotide is standard
            # or user wants only representative atoms
            if (
                self.only_representative
                or symbol in PURINES_STD + PYRIMIDINES_STD
            ):

                if symbol in PURINES:  # this includes modifications
                    rep_atom = "C4"

                elif symbol in PYRIMIDINES:  # this includes modifications
                    rep_atom = "C2"

            # if the dna/rna nucleotide is a modified nucleotide
            # rep atoms are all atoms in the modified nucleotide
            elif symbol in ALLOWED_DNA_MODS + ALLOWED_RNA_MODS:
                rep_atom = residue.child_list

            else:
                raise Exception(
                    f"Unknown DNA/RNA entity type for residue {symbol} "
                    f"in chain {residue.parent.id}."
                )

        elif residue.xtra.get("entityType") == "ion" and symbol in ION:

            rep_atom = symbol

        elif residue.xtra.get("entityType") == "ligand":

            # if the user wants only representative atoms
            # the first atom in the ligand is selected
            if self.only_representative:
                rep_atom = residue.child_list[0].get_name()
                warnings.warn(
                    f"""
                    Can not determine representative atom for ligand {symbol}
                    in chain {residue.parent.get_id()}
                    Setting representative atom to {rep_atom}.
                    """
                )

            # rep atoms are all atoms in the ligand
            else:
                rep_atom = residue.child_list

        # if the entity type is unknown
        # this could be a glycan modification
        else:

            # if the user wants only representative atoms
            # select the first atom as representative atom
            if self.only_representative:
                rep_atom = residue.child_list[0].get_name()
                warnings.warn(
                    f"""
                    Unknown entity type for residue {symbol}
                    in chain {residue.parent.id}.
                    It could be a glycan modification.
                    Setting representative atom to {rep_atom}.
                    """
                )

            # return all atoms in the residue
            else:
                rep_atom = residue.child_list

        if isinstance(rep_atom, list):
            rep_atom_list = [atom.id for atom in rep_atom]
        else:
            rep_atom_list = [rep_atom]

        # 2) return the requested quantity
        if quantity == "res_pos":
            if isinstance(rep_atom, list):
                # if rep_atom is a list, return residue positions for all
                res_pos_list = [residue.id[1]] * len(rep_atom_list)
                return res_pos_list, rep_atom_list
            else:
                return [residue.id[1]], rep_atom_list

        elif quantity == "coords":
            if isinstance(rep_atom, list):
                # if rep_atom is a list, return atom coordinates for all
                coords = [atom.coord for atom in rep_atom]
            else:
                # return coordinates of the representative atom
                coords = [residue[rep_atom].coord]
            return coords, rep_atom_list

        elif quantity == "plddt":
            if isinstance(rep_atom, list):
                # if rep_atom is a list, return atom pLDDT for all
                plddt = [atom.bfactor for atom in rep_atom]
            else:
                # return pLDDT of the representative atom
                plddt = [residue[rep_atom].bfactor]

            return plddt, rep_atom_list

        else:
            raise Exception(
                f"Specified quantity: {quantity} does not exist for {symbol}"
                f" in chain {residue.parent.id}.\n"
                f"Available quantities are: "
                "'res_pos', 'coords', 'plddt'."
            )

    def get_token_chain_res_ids(
        self,
        structure: Bio.PDB.Structure.Structure,
    ):
        """Get the token chain IDs and token residue IDs for all residues.

        This function extracts the chain IDs and residue IDs from the
        structure. \n
        The output is same as `token_chain_ids` and `token_res_ids` from the
        `.json` file provided by AF3 for the corresponding structure.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

            only_representative (bool):
                If True, only the representative atoms are considered for
                modified protein/dna/rna chains and ligands. \n
                (Default: False)

        Returns:

            token_chain_ids (list):
                List of chain IDs for each residue in the structure.

            token_res_ids (list):
                List of residue IDs for each residue in the structure.
                The residue ID is the position of the residue in the chain.
        """

        token_chain_ids = []
        token_res_ids = []

        for residue, chain_id in self.get_residues(structure):

            res_ids, _ = self.extract_perresidue_quantity(
                residue=residue,
                quantity="res_pos",
            )

            token_chain_ids.extend([chain_id]* len(res_ids))
            token_res_ids.extend(res_ids)

        return token_chain_ids, token_res_ids

    def get_rep_coordinates(
        self,
        structure: Bio.PDB.Structure.Structure,
    ):
        """Get the coordinates of representative atoms in all residues.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

            only_representative (bool):
                If True, only the representative atoms are considered for
                modified protein/dna/rna chains and ligands. \n
                (Default: False)

        Returns:

            coords_list (list):
                List containing the coordinates for representative atoms for
                each residue.
        """

        coords_list = []

        for residue, _chain_id in self.get_residues(structure):

            coords, _ = self.extract_perresidue_quantity(
                residue=residue,
                quantity="coords",
            )

            coords_list.extend(coords)

        return coords_list

    def get_rep_plddt(
        self,
        structure: Bio.PDB.Structure.Structure,
    ):
        """Get the pLDDT values for representative atoms in all residues.

        Args:

            structure (Bio.PDB.Structure.Structure):
                Biopython Structure object.

            only_representative (bool):
                If True, only the representative atoms are considered for
                modified protein/dna/rna chains and ligands. \n
                (Default: False)

        Returns:

            plddt_list (list):
                List containing the pLDDT values for representative atoms for
                each residue.
        """

        plddt_list = []

        for residue, _chain_id in self.get_residues(structure):

            plddt, _ = self.extract_perresidue_quantity(
                residue=residue,
                quantity="plddt",
            )

            plddt_list.extend(plddt)

        return plddt_list