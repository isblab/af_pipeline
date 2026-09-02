import os
import time
import Bio.PDB
from io import StringIO
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO
from Bio.Align import PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1
from typing import List, Dict, Mapping
from af_pipeline.constants.af_constants import FileFormat
from af_pipeline.parser.structure_parser import StructureParser
from af_pipeline.utils.misc_utils import get_res_range_from_key
from af_pipeline.tools.structure_tools import (
    ResidueSelect,
    save_structure_obj
)


def split_structure_by_chain(
    structure: Structure,
) -> Dict[str, Structure]:
    """ Split a Biopython Structure object into separate Structure objects for
    each chain.

    ## Arguments:

    - **structure (Structure)**:<br />
        Biopython Structure object to be split by chain.

    ## Returns:

    - **dict[str, Structure]**:<br />
        A dictionary mapping chain IDs to new Structure objects containing only
        that chain.
    """

    chain_structures = {}
    model = structure[0]  # Get the first model

    for chain in model:
        chain_id = chain.get_id()
        new_structure: Structure = Structure(f"{structure.get_id()}_{chain_id}")
        new_structure.add(Model(0))
        new_structure[0].add(chain.copy())
        # detach hetero residues from the new structure
        for res in list(new_structure[0][chain_id]):
            if res.id[0] != " ":  # Skip hetero residues
                new_structure[0][chain_id].detach_child(res.id)
        chain_structures[chain_id] = new_structure

    return chain_structures


def query_to_hit_mapping(
    query_aligned: str,
    template_aligned: str
) -> Mapping[int, int]:
    """ 0-based query index to hit index mapping.

    from https://github.com/google-deepmind/alphafold3

    ## Arguments:

    - **query_aligned (str)**:<br />
        Aligned query sequence.

    - **template_aligned (str)**:<br />
        Aligned template sequence.

    ## Returns:

    - **Mapping[int, int]**:<br />
        A mapping from 0-based query indices to 0-based hit indices.
    """

    query_to_hit_mapping_out = {}
    hit_index = 0
    query_index = 0

    for q_char, t_char in zip(query_aligned, template_aligned):

        # Gap inserted in the template
        if q_char == '-':
            hit_index += 1

        # Deleted residue in the template (would be a gap in the query).
        elif t_char == '-':
            query_index += 1

        # Normal aligned residue, in both query and template. Add to mapping.
        else:
            query_to_hit_mapping_out[query_index] = hit_index
            query_index += 1
            hit_index += 1

    return query_to_hit_mapping_out


def get_aligned_indices(
    query_seq: str,
    template_seq: str
) -> tuple[list[int], list[int]]:
    """ Align two sequences and map the indices

    ## Arguments:

    - **query_seq (str)**:<br />
        Query sequence to be aligned.

    - **template_seq (str)**:<br />
        Template sequence to be aligned.

    ## Returns:

    - **tuple**:<br />
        A tuple containing two lists:
        - List of query indices corresponding to the aligned template residues.
        - List of template indices corresponding to the aligned query residues.
    """

    # Perform pairwise alignment
    aligner = PairwiseAligner(scoring="blastp")
    alignments = aligner.align(seqA=query_seq, seqB=template_seq)
    alignment = next(iter(alignments))  # Take the best alignment
    query_aligned = alignment[1]
    template_aligned = alignment[0]

    # Map the aligned sequences
    aligned_mapping = query_to_hit_mapping(
        query_aligned=query_aligned,
        template_aligned=template_aligned
    )

    query_indices = []
    template_indices = []
    for template_index, query_index in aligned_mapping.items():
        query_indices.append(query_index)
        template_indices.append(template_index)

    return query_indices, template_indices


def get_custom_template_dict(
    input_seq: str,
    template_path: str,
    chain_id: str,
    residue_range: List[int] | None = None,
) -> Dict[str, str | list[int]]:
    """ Obtain dictionary for template settings for AlphaFold server JSON file

    For a given protein entity, return a dictionary:
    {
        "mmcif": mmcif_str,
        "queryIndices": query_indices,
        "templateIndices": template_indices,
    }

    ## Arguments:

    - **input_seq (str)**:<br />
        Input sequence for the query.

    - **template_path (str)**:<br />
        Path to the template file.

    - **chain_id (str)**:<br />
        Chain ID to be used from the template.

    - **residue_range (List[int] | None, optional):**:<br />
        Range of residues to be used from the template. If None, all residues are used.

    ## Returns:

    - **dict**:<br />
        Dictionary containing the template settings for the AlphaFold server JSON file.
    """

    if isinstance(residue_range, str):
        residue_range = get_res_range_from_key(residue_range)


    template_name = os.path.splitext(os.path.basename(template_path))[0]

    structure_parser = StructureParser(structure_file_path=template_path)
    structure = structure_parser.get_structure_obj()

    # collect metadata from the template mmCIF file
    mmcif_dict = MMCIF2Dict(template_path)
    headers_to_keep = [
        "_entry.id",
        "_entry.title",
        "_entry.deposition_date",
        "_pdbx_audit_revision_history.revision_date",
    ]
    filtered_metadata = {
        key: mmcif_dict[key] for key in headers_to_keep if key in mmcif_dict
    }

    # Make metadata if missing
    if "_pdbx_audit_revision_history.revision_date" not in filtered_metadata:
        filtered_metadata["_pdbx_audit_revision_history.revision_date"] = time.strftime(
            "%Y-%m-%d"
        )

    chain_structures = split_structure_by_chain(structure=structure)
    chain_path = os.path.join(
        os.path.dirname(template_path), f"{template_name}_{chain_id}.cif"
    )

    save_structure_obj(
        structure=chain_structures[chain_id],
        out_file=chain_path,
        res_select_obj=ResidueSelect({chain_id: residue_range}) if residue_range else None,
        save_type=FileFormat.CIF,
        preserve_header_footer=True
    )

    structure_parser = StructureParser(chain_path)
    structure = structure_parser.get_structure_obj()

    template_seq = ""
    for res, _ in structure_parser.get_residues(structure=structure):
        if res.id[0] != " ":  # Skip hetero residues
            continue
        # Get the first letter of the residue name
        template_seq += protein_letters_3to1.get(res.resname, "X")

    mmcif_dict = MMCIF2Dict(chain_path)
    mmcif_dict.update(filtered_metadata)

    io = MMCIFIO()
    stringio = StringIO()
    io.set_dict(mmcif_dict)
    io.save(stringio)

    mmcif_str = stringio.getvalue()

    query_indices, template_indices = get_aligned_indices(
        query_seq=input_seq,
        template_seq=template_seq,
    )

    return {
        "mmcif": mmcif_str,
        "queryIndices": query_indices,
        "templateIndices": template_indices,
    }


def pdb_to_mmcif(
    input_pdb: str,
    output_mmcif: str,
):
    """ Convert a .pdb file to a .cif file.

    ## Arguments:

    - **input_pdb (str)**:<br />
        Path to the input PDB file. Must have a .pdb extension.

    - **output_mmcif (str)**:<br />
        Path to the output mmCIF file. The directory will be created if it does not exist.
    """

    output_mmcif = os.path.abspath(output_mmcif)
    file_ext = os.path.splitext(input_pdb)[1].lower()
    if file_ext != ".pdb":
        raise ValueError(
            f"Input file must be a PDB file with .pdb extension, got {file_ext}"
        )

    os.makedirs(os.path.dirname(output_mmcif), exist_ok=True)

    structure = Bio.PDB.PDBParser().get_structure("structure", input_pdb)
    io = Bio.PDB.MMCIFIO()
    io.set_structure(structure)

    if not output_mmcif.lower().endswith(".cif"):
        output_mmcif += ".cif"

    io.save(output_mmcif)