"""
Parser for processing structure predictions
============================================

This module contains classes and methods for parsing the of structure predictions.
It is designed to handle the output from various sources, including:
- Predictions from the AlphaFold3[^af3], AlphaFold2[^af2] and ColabFold[^colabfold]
- Predictions from the AlphaFold Protein Structure Database[^afdb]

[^afdb]: Fleming, J et al. AlphaFold Protein Structure Database and 3D-Beacons: New Data and Capabilities. Journal of Molecular Biology, 168967 (2025). (https://alphafold.ebi.ac.uk/)

[^af3]: Abramson, J. et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024). (https://alphafoldserver.com/)

[^af2]: Jumper, J. et al. Highly Accurate Protein Structure Prediction with Alphafold. Nature 596, 583–589 (2021).

[^colabfold]: Mirdita, M. et al. ColabFold: making protein folding accessible to all. Nature Methods 19, 679–682 (2022). (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)

<hr>

## Organization

The parser module is organized into the following submodules:
- **initialize**: Contains the `Initialize` class, which stores the parsed attributes after reading the structure (PDB or mmCIF) and data (JSON or PKL) files using -
    - `af_pipeline.parser.structure_parser.StructureParser`
    - `af_pipeline.parser.data_parser.DataParser`

- **structure_parser**: Contains the `af_pipeline.parser.structure_parser.StructureParser` class, which is responsible for parsing the structure files (PDB or mmCIF) and extracting relevant information such as coordinates, pLDDT scores, and chain ids.

- **data_parser**: Contains the `af_pipeline.parser.data_parser.DataParser` class, which is responsible for parsing the data files (JSON or pickle) and extracting relevant information such as PAE matrices, contact probabilities, and token ids (for AlphaFold3).

```mermaid
---
config:
  class:
    hideEmptyMembersBox: true
---
classDiagram
    class Initialize {
        + str structure_file_path
        + str data_file_path
        + dict af_offset
        + dict rep_atom_dict
        + bool average_token_pae
        + bool average_token_plddt
        + str metric_level
        + bool use_fast_cif_parser
        + Bio.PDB.Structure.Structure structure
        + StructureParser structure_parser
        + DataParser data_parser
        + np.ndarray avg_pae
        + Dict[str, int] lengths_dict
        + RenumberResidues renumber
        - \_\_init__(self, data_file_path, structure_file_path, af_offset, rep_atom_dict, average_token_pae, average_token_plddt, metric_level, use_fast_cif_parser) None
        + set_attributes(self) None
        + get_attributes(self, metric_level) None
        + @staticmethod get_chain_lengths(token_chain_ids) Dict[str, int]$
        + @staticmethod get_idxs_to_keep(structure, rep_atom_dict) Dict[tuple, int]$
        + update_pae(self, token_res_ids, token_chain_ids) np.ndarray
        + update_contact_probs(self, token_chain_ids, token_res_ids) np.ndarray | None
        + @staticmethod get_min_pae(pae_matrix, lengths_dict, along_axis, return_type) np.ndarray | Dict[str, list] | list$
    }

    link Initialize "parser/initialize.html" "link to Initialize class documentation"

    Initialize *-- StructureParser

    link StructureParser "parser/structure_parser.html" "link to StructureParser class documentation"

    Initialize *-- DataParser

    link DataParser "parser/data_parser.html" "link to DataParser class documentation"

    Initialize *-- RenumberResidues

    link RenumberResidues "tools/structure_tools.html" "link to RenumberResidues class documentation"
```

```mermaid
classDiagram

    class StructureParser {
        + str structure_file_path
        + bool preserve_header_footer
        + bool use_fast_cif_parser
        - \_\_init__(self, structure_file_path, preserve_header_footer, use_fast_cif_parser) None
        + structure_type(self) str
        + parser(self) PDBParser | MMCIFParser | FastMMCIFParser
        + get_structure_obj(self) Structure
        + @staticmethod get_residues(structure) Generator[tuple[Residue, str], Any, None]$
        + @staticmethod get_atoms(structure) Generator[tuple[Atom, Residue, str], Any, None]$
        + @staticmethod sanity_check_quantities(quantities, requested_on)$
        + sanity_check_not_orphan(requested_on)
        + @staticmethod extract_peratom_quantities(atom, quantities) dict[str, Any]$
        + @staticmethod extract_perresidue_quantities(residue, quantities, rep_atom) dict[str, Any]$
        + @staticmethod get_rep_atom(residue) Bio.PDB.Atom.Atom$
        + @staticmethod get_token_atom_names(structure, rep_atom_dict, only_representative) list[str]$
        + @staticmethod get_token_chain_ids(structure, rep_atom_dict, only_representative) list[str]$
        + @staticmethod get_token_res_ids(structure, rep_atom_dict, only_representative) list[str]$
        + @staticmethod get_plddt(structure, per_atom, rep_atom_dict, average_token_plddt, only_representative) list[float]$
        + @staticmethod get_coordinates(structure, per_atom, rep_atom_dict, only_representative) list[np.ndarray]$
    }

    link StructureParser "parser/structure_parser.html" "link to StructureParser class documentation"
```

```mermaid
classDiagram
    class DataParser {
        + str data_file_path
        - \_\_init__(self, data_file_path) None
        + data_type(self) str
        + parser(self) Callable[[str], Dict | List]
        + get_data_dict(self) Dict
        + @staticmethod get_token_chain_ids(data) list | None$
        + @staticmethod get_token_res_ids(data) list | None$
        + @staticmethod get_pae(data) np.ndarray$
        + @staticmethod get_contact_probs_mat(data) np.ndarray | None$
        + @staticmethod get_atom_chain_ids(data) list | None$
        + @staticmethod get_atom_plddts(data) np.ndarray | None$
    }

    link DataParser "parser/data_parser.html" "link to DataParser class documentation"
```

## Usage

- User mainly creates the `Initialize` class instance during one of the analysis workflows (extracting and assessing rigid bodies or extracting interactions). Please refer to the [examples directory](https://github.com/isblab/af_pipeline/tree/main/examples) for sample scripts.
    - [`extract_interacting_patches.py`](https://github.com/isblab/af_pipeline/blob/main/examples/extract_interacting_patches.py#L118)
    - [`extract_rigid_bodies.py`](https://github.com/isblab/af_pipeline/blob/main/examples/extract_rigid_bodies.py#L170)

## Network

<body>
    <p>
        Double click on the node to go to the corresponding line in the source code.
    </p>
    <iframe src="../../docs/network_viz/network_parser.html" width="800" height="600" frameborder="0">
    </iframe>
</body>
"""