"""
[parser](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/parser)
============================================

- Parser for processing structure predictions
- This module contains classes and methods for parsing the of structure predictions. It is designed to handle the output from various sources, including:
    - Predictions from the AlphaFold3[^af3], AlphaFold2[^af2] and ColabFold[^colabfold]
    - Predictions from the AlphaFold Protein Structure Database[^afdb]

[^afdb]: Fleming, J et al. AlphaFold Protein Structure Database and 3D-Beacons: New Data and Capabilities. Journal of Molecular Biology, 168967 (2025). (https://alphafold.ebi.ac.uk/)

[^af3]: Abramson, J. et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024). (https://alphafoldserver.com/)

[^af2]: Jumper, J. et al. Highly Accurate Protein Structure Prediction with Alphafold. Nature 596, 583–589 (2021).

[^colabfold]: Mirdita, M. et al. ColabFold: making protein folding accessible to all. Nature Methods 19, 679–682 (2022). (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)

<hr>

## Usage

- User mainly creates the [`Initialize`](parser/initialize.html) class instance during one of the analysis workflows (extracting and assessing rigid bodies or extracting interactions). Please refer to the [examples directory](https://github.com/isblab/af_pipeline/tree/main/examples) for sample scripts.
    - [`extract_interacting_patches.py`](https://github.com/isblab/af_pipeline/blob/main/examples/extract_interacting_patches.py#L118)
    - [`extract_rigid_bodies.py`](https://github.com/isblab/af_pipeline/blob/main/examples/extract_rigid_bodies.py#L170)

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
    }

    class StructureParser {
        + str structure_file_path
        + bool preserve_header_footer
        + bool use_fast_cif_parser
    }

    class DataParser {
        + str data_file_path
    }

    link Initialize "parser/initialize.html" "link to Initialize class documentation"

    Initialize *-- StructureParser

    link StructureParser "parser/structure_parser.html" "link to StructureParser class documentation"

    Initialize *-- DataParser

    link DataParser "parser/data_parser.html" "link to DataParser class documentation"

    Initialize *-- RenumberResidues

    link RenumberResidues "tools/structure_tools.html" "link to RenumberResidues class documentation"
```

"""