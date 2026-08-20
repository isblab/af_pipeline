"""
[interaction](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/interaction)
=============================================================

- Module to handle interaction data for the predicted structure

- This module contains classes and methods for extracting confidently predicted
interactions between residue-pairs from AlphaFold predictions based on the
assessment metrics provided by AlphaFold: PAE and pLDDT.

<div align="center">
<img src="../assets/interaction_module.png" alt="interaction module" width="90%"/>
</div>

<hr>

## Usage

- Please refer to the [examples directory](https://github.com/isblab/af_pipeline/tree/main/examples)
  for sample scripts and config file.

- Use the following command to run the example script for extracting interacting patches:
  ```
  python extract_interacting_patches.py \\
      --i ./input/best_af_predictions.json \\
      --o ./output/interacting_patches \\
      --interaction_pae_cutoff 5.0 \\
      --plddt_cutoff 70.0 \\
      --contact_threshold 8.0
  ```

<hr>

## Organization

The module is organized into the following submodules:

- **interaction**: Contains the [`Interaction`](https://github.com/isblab/af_pipeline/blob/89b33286ff81b1f4075c65fd6d5dc88296d83c15/af_pipeline/interaction/interaction.py#L61) class, which can be used for extracting confidently predicted interactions between residue-pairs from AlphaFold predictions based on the assessment metrics provided by AlphaFold: PAE and pLDDT.

```mermaid
---
config:
  class:
    hideEmptyMembersBox: true
---
classDiagram

    class Interaction {
        + float contact_threshold
        + float plddt_cutoff
        + float pae_cutoff
        + Optional[float] plddt_cutoff_idr
        + Optional[list] idr_chains
        + Optional[bool] save_plot
        + Optional[bool] save_table
    }

    link Interaction "interaction/interaction.html#Interaction" "link to Interaction class documentation"

    Interaction ..> MatrixPatches

    link MatrixPatches "tools/matrix_patches.html#MatrixPatches" "link to MatrixPatches class documentation"

```

<hr>

## Prerequisites

- **best_af_predictions.json**: obtained by running [`rank_af_predictions.py`](https://github.com/isblab/af_pipeline/tree/main/examples/rank_af_predictions.py).
  This file contains information about the best prediction for each job.

- **Structure predictions**: obtained by running the prediction jobs on AlphaFold server
  or AlphaFold2 or ColabFold.

<hr>

## Workflows

- Workflow for extracting interacting patches from AF predictions:

```mermaid

graph TD
  A([data_file_path]) --> D[Initialize instance]
  click D "parser/initialize.html#Initialize" "Initialize" _blank
  B([structure_file_path]) --> D
  D --> E[Interaction instance]
  click E "interaction/interaction.html#Interaction" "Interaction" _blank
  E --> F[create_regions_of_interest] --> G[/regions_of_interest/]
  click F "interaction/interaction.html#Interaction.create_regions_of_interest" "create_regions_of_interest" _blank
  E --> H[save_ppair_interaction]
  G --> H
```

"""