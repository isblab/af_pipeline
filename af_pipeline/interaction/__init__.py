"""
[interaction](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/interaction)
=============================================================

- Module to handle interaction data for the predicted structure

- This module contains classes and methods for extracting confidently
interactions between residue-pairs from AlphaFold predictions based on the
assessment metrics provided by AlphaFold: PAE and pLDDT.

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
        - \_\_init__(self, contact_threshold, plddt_cutoff, pae_cutoff, **kwargs) None
        + check_is_set_up(self)
        + set_attributes_from(self, instance)
        + @staticmethod get_contact_map(coords1, coords2, contact_threshold) np.ndarray$
        + create_regions_of_interest(self) list
        + get_interaction_data(self, region_of_interest) tuple
        + apply_confidence_cutoffs(self, plddt1, plddt2, avg_pae) tuple
        + get_confident_interaction_map(self, region_of_interest) np.ndarray
        + get_interacting_patches(self, contact_map, region_of_interest) dict
        + save_ppair_interaction(self, region_of_interest, output_dir, save_plot, plot_type, p1_name, p2_name, concat_residues, contact_probability)
    }

    link Interaction "interaction/interaction.html#Interaction" "link to Interaction class documentation"

    Interaction ..> MatrixPatches

    link MatrixPatches "tools/matrix_patches.html#MatrixPatches" "link to MatrixPatches class documentation"

```

## Usage

- Please refer to the [examples directory](https://github.com/isblab/af_pipeline/tree/main/examples) for sample scripts and
  config file.

- Use the following command to run the example script for extracting interacting patches:
```
python extract_interacting_patches.py \\
    --i ./input/best_af_predictions.json \\
    --o ./output/interacting_patches \\
    --interaction_pae_cutoff 5.0 \\
    --plddt_cutoff 70.0 \\
    --contact_threshold 8.0
```

## Prerequisites

- **best_af_predictions.json**: obtained by running [`rank_af_predictions.py`](https://github.com/isblab/af_pipeline/tree/main/examples/rank_af_predictions.py).
  This file contains information about the best prediction for each job.

- **Structure predictions**: obtained by running the prediction jobs on AlphaFold server
  or AlphaFold2 or ColabFold.

## Workflows

- Workflow for extracting interacting patches from AF predictions:

```mermaid
sequenceDiagram
    autonumber
    participant User
    participant Initialize
    participant Interaction

    rect rgb(210, 250, 200)
    note over User, Interaction: Set up Interaction instance.
    User->>Initialize:
    User->>Interaction: parameters & instance of Initialize
    end

    rect rgb(240, 255, 255)
    note over User, Interaction: Extract and save interacting patches.
    User->>Interaction: create_regions_of_interest()
    User->>Interaction: save_ppair_interaction()
    end
```

"""