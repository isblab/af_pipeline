"""
[rigid_bodies](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/rigid_bodies)
=============================================================

- Module to extract confident region from AlphaFold predictions

- This module contains classes and methods for extracting **confidently predicted
  regions (rigid bodies)** from AlphaFold predictions based on the assessment
  metrics provided by AlphaFold: PAE and pLDDT.

- It uses graph-based community detection algorithms[^tristan] to identify
  communities of residues that are confidently predicted; a set of residues
  within each community satisfying a pLDDT cutoff constitute a rigid body.


<div align="center">
<img src="../assets/rb_extraction_module.png" alt="interaction module" width="90%"/>
</div>


[^tristan]: Tristan Croll, "Graph-based community clustering approach to extract protein domains from a predicted aligned error matrix": https://github.com/tristanic/pae_to_domains

<hr>

## Usage

- Please refer to the [examples directory](https://github.com/isblab/af_pipeline/tree/main/examples) for sample scripts and
  config file.

- Use the following command to run the example script for extracting rigid bodies:
  ```
  python extract_rigid_bodies.py \\
      -i ./input/best_af_predictions.json \\
      -o ./output/rigid_bodies \\
      --plddt_cutoff 70 \\
      --pae_cutoff 12 \\
      --pae_power 1 \\
      --resolution 0.5 \\
      --library networkx \\
      --min_res 10 \\
      --min_proteins 1 \\
      --apply_plddt_filter True
  ```

<hr>

## Organization

The module is organized into the following submodules:

- **rigid_bodies**: Contains the [`RigidBodies`](https://github.com/isblab/af_pipeline/blob/89b33286ff81b1f4075c65fd6d5dc88296d83c15/af_pipeline/rigid_bodies/rigid_bodies.py#L71) class, which can be used for extracting confidently predicted regions from AlphaFold predictions based on the assessment metrics provided by AlphaFold: PAE and pLDDT.

- **rigid_body_assessment**: Contains the [`RigidBodyAssessment`](https://github.com/isblab/af_pipeline/blob/89b33286ff81b1f4075c65fd6d5dc88296d83c15/af_pipeline/rigid_bodies/rigid_body_assessment.py#L1279) class, which can be used for assessing the quality of the extracted rigid bodies using various metrics such as average pLDDT, average PAE, and contact probabilities.

```mermaid
---
config:
  class:
    hideEmptyMembersBox: true
---
classDiagram

    class RigidBodies {
        + str library
        + float pae_cutoff
        + float plddt_cutoff
    }

    class RigidBodyAssessment {
        + Dict[str, List[Tuple[str, int]]] rb_dict
        + bool as_average
        + bool symmetric_pae
        + bool show_interface_residues_only
    }

    class RigidBodyChainAssessment {
        + bool as_average
        + bool show_interface_residues_only
    }

    class RigidBodyChainPairAssessment {
        + bool as_average
        + bool show_interface_residues_only
    }

    class PAEPatches {
        + dict num_to_idx
        + np.ndarray pae
        + dict lengths_dict
        + int rb_idx
    }

    class _Mask {
        + Dict[str, List[int]] rb_dict
        + np.ndarray pae
        + np.ndarray avg_pae
        + np.ndarray plddt_list
        + np.ndarray contact_map
        + Dict[str, int] lengths_dict
        + Dict[str, Dict[int, Dict[str, int]]] num_to_idx
        + Dict[int, Dict[str, int]] idx_to_num
    }

    link RigidBodyAssessment "rigid_bodies/rigid_body_assessment.html#RigidBodyAssessment" "link to RigidBodyAssessment class documentation"

    RigidBodies --> RigidBodyAssessment

    RigidBodyAssessment ..> RigidBodyChainAssessment

    link RigidBodyChainAssessment "rigid_bodies/rigid_body_assessment.html#RigidBodyChainAssessment" "link to RigidBodyChainAssessment class documentation"

    RigidBodyAssessment ..> RigidBodyChainPairAssessment

    link RigidBodyChainPairAssessment "rigid_bodies/rigid_body_assessment.html#RigidBodyChainPairAssessment" "link to RigidBodyChainPairAssessment class documentation"

    RigidBodyChainPairAssessment ..> _Mask
    RigidBodyChainAssessment ..> _Mask

    link RigidBodies "rigid_bodies/rigid_bodies.html#RigidBodies" "link to RigidBodies class documentation"

    RigidBodies --> PAEPatches

    link PAEPatches "rigid_bodies/rigid_bodies.html#PAEPatches" "link to PAEPatches class documentation"

    RigidBodies ..> pae_to_domains

    link pae_to_domains "pae_to_domains/pae_to_domains.html" "link to pae_to_domains module documentation"

    link _Mask "rigid_bodies/rigid_body_assessment.html" "link to _Mask class documentation"
```

<hr>

## Prerequisites

- **best_af_predictions.json**: obtained by running [`rank_af_predictions.py`](https://github.com/isblab/af_pipeline/tree/main/examples/rank_af_predictions.py).
  This file contains information about the best prediction for each job.

- **Structure predictions**: obtained by running the prediction jobs on AlphaFold server
  or AlphaFold2 or ColabFold.

<hr>

## Workflows

- Workflow for extracting rigid bodies from AlphaFold predictions:


```mermaid

graph TD
  A([best_af_predictions.json]) --> B[/input_dict/]
  B -- for each --> C[Initialize instance]
  click C "parser/initialize.html#Initialize" "Initialize class documentation" _blank
  C --> D[RigidBodies instance]
  click D "rigid_bodies/rigid_bodies.html#RigidBodies" "RigidBodies class documentation" _blank
  D --> E[extract_rigid_bodies]
  click E "rigid_bodies/rigid_bodies.html#RigidBodies.extract_rigid_bodies" "extract_rigid_bodies method documentation" _blank
  E --> F[/rigid bodies/]
  F --> G[save_rigid_bodies]
  click G "rigid_bodies/rigid_bodies.html#RigidBodies.save_rigid_bodies" "save_rigid_bodies method documentation" _blank
  G --> J(["Saved rigid bodies files (PDB/CIF and TXT/JSON)"])
  F --> H[assess_rigid_bodies]
  click H "rigid_bodies/rigid_body_assessment.html#RigidBodyAssessment.assess_rigid_bodies" "assess_rigid_bodies method documentation" _blank
  H --> K(["Assessment results (XLSX)"])
  F --> I[show_rigid_bodies_on_pae_matrix]
  click I "rigid_bodies/rigid_body_assessment.html#RigidBodyAssessment.show_rigid_bodies_on_pae_matrix" "show_rigid_bodies_on_pae_matrix method documentation" _blank
  I --> L(["Visualization of rigid bodies on PAE matrix (PNG)"])
```

"""