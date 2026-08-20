# AF-Pipeline

[![codecov](https://codecov.io/gh/isblab/af_pipeline/branch/main/graph/badge.svg)](https://codecov.io/gh/isblab/af_pipeline)

<img src="./docs/assets/af_pipeline_logo.png" alt="af_pipeline_logo" width="100%">

<hr>

## Overview & features

<img src="./docs/assets/af_pipeline_workflow1.png" alt="af_pipeline_logo" width="100%">

[AF-Pipeline](https://isblab.github.io/af_pipeline/) is a package to assist in
AlphaFold2[^af2] and AlphaFold3[^af3] related tasks. These include:
- Creating input files for prediction (for AlphaFold server[^afserver] or AlphaFold2[^af2] or Colabfold[^colabfold])
- Ranking the predictions based on confidence metrics
- Extracting **confidently predicted regions** from the predictions (based on [pae_to_domains](https://github.com/tristanic/pae_to_domains)[^pae_to_domains] package by Tristan Croll)
- Extracting **interacting residue patches** from the predictions

The workflow for the entire pipeline can be viewed [here](#workflow).

<hr>

## Installation

- Clone the repository
    ```bash
    git clone --recursive https://github.com/isblab/af_pipeline.git
    ```

  If cloned without `--recursive` flag, do the following.
  ```bash
  git submodule init
  git submodule update
  ```

- Run [`setup.py`](./setup.py)
  ```bash
  python setup.py
  ```

<hr>

## Tutorial

- Refer to the following Google colab notebook and slides (presented at Computational Biology Summer School 2026 at NCBS Bangalore)

  [![Colab notebook](https://img.shields.io/badge/google_colab-F9AB00?style=for-the-badge&logo=google-colab&logoColor=white&label)](https://colab.research.google.com/drive/18Xddb_Wm-_GVpCC1m0qA3p1SiTE2yGVF?usp=drive_link)
  [![Slides](https://img.shields.io/badge/slides-yellow?style=for-the-badge&logo=google-slides&logoColor=white&label)](https://github.com/isblab/af_pipeline/blob/main/tutorials/AF_pipeline_NCBS_summer_school_2026.pdf)

<hr>

## Documentation

- For a detailed explanation of different modules and their capabilites in AF-Pipeline, see [Documentation](https://isblab.github.io/af_pipeline/)

- See also:
  - [Changelog](changelog.md)
  - [Contributing](contributing.md)

<hr>

## Quick usage

- Refer to the scripts available in the [examples](./examples/) directory.
  Before running the scripts, make sure you unzip [`AF_predictions.zip`](./examples/input/AF_predictions.zip)
  ```bash
  unzip ./input/AF_predictions.zip -d input/
  ```

- Create `JSON` job files for AlphaFold server [^afserver].
  ```bash
  python create_af_jobs.py \\
      -i ./input/config.yaml \\
      -o ./output/af_input_jobs \\
      -p ./input/protein_sequences.fasta \\
      -n ./input/nucleic_acid_sequences.fasta \\
      -t AF3 # you can replace this with AF2 or ColabFold
  ```

- Upload the `JSON` files to AlphaFold server [^afserver] and download and extract the results.

- Rank the output predictions.
  ```bash
  python rank_af_predictions.py \\
      -i ./output/af_input_jobs.json \\
      -o ./output \\
      -t AF3 \\
      --pred_dirs ./input/AF_predictions/AF3
  ```

- Extract interacting patches from the predictions.
  ```bash
  python extract_interacting_patches.py \\
      --i ./input/best_af_predictions.json \\
      --o ./output/interacting_patches \\
      --interaction_pae_cutoff 5.0 \\
      --plddt_cutoff 70.0 \\
      --contact_threshold 8.0
  ```

- Extract confidently predicted regions from the predictions.
  ```bash
  python extract_rigid_bodies.py \\
      -i ./input/best_af_predictions.json \\
      -o ./output/rigid_bodies \\
      --apply_plddt_filter True
  ```

[^afserver]: AlphaFold Server. alphafoldserver.com. Available at: https://alphafoldserver.com/.

[^af3]: Abramson, J. et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024). (https://alphafoldserver.com/)

[^af2]: Jumper, J. et al. Highly Accurate Protein Structure Prediction with Alphafold. Nature 596, 583–589 (2021).

[^colabfold]: Mirdita, M. et al. ColabFold: making protein folding accessible to all. Nature Methods 19, 679–682 (2022). (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)

[^pae_to_domains]: pae_to_domains: Graph-based community clustering approach to extract protein domains from a predicted aligned error matrix. Tristan Croll (https://github.com/tristanic/pae_to_domains)

<hr>

## Workflow

```mermaid

graph TD
    A([Create job files for AlphaFold Server]) e1@--> B[[AlphaFoldServer]]
    click B "af_pipeline/af_input/alphafold3.html#AlphaFoldServer" "AlphaFoldServer" _blank
    B e2@--> C[/Input FASTA or JSON files for AlphaFold/] e3@--> D[\Submit jobs to AlphaFold server/]
    click D "https://alphafoldserver.com" "alphafoldserver" _blank
    D e4@--> E[/Output files from AlphaFold server/]
    E e16@--> R[[RankAF3JobSet]]
    click R "af_pipeline/rank_predictions/rank_af.html#RankAF3JobSet" "RankAF3JobSet" _blank
    R e17@--> S[/best prediction/]
    S e5@--> F[[Initialize]]
    click F "af_pipeline/parser/initialize.html#Initialize" "Initialize" _blank
    F e6@--> G[[RigidBodies]]
    click G "af_pipeline/rigid_bodies/rigid_bodies.html#RigidBodies" "RigidBodies" _blank
    F e7@--> H[[Interaction]]
    click H "af_pipeline/interaction/interaction.html#Interaction" "Interaction" _blank
    G e8@--> I[extract_rigid_bodies]
    click I "af_pipeline/rigid_bodies/rigid_bodies.html#RigidBodies.extract_rigid_bodies" "extract_rigid_bodies" _blank
    N e9@--> J[save_rigid_bodies]
    click J "af_pipeline/rigid_bodies/rigid_bodies.html#RigidBodies.save_rigid_bodies" "save_rigid_bodies" _blank
    N e10@--> K[assess_rigid_bodies]
    click K "af_pipeline/rigid_bodies/rigid_bodies.html#RigidBodies.assess_rigid_bodies" "assess_rigid_bodies" _blank
    H e11@--> M[save_ppair_interaction]
    click M "af_pipeline/interaction/interaction.html#Interaction.save_ppair_interaction" "save_ppair_interaction" _blank
    I e12@--> N(["Confidently predicted regions (rigid bodies) from AlphaFold predictions"])
    J e13@--> O[/"Output files for rigid bodies (txt/json and pdb)"/]
    K e14@--> P[/"Output files for rigid body assessment (xlsx)"/]
    M e15@--> Q[/"Output files for interacting patches (xlsx and plots)"/]
    classDef animate stroke-dasharray: 9,5,stroke-dashoffset: 900,animation: dash 25s linear infinite;
    class e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17 animate

```

<hr>

## Additional Information

**License**: [GPLv3](LICENSE.txt)

**Testable**: Yes
