# AF-Pipeline

[![codecov](https://codecov.io/gh/isblab/af_pipeline/branch/main/graph/badge.svg)](https://codecov.io/gh/isblab/af_pipeline)

<img src="./docs/assets/af_pipeline_logo.png" alt="af_pipeline_logo" width="100%">

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

## Overview

[AF-Pipeline](https://isblab.github.io/af_pipeline/) is a package to assist in
AlphaFold2 and AlphaFold3 related tasks. These include:
- Creating input files for prediction (for AlphaFold server or AlphaFold2 or Colabfold)
- Ranking the predictions based on confidence metrics
- Extracting confidently predicted regions from the predictions
- Extracting interacting regions from the predictions

The workflow for the entire pipeline can be viewed [here](#workflow).

See also:
- [Changelog](changelog.md)
- [Contributing](contributing.md)
- [Documentation](https://isblab.github.io/af_pipeline/)

## Usage

- Example scripts are available in the [examples](./examples/) directory.

## Workflow

```mermaid

graph TD
    A([Create job files for AlphaFold Server]) e1@--> B[[AlphaFoldServer]]
    click B "af_pipeline/af_input/alphafold3.html#AlphaFoldServer" "AlphaFoldServer" _blank
    B e2@--> C[/Input FASTA or JSON files for AlphaFold/] e3@--> D([Submit jobs to AlphaFold server])
    click D "https://alphafoldserver.com" "alphafoldserver" _blank
    D e4@--> E[/Output files from AlphaFold server/]
    E e5@--> F[[Initialize]]
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
    class e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15 animate

```

## Additional Information

**License**: GPLv3

**Testable**: Yes