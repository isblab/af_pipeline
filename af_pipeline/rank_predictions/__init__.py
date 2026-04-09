"""
[rank_predictions](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/rank_predictions)
=====================================

- Methods to rank AlphaFold predictions.
- This module provides functionality to rank AlphaFold3[^af3] predictions based on various metrics.
- Currently supports ranking of AlphaFold3 predictions from AlphaFold server.
- Note the following directory structure for AlphaFold3 predictions:
```
master_directory
    └── job_cycle/
        └── job_set/
            └── model_seed/
                ├── summary_confidences_0.json
                ├── summary_confidences_1.json
                ├── summary_confidences_2.json
                ├── summary_confidences_3.json
                ├── summary_confidences_4.json
                ├── job_request.json
                ├── full_data_0.json
                ├── full_data_1.json
                ├── full_data_2.json
                ├── full_data_3.json
                ├── full_data_4.json
                ├── model_0.cif
                ├── model_1.cif
                ├── model_2.cif
                ├── model_3.cif
                └── model_4.cif
```
- The ranking is done at the job set level, and the best model is selected from
  all the jobs in the job set.
- The best model is selected based on the following metrics in order of priority:
    1. Ranking Score
    2. Model Index (lower is better)
    3. ipTM
    4. pTM
    5. Fraction Disordered
- The best way to use this module is to provide the master directory containing all job cycles.
- To rank a subset of predictions, one can also provide a specific job cycle directory or job set directory.

[^af3]: Abramson, J. et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024). (https://alphafoldserver.com/)

## Organization

The rank_predictions module is organized into the following submodules:
- **rank_af**: Contains the `RankAF3JobSet` class, which is repsonsible for ranking
  the AlphaFold3 predictions for a given job set and selecting the best model based
  on the above-mentioned metrics.

```mermaid
---
config:
  class:
    hideEmptyMembersBox: true
---
classDiagram
    class RankAF3JobSet {
        + str job_set_dir
        + str job_set_name
        + str job_cycle_name
        + int job_set_id
        + bool try_af_offset_from_path
    }
```

## Prerequisites

- **af_input_jobs.json**: obtained by running [`create_af_jobs.py`](https://github.com/isblab/af_pipeline/tree/main/examples/create_af_jobs.py).
   This file contains the input specifications for the prediction jobs and is
   used to match the predictions with the corresponding input specifications.

- **AF3 predictions**: obtained by running the prediction jobs on AlphaFold server.
  The predictions are expected to be organized in the directory structure mentioned above.

## Usage

- Please refer to the [examples directory](https://github.com/isblab/af_pipeline/tree/main/examples) for sample scripts.

- Use the following command to run the example script:
```
python rank_af_predictions.py \\
    -i ./output/af_input_jobs.json \\
    -o ./output \\
    --pred_dirs ./input/AF_predictions/AF3
```

## Workflows

- Workflow to rank AF3 predictions:

```mermaid

graph TD
  A[pred_dir] --> B[get_job_set_dirs]
  B[get_job_set_dirs] --> C[/job_set_dirs/]
  C -- for each --> D[RankAF3JobSet instance]
  click D "rank_predictions/rank_af.html#RankAF3JobSet" "RankAF3JobSet" _blank
  D --> E[add_job_set_id]
  click E "rank_predictions/rank_af.html#RankAF3JobSet.add_job_set_id" "add_job_set_id" _blank
  E --> F[extract_af3_best_pred_data]
  click F "rank_predictions/rank_af.html#RankAF3JobSet.extract_af3_best_pred_data" "extract_af3_best_pred_data" _blank
  F -- update --> G[/best_predictions/]
  G --> I([best_af_predictions.json])

```

"""