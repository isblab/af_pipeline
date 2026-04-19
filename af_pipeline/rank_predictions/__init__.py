"""
[rank_predictions](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/rank_predictions)
=====================================

- Methods to rank AlphaFold predictions.
- This module provides functionality to rank the structure predictions based on various metrics.
- Currently supports:
  - Ranking of AlphaFold3[^af3] predictions from AlphaFold server[^afserver].
  - Ranking of AlphaFold2[^af2]/ColabFold[^colabfold] predictions.
- Note the following directory structure for AlphaFold3 predictions:
```
predictions/
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
- Note the following directory structure for ColabFold predictions:
```
predictions/
└── job_set/
    ├── job_set_name_scores_rank_001_alphafold2_multimer_v3_model_3_seed_000.json
    ├── job_set_name_scores_rank_002_alphafold2_multimer_v3_model_5_seed_000.json
    ├── job_set_name_scores_rank_003_alphafold2_multimer_v3_model_2_seed_000.json
    ├── job_set_name_scores_rank_004_alphafold2_multimer_v3_model_1_seed_000.json
    ├── job_set_name_scores_rank_005_alphafold2_multimer_v3_model_4_seed_000.json
    ├── job_set_name_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_000.pdb
    ├── job_set_name_unrelaxed_rank_002_alphafold2_multimer_v3_model_5_seed_000.pdb
    ├── job_set_name_unrelaxed_rank_003_alphafold2_multimer_v3_model_2_seed_000.pdb
    ├── job_set_name_unrelaxed_rank_004_alphafold2_multimer_v3_model_1_seed_000.pdb
    ├── job_set_name_unrelaxed_rank_005_alphafold2_multimer_v3_model_4_seed_000.pdb
    └── config.json
```
- Note the following directory structure for AlphaFold2 predictions:
```
predictions/
└── job_set/
    ├── result_model_1_multimer_v3_pred_0.pkl/json
    ├── result_model_1_multimer_v3_pred_1.pkl/json
    ├── result_model_1_multimer_v3_pred_2.pkl/json
    ├── result_model_1_multimer_v3_pred_3.pkl/json
    ├── result_model_1_multimer_v3_pred_4.pkl/json
    ├── unrelaxed_model_1_multimer_v3_pred_0.pdb/cif
    ├── unrelaxed_model_1_multimer_v3_pred_1.pdb/cif
    ├── unrelaxed_model_1_multimer_v3_pred_2.pdb/cif
    ├── unrelaxed_model_1_multimer_v3_pred_3.pdb/cif
    ├── unrelaxed_model_1_multimer_v3_pred_4.pdb/cif
    └── ranking_debug.json
```
- The ranking is done at the job set level, and the best model is selected from
  all the jobs in the job set.
- The best model is selected based on the following metrics in order of priority:
    1. Ranking Score
    2. Model Index (lower is better)
    3. ipTM
    4. pTM
    5. Fraction Disordered
- The best way to use this module is to provide the path to the `predictions` directory containing all the predictions.
- To rank a subset of predictions, one can also provide a specific job set directory.

[^afserver]: AlphaFold Server. alphafoldserver.com. Available at: https://alphafoldserver.com/.

[^af3]: Abramson, J. et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024). (https://alphafoldserver.com/)

[^af2]: Jumper, J. et al. Highly Accurate Protein Structure Prediction with Alphafold. Nature 596, 583–589 (2021).

[^colabfold]: Mirdita, M. et al. ColabFold: making protein folding accessible to all. Nature Methods 19, 679–682 (2022). (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)

## Organization

The module is organized into the following submodule:
- **rank_af**: Contains the `RankAF3JobSet` and `RankAF2JobSet` classes, that are
  repsonsible for ranking the predictions for a given job set and
  selecting the best model.

```mermaid
---
config:
  class:
    hideEmptyMembersBox: true
---
classDiagram
    class RankAF2JobSet {
        + str job_set_dir
        + str job_set_name
        + int job_set_id
        + list af_input_jobs
        + str data_format
    }
    class RankAF3JobSet {
        + str job_set_dir
        + str job_set_name
        + int job_set_id
        + list af_input_jobs
    }

    link RankAF2JobSet "rank_predictions/rank_af.html#RankAF2JobSet" "link to RankAF2JobSet class documentation"
    link RankAF3JobSet "rank_predictions/rank_af.html#RankAF3JobSet" "link to RankAF3JobSet class documentation"
```

## Prerequisites

- **af_input_jobs.json**: obtained by running [`create_af_jobs.py`](https://github.com/isblab/af_pipeline/tree/main/examples/create_af_jobs.py).
   This file contains the input specifications for the prediction jobs and is
   used to match the predictions with the corresponding input specifications.

- **AF3 predictions**: obtained by running the prediction jobs on AlphaFold server[^afserver].
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
  click B "rank_predictions/rank_af.html#get_job_set_dirs" "get_job_set_dirs" _blank
  C -- for each --> D[RankAF3JobSet instance]
  click D "rank_predictions/rank_af.html#RankAF3JobSet" "RankAF3JobSet" _blank
  D --> E[assign_job_set_id]
  click E "rank_predictions/rank_af.html#assign_job_set_id" "assign_job_set_id" _blank
  E --> F[extract_af3_best_pred_data]
  click F "rank_predictions/rank_af.html#RankAF3JobSet.extract_af3_best_pred_data" "extract_af3_best_pred_data" _blank
  F -- update --> G[/best_predictions/]
  G --> I([best_af_predictions.json])

```


- Workflow to rank AF2 or ColabFold predictions:

```mermaid

graph TD
  A[pred_dir] --> B[get_job_set_dirs]
  B[get_job_set_dirs] --> C[/job_set_dirs/]
  click B "rank_predictions/rank_af.html#get_job_set_dirs" "get_job_set_dirs" _blank
  C -- for each --> D[RankAF2JobSet instance]
  click D "rank_predictions/rank_af.html#RankAF2JobSet" "RankAF2JobSet" _blank
  D --> E[assign_job_set_id]
  click E "rank_predictions/rank_af.html#assign_job_set_id" "assign_job_set_id" _blank
  E --> F[extract_af2_best_pred_data]
  E --> J[extract_colabfold_best_pred_data]
  click F "rank_predictions/rank_af.html#RankAF2JobSet.extract_af2_best_pred_data" "extract_af2_best_pred_data" _blank
  click J "rank_predictions/rank_af.html#RankAF2JobSet.extract_colabfold_best_pred_data" "extract_colabfold_best_pred_data" _blank
  F -- update --> G[/best_predictions/]
  J -- update --> G[/best_predictions/]
  G --> I([best_af_predictions.json])

```

"""