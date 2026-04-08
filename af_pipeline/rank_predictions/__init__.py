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
sequenceDiagram
    autonumber
    participant User
    participant RankAF3JobSet

    rect rgb(210, 250, 200)
    User->>RankAF3JobSet: job_set_dir
    User->>RankAF3JobSet: add_job_set_id()
    create participant best_pred_info
    User->>best_pred_info: extract_af3_best_pred_data()
    best_pred_info->>best_predictions.json: write_json()
    end

```

"""