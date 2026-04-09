```mermaid
---
title: rank_af.py
---
classDiagram
    class RankAF3JobSet {
        + str job_set_dir
        + str job_set_name
        + str job_cycle_name
        + int job_set_id
        + bool try_af_offset_from_path
        - \_\_init__(self, job_set_dir, try_af_offset_from_path) None
        + extract_af3_best_pred_data(self, af_input_jobs) list
        + rank_seeds(self) tuple
        + get_af3_model_metrics_per_seed(self) dict | None
        + @staticmethod parse_af3_summary_confidences(af3_summary_confidence_file) dict$
        + @staticmethod get_seed_from_job_request(af3_job_request_file) int$
        + @staticmethod get_data_path_from_structure_path(structure_path) str$
        + extract_af_offset_from_af_input_jobs(self, af_input_jobs) dict
        + @staticmethod extract_af_offset_from_path(structure_path) dict$
        + add_job_set_id(self, af_input_jobs, soft_match) int
    }
```
