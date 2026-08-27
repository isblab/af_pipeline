"""
[colabfold](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/af_input/colabfold.py)
============================
- Create input `FASTA` files for ColabFold jobs.
- For ColabFold, only `proteinChain` entities are supported.
- Since most of the code is similar to AlphaFold2 input file creator,
  this class inherits from `af_pipeline.af_input.alphafold2.AlphaFold2`.
"""

import os
import warnings
from af_pipeline.af_input.alphafold2 import AF2Config
from typing import Any, Dict, overload
from af_pipeline.constants.af_constants import (
    ConfigYaml,
    AFInputJobFields,
)
from af_pipeline.utils.file_utils import write_json

class ColabConfig(AF2Config):
    """Class to create FASTA files for ColabFold jobs."""

    def __init__(
        self,
        config_dict: Dict[str, Any],
        protein_sequences: Dict[str, str],
    ):

        super().__init__(
            config_dict=config_dict,
            protein_sequences=protein_sequences,
        )

    @overload
    def write_job_files(self, output_dir):
        ...

    def write_job_files(
        self,
        output_dir: str,
    ):
        """Convert the input information into the format required by ColabFold
        and write the job files to the output directory.

        Arguments:

        - **output_dir (str)**:<br /> Output directory to save the job files.
        """

        for job_set_idx, job_set_info in enumerate(self.input_job_sets):

            sequences_to_add, job_set_name = self.generate_job_entities(
                job_info=job_set_info
            )
            if len(sequences_to_add) == 0:
                warnings.warn(f"""
                    No valid entities found for job.
                    Skipping job file creation for this job."""
                )
                continue
            os.makedirs(output_dir, exist_ok=True)

            fasta_dict = {job_set_name: ":\n".join(list(sequences_to_add.values()))}

            AF2Config.write_to_fasta(
                fasta_dict=fasta_dict,
                file_name=job_set_name,
                output_dir=os.path.join(output_dir, job_set_name),
            )

            self.input_job_sets[job_set_idx] = {
                AFInputJobFields.JOB_SET_NAME: job_set_name,
                AFInputJobFields.AF_OFFSET: self.job_set_af_offset,
            }

            write_json(
                file_path=os.path.join(os.path.dirname(output_dir), "af_input_jobs.json"),
                data={
                    ConfigYaml.PROTEIN_UNIPROT_MAP: self.entities_map,
                    ConfigYaml.AF_INPUT_JOBS: self.input_job_sets,
                },
            )