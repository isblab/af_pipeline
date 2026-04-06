"""
[colabfold](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/af_input/colabfold.py)
============================
- Create input `FASTA` files for ColabFold jobs.
- For ColabFold, only `proteinChain` entities are supported.
- Since most of the code is similar to AlphaFold2 input file creator,
  this class inherits from `af_pipeline.af_input.alphafold2.AlphaFold2`.
"""

from af_pipeline.af_input.alphafold2 import AlphaFold2
from typing import Any, Dict, List, Tuple

class ColabFold(AlphaFold2):
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

    def create_colabfold_job_cycles(
        self,
    ) -> Dict[str, List[Tuple[Dict[str, str], str]]]:
        """Create job cycles for ColabFold

        Convert the input information into the format required by
        the ColabFold.

        Each job within a cycle is a tuple -> (`fasta_dict`, `job_name`)<br />
        where, `fasta_dict` = `{job_name: all_sequence_str}`.

        `all_sequence_str` is concatenation of all sequences in the job
        joined by ":".
        """

        self.job_cycles = {}

        for job_cycle, jobs_info in self.input_dict.items():

            job_list = []

            for job_info in jobs_info:
                sequences_to_add, job_name = self.generate_job_entities(
                    job_info=job_info
                )

                fasta_dict = {job_name: ":\n".join(list(sequences_to_add.values()))}

                job_list.append((fasta_dict, job_name))

            self.job_cycles[job_cycle] = job_list