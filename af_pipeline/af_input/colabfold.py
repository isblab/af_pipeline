"""
ColabFold input file creator
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

    input_dict: Dict[str, List[Dict[str, Any]]]
    """Dictionary with:<br />

    - `key` -> `job_cycle_id` <br />
      Unique string identifier for the job cycle.<br />

    - `val` -> `job_sets_list` <br />
      List of `AFJobSet.job_set_info`s, each of which specifies
      the entities, model seeds, job name, etc."""

    protein_sequences: Dict[str, str] | None
    """Dictionary with:<br />

    - `key` -> `identifier` <br />
       Usually `uniprot_id` in case of `proteinChain` entities.<br />
       `identifier != entity_name` necessitates `entities_map`.<br />

    - `val` -> `sequence` <br />
      Amino acid sequence of the protein chain.
    """

    entities_map: Dict[str, str]
    """Dictionary with:<br />

    - `key` -> `entity_name` <br />

    - `val` -> `identifier` <br />
      `identifier` is usually `uniprot_id` in case of `proteinChain` entities."""

    def __init__(
        self,
        input_dict: Dict[str, List[Dict[str, Any]]],
        protein_sequences: Dict[str, str],
        entities_map: Dict[str, str] = {},
    ):

        self.entities_map = entities_map
        self.protein_sequences = protein_sequences
        self.input_dict = input_dict

        super().__init__(
            input_dict=input_dict,
            protein_sequences=protein_sequences,
            entities_map=entities_map,
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

        Returns:

        - **job_cycles (dict)**:<br />
            Dictionary with:<br />

            - `key` -> `job_cycle_id` <br />
                Unique string identifier for the job cycle.<br />

            - `val` -> `job_list` <br />
        """

        job_cycles = {}

        for job_cycle, jobs_info in self.input_dict.items():

            job_list = []

            for job_info in jobs_info:
                sequences_to_add, job_name = self.generate_job_entities(
                    job_info=job_info
                )

                fasta_dict = {job_name: ":\n".join(list(sequences_to_add.values()))}

                job_list.append((fasta_dict, job_name))

            job_cycles[job_cycle] = job_list

        return job_cycles