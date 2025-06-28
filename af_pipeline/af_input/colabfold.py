
from af_pipeline.af_input.alphafold2 import AlphaFold2
from typing import Any, Dict, List, Tuple

class ColabFold(AlphaFold2):
    """Class to handle the creation of ColabFold input files

    Attributes:

        input_yml (dict):
            Input dictionary containing job cycles and jobs.
            Usually loaded from a YAML file.

        protein_sequences (dict):
            Dictionary containing protein sequences.
            Format: {header: sequence}

        entities_map (dict):
            Mapping of entity headers to their corresponding sequences.
            Format: {protein: header}
    """

    def __init__(
        self,
        input_yml: Dict[str, List[Dict[str, Any]]],
        protein_sequences: Dict[str, str],
        entities_map: Dict[str, str] = {},
    ):
        """ Initialize the ColabFold class.

        Args:

            input_yml (dict):
                Input dictionary containing job cycles and jobs.
                Usually loaded from a YAML file.

            protein_sequences (dict):
                Dictionary containing protein sequences.
                Format: {header: sequence}

            entities_map (dict):
                Mapping of entity headers to their corresponding sequences.
                Format: {protein: header}
                Defaults to {}.
        """

        self.entities_map = entities_map
        self.protein_sequences = protein_sequences
        self.input_yml = input_yml

        super().__init__(
            input_yml=input_yml,
            protein_sequences=protein_sequences,
            entities_map=entities_map,
        )

    def create_colabfold_job_cycles(
        self,
    ) -> Dict[str, List[Tuple[Dict[str, str], str]]]:
        """Create job cycles for ColabFold

        each job cycle is a list of jobs. \n
        each job is a tuple of `sequences_to_add` and `job_name`. \n
        `sequences_to_add` is a dictionary of fasta sequences {header: sequence} \n

        Returns:

            job_cycles (dict):
                Dictionary of job cycles {job_cycle: job_list}
        """

        job_cycles = {}

        for job_cycle, jobs_info in self.input_yml.items():

            job_list = []

            for job_info in jobs_info:
                sequences_to_add, job_name = self.generate_job_entities(
                    job_info=job_info
                )

                fasta_dict = {job_name: ":\n".join(list(sequences_to_add.values()))}

                job_list.append((fasta_dict, job_name))

            job_cycles[job_cycle] = job_list

        return job_cycles