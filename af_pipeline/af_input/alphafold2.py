"""
AlphaFold2 input file creator
=============================
- Create input `FASTA` files for AlphaFold2 jobs.
- For AlphaFold2, only `proteinChain` entities are supported.
"""

import os
import warnings
from collections import defaultdict
from typing import Any, Dict, List, Tuple

from af_pipeline.constants.af_constants import RES_RANGE_SEP

class AlphaFold2:
    """Class to create FASTA files for AlphaFold2 jobs."""

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

    def create_af2_job_cycles(
        self
    ) -> Dict[str, List[Tuple[Dict[str, str], str]]]:
        """Create job cycles for AlphaFold2.

        Convert the input information into the format required by
        the AlphaFold2.

        Each job within a cycle is a tuple ->
        (`sequences_to_add`, `job_name`)<br />
        where, `sequences_to_add` = `{identifier: sequence}`

        Returns:

        - **job_cycles (dict)**:<br />
            Dictionary with:<br />

            - `key` -> `job_cycle_id` <br />
                Unique string identifier for the job cycle.<br />

            - `val` -> `job_list` <br />
        """

        job_cycles = {}

        for job_cycle_id, jobs_info in self.input_dict.items():

            job_list = []

            for job_info in jobs_info:
                sequences_to_add, job_name = self.generate_job_entities(
                    job_info=job_info
                )
                job_list.append((sequences_to_add, job_name))

            job_cycles[job_cycle_id] = job_list

        return job_cycles

    @staticmethod
    def write_to_fasta(
        fasta_dict: Dict[str, str],
        file_name: str,
        output_dir: str,
    ):
        """Write the sequences to a `FASTA` file

        Arguments:

        - **fasta_dict (dict)**:<br />
            Dictionary with:<br />

            - `key` -> `identifier` <br />

            - `val` -> `sequence` <br />

        - **file_name (str)**:<br /> Name of the output file

        - **output_dir (str)**:<br /> Directory to save the output fasta file.
        """

        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"{file_name}.fasta")

        with open(save_path, "w") as f:
            for identifier, sequence in fasta_dict.items():
                f.write(f">{identifier}\n{sequence}\n")

        print(f"\nFasta file written to {save_path}")

    @staticmethod
    def write_job_files(
        job_cycles: Dict[str, List[Tuple[Dict[str, str], str]]],
        output_dir: str,
    ):
        """Write the generated job files to the output directory.

        Arguments:

        - **job_cycles (dict)**:<br />
            Dictionary with:<br />

            - `key` -> `job_cycle_id` <br />

            - `val` -> `job_list` <br />

        - **output_dir (str)**:<br /> Output directory to save the job files.
        """

        for job_cycle, job_list in job_cycles.items():

            os.makedirs(os.path.join(output_dir, job_cycle), exist_ok=True)

            for fasta_dict, job_name in job_list:

                AlphaFold2.write_to_fasta(
                    fasta_dict=fasta_dict,
                    file_name=job_name,
                    output_dir=os.path.join(output_dir, job_cycle),
                )

        print("\nAll job files written to", output_dir)

    def generate_job_entities(
        self,
        job_info: Dict[str, Any],
    ) -> Tuple[Dict[str, str], str]:
        """Generate job entities.

        Job entities are the collection of entities within a job.<br />
        Each entity is a `proteinChain` with a identifier and sequence. \n

        Arguments:

        - **job_info (dict)**:<br /> Job information (`name`, `range`, `count`, `type`)

        Returns:

        - **(tuple)**:<br /> `(sequences_to_add, job_name)`
        """

        # get the job name if provided
        job_name = job_info.get("job_set_name", None)

        # get the information for each proteinChain
        identifiers = self.get_entity_info(job_info, "name", None)
        ranges = self.get_entity_info(job_info, "range", None)
        counts = self.get_entity_info(job_info, "count", 1)

        sequences = self.get_entity_sequences(ranges=ranges, identifiers=identifiers)

        job_dict = {
            "job_name": job_name,
            "entities": [],
        }

        for entity_count, (identifier, sequence, range_, count_) in enumerate(
            zip(identifiers, sequences, ranges, counts)
        ):
            for count_ in range(1, count_ + 1):
                job_dict["entities"].append(
                    {
                        "identifier": identifier,
                        "sequence": sequence,
                        "range": range_ if range_ else [1, len(sequence)],
                        "count": count_,
                    }
                )

        # generate job name if not provided
        if not job_name:
            job_name = self.generate_job_name(job_dict)

        # create fasta dictionary for each job {header: sequence}
        sequences_to_add = {}

        for entity in job_dict["entities"]:
            for entity_count in range(1, entity["count"] + 1):
                identifier = entity["identifier"]
                sequence = entity["sequence"]
                start, end = entity["range"]

                sequences_to_add[
                    f"{identifier}_{entity_count}_{start}{RES_RANGE_SEP}{end}"
                ] = sequence

        # warn if any of the entities is not a proteinChain
        self.warning_not_protien(job_info, job_name)

        return (sequences_to_add, job_name)

    def get_entity_info(
        self,
        job_info: Dict[str, Any],
        info_type: str,
        default_val: Any
    ) -> List[Dict[str, Any]]:
        """Get the entity information.

        Arguments:

        - **job_info (dict)**:<br />
            Job information (`name`, `range`, `count`, `type`)

        - **info_type (str)**:<br />
            Information type to get (`name`, `range`, `count`, `type`)

        - **default_val (Any)**:<br />
            Default value of the `info_type` if not found

        Returns:

        - **(list)**:<br /> List of entity information for the given type
        """

        return [
            entity.get(info_type, default_val)
            for entity in job_info["entities"]
            if entity["type"] == "proteinChain"
        ]

    def get_entity_sequences(
        self,
        ranges: List[Tuple[int, int]],
        identifiers: List[str],
    ) -> List[str]:
        """Get the entity sequences for given entity identifiers.

        Try to get the sequence from the `protein_sequences`.<br />
        - First, try using the provided identifier.<br />
        - If not found, try using the entities_map to get the identifier.<br />
        - If still not found, raise an exception.

        Depending on the provided range, slice the sequence.

        Arguments:

        - **ranges (list)**:<br /> `[start, end]` of the entities.

        - **identifiers (list)**:<br /> List of entity identifiers.

        Returns:

        - **sequences (list)**:<br />
            List of entity sequences to be used as input for the prediction.
        """

        sequences = []

        for identifier in identifiers:
            try:
                sequences.append(self.protein_sequences[identifier])
            except KeyError:
                try:
                    sequences.append(
                        self.protein_sequences[self.entities_map[identifier]]
                    )
                except KeyError:
                    raise Exception(
                        f"Could not find the entity sequence for {identifier}"
                    )

        for i, range_ in enumerate(ranges):
            if range_:
                start, end = range_
                sequences[i] = sequences[i][start - 1 : end]

        return sequences

    def generate_job_name(
        self,
        job_dict: Dict[str, Any],
    ) -> str:
        """Generate job name (if not provided).

        Arguments:

        - **job_dict (dict)**:<br /> Job dictionary.

        Returns:
        - **job_name (str)**:<br /> Job name.
        """

        job_name = ""

        fragments = defaultdict(list)

        for entity in job_dict["entities"]:
            identifier = entity["identifier"]
            start, end = entity["range"]
            count = entity["count"]

            fragments[f"{identifier}_{start}{RES_RANGE_SEP}{end}"].append(count)

        fragments = {k: max(v) for k, v in fragments.items()}

        for identifier, count in fragments.items():
            identifier_, range_ = identifier.split("_")
            job_name += f"{identifier_}_{count}_{range_}_"

        job_name = job_name[:-1] if job_name[-1] == "_" else job_name

        return job_name

    @staticmethod
    def warning_not_protien(
        job_info: Dict[str, Any],
        job_name: str
    ):
        """Warn if the entity is not a protein.

        AF2/ ColabFold only supports proteinChain entities.<br />
        Will skip the entities which are not proteins.

        Arguments:

        - **job_info (dict)**:<br /> Job information.

        - **job_name (str)**:<br /> Job name.
        """

        if any([
                entity_type != "proteinChain"
                for entity_type in [
                    entity["type"] for entity in job_info["entities"]
                ]
        ]):
            warnings.warn(
                f"""

                AF2/ ColabFold only supports proteinChain entities.
                Will skip the entities which are not proteins.
                {job_name} will be created with only proteinChain entities.
                """
            )