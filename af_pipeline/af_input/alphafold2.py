import os
import warnings
from collections import defaultdict
from typing import Any, Dict, List, Tuple


class AlphaFold2:
    """Class to handle the creation of AlphaFold2 input files

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
        """ Initialize the AlphaFold2 class

        Args:

            input_yml (Dict[str, List[Dict[str, Any]]]):
                Input dictionary containing job cycles and jobs.

            protein_sequences (Dict[str, str]):
                Dictionary containing protein sequences.
                Format: {header: sequence}

            entities_map (Dict[str, str], optional):
                Mapping of entity headers to their corresponding sequences.
                Format: {protein: header}
                Defaults to an empty dictionary.
        """

        self.entities_map = entities_map
        self.protein_sequences = protein_sequences
        self.input_yml = input_yml

    def create_af2_job_cycles(
        self
    ) -> Dict[str, List[Tuple[Dict[str, str], str]]]:
        """Create job cycles for AlphaFold2

        Each job cycle is a list of jobs. \n
        Each job is a tuple of `sequences_to_add` and `job_name`. \n
        `sequences_to_add` is a dictionary of fasta sequences 
        {header: sequence} \n

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
                job_list.append((sequences_to_add, job_name))

            job_cycles[job_cycle] = job_list

        return job_cycles

    def write_to_fasta(
        self,
        fasta_dict: Dict[str, str],
        file_name: str,
        output_dir: str = "./output/af_input",
    ):
        """Write the fasta sequences to a file

        Args:

            fasta_dict (dict):
                Dictionary of fasta sequences {header: sequence}

            file_name (str):
                Name of the output file

            output_dir (str, optional):
                Directory to save the output fasta file.
                Defaults to "./output/af_input".
        """

        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"{file_name}.fasta")

        with open(save_path, "w") as f:
            for header, sequence in fasta_dict.items():
                f.write(f">{header}\n{sequence}\n")

        print(f"\nFasta file written to {save_path}")

    def write_job_files(
        self,
        job_cycles: Dict[str, List[Tuple[Dict[str, str], str]]],
        output_dir: str = "./output/af_input",
    ):
        """Write job files to the output directory

        Args:

            job_cycles (dict):
                Dictionary of job cycles {job_cycle: job_list}

            output_dir (str, optional):
                Output directory to save the job files.
                Defaults to "./output/af_input".
        """

        for job_cycle, job_list in job_cycles.items():

            os.makedirs(os.path.join(output_dir, job_cycle), exist_ok=True)

            for fasta_dict, job_name in job_list:

                self.write_to_fasta(
                    fasta_dict=fasta_dict,
                    file_name=job_name,
                    output_dir=os.path.join(output_dir, job_cycle),
                )

        print("\nAll job files written to", output_dir)

    def generate_job_entities(
        self,
        job_info: Dict[str, Any],
    ) -> Tuple[Dict[str, str], str]:
        """Generate job entities

        job entities are the collection of entities within a job. \n
        Each entity is a proteinChain with a header and sequence. \n

        Args:

            job_info (dict):
                job information (name, range, count, type)

        Returns:

            Tuple[Dict[str, str], str]: sequences_to_add, job_name
        """

        # get the job name if provided
        job_name = job_info.get("name", None)

        # get the information for each proteinChain
        headers = self.get_entity_info(job_info, "name", None)
        ranges = self.get_entity_info(job_info, "range", None)
        counts = self.get_entity_info(job_info, "count", 1)

        sequences = self.get_entity_sequences(ranges=ranges, headers=headers)

        job_dict = {
            "job_name": job_name,
            "entities": [],
        }

        for entity_count, (header, sequence, range_, count_) in enumerate(
            zip(headers, sequences, ranges, counts)
        ):
            for count_ in range(1, count_ + 1):
                job_dict["entities"].append(
                    {
                        "header": header,
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
                header = entity["header"]
                sequence = entity["sequence"]
                start, end = entity["range"]

                sequences_to_add[
                    f"{header}_{entity_count}_{start}to{end}"
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
        """Get the entity information

        Get the required information for each entity in the job

        Args:
            job_info (dict): job information (name, range, count, type)
            info_type (str): type of information to get (name, range, count, type)
            default_val (Any): default value if not found

        Returns:
            List[Dict[str, Any]]: list of entity information for the given type
        """

        return [
            entity.get(info_type, default_val)
            for entity in job_info["entities"]
            if entity["type"] == "proteinChain"
        ]

    def get_entity_sequences(
        self,
        ranges: List[Tuple[int, int]],
        headers: List[str],
    ) -> List[str]:
        """Get the entity sequences

        First try to get the sequence from the protein_sequences dictionary. \n
        If not found, try to get the sequence from the proteins dictionary. \n
        If not found, raise an exception.

        If a range is provided, get the sequence within the range.

        Args:

            ranges (list):
                [start, end] of the entities

            headers (list):
                fasta headers

        Returns:

            sequences (list):
                list of entity sequences
        """

        sequences = []

        for header in headers:
            try:
                sequences.append(self.protein_sequences[header])
            except KeyError:
                try:
                    sequences.append(
                        self.protein_sequences[self.entities_map[header]]
                    )
                except KeyError:
                    raise Exception(
                        f"Could not find the entity sequence for {header}"
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
        """Generate job name (if not provided)

        see :py:mod:`AFInput.generate_job_entities` for the job dictionary format.

        Args:

            job_dict (dict):
                job dictionary

        Returns:

            job_name (str):
                job name
        """

        job_name = ""

        fragments = defaultdict(list)

        for entity in job_dict["entities"]:
            header = entity["header"]
            start, end = entity["range"]
            count = entity["count"]

            fragments[f"{header}_{start}to{end}"].append(count)

        fragments = {k: max(v) for k, v in fragments.items()}

        for header, count in fragments.items():
            header_, range_ = header.split("_")
            job_name += f"{header_}_{count}_{range_}_"

        job_name = job_name[:-1] if job_name[-1] == "_" else job_name

        return job_name

    @staticmethod
    def warning_not_protien(
        job_info: Dict[str, Any],
        job_name: str
    ):
        """Warn if entity is not a protein

        AF2/ ColabFold only supports proteinChain entities. \n
        Will skip the entities which are not proteins. \n

        Args:

            job_info (dict):
                job information

            job_name (str):
                job name
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