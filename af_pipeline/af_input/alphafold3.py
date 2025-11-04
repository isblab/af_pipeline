"""
AlphaFold server input file creator
===================================
- Create input `JSON` files for AlphaFold server (https://alphafoldserver.com)
- Note the following hierarchy:
```
    job_cycle
    └── job_sets_list or job_list
        └── job (seeded from job_set)
            └── entities (list of AFSequence)
```
- A `job_cycle` contains a list of job sets.
- Each `job_set` is converted to multiple `jobs` based on the `model_seeds`.
- Each `job` contains a list of entities.
- An entity can be of one of the following types:
    - `proteinChain`
    - `dnaSequence`
    - `rnaSequence`
    - `ligand`
    - `ion`
- Each entity is an instance of :py:class:`AFSequence`.
"""

import os
import json
import random
import warnings
from typing import List, Dict, Any, Final, Tuple

from af_pipeline.tools.misc_tools import chain_id_gen
from af_pipeline.constants.af_constants import (
    PTM, DNA_MOD, RNA_MOD, LIGAND, ION, ENTITY_TYPES, MAX_TEMPLATE_DATE,
    JOB_LIMIT_PER_JSON
)

SEED_MULTIPLIER: Final[int] = 100
random.seed(47)  # for reproducibility in model seed generation

class AlphaFoldServer:
    """Class to handle the creation of AlphaFold server input JSON files"""

    input_dict: Dict[str, List[Dict[str, Any]]]
    """`{cycle_id : job_sets_list}`.\n
    `cycle_id` is a string identifier for the job cycle.\n
    `job_sets_list` is a list of `job_set`."""

    protein_sequences: Dict[str, str] | None
    """`{uniprot_id: sequence}`."""

    nucleic_acid_sequences: Dict[str, str] | None
    """`{nucleic_acid_id: sequence}`."""

    entities_map: Dict[str, str]
    """`{entity_name: identifier}`.
    `identifier` can be `uniprot_id` or `nucleic_acid_id`."""

    def __init__(
        self,
        input_dict: Dict[str, List[Dict[str, Any]]],
        protein_sequences: Dict[str, str] | None = None,
        nucleic_acid_sequences: Dict[str, str] | None = None,
        entities_map: Dict[str, str] = {},
    ):

        self.entities_map = entities_map
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.input_dict = input_dict

    def create_af3_job_cycles(self) -> Dict[str, List[Dict[str, Any]]]:
        """Create job cycles

        Each `job_cycle` is a list of jobs with each job being a dictionary.\n
        See :py:mod:`AFCycle.seed_jobs` or :py:mod:`AFJobSet.create_job_set`
        for the job dictionary format

        Returns:
            `job_cycles (dict)`:`{job_cycle: job_list}`.
        """

        job_cycles = {}
        job_set_names = {}
        af_offsets = {}

        for job_cycle_id, job_sets_list in self.input_dict.items():

            print("Creating job cycle", job_cycle_id, "\n")

            af_cycle = AFCycle(
                job_sets_list=job_sets_list,
                protein_sequences=self.protein_sequences,
                nucleic_acid_sequences=self.nucleic_acid_sequences,
                entities_map=self.entities_map,
            )
            af_cycle.update_cycle()
            job_cycles[job_cycle_id] = af_cycle.job_list
            job_set_names[job_cycle_id] = af_cycle.job_set_names
            af_offsets[job_cycle_id] = af_cycle.job_set_af_offsets

        return job_cycles, job_set_names, af_offsets

    def write_to_json(
        self,
        sets_of_n_jobs: List[List[Dict[str, Any]]],
        file_name: str,
        output_dir: str = "./output/af_input",
    ):
        """Write the sets of n jobs to json files

        Args:
            sets_of_n_jobs (list):
                List of lists, each list containing n jobs
            file_name (str):
                Name of the file
            output_dir (str, optional):
                Directory to save the job_files.
                Defaults to "./output/af_input".
        """

        os.makedirs(output_dir, exist_ok=True)
        for i, job in enumerate(sets_of_n_jobs):

            save_path = os.path.join(output_dir, f"{file_name}_set_{i}.json")

            with open(save_path, "w") as f:
                json.dump(job, f, indent=4)

            print(f"{len(job)} jobs written for {file_name}_set_{i}")

    def write_job_files(
        self,
        job_cycles: Dict[str, List[Dict[str, Any]]],
        output_dir: str = "./output/af_input",
        num_jobs_per_file: int = 20,
    ):
        """Write job files to the output directory

        Args:
            job_cycles (dict): `{job_cycle: job_list}`
            output_dir (str, optional): Directory to save the job_files.
            num_jobs_per_file (int, optional): Number of jobs per file.
        """

        assert JOB_LIMIT_PER_JSON + 1 > num_jobs_per_file > 0; (
            "Number of jobs per file must be within 1 and 100"
        )

        for job_cycle, job_list in job_cycles.items():

            sets_of_n_jobs = [
                job_list[i : i + num_jobs_per_file]
                for i in range(0, len(job_list), num_jobs_per_file)
            ]
            os.makedirs(output_dir, exist_ok=True)

            self.write_to_json(
                sets_of_n_jobs=sets_of_n_jobs,
                file_name=job_cycle,
                output_dir=os.path.join(output_dir, job_cycle),
            )

class AFCycle:

    job_sets_list: List[Dict[str, Any]]
    """List of dictionaries containing job information for each job set.
    Each job set can correspond to multiple jobs based on `model_seeds`.
    All jobs corresponding to a job set are same except for the model seed."""

    protein_sequences: Dict[str, str] | None
    """`{uniprot_id: sequence}`."""

    nucleic_acid_sequences: Dict[str, str] | None
    """`{nucleic_acid_id: sequence}`."""

    entities_map: Dict[str, str]
    """`{entity_name: identifier}`.
    `identifier` can be `uniprot_id` or `nucleic_acid_id`."""

    job_list: List[Dict[str, Any]]
    """List of jobs created from the job information in `job_sets_list`.
    In this, the jobs are in their final form."""

    job_set_names: List[str]
    """List of job set names in the cycle."""

    job_set_af_offsets: List[Dict[str, List[int]]]
    """List of offsets for each chain within each job set in the cycle.
    This is same as range of each entity in the job set with only difference
    being that the chain ids are mapped to the range instead of `entity_name`."""

    def __init__(
        self,
        job_sets_list: List[Dict[str, Any]],
        protein_sequences: Dict[str, str] | None = None,
        nucleic_acid_sequences: Dict[str, str] | None = None,
        entities_map: Dict[str, str] = {},
    ):

        self.job_sets_list = job_sets_list  # all jobs within the cycle
        self.entities_map = entities_map
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.job_list = []
        self.job_set_names = []
        self.job_set_af_offsets = []  # offset for each job set in the cycle

    def update_cycle(self):
        """Update the cycle with the jobs

        For each job in `job_sets_list`, creates an `AFJobSet` instance and uses
        it to create a job dictionary. The job dictionary is then seeded using
        `seed_jobs` to create multiple jobs per set based on `model_seeds`.
        """

        for job_set_info in self.job_sets_list:
            af_job_set = AFJobSet(
                job_set_info=job_set_info,
                protein_sequences=self.protein_sequences,
                nucleic_acid_sequences=self.nucleic_acid_sequences,
                entities_map=self.entities_map,
            )

            job_set_dict = af_job_set.create_job_set()
            self.job_set_names.append(af_job_set.job_set_name)
            self.job_set_af_offsets.append(af_job_set.job_set_af_offset)
            self.seed_jobs(job_set_dict)

    def seed_jobs(
        self,
        job_set_dict: Dict[str, Any],
    ):
        """Create a job for each model seed in the job set

        job dictionary in the following format:
        ```
        {
            "name": "job_set_name",
            "modelSeeds": [1, 2],
            "sequences": [... ]
        }
        ```
        will lead to -->
        ```
        {
            "name": "job_set_name",
            "modelSeeds": [1],
            "sequences": [... ]
        },
        {
            "name": "job_set_name",
            "modelSeeds": [2],
            "sequences": [... ]
        }
        ```
        Args:
            job_set_dict (dict): job dictionary in the format mentioned above
        """

        # this will lead to AF3 server deciding the seed
        if len(job_set_dict["modelSeeds"]) == 0:
            self.job_list.append(job_set_dict)

        else:
            for seed in job_set_dict["modelSeeds"]:
                job_copy = job_set_dict.copy()
                job_copy["modelSeeds"] = [seed]
                job_copy["name"] = f"{job_set_dict['name']}_{seed}"
                self.job_list.append(job_copy)


class AFJobSet:
    """AlphaFold job constructor \n

    A job is a dictionary with the following keys
    - name
    - modelSeeds
    - sequences
    """

    job_set_info: Dict[str, Any]
    """Dictionary containing job attributes such as name, modelSeeds, and entities."""

    protein_sequences: Dict[str, str] | None
    """`{uniprot_id: sequence}`."""

    nucleic_acid_sequences: Dict[str, str] | None
    """`{nucleic_acid_id: sequence}`."""

    entities_map: Dict[str, str]
    """`{entity_name: identifier}`.
    `identifier` can be `uniprot_id` or `nucleic_acid_id`."""

    job_set_name: str | None
    """Name of the job.
    If not provided, it will be generated using `generate_job_set_name`."""

    model_seeds: List[int]
    """List of `model_seeds` for the job.
    If not provided, it will be generated using `generate_seeds`."""

    af_sequences: List[Dict[str, Any]]
    """List of `AFSequence` dictionaries each representing an entity in the job."""

    name_fragments: List[str]
    """List of name fragments derived from the entities in the job."""

    job_set_af_offset: Dict[str, List[int]]
    """Dictionary mapping chain IDs to their start and end positions
    within the job set."""

    def __init__(
        self,
        job_set_info: Dict[str, Any],
        protein_sequences: Dict[str, str] | None = None,
        nucleic_acid_sequences: Dict[str, str] | None = None,
        entities_map: Dict[str, str] = {},
    ):

        self.job_set_info = job_set_info
        self.entities_map = entities_map
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.job_set_name = None
        self.model_seeds = []
        self.af_sequences = []
        self.name_fragments = []
        self.job_set_af_offset = {}

    def create_job_set(self) -> Dict[str, Any]:
        """Create a job from the job info

        Returns:
            `job_set_dict (dict)`:
                dictionary in the following format-
                ```
                {
                    "name": "job_set_name",
                    "modelSeeds": [1, 2],
                    "sequences": [... ]
                }
                ```
        """

        self.update_job_set_name()
        self.update_model_seeds()
        self.update_af_sequences()

        if self.job_set_name is None:
            self.generate_job_set_name()

        job_set_dict = {
            "name": self.job_set_name,
            "modelSeeds": self.model_seeds,
            "sequences": self.af_sequences,
        }

        return job_set_dict

    def update_job_set_name(self):
        """Create a job from the job set info

        Tries to get the job name from the `job_set_info`.
        """

        self.job_set_name = self.job_set_info.get("name")

    def update_model_seeds(self):
        """Update the `model_seeds`

        - If modelSeeds is an integer, generate that many seeds
        - If modelSeeds is a list, use those seeds
        - If modelSeeds is not provided, empty list (auto seed by AF3)
        """

        model_seeds = self.job_set_info.get("modelSeeds")

        if "modelSeeds" in self.job_set_info:

            if isinstance(model_seeds, int):
                self.model_seeds = self.generate_seeds(num_seeds=model_seeds)

            elif isinstance(model_seeds, list):
                self.model_seeds = model_seeds

            else:
                raise Exception("modelSeeds must be an integer or a list")

    def update_af_sequences(self):
        """Update the AF sequences

        - For each entity, create an `AFSequence` instance
        - Get the name fragment for each entity
            (used in job name if `job_name` is not provided)
        """

        chainGen = chain_id_gen()

        # add af_sequence for each entity
        for entity_info in self.job_set_info["entities"]:
            af_sequence = AFSequence(
                entity_info=entity_info,
                protein_sequences=self.protein_sequences,
                nucleic_acid_sequences=self.nucleic_acid_sequences,
                entities_map=self.entities_map,
            )
            af_sequence_dict = af_sequence.create_af_sequence()
            self.af_sequences.append(af_sequence_dict)

            for count in range(af_sequence.count):
                entity_chain_id = next(chainGen)
                self.job_set_af_offset[entity_chain_id] = [
                    af_sequence.start, af_sequence.end
                ]
            self.name_fragments.append(af_sequence.get_name_fragment())

    def generate_job_set_name(self):
        """Generate a job name

        Currently, the job name is generated by concatenating the name fragments
        of each entity in the job, separated by underscores.

        NOTE: Two jobs that only differ by modifications or glycans can not be
        distinguished by this naming scheme.
        """

        job_set_name = "_".join(self.name_fragments)
        self.job_set_name = job_set_name

    def generate_seeds(self, num_seeds: int) -> List[int]:
        """Generate `model_seeds`"""

        model_seeds = random.sample(range(1, SEED_MULTIPLIER * num_seeds), num_seeds)

        return model_seeds


class Entity:
    """Entity constructor in the AlphaFold job \n
    an entity can be a proteinChain, dnaSequence, rnaSequence, ligand or ion \n
    See :py:mod:`AfSequence.create_af_sequence` to check attributes for each entity type
    """

    entity_info: Dict[str, Any]
    """Dictionary containing entity attributes such as name, type, count,
    range, glycans, and modifications."""

    protein_sequences: Dict[str, str] | None
    """`{uniprot_id: sequence}`."""

    nucleic_acid_sequences: Dict[str, str] | None
    """`{nucleic_acid_id: sequence}`."""

    entities_map: Dict[str, str]
    """`{entity_name: identifier}`.
    `identifier` can be `uniprot_id` or `nucleic_acid_id`."""

    entity_name: str
    """Name of the entity."""

    entity_type: str
    """Type of the entity.
    (e.g., proteinChain, dnaSequence, rnaSequence, ligand, ion)."""

    entity_count: int
    """Count or copy number of the entity (default: 1)."""

    real_sequence: str
    """Real sequence of the entity.
    This is the amino acid or nucleic acid sequence."""

    start: int
    """Start position of the entity in the sequence (default: 1)."""

    end: int | None
    """End position of the entity in the sequence"""

    glycans: List[Dict[str, Any]] | None
    """List of glycans associated with the entity, if applicable."""

    modifications: List[Dict[str, Any]] | None
    """List of modifications associated with the entity, if applicable."""

    template_settings: Dict[str, Any]
    """Template settings for the entity, used for proteinChain entities."""

    def __init__(
        self,
        entity_info: Dict[str, Any],
        protein_sequences: Dict[str, str] | None = None,
        nucleic_acid_sequences: Dict[str, str] | None = None,
        entities_map: Dict[str, str] = {},
    ):

        self.entities_map = entities_map
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences

        self.entity_info = entity_info
        self.sanity_check_entity_type(entity_type=entity_info["type"])
        self.entity_type = entity_info["type"]
        self.entity_name = entity_info["name"]

        self.entity_count = 1
        self.real_sequence = ""
        self.start = 1
        self.end = None
        self.glycans = []
        self.modifications = []

        self.fill_up_entity()
        self.sanity_check_glycans()
        self.sanity_check_modifications()
        self.sanity_check_small_molecule(
            entity_type=self.entity_type,
            entity_name=self.entity_name
        )

    def get_template_settings(self):
        """Get the template settings for the entity

        - For proteinChain, get the settings from the entity_info dictionary.
        - For dnaSequence, rnaSequence, ligand and ion, return {}.
        - For proteinChain, `useStructureTemplate` by default is True.
        - For proteinChain, `maxTemplateDate` by default is 2021-09-30.

        Returns:
            `template_dict (dict)`:
                Template settings for the entity as dictionary:
                ```
                {
                    "maxTemplateDate": "YYYY-MM-DD",
                    "useStructureTemplate": True,
                }
                ```
        """

        template_dict = {}

        if self.entity_type != "proteinChain":
            return template_dict

        max_template_date = self.entity_info.get(
            "maxTemplateDate", MAX_TEMPLATE_DATE
        )
        use_structure_template = self.entity_info.get(
            "useStructureTemplate", True
        )

        assert isinstance(use_structure_template, bool), \
            "useStructureTemplate must be a boolean value."

        template_dict_configs = {
            True: {
                "useStructureTemplate": True,
                "maxTemplateDate": max_template_date,
            },
            False: {
                "useStructureTemplate": False,
            },
        }

        template_dict = template_dict_configs[use_structure_template]

        if (
            "maxTemplateDate" in self.entity_info
            and use_structure_template is False
        ):
            warnings.warn(
                f"maxTemplateDate is provided for {self.entity_name} \
                but useStructureTemplate is False. \
                Ignoring maxTemplateDate."
            )

        return template_dict

    def get_entity_count(self)-> int:
        """Get the count of the entity

        Returns:
            `entity_count (int)`:
                Count or copy number of the entity (default: 1)
        """

        entity_count = self.entity_info.get("count", 1)

        return entity_count

    def get_real_sequence(self)-> str:
        """Get the real sequence of the entity

        - For proteinChain, get the sequence from the `protein_sequences`
        - For dnaSequence and rnaSequence, get the sequence from
            the `nucleic_acid_sequences`

        Returns:
            `real_sequence (str)`:
                Amino acid or nucleic acid sequence of the entity
        """

        real_sequence = ""

        if self.entity_type == "proteinChain":

            try: # try with uniprot id as a key
                uniprot_id = self.entities_map[self.entity_name]
                real_sequence = self.protein_sequences[uniprot_id]

            except KeyError:

                try: # try with entity name as a key
                    real_sequence = self.protein_sequences[self.entity_name]

                except KeyError:
                    raise Exception(
                        f"Could not find the entity sequence for {self.entity_name}"
                    )

        elif (
            self.entity_type in ["dnaSequence", "rnaSequence"]
            and self.nucleic_acid_sequences is not None
        ):

            try: # try with nucleic acid id as a key
                nucleic_acid_id = self.entities_map[self.entity_name]
                real_sequence = self.nucleic_acid_sequences[nucleic_acid_id]

            except KeyError:

                try: # try with entity name as a key
                    real_sequence = self.nucleic_acid_sequences[self.entity_name]

                except KeyError:
                    raise Exception(
                        f"Could not find the entity sequence for {self.entity_name}"
                    )

        return real_sequence

    def get_entity_range(self)-> Tuple[int, int]:
        """Get the range of the entity

        what part of the sequence to use? (defined by start and end)
        - If range is provided, use that
        - If no range is provided, use the full sequence
        - If no sequence is found (e.g. ligand or ion), use a range of [1, 1]

        Returns:
            `tuple (start, end)`: start and end of the entity
        """

        if "range" in self.entity_info:

            assert (
                len(self.entity_info["range"]) == 2
            ), "Invalid range; must be a list of two integers (start and end)"

            start, end = self.entity_info["range"]

        else: # use full sequence or [1, 1] for small molecules
            start, end = 1, max(1, len(self.real_sequence))

        return start, end

    def get_glycans(self) -> List[Dict[str, Any]]:
        """Get the glycans of the protein chains

        - If glycans are provided, use those else return []
        - For proteinChain, get the glycans from the entity_info dictionary

        Returns:
            `glycans (list)`:
                List of glycans associated with the entity, if applicable.
        """

        glycans = []

        if self.entity_type == "proteinChain" and "glycans" in self.entity_info:
            glycans = self.entity_info["glycans"]
            glycans = [
                {
                    "residues": glycan[0],
                    "position": glycan[1] - self.start + 1,
                }
                for glycan in glycans
            ]

        return glycans

    def get_modifications(self)-> List[Dict[str, Any]]:
        """Get the modifications of the entity

        - If `modifications` are provided, use those else empty list

            - For proteinChain, get the `modifications` from the
                `entity_info` dictionary (`ptmType`, `ptmPosition`)

            - For dnaSequence and rnaSequence, get the `modifications` from the
                `entity_info` dictionary (`modificationType`, `basePosition`)

        Returns:
            `modifications (list)`:
                List of modifications associated with the entity, if applicable.
        """

        modifications = self.entity_info.get("modifications", [])

        if len(modifications) == 0:
            return modifications

        modification_keys = {
            "proteinChain": ["ptmType", "ptmPosition"],
            "dnaSequence": ["modificationType", "basePosition"],
            "rnaSequence": ["modificationType", "basePosition"],
        }

        if self.entity_type not in modification_keys:
            raise Exception(
                "Modifications are not supported for this entity type"
            )

        modifications = [
            {
                modification_keys[self.entity_type][0]: mod[0],
                modification_keys[self.entity_type][1]: mod[1] - self.start + 1
            }
            for mod in modifications
        ]

        return modifications

    @staticmethod
    def sanity_check_entity_type(entity_type):
        """Sanity check the entity type

        Allowed entity types:
            `proteinChain`, `dnaSequence`, `rnaSequence`, `ligand`, `ion`
        """

        if entity_type not in ENTITY_TYPES:
            raise Exception(f"Invalid entity type {entity_type}")

    @staticmethod
    def sanity_check_small_molecule(entity_type, entity_name):
        """Sanity check the small molecules"""

        if (entity_type == "ligand" and entity_name not in LIGAND) or (
            entity_type == "ion" and entity_name not in ION
        ):
            raise Exception(f"Invalid small molecule {entity_name}")

    def sanity_check_glycans(self):
        """Sanity check the glycans

        - check if the glycosylation position is valid
            (should be within the provided sequence)
        - glycans are only supported for proteinChain, raise exception otherwise
        """

        if self.entity_type == "proteinChain" and len(self.glycans) > 0:

            # check if the glycosylation position is valid
            for glycan in self.glycans:
                glyc_pos = glycan["position"]

                if glyc_pos < 1 or glyc_pos > len(self.real_sequence):
                    raise Exception(
                        f"Invalid glycan position at {glyc_pos} \
                        in {self.entity_name}"
                    )

        if self.entity_type != "proteinChain" and len(self.glycans) > 0:
            raise Exception(
                """

                Glycosylation is not supported for this entity type.
                """
            )

    def sanity_check_modifications(self):
        """Sanity check the modifications

        - check if the modification type is valid
            (should be in the allowed modifications)
        - check if the modification position is valid
            (should be within the provided sequence)
        - modifications are only supported for
            `proteinChain`, `dnaSequence`, `rnaSequence`;
            raise exception otherwise
        """

        if (
            self.entity_type not in [
                "proteinChain", "dnaSequence", "rnaSequence"
            ]
            and len(self.modifications) > 0
        ):
            raise Exception(
                """

                Modifications are not supported for this entity type
                """
            )

        # check if the modification type is valid
        if self.entity_type == "proteinChain":

            if not all([
                mod["ptmType"] in PTM
                for mod in self.modifications
            ]):
                raise Exception("Invalid modification type")

        elif self.entity_type == "dnaSequence":

            if not all([
                mod["modificationType"] in DNA_MOD
                for mod in self.modifications
            ]):
                raise Exception("Invalid modification type")

        elif self.entity_type == "rnaSequence":

            if not all([
                mod["modificationType"] in RNA_MOD
                for mod in self.modifications]
            ):
                raise Exception("Invalid modification type")

        # check if the modification position is valid
        for mod in self.modifications:

            mod_pos = (
                mod["ptmPosition"]
                if self.entity_type == "proteinChain"
                else mod["basePosition"]
            )

            if mod_pos < 1 or mod_pos > len(self.real_sequence):
                raise Exception(
                    f"""

                    Invalid modification at {mod_pos} in {self.entity_name}.
                    """
                )

    def fill_up_entity(self):
        """Fill up the entity with the required information"""

        self.entity_count = self.get_entity_count()
        self.real_sequence = self.get_real_sequence()
        self.start, self.end = self.get_entity_range()
        self.glycans = self.get_glycans()
        self.modifications = self.get_modifications()
        self.template_settings = self.get_template_settings()


class AFSequence(Entity):
    """AlphaFold sequence constructor \n

    This class inherits from `Entity` and is used to create a sequence
    for the AlphaFold job. It initializes the sequence based on the entity
    information provided in the `entity_info` dictionary. \n

    A sequence is an entity ready to be used in the AlphaFold job \n
    'sequences' key in AF job holds a list of sequences \n
    each sequence is a dictionary with following keys:
    - for `proteinChain`:
        1. sequence
        2. glycans
        3. modifications
        4. count
    - for `dnaSequence` or `rnaSequence`:
        1. sequence
        2. modifications
        3. count
    - for `ligand` or `ion`:
        1. ligand or ion identifier
        2. count
    """

    entity_info: Dict[str, Any]
    """Dictionary containing entity attributes such as name, type, count,
    range, glycans, and modifications."""

    protein_sequences: Dict[str, str] | None
    """`{uniprot_id: sequence}`."""

    nucleic_acid_sequences: Dict[str, str] | None
    """`{nucleic_acid_id: sequence}`."""

    entities_map: Dict[str, str]
    """`{entity_name: identifier}`.
    `identifier` can be `uniprot_id` or `nucleic_acid_id`."""

    name: str
    """Name of the entity."""

    type: str
    """Type of the entity.
    (e.g., proteinChain, dnaSequence, rnaSequence, ligand, ion)."""

    count: int
    """Count or copy number of the entity (default: 1)."""

    real_sequence: str
    """Real sequence of the entity.
    This is the amino acid or nucleic acid sequence."""

    glycans: List[Dict[str, Any]] | None
    """List of glycans associated with the entity, if applicable."""

    modifications: List[Dict[str, Any]] | None
    """List of modifications associated with the entity, if applicable."""

    template_settings: Dict[str, Any]
    """Template settings for the entity, used for proteinChain entities."""

    start: int
    """Start position of the entity in the sequence (default: 1)."""

    end: int | None
    """End position of the entity in the sequence"""

    def __init__(
        self,
        entity_info: Dict[str, Any],
        protein_sequences: Dict[str, str] | None = None,
        nucleic_acid_sequences: Dict[str, str] | None = None,
        entities_map: Dict[str, str] = {},
    ):

        super().__init__(
            entity_info=entity_info,
            protein_sequences=protein_sequences,
            nucleic_acid_sequences=nucleic_acid_sequences,
            entities_map=entities_map,
        )
        self.name = self.entity_name
        self.type = self.entity_type
        self.count = self.entity_count
        self.real_sequence = self.update_real_sequence()

    def create_af_sequence(self)-> Dict[str, Any]:
        """Create an AF sequence dictionary

        Returns:
            `af_sequence_dict (dict)`:
            AF sequence dictionary in the following format:
            - for `proteinChain`:
            ```
                {
                    "proteinChain": {
                        "sequence": "AAAA",
                        "glycans": [... ],
                        "modifications": [... ],
                        "count": 1
                    }
                }
            ```
            - for `dnaSequence` or `rnaSequence`:
            ```
                {
                    "dnaSequence"("rnaSequence"): {
                        "sequence": "ATCG",
                        "modifications": [... ],
                        "count": 1
                }
            ```
            - for `ligand` or `ion`:
            ```
                {
                    "ligand"("ion"): {
                        "ligand": "ATP",
                        "count": 1
                    }
                }
            ```
        """

        af_sequence_dict = {}

        if self.type == "proteinChain":
            af_sequence_dict = {
                self.type: {
                    "sequence": self.real_sequence,
                    "glycans": self.glycans,
                    "modifications": self.modifications,
                    "count": self.count,
                }
            }

            af_sequence_dict[self.type].update(self.template_settings)

        elif self.type in ["dnaSequence", "rnaSequence"]:
            af_sequence_dict = {
                self.type: {
                    "sequence": self.real_sequence,
                    "modifications": self.modifications,
                    "count": self.count,
                }
            }

        elif self.type in ["ligand", "ion"]:
            af_sequence_dict = {
                self.type: {self.type: self.name, "count": self.count}
            }

        return af_sequence_dict

    def update_real_sequence(self)-> str:
        """Update the real sequence of the entity

        A real sequence is:
        - amino acid sequence for proteinChain
        - and nucleic acid sequence for dnaSequence and rnaSequence

        If a range [start, end] is provided, slice the sequence accordingly.

        Returns:
            `real_sequence (str)`:
                Amino acid or nucleic acid sequence of the entity
        """

        real_sequence = self.real_sequence
        start, end = self.start, self.end

        if self.type in ["proteinChain", "dnaSequence", "rnaSequence"]:
            real_sequence = real_sequence[start - 1 : end]

        return real_sequence

    def get_name_fragment(self)-> str:
        """Get the name fragments of the entity

        Name format: `"name_count_start-end"`

        Example:

            For an entity `"protA"` (`count==2`, `start==1`, `end==100`),
            the `name_fragment` will be `"protA_2_1-100"`.

        Returns:
            `name_fragment (str)`: Name fragment of the entity
        """

        return f"{self.name}_{self.count}_{self.start}-{self.end}"
