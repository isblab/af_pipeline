"""
[alphafold3](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/af_input/alphafold3.py)
===================================

- Create input `JSON` files for AlphaFold server (https://alphafoldserver.com)

- Keep in mind the following hierarchy:

```mermaid
graph LR
    job_set -->|contains| job
    job -->|contains| entity
```

- An entity can be of one of the following types:
    - `proteinChain`
    - `dnaSequence`
    - `rnaSequence`
    - `ligand`
    - `ion`
- Each entity in the `job` is an instance of :py:class:`af_pipeline.af_input.alphafold3.AFSequence`.

"""

import os
import json
import warnings
from typing import List, Dict, Any, Tuple

from af_pipeline.utils.file_utils import write_json
from af_pipeline.utils.misc_utils import (
    chain_id_gen,
    generate_seeds
)
from af_pipeline.constants.af_constants import (
    PTM, DNA_MOD, RNA_MOD, LIGAND, ION, ENTITY_TYPES, ALLOWED_MSA_TYPES,
    MAX_TEMPLATE_DATE, JOB_LIMIT_PER_JSON, AF_JOB_FILE,
    ConfigYaml, MiscStrEnum, NucleicAcidStrand, MSAType, MSAFields,
)
from af_pipeline.af_input.msa import parse_a3m, MMSeqs2
from af_pipeline.constants import af_constants
from af_pipeline.constants.af_constants import (
    NucleicAcidModificationFields,
    ProteinModificationFields,
    GlycanModificationFields,
    AFServerSequenceFields,
    AFInputEntityFields,
    AFInputJobFields,
    EntityType,
    FileFormat,
)

class AFServerConfig:
    """Class to help create JSON files for AlphaFold server jobs."""

    input_dict: Dict[str, List[Dict[str, Any]]]
    """List of job sets, each of which specifies the entities, model seeds, job name, etc."""

    protein_sequences: Dict[str, str] | None
    """Dictionary with:<br />

    - `key` -> `fasta_header` <br />
       Usually `uniprot_id` in case of `proteinChain` entities.<br />
       If `fasta_header != entity_name`, `entities_map` should be provided.<br />

    - `val` -> `sequence` <br />
      Amino acid sequence of the protein chain.
    """

    nucleic_acid_sequences: Dict[str, str] | None
    """Dictionary with:<br />

    - `key` -> `fasta_header` <br />
      If `fasta_header != entity_name`, `entities_map` should be provided.<br />

    - `val` -> `sequence` <br />
      Nucleotide sequence of the nucleic acid.
    """

    entities_map: Dict[str, str]
    """Dictionary with:<br />

    - `key` -> `entity_name` <br />

    - `val` -> `fasta_header` <br />
      `fasta_header` is usually `uniprot_id` in case of `proteinChain` entities."""

    set_seed: int
    """Seed used to generate `model_seeds` if not provided."""

    def __init__(
        self,
        config_dict: Dict[str, Any],
        protein_sequences: Dict[str, str] | None = None,
        nucleic_acid_sequences: Dict[str, str] | None = None,
        set_seed: int = 47,
    ):

        self.config_dict = config_dict
        self.entities_map = config_dict.get(ConfigYaml.PROTEIN_UNIPROT_MAP, {}).copy()
        self.input_job_sets = config_dict.get(ConfigYaml.AF_INPUT_JOBS, []).copy()
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.set_seed = set_seed

    def write_job_files(
        self,
        output_dir: str,
        num_jobs_per_file: int = 20,
        indent: int = 4,
    ) -> None:
        """Convert the input information into the dictionary format required by
        the AlphaFold server and write job files to the output directory.

        For each set in `input_job_sets`, split the jobs into `sets_of_n_jobs`
        (depending on `num_jobs_per_file`) and write `JSON` files in
        set-specific directory.

        The files are stored as:\n
            {output_dir}/{job_set_id}/{job_set_id}_set_{i}.json
        where, `i` denotes the index.

        There is a upper limit of 100 jobs per `JSON` file imposed by
        AlphaFold server.

        <a href="https://github.com/google-deepmind/alphafold/blob/main/server/example.json">
        Check example JSON file.</a>

        Arguments:

        - **output_dir (str)**:<br />
            Directory to save the `JSON` files.

        - **num_jobs_per_file (int, optional)**:<br />
            Number of jobs per file.
        """

        assert JOB_LIMIT_PER_JSON + 1 > num_jobs_per_file > 0; (
            "Number of jobs per file must be within 1 and 100"
        )

        for job_set_idx, job_set_info in enumerate(self.input_job_sets):
            job_set = AFJobSet(
                job_set_info=job_set_info,
                protein_sequences=self.protein_sequences,
                nucleic_acid_sequences=self.nucleic_acid_sequences,
                entities_map=self.entities_map,
                set_seed=self.set_seed,
            )
            job_set_dict = job_set.create_job_set()
            job_list = self.seed_jobs(job_set_dict=job_set_dict)

            sets_of_n_jobs = [
                job_list[i:i + num_jobs_per_file]
                for i in range(0, len(job_list), num_jobs_per_file)
            ]
            os.makedirs(output_dir, exist_ok=True)

            AFServerConfig.write_to_json(
                sets_of_n_jobs=sets_of_n_jobs,
                file_name=job_set.job_set_name,
                output_dir=os.path.join(output_dir, job_set.job_set_name),
            )

            self.input_job_sets[job_set_idx] = {
                AFInputJobFields.JOB_SET_NAME: job_set.job_set_name,
                AFInputJobFields.MODEL_SEEDS: job_set.model_seeds,
                **{
                    k:v for k,v in self.input_job_sets[job_set_idx].items()
                    if k not in [
                        AFInputJobFields.JOB_SET_NAME,
                        AFInputJobFields.MODEL_SEEDS,
                        AFInputJobFields.AF_OFFSET,
                    ]
                },
                AFInputJobFields.AF_OFFSET: job_set.job_set_af_offset,
            }

            write_json(
                file_path=os.path.join(os.path.dirname(output_dir), "af_input_jobs.json"),
                data={
                    ConfigYaml.PROTEIN_UNIPROT_MAP: self.entities_map,
                    ConfigYaml.AF_INPUT_JOBS: self.input_job_sets,
                },
                indent=indent,
            )

    @staticmethod
    def write_to_json(
        sets_of_n_jobs: List[List[Dict[str, Any]]],
        file_name: str,
        output_dir: str,
    ):
        """Write the sets of "n" jobs to `JSON` files.

        Arguments:

        - **sets_of_n_jobs (list)**:<br />
            List of lists, each list containing n jobs.

        - **file_name (str)**:<br />
            Name of the file.

        - **output_dir (str, optional)**:<br />
            Directory to save the `JSON` files.
        """

        os.makedirs(output_dir, exist_ok=True)
        for i, job in enumerate(sets_of_n_jobs):

            f_name = AF_JOB_FILE.substitute(fname=file_name, set_idx=i)

            save_path = os.path.join(
                output_dir,
                f"{f_name}.{FileFormat.JSON}",
            )

            with open(save_path, "w") as f:
                json.dump(job, f, indent=4)

            print(f"{len(job)} jobs written for {f_name}")

    def seed_jobs(
        self,
        job_set_dict: Dict[str, Any],
    ):
        """Create a job for each model seed in the job set

        Job set dictionary in the following format:
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

        Arguments:

        - **job_set_dict (dict)**:<br />
            Job set dictionary in the format mentioned above.
        """

        job_list = []
        # this will lead to AF3 server deciding the seed
        if len(job_set_dict[AFInputJobFields.MODEL_SEEDS]) == 0:
            job_list.append(job_set_dict)

        else:
            for seed in job_set_dict[AFInputJobFields.MODEL_SEEDS]:
                job_copy = job_set_dict.copy()
                job_copy[AFInputJobFields.MODEL_SEEDS] = [seed]
                job_copy[AFInputJobFields.NAME] = f"{job_set_dict[AFInputJobFields.NAME]}_{seed}"
                job_list.append(job_copy)

        return job_list

class AFJobSet:
    """AlphaFold job set constructor"""

    job_set_info: Dict[str, Any]
    """Dictionary containing job attributes such as name, modelSeeds, and entities."""

    protein_sequences: Dict[str, str] | None
    """Dictionary with:<br />

    - `key` -> `identifier` <br />
       Usually `uniprot_id` in case of `proteinChain` entities.<br />
       `identifier != entity_name` necessitates `entities_map`.<br />

    - `val` -> `sequence` <br />
      Amino acid sequence of the protein chain.
    """

    nucleic_acid_sequences: Dict[str, str] | None
    """Dictionary with:<br />

    - `key` -> `identifier` <br />
      `identifier != entity_name` necessitates `entities_map`.<br />

    - `val` -> `sequence` <br />
      Nucleotide sequence of the nucleic acid.
    """

    entities_map: Dict[str, str]
    """Dictionary with:<br />

    - `key` -> `entity_name` <br />

    - `val` -> `identifier` <br />
      `identifier` is usually `uniprot_id` in case of `proteinChain` entities."""

    set_seed: int
    """Seed used to generate `model_seeds` if not provided."""

    job_set_name: str | None
    """Name of the job set.
    If not provided, it will be generated using `generate_job_set_name`."""

    model_seeds: List[int]
    """List of `model_seeds` for the job.
    If not provided, it will be generated using `generate_seeds`."""

    af_sequences: List[Dict[str, Any]]
    """List of dictionaries each corresponding to an entity in the job.<br />
    These are created using :py:meth:`AFSequence.create_af_sequence`."""

    name_fragments: List[str]
    """List of name fragments derived from the entities in the job.<br />
    Each fragment is generated using :py:meth:`AFSequence.get_name_fragment`."""

    job_set_af_offset: Dict[str, List[int]]
    """Dictionary mapping chain ids to their start and end positions as defined
    in the job set."""

    def __init__(
        self,
        job_set_info: Dict[str, Any],
        protein_sequences: Dict[str, str] | None = None,
        nucleic_acid_sequences: Dict[str, str] | None = None,
        entities_map: Dict[str, str] = {},
        set_seed: int = 47,
    ):

        self.job_set_info = job_set_info
        self.msa_type = job_set_info.get(AFInputJobFields.MSA, {}).get(MSAFields.TYPE)
        self.msa_path = job_set_info.get(AFInputJobFields.MSA, {}).get(MSAFields.PATH)
        self.entities_map = entities_map
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.set_seed = set_seed
        self.job_set_name = None
        self.model_seeds = []
        self.af_sequences = []
        self.name_fragments = []
        self.job_set_af_offset = {}
        self.sanity_check_msa_fields()
        self.a3m_lines = self.prepare_msa(use_env=True, use_filter=True)

    def prepare_msa(
        self,
        use_env: bool = False,
        use_filter: bool = False,
        overwrite: bool = False,
    ) -> Dict[str, str]:
        """ Prepare MSA information for the given AF3 job.
        Also see :py:class:`MMSeqs2` for more details on MSA generation.

        ## Arguments:

        - **use_env (bool, optional):**:<br />
            Whether to use sequences from BFD in the MSA.
            Defaults to False.

        - **use_filter (bool, optional):**:<br />
            Whether to apply filtering to the MSA. Defaults to False.

        - **overwrite (bool, optional):**:<br />
            Whether to overwrite existing MSA results. Defaults to False.

        ## Returns:

        - **dict**:<br />
            Dictionary with:<br />
            - `header`: `a3m_content` <br />
        """

        a3m_lines = {}

        if isinstance(self.msa_path, str) is False:
            return a3m_lines

        use_pairing = self.msa_type == MSAType.PAIRED

        if os.path.exists(self.msa_path) is False:

            mmseqs2 = MMSeqs2(
                sequences=self.prepare_msa_sequences(),
                targz_file=f"{os.path.dirname(self.msa_path)}.tar.gz",
                use_env=use_env,
                use_filter=use_filter,
                use_pairing=use_pairing,
            )
            mmseqs2.main(overwrite=overwrite)

            resultdir = os.path.dirname(self.msa_path)

            if use_pairing is False:
                a3m_lines = parse_a3m(
                    a3m_file=os.path.join(resultdir, "uniref.a3m")
                )
                if use_env:
                    env_a3m = parse_a3m(
                        a3m_file=os.path.join(
                            resultdir, "bfd.mgnify30.metaeuk30.smag30.a3m"
                        )
                    )
                    a3m_lines.update(env_a3m)
            else:
                a3m_lines = parse_a3m(
                    a3m_file=os.path.join(resultdir, "pair.a3m")
                )

        else:

            a3m_lines = parse_a3m(a3m_file=self.msa_path)

        return a3m_lines

    def prepare_msa_sequences(self) -> Dict[str, str]:
        """ Prepare input for MMSeqs2 API from the entity sequences.

        ## Returns:

        - **dict**:<br />
            Dictionary with:<br />
            - `header`: `sequence`
        """

        sequences = {}

        for entity_info in self.job_set_info[AFInputJobFields.ENTITIES]:

            entity = Entity(
                entity_info=entity_info,
                protein_sequences=self.protein_sequences,
                nucleic_acid_sequences=self.nucleic_acid_sequences,
                entities_map=self.entities_map,
            )

            if entity.entity_type == EntityType.PROTEIN_CHAIN:

                header = (
                    entity.entity_name + "_" +
                    str(entity.start) + af_constants.RES_RANGE_SEP + str(entity.end)
                )

                sequences[header] = (
                    entity.real_sequence[entity.start - 1 : entity.end]
                )


        return sequences

    def create_job_set(self) -> Dict[str, Any]:
        """Create a job from the job info

        Returns:
        - **job_set_dict (dict)**:<br />
            Dictionary with following `key`:`val` pair:<br />
            - `name`:`job_set_name` <br />
            - `modelSeeds`:`model_seeds` <br />
            - `sequences`:`af_sequences` <br />
        """

        self.update_model_seeds()
        self.update_nucleic_acid_entities()
        self.update_af_sequences()
        self.update_job_set_name()

        job_set_dict = {
            AFInputJobFields.NAME: self.job_set_name,
            AFInputJobFields.MODEL_SEEDS: self.model_seeds,
            AFInputJobFields.SEQUENCES: self.af_sequences,
        }

        return job_set_dict

    def update_job_set_name(self):
        """Update the `job_set_name`.

        Tries to get the job name from the `job_set_info`.
        """

        self.job_set_name = self.job_set_info.get(AFInputJobFields.JOB_SET_NAME)
        if self.job_set_name is None:
            self.generate_job_set_name()

    def update_nucleic_acid_entities(self):
        """Update the nucleic acid entities in the job set.

        For each nucleic acid entity, create a new entity for its reverse complement
        and add it to the job set.
        """

        entities_to_add = []

        for idx, entity_info in enumerate(self.job_set_info[AFInputJobFields.ENTITIES]):

            if entity_info[AFInputEntityFields.TYPE] in [
                EntityType.DNA_SEQUENCE,
                EntityType.RNA_SEQUENCE,
            ]:

                strand = entity_info.get(
                    AFInputEntityFields.STRAND, NucleicAcidStrand.SINGLE
                )

                if strand == NucleicAcidStrand.SINGLE:
                    continue

                reverse_entity_info = entity_info.copy()
                reverse_entity_info[AFInputEntityFields.NAME] = (
                    f"{entity_info[AFInputEntityFields.NAME]}{MiscStrEnum.REVERSE_COMPLEMENT}"
                )
                reverse_entity_info[AFInputEntityFields.STRAND] = NucleicAcidStrand.SINGLE

                # switch the entity info strand to single
                self.job_set_info[AFInputJobFields.ENTITIES][idx][AFInputEntityFields.STRAND] = NucleicAcidStrand.SINGLE

                entities_to_add.append([idx, reverse_entity_info])

        # reverse sort entities to avoid index shift when adding new entities
        entities_to_add.sort(key=lambda x: x[0], reverse=True)
        for idx, entity_info in entities_to_add:
            self.job_set_info[AFInputJobFields.ENTITIES].insert(idx + 1, entity_info)

    def update_model_seeds(self):
        """Update the `model_seeds`.

        - If "modelSeeds" is an integer, generate that many seeds.
        - If "modelSeeds" is a list, use those seeds.
        - If "modelSeeds" is not provided, return an empty list (auto seed by AF3).
        """

        model_seeds = self.job_set_info.get(AFInputJobFields.MODEL_SEEDS)

        if AFInputJobFields.MODEL_SEEDS in self.job_set_info:

            if isinstance(model_seeds, int):
                self.model_seeds = generate_seeds(
                    num_seeds=model_seeds,
                    set_seed=self.set_seed,
                )

            elif isinstance(model_seeds, list):
                self.model_seeds = model_seeds

            else:
                raise Exception("modelSeeds must be an integer or a list")

    def update_af_sequences(self):
        """Update the AF sequences.

        - For each entity, create an `AFSequence` instance
        - Get the name fragment for each entity
            (used in job name if `job_name` is not provided)
        """

        chainGen = chain_id_gen()

        # add af_sequence for each entity
        for entity_info in self.job_set_info[AFInputJobFields.ENTITIES]:

            af_sequence = AFSequence(
                entity_info=entity_info,
                protein_sequences=self.protein_sequences,
                nucleic_acid_sequences=self.nucleic_acid_sequences,
                entities_map=self.entities_map,
                a3m_lines=self.a3m_lines,
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
        """Generate `job_set_name`.

        Currently, the `job_set_name` is generated by concatenating the name fragments
        of each entity in the job, separated by underscores.

        > [!WARNING]
        > Jobs that only differ by modifications or glycans can not be
        distinguished by this naming scheme.
        """

        job_set_name = "_".join(self.name_fragments)
        self.job_set_name = job_set_name

        if len(self.job_set_name) > 100:
            warnings.warn(
                f"Job name {self.job_set_name} is too long (>100 characters). \
                Consider providing a custom job name for job set. \n \
                {self.job_set_info}"
            )

    def sanity_check_msa_fields(self):

        if self.msa_type is not None:
            assert self.msa_type in ALLOWED_MSA_TYPES, (
                f"MSA type {self.msa_type} is not allowed."
                f" Allowed types: {ALLOWED_MSA_TYPES}"
            )

        if self.msa_type in ALLOWED_MSA_TYPES:
            assert self.msa_path is not None, (
                f"MSA type {self.msa_type} requires a path to the a3m MSA file."
            )


class Entity:
    """Entity constructor in the AlphaFold job

    An entity can be a `proteinChain`, `dnaSequence`, `rnaSequence`, `ligand` or `ion` \n
    See :py:mod:`AFSequence.create_af_sequence` to check attributes for each entity type
    """

    entity_info: Dict[str, Any]
    """Dictionary containing entity attributes such as name, type, count,
    range, glycans, and modifications."""

    protein_sequences: Dict[str, str] | None
    """Dictionary with:<br />

    - `key` -> `identifier` <br />
       Usually `uniprot_id` in case of `proteinChain` entities.<br />
       `identifier != entity_name` necessitates `entities_map`.<br />

    - `val` -> `sequence` <br />
      Amino acid sequence of the protein chain.
    """

    nucleic_acid_sequences: Dict[str, str] | None
    """Dictionary with:<br />

    - `key` -> `identifier` <br />
      `identifier != entity_name` necessitates `entities_map`.<br />

    - `val` -> `sequence` <br />
      Nucleotide sequence of the nucleic acid.
    """

    entities_map: Dict[str, str]
    """Dictionary with:<br />

    - `key` -> `entity_name` <br />

    - `val` -> `identifier` <br />
      `identifier` is usually `uniprot_id` in case of `proteinChain` entities."""

    entity_name: str
    """Name of the entity."""

    entity_type: str
    """Type of the entity. It can be one of the following:<br />
    `proteinChain`, `dnaSequence`, `rnaSequence`, `ligand`, `ion`."""

    entity_count: int
    """Count or copy number of the entity."""

    real_sequence: str
    """Real sequence of the entity.
    Refers to the amino acid or nucleic acid sequence."""

    start: int
    """Start position of the entity in the sequence."""

    end: int | None
    """End position of the entity in the sequence."""

    glycans: List[Dict[str, Any]] | None
    """List of glycans associated with the entity, if applicable."""

    modifications: List[Dict[str, Any]] | None
    """List of modifications associated with the entity, if applicable."""

    template_settings: Dict[str, Any]
    """Template settings for the entity, used for `proteinChain` entities."""

    def __init__(
        self,
        entity_info: Dict[str, Any],
        protein_sequences: Dict[str, str] | None = None,
        nucleic_acid_sequences: Dict[str, str] | None = None,
        entities_map: Dict[str, str] = {},
        a3m_lines: Dict[str, str] = {},
    ):

        self.entities_map = entities_map
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.a3m_lines = a3m_lines

        self.entity_info = entity_info
        self.sanity_check_entity_type(
            entity_type=entity_info[AFInputEntityFields.TYPE]
        )
        self.entity_type = entity_info[AFInputEntityFields.TYPE]
        self.entity_name = entity_info[AFInputEntityFields.NAME]
        self.sanity_check_entity_name()

        self.entity_count = 1
        self.real_sequence = ""
        self.start = 1
        self.end = None
        self.glycans = []
        self.modifications = []

        self.update_entity()
        self.sanity_check_glycans()
        self.sanity_check_modifications()
        self.sanity_check_nucleic_acid_strand()
        self.sanity_check_small_molecule(
            entity_type=self.entity_type,
            entity_name=self.entity_name
        )

    def get_msa_settings(self) -> Dict[str, str]:
        """ Get MSA settings for the entity

        - For proteinChain, get the unpaired MSA from the `a3m_lines` dictionary.
        - For dnaSequence, rnaSequence, ligand and ion, return {}.

        > [!NOTE]
        > Even if the msa_type is `pairedMsa`, the `unpairedMsa` field is filled.
        > This is the recommended approach.
        >
        > See https://alphafoldserver.com/faq#what-structure-templates-and-msa-are-used-by-alphafold-server-can-i-customize-these
        > See https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md#msa-pairing

        ## Returns:

        - **dict**:<br />
            Dictionary with:<br />
            - `pairedMsa`:`""` <br />
            - `unpairedMsa`:`a3m_content` <br />
        """

        if self.entity_type != EntityType.PROTEIN_CHAIN:
            return {}

        header = (
            self.entity_name + "_" +
            str(self.start) + af_constants.RES_RANGE_SEP + str(self.end)
        )

        unpaired_msa = "".join(self.a3m_lines.get(header, []))
        if len(unpaired_msa) == 0:
            return {}

        return {
            AFInputEntityFields.PAIRED_MSA: "",
            AFInputEntityFields.UNPAIRED_MSA: unpaired_msa,
        }

    def get_template_settings(self):
        """Get the template settings for the entity

        - For proteinChain, get the settings from the entity_info dictionary.
        - For dnaSequence, rnaSequence, ligand and ion, return {}.
        - For proteinChain, `useStructureTemplate` by default is True.
        - For proteinChain, `maxTemplateDate` by default is 2021-09-30.

        Returns:
        - **template_dict (dict)**:<br />
            Dictionary with following `key`:`val` pairs:<br />
            - `maxTemplateDate`:`YYYY-MM-DD` <br />
            - `useStructureTemplate`:`True` or `False` <br />
        """

        template_dict = {}

        if self.entity_type != EntityType.PROTEIN_CHAIN:
            return template_dict

        max_template_date = self.entity_info.get(
            AFInputEntityFields.MAX_TEMPLATE_DATE, MAX_TEMPLATE_DATE
        )
        use_structure_template = self.entity_info.get(
            AFInputEntityFields.USE_STRUCTURE_TEMPLATE, True
        )

        assert isinstance(use_structure_template, bool), \
            "useStructureTemplate must be a boolean value."

        template_dict_configs = {
            True: {
                AFInputEntityFields.USE_STRUCTURE_TEMPLATE: True,
                AFInputEntityFields.MAX_TEMPLATE_DATE: max_template_date,
            },
            False: {
                AFInputEntityFields.USE_STRUCTURE_TEMPLATE: False,
            },
        }

        template_dict = template_dict_configs[use_structure_template]

        if (
            AFInputEntityFields.MAX_TEMPLATE_DATE in self.entity_info
            and use_structure_template is False
        ):
            warnings.warn(
                f"maxTemplateDate is provided for {self.entity_name} \
                but useStructureTemplate is False. \
                Ignoring maxTemplateDate."
            )

        return template_dict

    def get_entity_count(self)-> int:
        """Get the count of the entity.

        Returns:
        - **entity_count (int)**:<br />
            Count or copy number of the entity (default: 1).
        """

        entity_count = self.entity_info.get(AFInputEntityFields.COUNT, 1)

        return entity_count

    def get_real_sequence(self)-> str:
        """Get the real sequence of the entity.

        - For `proteinChain`, get the sequence from the `protein_sequences`.
        - For `dnaSequence` and `rnaSequence`, get the sequence from
            the `nucleic_acid_sequences`.

        Returns:
        - **real_sequence (str)**:<br />
            Amino acid or nucleic acid sequence of the entity.
        """

        real_sequence = ""

        if (
            self.entity_type == EntityType.PROTEIN_CHAIN
            and self.protein_sequences is not None
        ):

            try: # try with uniprot id as a key
                uniprot_id = self.entities_map[self.entity_name]
                real_sequence = self.protein_sequences[uniprot_id]

            except KeyError:

                try: # try with entity name as a key
                    real_sequence = self.protein_sequences[self.entity_name]

                except KeyError:
                    raise Exception(
                        f"Could not find the entity sequence for {self.entity_name}."
                    )

        elif (
            self.entity_type in [EntityType.DNA_SEQUENCE, EntityType.RNA_SEQUENCE]
            and self.nucleic_acid_sequences is not None
        ):

            entity_name = self.entity_name.replace(MiscStrEnum.REVERSE_COMPLEMENT, "")

            try: # try with nucleic acid id as a key
                nucleic_acid_id = self.entities_map[entity_name]
                real_sequence = self.nucleic_acid_sequences[nucleic_acid_id]

            except KeyError:

                try: # try with entity name as a key
                    real_sequence = self.nucleic_acid_sequences[entity_name]

                except KeyError:
                    raise Exception(
                        f"Could not find the entity sequence for {self.entity_name}."
                    )

            if MiscStrEnum.REVERSE_COMPLEMENT in self.entity_name:
                from Bio.Seq import Seq
                if self.entity_type == EntityType.DNA_SEQUENCE:
                    reverse_comp = str(Seq(real_sequence).reverse_complement())
                elif self.entity_type == EntityType.RNA_SEQUENCE:
                    reverse_comp = str(Seq(real_sequence).reverse_complement_rna())
                # NOTE: the reverse_comp is already in the 5' to 3' direction,
                # so no need to reverse it again
                real_sequence = reverse_comp

        ivalid_case = (
            self.entity_type in [
                EntityType.PROTEIN_CHAIN,
                EntityType.DNA_SEQUENCE,
                EntityType.RNA_SEQUENCE
            ]
            and real_sequence == ""
        )

        if ivalid_case:
            raise Exception(
                f"Could not find the entity sequence for {self.entity_name}."
            )

        return real_sequence

    def get_entity_range(self)-> Tuple[int, int]:
        """Get the range of the entity.

        Region of the sequence to be used in the job? (defined by start and end)
        - If range is provided, slice the full sequence.
        - If no range is provided, use the full sequence.
        - If no sequence is found (e.g. ligand or ion), use `(1, 1)`.

        Returns:
        - **(tuple)**:<br />
            `start` and `end` of the entity.
        """

        if AFInputEntityFields.RANGE in self.entity_info:

            assert (
                len(self.entity_info[AFInputEntityFields.RANGE]) == 2
            ), "Invalid range; must be a list of two integers (start and end)"

            start, end = self.entity_info[AFInputEntityFields.RANGE]

        else: # use full sequence or (1, 1) for small molecules
            start, end = 1, max(1, len(self.real_sequence))

        return start, end

    def get_glycans(self) -> List[Dict[str, Any]]:
        """Get the glycans for `proteinChain` entity.

        Get the glycans from the `entity_info`.

        Returns:
        - **glycans (list)**:<br />
            List of glycans associated with the entity.
        """

        glycans = []

        if (
            self.entity_type == EntityType.PROTEIN_CHAIN and
            AFInputEntityFields.GLYCANS in self.entity_info
        ):
            glycans = self.entity_info[AFInputEntityFields.GLYCANS]
            glycans = [
                {
                    GlycanModificationFields.RESIDUES: glycan[0],
                    GlycanModificationFields.POSITION: glycan[1] - self.start + 1,
                }
                for glycan in glycans
            ]

        return glycans

    def get_modifications(self)-> List[Dict[str, Any]]:
        """Get the modifications of the entity.

        - For proteinChain, get the `modifications` from the `entity_info`<br />
          (`ptmType`, `ptmPosition`)

        - For dnaSequence and rnaSequence, get the `modifications` from the `entity_info`<br />
          (`modificationType`, `basePosition`)

        Returns:
        - **modifications (list)**:<br />
            List of modifications associated with the entity.
        """

        modifications = self.entity_info.get("modifications", [])

        if len(modifications) == 0:
            return modifications

        modification_keys = {
            EntityType.PROTEIN_CHAIN: [
                ProteinModificationFields.PTM_TYPE,
                ProteinModificationFields.PTM_POSITION
            ],
            EntityType.DNA_SEQUENCE: [
                NucleicAcidModificationFields.MODIFICATION_TYPE,
                NucleicAcidModificationFields.BASE_POSITION
            ],
            EntityType.RNA_SEQUENCE: [
                NucleicAcidModificationFields.MODIFICATION_TYPE,
                NucleicAcidModificationFields.BASE_POSITION
            ],
        }

        if self.entity_type not in modification_keys:
            raise Exception(
                "Modifications are not supported for this entity type."
            )

        modifications = [
            {
                modification_keys[self.entity_type][0]: mod[0],
                modification_keys[self.entity_type][1]: mod[1] - self.start + 1
            }
            for mod in modifications
        ]

        return modifications

    def sanity_check_entity_name(self):

        if (self.entity_type in [
            EntityType.PROTEIN_CHAIN,
            EntityType.DNA_SEQUENCE,
            EntityType.RNA_SEQUENCE,
        ]):

            assert "_" not in self.entity_name, (
                f"""Underscore is not allowed in entity name for {self.entity_type}
                entities. Found in {self.entity_name}."""
            )

    @staticmethod
    def sanity_check_entity_type(entity_type):
        """Sanity check the entity type.

        Allowed entity types:
            `proteinChain`, `dnaSequence`, `rnaSequence`, `ligand`, `ion`
        """

        if entity_type not in ENTITY_TYPES:
            raise Exception(f"Invalid entity type {entity_type}")

    @staticmethod
    def sanity_check_small_molecule(entity_type, entity_name):
        """Sanity check the small molecules."""

        if (entity_type == EntityType.LIGAND and entity_name not in LIGAND) or (
            entity_type == EntityType.ION and entity_name not in ION
        ):
            raise Exception(f"Invalid small molecule {entity_name}.")

    def sanity_check_glycans(self):
        """Sanity check the `glycans`.

        - Check if the glycosylation position is valid (i.e. within the
          provided sequence) or not.
        - `glycans` are only supported for `proteinChain`, raise exception otherwise
        """

        if self.entity_type == EntityType.PROTEIN_CHAIN and len(self.glycans) > 0:

            # check if the glycosylation position is valid
            for glycan in self.glycans:
                glyc_pos = glycan[GlycanModificationFields.POSITION]

                if glyc_pos < 1 or glyc_pos > len(self.real_sequence):
                    raise Exception(
                        f"Invalid glycan position at {glyc_pos} \
                        in {self.entity_name}"
                    )

        if self.entity_type != EntityType.PROTEIN_CHAIN and len(self.glycans) > 0:
            raise Exception(
                """

                Glycosylation is not supported for this entity type.
                """
            )

    def sanity_check_nucleic_acid_strand(self):
        """Sanity check the `strand` for nucleic acid entities.

        Allowed strands:
            `single`, `double`

        - For double-stranded nucleic acids, modifications are currently not
          supported. Provide the modifications for each strand separately.
        """

        if self.entity_type in [EntityType.DNA_SEQUENCE, EntityType.RNA_SEQUENCE]:

            strand = self.entity_info.get(
                AFInputEntityFields.STRAND, NucleicAcidStrand.SINGLE
            )

            if strand not in [NucleicAcidStrand.SINGLE, NucleicAcidStrand.DOUBLE]:
                raise Exception(f"Invalid strand {strand} for {self.entity_name}.")

            if len(self.modifications) > 0 and strand == NucleicAcidStrand.DOUBLE:
                raise Exception(
                    f"Modifications are not supported for double-stranded nucleic acids. \
                    Please provide the modifications for each strand separately."
                )

    def sanity_check_modifications(self):
        """Sanity check the `modifications`.

        - check if the modification type is valid or not.
            (should be in the allowed modifications)
        - check if the modification position is valid or not.
            (should be within the provided sequence)
        - `modifications` are only supported for
            `proteinChain`, `dnaSequence`, `rnaSequence`;
            raise exception otherwise.
        """

        if (
            self.entity_type not in [
                EntityType.PROTEIN_CHAIN,
                EntityType.DNA_SEQUENCE,
                EntityType.RNA_SEQUENCE,
            ]
            and len(self.modifications) > 0
        ):
            raise Exception(
                """

                Modifications are not supported for this entity type.
                """
            )

        # check if the modification type is valid
        if self.entity_type == EntityType.PROTEIN_CHAIN:

            if not all([
                mod[ProteinModificationFields.PTM_TYPE] in PTM
                for mod in self.modifications
            ]):
                raise Exception("Invalid modification type.")

        elif self.entity_type == EntityType.DNA_SEQUENCE:

            if not all([
                mod[NucleicAcidModificationFields.MODIFICATION_TYPE] in DNA_MOD
                for mod in self.modifications
            ]):
                raise Exception("Invalid modification type.")

        elif self.entity_type == EntityType.RNA_SEQUENCE:

            if not all([
                mod[NucleicAcidModificationFields.MODIFICATION_TYPE] in RNA_MOD
                for mod in self.modifications]
            ):
                raise Exception("Invalid modification type.")

        # check if the modification position is valid
        for mod in self.modifications:

            mod_pos = (
                mod[ProteinModificationFields.PTM_POSITION]
                if self.entity_type == EntityType.PROTEIN_CHAIN
                else mod[NucleicAcidModificationFields.BASE_POSITION]
            )

            if mod_pos < 1 or mod_pos > len(self.real_sequence):
                raise Exception(
                    f"""

                    Invalid modification at {mod_pos} in {self.entity_name}.
                    """
                )

    def update_entity(self):
        """Fill up the entity with the information.

        This method updates the following attributes of the entity:
        - `entity_count`
        - `real_sequence`
        - `start`
        - `end`
        - `glycans`
        - `modifications`
        - `template_settings`
        """

        self.entity_count = self.get_entity_count()
        self.real_sequence = self.get_real_sequence()
        self.start, self.end = self.get_entity_range()
        self.glycans = self.get_glycans()
        self.modifications = self.get_modifications()
        self.template_settings = self.get_template_settings()
        self.msa_settings = self.get_msa_settings()


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
    """Dictionary with:<br />

    - `key` -> `identifier` <br />
       Usually `uniprot_id` in case of `proteinChain` entities.<br />
       `identifier != entity_name` necessitates `entities_map`.<br />

    - `val` -> `sequence` <br />
      Amino acid sequence of the protein chain.
    """

    nucleic_acid_sequences: Dict[str, str] | None
    """Dictionary with:<br />

    - `key` -> `identifier` <br />
      `identifier != entity_name` necessitates `entities_map`.<br />

    - `val` -> `sequence` <br />
      Nucleotide sequence of the nucleic acid.
    """

    entities_map: Dict[str, str]
    """Dictionary with:<br />

    - `key` -> `entity_name` <br />

    - `val` -> `identifier` <br />
      `identifier` is usually `uniprot_id` in case of `proteinChain` entities."""

    name: str
    """Name of the entity."""

    type: str
    """Type of the entity. It can be one of the following:<br />
    `proteinChain`, `dnaSequence`, `rnaSequence`, `ligand`, `ion`."""

    count: int
    """Count or copy number of the entity."""

    real_sequence: str
    """Real sequence of the entity.
    Refers to the amino acid or nucleic acid sequence."""

    start: int
    """Start position of the entity in the sequence."""

    end: int | None
    """End position of the entity in the sequence"""

    glycans: List[Dict[str, Any]] | None
    """List of glycans associated with the entity, if applicable."""

    modifications: List[Dict[str, Any]] | None
    """List of modifications associated with the entity, if applicable."""

    template_settings: Dict[str, Any]
    """Template settings for the entity, used for `proteinChain` entities."""

    def __init__(
        self,
        entity_info: Dict[str, Any],
        protein_sequences: Dict[str, str] | None = None,
        nucleic_acid_sequences: Dict[str, str] | None = None,
        entities_map: Dict[str, str] = {},
        a3m_lines: Dict[str, str] = {},
    ):

        super().__init__(
            entity_info=entity_info,
            protein_sequences=protein_sequences,
            nucleic_acid_sequences=nucleic_acid_sequences,
            entities_map=entities_map,
            a3m_lines=a3m_lines,
        )
        self.name = self.entity_name
        self.type = self.entity_type
        self.count = self.entity_count
        self.real_sequence = self.update_real_sequence()

    def create_af_sequence(self)-> Dict[str, Any]:
        """Create an AF sequence dictionary.

        Returns:
        - **af_sequence_dict (dict)**:<br />
            Dictionary in the following format:

        ```python
        # for proteinChain
        {"proteinChain": {"sequence": "AAAA",
                            "glycans": [... ],
                            "modifications": [... ],
                            "count": 1}}

        # for dnaSequence or rnaSequence
        {"dnaSequence"|"rnaSequence": {"sequence": "ACGA",
                                       "modifications": [... ],
                                       "count": 1}}

        # for ligand
        {"ligand": {"ligand": "ATP", "count": 1}}

        # for ion
        {"ion": {"ion": "MG", "count": 1}}
        ```
        """

        af_sequence_dict = {}

        if self.type == EntityType.PROTEIN_CHAIN:
            af_sequence_dict = {
                self.type: {
                    AFServerSequenceFields.SEQUENCE: self.real_sequence,
                    AFServerSequenceFields.GLYCANS: self.glycans,
                    AFServerSequenceFields.MODIFICATIONS: self.modifications,
                    AFServerSequenceFields.COUNT: self.count,
                }
            }

            af_sequence_dict[self.type].update(self.template_settings)
            af_sequence_dict[self.type].update(self.msa_settings)

        elif self.type in [EntityType.DNA_SEQUENCE, EntityType.RNA_SEQUENCE]:
            af_sequence_dict = {
                self.type: {
                    AFServerSequenceFields.SEQUENCE: self.real_sequence,
                    AFServerSequenceFields.MODIFICATIONS: self.modifications,
                    AFServerSequenceFields.COUNT: self.count,
                }
            }

        elif self.type in [EntityType.LIGAND, EntityType.ION]:
            af_sequence_dict = {
                self.type: {
                    self.type: self.name,
                    AFServerSequenceFields.COUNT: self.count
                }
            }

        return af_sequence_dict

    def update_real_sequence(self)-> str:
        """Update the real sequence of the entity.

        A real sequence is:
        - Amino acid sequence for `proteinChain`.
        - Nucleic acid sequence for `dnaSequence` and `rnaSequence`.

        If a `range` (i.e. [`start`, `end`]) is provided, slice the sequence accordingly.

        Returns:
        - **real_sequence (str)**:<br />
            Amino acid or nucleic acid sequence of the entity.
        """

        real_sequence = self.real_sequence
        start, end = self.start, self.end

        if self.type in [
            EntityType.PROTEIN_CHAIN,
            EntityType.DNA_SEQUENCE,
            EntityType.RNA_SEQUENCE
        ]:
            real_sequence = real_sequence[start - 1 : end]

        return real_sequence

    def get_name_fragment(self)-> str:
        """Get the name fragments of the entity.

        Name format: `"name_count_start-end"`

        Example:

        - For an entity with following attributes:\n
        ```
            name = "protA"
            count = 2
            start = 1
            end = 100
        ```

        the `name_fragment` will be `"protA_2_1-100"`.

        Returns:
        - **name_fragment (str)**:<br />
            Name fragment of the entity.
        """

        return f"{self.name}_{self.count}_{self.start}{af_constants.RES_RANGE_SEP}{self.end}"
