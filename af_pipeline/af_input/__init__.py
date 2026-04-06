"""
[af_input](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/af_input)
===================================

- Module for generating input files for AlphaFold3, AlphaFold2 and ColabFold predictions.

- This module contains classes and methods to generate input files for
AlphaFold3[^af3], AlphaFold2[^af2] and ColabFold predictions[^colabfold].

[^af3]: Abramson, J. et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024). (https://alphafoldserver.com/)

[^af2]: Jumper, J. et al. Highly Accurate Protein Structure Prediction with Alphafold. Nature 596, 583–589 (2021).

[^colabfold]: Mirdita, M. et al. ColabFold: making protein folding accessible to all. Nature Methods 19, 679–682 (2022). (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)

- The user provides an input config file that contains details about the jobs to be run. The format of the config file is as follows.
Check the [examples directory](https://github.com/isblab/af_pipeline/blob/main/examples/input/config.yaml) for a sample config file.

```yaml

af_input_jobs:

  job_cycle_1:

    - job_set_name: "job_set_a"
      modelSeeds: [1, 2] # -> will lead to 2 jobs
      entities: ...

    - job_set_name: "job_set_b"
      modelSeeds: 3 # -> will lead to 3 jobs
      entities: ...

  job_cycle_2:

    - modelSeeds: [4, 5, 6] # -> will lead to 3 jobs
      entities: ...

```

- `af_input_jobs` contains the input specifications for prediction jobs. It is organized into multiple job cycles.
- Each `job_cycle` is a group of related predictions on one or more systems. A system corresponds to a single set of input sequences.
- Each `job_cycle` contains a list of `job_set` instances (in config.yaml).
- Each `job_set` is a set of predictions (jobs) on a single system, with the constituent jobs varying only in the model seed.
- User only specifies the `job_set` and the `modelSeeds`. Each `job_set` is converted to a list of `job` instances based on `modelSeeds`.

<hr>

## Organization

- Keep in mind the following structure while using or adding new classes or methods to
  any of the submodules.

```mermaid
graph LR

    config.yaml -->|contains| af_input_jobs
    af_input_jobs -->|contains| job_cycle
    job_cycle -->|contains| job_set
    job_set -->|contains| job
    job -->|contains| entity

```

- Each `job` or `job_set` contains a list of entities that can be of one of the following types:

    - [`proteinChain`](https://github.com/google-deepmind/alphafold/blob/main/server/README.md#protein-chains)
    - [`dnaSequence`](https://github.com/google-deepmind/alphafold/blob/main/server/README.md#dna-chains)
    - [`rnaSequence`](https://github.com/google-deepmind/alphafold/blob/main/server/README.md#rna-chains)
    - [`ligand`](https://github.com/google-deepmind/alphafold/blob/main/server/README.md#ligands)
    - [`ion`](https://github.com/google-deepmind/alphafold/blob/main/server/README.md#ions)

</newline>

- Each entity in the `job` is an instance of :py:class:`af_pipeline.af_input.alphafold3.AFSequence`.


- In all three cases (AlphaFoldServer, AlphaFold2, ColabFold), the input data is provided as a dictionary, which is
  stored in a `YAML` file. See the [examples directory](https://github.com/isblab/af_pipeline/tree/main/examples)
  for sample input file.

- The user can create job cycles using the method for creating job cycles in the corresponding class:

  - **AlphaFoldServer**: `af_pipeline.af_input.alphafold3.AlphaFoldServer.create_af3_job_cycles`
  - **AlphaFold2**: `af_pipeline.af_input.alphafold2.AlphaFold2.create_af2_job_cycles`
  - **ColabFold**: `af_pipeline.af_input.colabfold.ColabFold.create_colabfold_job_cycles`

</newline>

- The method `write_job_files` corresponding to each submodule is used to
  save the generated files in appropriate format:

  - **AlphaFoldServer**: `af_pipeline.af_input.alphafold3.AlphaFoldServer.write_job_files`
  - **AlphaFold2**: `af_pipeline.af_input.alphafold2.AlphaFold2.write_job_files`
  - **ColabFold**: `af_pipeline.af_input.colabfold.ColabFold.write_job_files`

<hr>

## Usage

- Please refer to the [examples directory](https://github.com/isblab/af_pipeline/tree/main/examples) for sample scripts and
  config file.

- Use the following command to run the example script:
```
python create_af_jobs.py \\
    -i ./input/config.yaml \\
    -o ./output/af_input_jobs \\
    -p ./input/protein_sequences.fasta \\
    -n ./input/nucleic_acid_sequences.fasta
```

> [!TIP]
> You can use `FetchSequences` class from `IMP_Toolbox` to get the protein
sequences fasta file from a list of UniProt IDs.

<hr>

## Workflows

- Workflow for creating AlphaFold3 job `JSON` files:

```mermaid
sequenceDiagram
    autonumber
    participant User
    participant FetchSequences
    participant AlphaFold3

    rect rgb(210, 250, 200)
    User->>FetchSequences: uniprot_ids
    User->>FetchSequences: query_uniprot_for_sequences()
    FetchSequences->>User: protein_sequences (fasta format)
    end

    rect rgb(240, 255, 255)
    note over User, AlphaFold3: Create AlphaFold3 instance.
    User->>AlphaFold3: config_dict
    User->>AlphaFold3: protein_sequences (dict)
    note over User, AlphaFold3: The following parameters are optional.
    User->>AlphaFold3: nucleic_acid_sequences (dict)
    end

    rect rgb(245, 255, 200)
    note over User, AlphaFold3: Create AlphaFold3 job cycles.
    User->>output(.json): create_af3_job_cycles() & write_job_files()
    end
```

- Workflow for creating AlphaFold2/ColabFold job `FASTA` files:
```mermaid
sequenceDiagram
    autonumber
    participant User
    participant FetchSequences
    participant AlphaFold2/ColabFold

    rect rgb(210, 250, 200)
    User->>FetchSequences: uniprot_ids
    User->>FetchSequences: query_uniprot_for_sequences()
    FetchSequences->>User: protein_sequences (fasta format)
    end

    rect rgb(240, 255, 255)
    note over User, AlphaFold2/ColabFold: Create AlphaFold2/ColabFold instance.
    User->>AlphaFold2/ColabFold: config_dict
    User->>AlphaFold2/ColabFold: protein_sequences (dict)
    end

    rect rgb(245, 255, 200)
    note over User, AlphaFold2/ColabFold: Create AlphaFold2/ColabFold job cycles.
    User->>output(.fasta): create_af2_job_cycles() / create_colabfold_job_cycles() & write_job_cycles()
    end
```

- You can customize the above workflows as per the requirements.

<hr>

## FAQ

**Q. Can the same config file be used for AlphaFold3, AlphaFold2 and ColabFold?**
<details>
<summary>Yes, the same config file can be used for all three methods.</summary>
However, if the config file has non-protein entities, then they will be ignored
for AlphaFold2 and ColabFold.
Similarly, if the config file has `modelSeeds` attribute, then it will be ignored
for AlphaFold2 and ColabFold.

</details>


**Q. Is it necessary to specify all the attributes for each entity in the config file?**
<details>
<summary>No, only the `name` and `type` are mandatory for each entity.</summary>
Other attributes have default values as follows:
```
count: 1
range: [1, len(sequence)]
useStructureTemplate: true
maxTemplateDate: 2021-09-30
glycans: []
modifications: []
af_offset: {}
```

Some of the attributes can be retrieved in the config file while creating the
input `JSON` or `FASTA` files within the script.
For e.g., `af_offset` or `job_name` need not be specified in the config file.
They can be generated within the script and added to the config file as shown in
the [example script](https://github.com/isblab/af_pipeline/tree/main/examples/create_af_jobs.py).

> [!CAUTION]
> Doing the above will modify the original config file.
> This modification does not preserve comments in the config file.

> [!NOTE]
> If `job_name` is not specified in the config file, then it will be generated
> automatically using the information provided in the job set.
> It is a combination of `entity_name`, `range` and `count` for each entity and
> the `model_seed`.
> See :py:meth:`af_pipeline.af_input.alphafold3.AlphaFoldServer.AFJobSet.generate_job_set_name`.

> [!WARNING]
> - The current naming scheme can not distinguish between two job sets which differ
> only due to modifications or glycans.
>
> - Alphafold server has a limit of 100 characters for job names.
> Nonetheless, if such a case occurs, a warning is provided to provide a custom
> name.
</details>


**Q. Is the `YAML` format of the config file necessary?**
<details>
<summary>No, it is not necessary.</summary>
As long as the input data is provided as a dictionary
in the required format, it can be read from any file format (e.g., `JSON`, `TOML`,
etc.) or directly created in the wrapper script.
</details>


**Q. How can I specify model seeds?**
<details>
<summary>You can specify model seeds in the config file in the "modelSeeds" attribute.</summary>
Each job set within a job cycle has a "modelSeeds" attribute to specify the model seeds.

You can provide a single integer value (to directly denote the number of jobs) or
a list of integer values (to specifically denote the seed values).
In the former case, the seeds are generated using
:py:meth:`af_pipeline.utils.misc_utils.generate_seeds`.

> [!NOTE]
> model seeds are only applicable for AlphaFold3 predictions.
</details>

## Network

<body>
    <p>
        Double click on the node to go to the corresponding line in the source code.
    </p>
    <iframe src="../../docs/network_viz/network_af_input.html" width="800" height="600" frameborder="0">
    </iframe>
</body>

"""
