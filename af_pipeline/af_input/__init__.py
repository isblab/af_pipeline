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

  - job_set_name: "job_set_a"
    modelSeeds: [1, 2] # -> will lead to 2 jobs with seeds 1 and 2
    entities: ...

  - job_set_name: "job_set_b"
    modelSeeds: 3 # -> will lead to 3 jobs with randomly generated seeds
    entities: ...

  - modelSeeds: [4, 5, 6] # -> will lead to 3 jobs
    entities: ...

```

- `af_input_jobs` contains the input specifications for prediction jobs. It is organized into multiple job sets.
- Each `job_set` is a set of predictions (jobs) with the same input specifications except for the model seed.
- Each `job_set` is converted to a list of `job` instances based on `modelSeeds`.

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
      -n ./input/nucleic_acid_sequences.fasta \\
      -t AF3 # replace with AF2 or ColabFold
  ```

> [!TIP]
> You can use [`sequence`](https://github.com/isblab/IMP_Toolbox/blob/main/IMP_Toolbox/sequence/sequence.py)
> module from [`IMP_Toolbox`](https://github.com/isblab/IMP_Toolbox) to get the protein sequences fasta file 
> from a list of UniProt IDs.

<hr>

## Organization

- Keep in mind the following structure while using or adding new classes or methods to
  any of the submodules.

```mermaid
graph LR

    config.yaml -->|contains| af_input_jobs
    af_input_jobs -->|contains| job_set
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

- The method `write_job_files` in each submodule is used to save the generated files in appropriate format:

  - **AlphaFoldServer**: `af_pipeline.af_input.alphafold3.AlphaFoldServer.write_job_files`
  - **AlphaFold2**: `af_pipeline.af_input.alphafold2.AlphaFold2.write_job_files`
  - **ColabFold**: `af_pipeline.af_input.colabfold.ColabFold.write_job_files`

<hr>

## Workflows

- Workflow for creating AlphaFold3 job `JSON` files or AlphaFold2/ColabFold job `FASTA` files:

```mermaid

graph TD
  A([config.yaml]) --> B[/config_dict/]
  C([protein_sequences.fasta]) --> D[/protein_sequences/]
  E([nucleotide_sequences.fasta]) --> F[/nucleic_acid_sequences/]
  B --> G[AlphaFoldServer instance]
  click G "af_input/alphafold3.html#AlphaFoldServer" "AlphaFoldServer" _blank
  D --> G
  F -- optional --> G
  G --> H[write_job_files]
  click H "af_input/alphafold3.html#AlphaFoldServer.write_job_files" "write_job_files" _blank
  H --> J([JSON files for AF3 jobs or FASTA files for AF2/ColabFold jobs])
```

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
For e.g., `af_offset` or `job_set_name` need not be specified in the config file.
The accompanying `af_input_jobs.json` file will have these attributes filled in
based on the information provided in the config file and the job set.

> [!NOTE]
> If `job_set_name` is not specified in the config file, then it will be generated
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
Each job set has a "modelSeeds" attribute to specify the model seeds.

You can provide a single integer value (to directly denote the number of jobs) or
a list of integer values (to specifically denote the seed values).
In the former case, the seeds are generated using
:py:meth:`af_pipeline.utils.misc_utils.generate_seeds`.

> [!NOTE]
> model seeds are only applicable for AlphaFold3 predictions.
</details>

"""
