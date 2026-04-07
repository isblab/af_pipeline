# AF-Pipeline

<!-- ![Coverage Status](./tests/badges/coverage-badge.svg?raw=true) -->

[![codecov](https://codecov.io/gh/isblab/af_pipeline/branch/main/graph/badge.svg)](https://codecov.io/gh/isblab/af_pipeline)

![af_pipeline_logo](./docs/assets/af_pipeline_logo.png)

A module to assist in creating input files and analyzing and assessing output
predictions from AlphaFold2 and AlphaFold3.

## Installation

- Clone the repository
    ```bash
    git clone --recursive https://github.com/isblab/af_pipeline.git
    ```

  If cloned without `--recursive` flag, do the following.
  ```bash
  git submodule init
  git submodule update
  ```

- Run `setup.py`
  ```bash
  python setup.py
  ```

## Usage
- Refer to the [documentation](./docs/).
- Example scripts are available in the [examples](./examples/) directory.

<!-- 
## Documentation

 - Currently, only network visualizations are available in [docs](./docs/) directory.

> [!NOTE]
> *Assuming your working directory is `/path/to/af_pipeline`*
>
> For now, generate documentation using `pdoc` manually as follows:
> ```bash
> pdoc -t ./docs/template --search -e af_pipeline=https://github.com/isblab/af_pipeline/tree/main/af_pipeline/ --mermaid af_pipeline
> ```
>
> The generated documentation will be available at `http://localhost:8080`
>
> To generate the network visualizations, run:
> ```bash
> python docs/network_viz.py -d af_pipeline
> ```

## Testing
- Refer to the [tests](./tests/) directory.

> [!NOTE]
> *Assuming your working directory is `/path/to/af_pipeline`*
>
> For now, run tests using `pytest` manually as follows:
> ```bash
> pytest --cov=af_pipeline --doctest-modules --cov-context=test
> ```
>
> The coverage badge was generated locally using `genbadge` as follows.
> ```bash
> coverage xml -o tests/reports/coverage.xml
> genbadge coverage -i tests/reports/coverage.xml -l -o tests/badges/coverage-badge.svg
> ``` -->

## Information

__Author(s):__ Omkar Golatkar

__Date:__ 21st March, 2026

__License:__ GPLv3

__Testable:__ Yes

__Parallelizeable:__

__Publications:__
