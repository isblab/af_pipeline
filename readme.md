# AF Pipeline

![Coverage Status](./tests/badges/coverage-badge.svg?raw=true)

A module to assist in creating input files and processing and analyzing output
predictions from AlphaFold predictions.

## Installation
- Clone the repository
- Install the required packages using pip:
    ```bash
    pip install -r requirements.txt
    ```
- Add the PYTHONPATH environment variable to include the af_pipeline directory:
    ```bash
    export PYTHONPATH=$PYTHONPATH:/path/to/af_pipeline
    ```

## Usage
- Refer to the documentation.
- Example scripts are available in the [examples](./examples/) directory.

> [!NOTE]
> *Assuming your working directory is `/path/to/af_pipeline`*
>
> For now, generate documentation using `pdoc` manually as follows:
> ```bash
> pdoc --mermaid af_pipeline
> ```
>
> The generated documentation will be available at `http://localhost:8080`
>
> To generate the network visualizations, run:
> ```bash
> python docs/network_viz.py
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
> ```

## Information

__Author(s):__

__Date:__

__License:__

__Testable:__

__Parallelizeable:__

__Publications:__