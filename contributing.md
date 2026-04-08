
# Contributing

## Code and documentation

- Please raise an [issue](https://github.com/isblab/af_pipeline/issues) regarding
  the bug you found or feature you would like to implement. If you would like to
  work on the issue see follow the steps below.

- [Fork](https://github.com/isblab/af_pipeline/fork) the af_pipeline repository.

- Clone the forked repository.
    ```bash
    git clone --recursive https://github.com/<USER_NAME>/af_pipeline.git
    ```

- Run [`setup.py`](./setup.py)
  ```bash
  python setup.py
  ```

- Make changes to the code and documention. Add tests if necessary.

  > [!IMPORTANT]
  > Before making any change in the code, use the network diagram to see how the
  > change might affect the existing code.

- Do the ["Documentation and tests"](#documentation-and-tests-check) check specified below.

- Create a [pull request](https://github.com/isblab/af_pipeline/pulls) and don't
  forget to mention the addressed issue.

## Documentation and tests

> [!NOTE]
> *Assuming your working directory is `/path/to/af_pipeline`*

-  Generate the documentation using `pdoc` manually as follows:
   ```bash
   python ./docs/network_viz.py # to update the network visualizations
   pdoc -t ./docs/template --search -e af_pipeline=https://github.com/isblab/af_pipeline/tree/main/af_pipeline/ --mermaid af_pipeline -o docs_/
   ```

- Alternatively, use [`generate_docs.py`](./docs/generate_docs.py) as:
  ```bash
  python ./docs/generate_docs.py
  ```

- The generated documentation will be available inside `docs_` directory.

- Refer to the [tests](./tests/) directory.

- Run tests using `pytest` as follows:
  ```bash
  pytest --cov=af_pipeline --doctest-modules --cov-context=test --cov-report=html
  ```

- Make sure that none of the tests fail and the documentation is generated successfully.