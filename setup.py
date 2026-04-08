import os
import json
import logging
import textwrap
import subprocess

## Prerequisites
# conda or miniconda

SUPER_DIR = os.path.dirname(os.path.abspath(__file__))
SETUP_CONFIG_PATH = os.path.join(SUPER_DIR, "config_setup.json")
logging.basicConfig(
    filename=os.path.join(SUPER_DIR, 'setup.log'),
    filemode="a",
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)

BASHRC_PATH = os.path.expanduser("~/.bashrc")

PYTHON_VERSION = "3.12"
CONDA_ENV_NAME = "af_pipeline"

SUCCESS_MSG = f"""
Full setup completed successfully!
Please restart your terminal or run 'source ~/.bashrc' to apply changes.
Don't forget to activate the conda environment using:

    conda activate {CONDA_ENV_NAME}
"""

FAIL_MSG = f"""
Setup completed, but some things encountered errors.
Please check setup.log for details.
"""

config = {
    "SUPER_DIR": SUPER_DIR,
    "PYTHON_VERSION": PYTHON_VERSION,
    "CONDA_ENV_NAME": CONDA_ENV_NAME,
}

def get_conda_environments() -> tuple[dict[str, dict[str, str]], str | None]:
    """
    Runs 'conda info --envs' and parses the output to return a list of
    conda environments (name and path).
    Courtesy of Google Gemini
    """
    try:
        # Run the command and capture the output
        result = subprocess.run(
            ['conda', 'info', '--envs'],
            capture_output=True,
            text=True,
            check=True,
            shell=False # Set to True on Windows if needed, but usually False is better
        )
        output = result.stdout

        environments = {}
        # Split the output by lines and skip the first line (header)
        lines = output.strip().split('\n')
        for line in lines:
            if line.startswith('#'):
                continue # Skip header lines
            parts = line.split()
            if parts:
                name = parts[0]
                path = parts[-1]
                # An asterisk indicates the active environment
                is_active = '*' in parts
                environments[name] = {'path': path, 'active': is_active}

        return environments, None

    except FileNotFoundError:
        error = "Error: 'conda' command not found. Make sure conda is installed and set up."
        return [], error

    except subprocess.CalledProcessError as e:
        error = f"Error running conda command: {e.stderr}"
        return [], error

def add_statement_to_bashrc(statement: str):
    with open(BASHRC_PATH, 'r') as f:
        bashrc_content = f.readlines()
        bashrc_content = [line.strip() for line in bashrc_content if not line.strip().startswith("#")]
        if statement not in bashrc_content:
            # print(f"Adding {statement} to ~/.bashrc...")
            with open(BASHRC_PATH, 'a') as f_append:
                f_append.write(f'\n# Added by {os.path.abspath(__file__)}\n{statement}\n')

def setup_conda_environment(
    environment_name: str,
    python_version: str,
    logger: logging.Logger,
) -> tuple[list[str], str | None]:

    error = None
    logger.info(f"Creating conda environment '{environment_name}'.")

    existing_envs, error = get_conda_environments()

    if error:
        logger.warning(f"Error getting existing conda environments: {error}")

    if environment_name in existing_envs:
        logger.info(f"Conda environment '{environment_name}' already exists.")

    else:
        logger.info(f"Conda environment '{environment_name}' does not exist. Creating ...")
        result = subprocess.run(
            ["conda", "create", "-n", environment_name, f"python={python_version}", "-y"],
            check=False,
            text=True,
            capture_output=True
        )

        if result.returncode == 0:
            logger.info(f"Successfully created conda environment '{environment_name}'.")

        else:
            logger.error(f"Error creating conda environment '{environment_name}': {result.stderr.replace('\n', ' ')}")
            error = result.stderr

    conda_prefix = ["conda", "run", "-n", environment_name]

    return conda_prefix, error

def is_safe_to_delete(path_to_delete: str) -> tuple[bool, ValueError | None]:
    if path_to_delete in ["", "*", "/"]:
        var_name = f"{path_to_delete=}".split('=')[0]
        return False, ValueError(
            f"{var_name} is set to '{path_to_delete}'. \
            Unless you want to delete everything in your system, \
            please change the {var_name} variable in the script to a safe path."
        )

    else:
        return True, None

if __name__ == "__main__":

    logger = logging.getLogger(__name__)
    logger.info("="*30 + " Starting setup.py " + "="*30)
    things_to_install = ["python_packages", "pae_to_domains"]
    error_occurred = {k: False for k in things_to_install}

    conda_prefix, error = setup_conda_environment(
        environment_name=CONDA_ENV_NAME,
        python_version=PYTHON_VERSION,
        logger=logger,
    )

    if error:
        logger.error(f"Error setting up conda environment: {error}")
        error_occurred["conda_setup"] = True
    else:
        logger.info(f"Conda environment '{CONDA_ENV_NAME}' is set up.")

    # Install required python packages
    result = subprocess.run(
        conda_prefix +
        ["python", "-m", "pip", "install", "-r", "requirements.txt"],
        check=False,
        text=True,
        capture_output=True
    )
    if result.returncode == 0:
        logger.info("Successfully installed required python packages.")
    else:
        error_occurred["python_packages"] = True
        logger.error(f"Error installing required python packages: {result.stderr.replace("\n", " ")}")

    # verify if pae_to_domains submodule is set up correctly, set up if not
    if not os.path.exists(os.path.join(SUPER_DIR, "af_pipeline/pae_to_domains/pae_to_domains.py")):
        logger.info("submodule pae_to_domains is not initialized. Initializing git submodules ...")
        results = subprocess.run(
            ["git", "submodule", "update", "--init", "--recursive"],
            check=True,
            text=True,
            capture_output=True,
        )
        if results.returncode == 0:
            logger.info("Successfully initialized git submodules.")
        else:
            error_occurred["pae_to_domains"] = True
            logger.error(f"Error initializing git submodules: {results.stderr.replace('\n', ' ')}")
    else:
        logger.info("pae_to_domains submodule already initialized. Skipping.")

    # add af_pipeline to bashrc
    export_statement = f'export PYTHONPATH="$PYTHONPATH:{SUPER_DIR}"'
    add_statement_to_bashrc(export_statement)
    logger.info(f"Added '{export_statement}' to ~/.bashrc.")

    logger.info("-"*79)

    if all(error is False for error in error_occurred.values()):
        print(textwrap.dedent(SUCCESS_MSG))

    else:
        print(textwrap.dedent(FAIL_MSG))
        for thing, error in error_occurred.items():
            if error is False:
                continue
            print(f"- Error occured in {thing}.")

    with open(SETUP_CONFIG_PATH, 'w') as f:
        json.dump(config, f, indent=4)