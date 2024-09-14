# Documentation for The Adaptable Input Output System

# Generate User Guide in html with Sphinx

This user guide is hosted in readthedocs: https://adios2.readthedocs.io/en/latest/

To generate the User Guide under docs/user_guide/build/html format from the Sphinx source files, using pip:

```bash
$ cd ADIOS2/docs
docs$ python3 -m venv .
docs$ . bin/activate
docs$ pip3 install -r requirements.txt
user_guide$ cd user_guide
user_guide$ make html
```

Or using conda:

```bash
$ cd ADIOS2/docs
docs$ conda env create -f environment.yml
docs$ conda activate adios2-python-docs
user_guide$ cd user_guide
user_guide$ make html
```

# Updating dependencies

Read the Docs uses only the environment.yml file, so if you make changes to dependencies in requirements.txt, they will not take effect in RTD.
The requirements.txt is provided only for building locally for those who don't have conda installed.
If you make changes to the conda environment, you should update requirements.txt as well, either manually, or by running the following command (with the conda adios2-python-docs environment activated):

```bash
docs$ pip list --format=freeze > requirements.txt
```
