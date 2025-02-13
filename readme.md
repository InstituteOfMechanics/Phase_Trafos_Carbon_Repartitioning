# Modelling and finite element simulation of martensite and bainite phase transformations during quenching under consideration of carbon repartitioning

## Publication

This repository contains data and source code accompanying the publication

Tim Furlan, Markus Schewe, Philipp Scherm, Philipp Retzl, Ernst Kozeschnik and Andreas Menzel:
Modelling and finite element simulation of martensite and bainite phase transformations during quenching under consideration of carbon repartitioning (2025).
Mechanics of Materials, 105275.
DOI: [10.1016/j.mechmat.2025.105275](https://doi.org/10.1016/j.mechmat.2025.105275)

If you use the provided model in your project (or build upon it), please cite the journal article:

    @article{furlan2025,
	title = {Modelling and finite element simulation of martensite and bainite phase transformations during quenching under consideration of carbon repartitioning},
	journal = {Mechanics of Materials},
	pages = {105275},
	year = {2025},
	issn = {0167-6636},
	doi = {https://doi.org/10.1016/j.mechmat.2025.105275},
	author = {Tim Furlan and Markus Schewe and Philipp Scherm and Philipp Retzl and Ernst Kozeschnik and Andreas Menzel},
    }


## Repository structure

The repository contains the following directories:

- ttt-fit

  contains data and Python source code for the fitting of the JMAK model parameters with different strategies.

- simulations

  contains all source files for the simulations in Abaqus. This includes Abaqus python script files for model creation and result export, Fortran user subroutines for the material model, and shell scripts to run the complete simulation pipeline on Linux, as well as selected simulation results in the vtk format.

## Python venv

Python has to be installed on your system to execute any Python scripts in this repository. A Python venv can be created from the requirements.txt file to install all Python packages that the scripts depend on. To create the venv, open a terminal in the source directory of the repository and use the following commands

```
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Activate the new environment at any time by typing

```
source .venv/bin/activate
```

## License

All source code in this repository is made available under the GPLv3 license. Some of the scripts generate figures that have been published the the article mentioned above, and are subject to the corresponding copyright conditions.


## Contributions

This repository accompanies a journal article, therefore no external contributions are accepted. You are however welcome to fork the repository to develop your own simulations based on the source code under the terms of the GPLv3 license.


