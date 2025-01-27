# Modelling and finite element simulation of martensite and bainite phase transformations during quenching under consideration of carbon repartitioning

This repository contains the source code to reproduce the examples published in [@furlan2025]. 


## Repository structure

The repository contains the following directories:

- ttt-fit

  contains data and Python source code for the fitting of the JMAK model parameters with different strategies.

- simulations

  contains all source files for the simulations in Abaqus. This includes Abaqus python script files for model creation and result export, Fortran user subroutines for the material model, and shell scripts to run the complete simulation pipeline on Linux.

- vtk-results

  contains simulation results for the bearing ring bvp, exported to the open vtk format using [Paraqus](https://www.github.com/InstituteOfMechanics/Paraqus).


## Python venv

Python has to be installed on your system to execute any Python scripts in this repository. A Python venv can be created from the requirements.txt file to install all Python packages that the scripts depend on. To create the venv, open a Terminal in the source directory of the repository and use the following commands

```
python -m pip install -r requirements.txt
```

Activate the new environment by

```
source .venv/bin/activate
```

## License

All source code in this repository is made available under the GPLv3 license. Some of the scripts generate figures that have been published in [@furlan2025], and are subject to the corresponding copyright conditions.


## Contributions

This repository accompanies the publication [@furlan2025], therefore no external contributions are accepted. You are however welcome to fork the repository to develop your own simulations based on the source code.


## References
