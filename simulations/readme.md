# Abaqus simulations

This directory contains the source code for all Abaqus simulations presented in [@furlan2025].

## Requirements

TODO

Make sure to set the paths in the top section of all shell scripts (``pipeline_XXX.sh``) so that they match your installations.


### Abaqus

The simulations are performed with Abaqus, so you need to have it installed.

### OneAPI

The fortran subroutines use the Intel OneAPI library, so you need to install it (base and HPC toolkits).

### Python

The post-processing and plotting is performed in Python. The easiest way is to install a venv from the requirements.txt in the repository base directory (see the readme there). Alternatively, you can install the packages yourself.

### LaTeX

Some of the plots use LaTeX to create nice captions.


## Reproduction of isothermal bainite transformations

These simulations are not presented in the paper, but were performed to validate the implementation of the JMAK model. The simulations contain a single hex element with prescribed, constant temperature above martensite start. Only the model for 100Cr6, i.e. with constant austenite phase carbon content, is considered. The times until 1% and 99% transformation to bainite are achieved are interpolated from the simulation results, and compared to the reference values from the underlying ITT diagram in a text output.

The shell script ``pipeline_homogeneous_isothermal.sh`` automates the model creation, simulations, export of results, and analysis. The following steps are performed:

  - activate the python venv
  - load ``oneapi`` so that the Abaqus user subroutines can use the library
  - create a subdirectory ``simulations_homogeneous_isothermal``, which will hold all files created
  - copy the environment file to the subdirectory, so that it will be used for compiler configuration when the subroutines are compiled
  - execute the script ``cae_create_model_homogeneous.py`` in Abaqus python to generate input files for temperatures of 320, 350, and 380°C.
  - for each input file:
    - perform the corresponding simulation
    - execute the script ``cae_export_results_homogeneous.py`` in Abaqus python to store selected output in a numpy .npz file
  - execute the script ``plot_results_homogeneous_isothermal.py`` in Python outside of Abaqus to print the comparison of transformation times and generate a plot of the transformations over time, which is saved as ``plots/results_homogeneous_isothermal.png``



## Comparison of thermal conditions

Four simulations are performed at a constant ambient temperature of 300°C. The following cases are evaluated:

| Material     | Temperature condition                                          |
|--------------|----------------------------------------------------------------|
| 100Cr6       | Prescribed constant temperature                                |
| 100Cr6       | Temperature change due to latent heat and convectional cooling |
| 100CrMnSi6-4 | Prescribed constant temperature                                |
| 100CrMnSi6-4 | Temperature change due to latent heat and convectional cooling |


The shell script ``pipeline_homogeneous_no_martensite.sh`` automates the model creation, simulations, export of results, and analysis. The following steps are performed:

  - activate the python venv
  - load ``oneapi`` so that the Abaqus user subroutines can use the library
  - create a subdirectory ``simulations_homogeneous_no_martensite``, which will hold all files created
  - copy the environment file to the subdirectory, so that it will be used for compiler configuration when the subroutines are compiled
  - execute the script ``cae_create_model_homogeneous.py`` in Abaqus python to generate input files for the four test cases
  - for each input file:
    - perform the corresponding simulation
    - execute the script ``cae_export_results_homogeneous.py`` in Abaqus python to store selected output in a numpy .npz file
  - execute the script ``plot_results_homogeneous_no_martensite.py`` in Python outside of Abaqus to create the plot of temperature, bainite fraction, carbon content, and transformation factor, which is saved under ``plots/results_homogeneous_no_martensite.png``.


## Homogeneous cooling

Simulations are performed for the homogeneous cooling of a single hex element from 500°C to 60°C. 100Cr6 and 100CrMnSi6-4 are considered.

The shell script ``pipeline_homogeneous_cooling.sh`` automates the model creation, simulations, export of results, and analysis. The following steps are performed:

  - activate the python venv
  - load ``oneapi`` so that the Abaqus user subroutines can use the library
  - create a subdirectory ``simulations_homogeneous_cooling``, which will hold all files created
  - copy the environment file to the subdirectory, so that it will be used for compiler configuration when the subroutines are compiled
  - execute the script ``cae_create_model_homogeneous.py`` in Abaqus python to generate input files for the two test cases
  - for each input file:
    - perform the corresponding simulation
    - execute the script ``cae_export_results_homogeneous.py`` in Abaqus python to store selected output in a numpy .npz file
  - execute the script ``plot_results_homogeneous_cooling.py`` in Python outside of Abaqus to create the plots of temperature, volume fractions, and carbon content (``plots/results_homogeneous_cooling.png``), as well as strain over time (``plots/results_homogeneous_cooling_strain.png``).


## Quenching of bearing ring

The quenching of a bearing ring in water is simulated for the following cases:

| Material     | Quenching bath |
|--------------|----------------|
| 100Cr6       | Water, 20°C    |
| 100Cr6       | Water. 60°C    |
| 100CrMnSi6-4 | Water, 20°C    |
| 100CrMnSi6-4 | Water, 60°C    |

The shell script ``pipeline_bearing_race.sh`` automates the model creation, simulations, export of results, and analysis. The following steps are performed:

  - activate the python venv
  - load ``oneapi`` so that the Abaqus user subroutines can use the library
  - create a subdirectory ``simulations_bearing_race``, which will hold all files created
  - copy the environment file to the subdirectory, so that it will be used for compiler configuration when the subroutines are compiled
  - execute the script ``cae_create_model_bearing_race.py`` in Abaqus python to generate input files for the four test cases
  - for each input file:
    - perform the corresponding simulation
    - execute the script ``cae_export_results_bearing_race.py`` in Abaqus python to store selected output in a numpy .npz file
    - execute the script ``cae_export_results_bearing_race_contours.py`` in Abaqus python to store displacement output related to the resulting contour in a numpy .npz file
    - execute the script ``cae_export_vtk_bearing_race.py`` in Abaqus python to export vtk files from the simulation results
  - execute the script ``plot_results_bearing_race.py`` in Python outside of Abaqus to create the plots of the evolution at the selected points P1 and P2 (``plots/...`` and ``plots/...``)
  - execute the script ``plot_results_bearing_race_contours.py`` in Python outside of Abaqus to create the plots of the inner bore contours after quenching (``plots/...``)

TODO: Paraview



