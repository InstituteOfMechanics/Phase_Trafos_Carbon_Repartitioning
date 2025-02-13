# Fitting of time-temperature-transformation (TTT) diagrams

This folder contains the source data and Python script for the fitting of JMAK model parameters based on digitised TTT diagrams. The reference data is stored in the file ``ttt-diagram_100cr6_kaymak.csv`` and based on the diagram published by [Kaymak (2007)](http://dx.doi.org/10.25673/4872). 


## Requirements

The processing and plotting is performed in Python. The easiest way is to install a venv from the requirements.txt in the repository base directory (see the readme there). Alternatively, you can install the packages yourself.


## Running the scripts

By executing the python scrip ``fit_jmak_parameters_100cr6.py``, all variants for the fitting of the parameters are performed, and the TTT diagrams published in the [journal article](https://doi.org/10.1016/j.mechmat.2025.105275) are recreated. The final polynomial coefficients are also printed.
