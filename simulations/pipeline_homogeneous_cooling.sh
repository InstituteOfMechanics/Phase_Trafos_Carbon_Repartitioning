#!/usr/bin/env bash

# Main script for the homogeneous cooling example.
#
# Executing this script will create the different models, perform the
# simulations, export results, and generate the plots.

VENV_PATH="../.venv/bin/activate"
ONEAPI_SCRIPT="/opt/intel/oneapi/setvars.sh"
ABAQUS_BIN="/opt/abaqus/2023/SIMULIA/Commands/abaqus"
UMAT_DIR="UMAT"
UMAT_MAIN_FILE="umat_main.f"

SIMULATION_DIR="simulations_homogeneous_cooling"

source $VENV_PATH
source $ONEAPI_SCRIPT

# create the simulation directory if needed
mkdir $SIMULATION_DIR -p

# copy the abaqus environment file to the simulation folder
cp abaqus_v6.env $SIMULATION_DIR

# create the input for the jobs
$ABAQUS_BIN cae nogui=cae_create_model_homogeneous.py -- --dir $SIMULATION_DIR --initial_temperature 500 --htc 0.00001 --material_flag 0
$ABAQUS_BIN cae nogui=cae_create_model_homogeneous.py -- --dir $SIMULATION_DIR --initial_temperature 500 --htc 0.00001 --material_flag 1

# go to the simulation directory and execute all simulations
cd $SIMULATION_DIR

for FILE_NAME in homogeneous*.inp; do
	# extract job name without extension
	JOB_NAME="${FILE_NAME%.inp}"
	
	echo "============================================================="
	echo "Processing job ${JOB_NAME}"
	echo "============================================================="
	
    $ABAQUS_BIN job=$JOB_NAME user=../$UMAT_DIR/$UMAT_MAIN_FILE ask_delete=OFF -interactive
    $ABAQUS_BIN cae nogui="../cae_export_results_homogeneous.py" -- --job_name $JOB_NAME
done

# create the plots from the exported results
cd ..
python "plot_results_homogeneous_cooling.py"
