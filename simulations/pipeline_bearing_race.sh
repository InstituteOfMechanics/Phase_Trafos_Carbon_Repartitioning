#!/usr/bin/env bash

# Main script for the bearing race example
#
# Executing this script will create the different models, perform the
# simulations, export results, and generate the plots.

VENV_PATH="../.venv/bin/activate"
ONEAPI_SCRIPT="/opt/intel/oneapi/setvars.sh"
ABAQUS_BIN="/opt/abaqus/2023/SIMULIA/Commands/abaqus"
UMAT_DIR="UMAT"
UMAT_MAIN_FILE="umat_main.f"

SIMULATION_DIR="simulations_bearing_race"
VTK_DIR="vtk-results"

source $VENV_PATH
MESH_SIZE=1

source $ONEAPI_SCRIPT

# create the simulation directory if needed
mkdir $SIMULATION_DIR -p

# copy the abaqus environment file to the simulation folder
cp abaqus_v6.env $SIMULATION_DIR

# create the input for the jobs
$ABAQUS_BIN cae nogui=cae_create_model_bearing_race.py -- --dir $SIMULATION_DIR --material_flag 0 --quench_temperature 20 --mesh_size $MESH_SIZE
$ABAQUS_BIN cae nogui=cae_create_model_bearing_race.py -- --dir $SIMULATION_DIR --material_flag 1 --quench_temperature 20 --mesh_size $MESH_SIZE
$ABAQUS_BIN cae nogui=cae_create_model_bearing_race.py -- --dir $SIMULATION_DIR --material_flag 0 --quench_temperature 60 --mesh_size $MESH_SIZE
$ABAQUS_BIN cae nogui=cae_create_model_bearing_race.py -- --dir $SIMULATION_DIR --material_flag 1 --quench_temperature 60 --mesh_size $MESH_SIZE

# go to the simulation directory and execute all simulations
cd $SIMULATION_DIR

for FILE_NAME in bearing_race*.inp; do
	# extract job name without extension
	JOB_NAME="${FILE_NAME%.inp}"
	
	echo "============================================================="
	echo "Processing job ${JOB_NAME}"
	echo "============================================================="
	
    $ABAQUS_BIN job=$JOB_NAME user=../$UMAT_DIR/$UMAT_MAIN_FILE cpus=4 ask_delete=OFF -interactive
    $ABAQUS_BIN cae nogui="../cae_export_results_bearing_race.py" -- --job_name $JOB_NAME
    $ABAQUS_BIN cae nogui="../cae_export_results_bearing_race_contours.py" -- --job_name $JOB_NAME
    $ABAQUS_BIN cae nogui="../cae_export_vtk_bearing_race.py" -- --job_name $JOB_NAME --output_dir $VTK_DIR
done

cd ..
python plot_results_bearing_race.py
python plot_results_bearing_race_contours.py
