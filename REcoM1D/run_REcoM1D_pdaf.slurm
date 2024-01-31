#!/bin/bash 
#SBATCH --job-name=pdaf
#SBATCH -p mpp
#SBATCH --qos=30min
#SBATCH --account=computing.computing
#SBATCH --ntasks=512
#SBATCH --cpus-per-task=2
#SBATCH --time=00:15:00



set -x

ulimit -s unlimited

# determine JOBID
JOBID=`echo $SLURM_JOB_ID |cut -d"." -f1`

# source environment file
source env.sh 

#----------------------------------------------------------------------------------------------
#
# script to run REcoM1D along a selected track/trajectory
# 
#----------------------------------------------------------------------------------------------
export PDAF_OUT_DIR=/albedo/work/user/nmamnun/nuarctic/REcoM1D/da_outputs
rm $PDAF_OUT_DIR/*

# manage time file infos consistently to mesh and to time.recom
cd config/scripts/
./REcoM_time.sh
cd ../..
# create symbolic link to the data
# mesh
cd grid
rm *.nc
ln -s ../../data/MESH/REcoM1D_mesh_v2.nc
# forcing
cd ../forcing
rm *.nc
ln -s ../../data/REcoM_forcing_data/REcoM1D_forcing_v2.nc
# atm deposition
cd ../data
rm *.nc
ln -s ../../data/atm_deposition/atm_deposition_v2.nc

# tracer initialisation
ln -s ../../data/initialization/tracer_initialization.nc

# perform computation 
cd ../bin
# compile path
. set_path_REcoM.sh

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# run REcoM1D
date
# ./REcoM1D.x > "REcoM1D.out"
srun ./REcoM1D.x -dim_ens 512
date

cd ..
