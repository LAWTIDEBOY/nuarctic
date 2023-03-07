#!/bin/bash 
#SBATCH --job-name=simu1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH -o slurm-out.out
#SBATCH -e slurm-err.out

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

# run REcoM1D
date
./REcoM1D.x > "REcoM1D.out"
date

cd ..
