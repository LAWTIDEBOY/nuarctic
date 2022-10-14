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

#----------------------------------------------------------------------------------------------
#
# script to run REcoM1D along a selected track/trajectory
# 
#----------------------------------------------------------------------------------------------
# create symbolic link to the data
# mesh
cd grid
rm *.nc
ln -s ../../data/MESH/REcoM1D_mesh.nc
# forcing
cd ../forcing
rm *.nc
ln -s ../../data/REcoM_forcing_data/REcoM1D_forcing.nc
# atm deposition
cd ../data
rm *.nc
ln -s ../../data/atm_deposition/atm_deposition.nc

# perform computation 
cd ../bin
# compile path
. set_path_REcoM.sh

# run REcoM1D
date
./REcoM1D.x > "REcoM1D.out"
date

cd ..
