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
cd bin
# compile path
. set_path_REcoM.sh

# run REcoM1D
date
srun ./REcoM1D.x > "REcoM1D.out"
date

cd ..
