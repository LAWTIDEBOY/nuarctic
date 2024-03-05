#!/bin/bash 
#SBATCH --job-name=pdaf
#SBATCH -p mpp
#SBATCH --qos=30min
#SBATCH --account=computing.computing
#SBATCH --ntasks=512
#SBATCH --cpus-per-task=1
#SBATCH --time=00:29:00


set -x

ulimit -s unlimited

# determine JOBID
JOBID=`echo $SLURM_JOB_ID |cut -d"." -f1`

# source environment file
source env.sh 

export RECOM_MAIN_DIR="$PWD"
export PDAF_OUT_DIR=$RECOM_MAIN_DIR/da_outputs

if [ ! -d "$PDAF_OUT_DIR" ]; then
    echo "Directory $PDAF_OUT_DIR does not exist. Creating now."
    mkdir -p $PDAF_OUT_DIR
else
    echo "Directory $PDAF_OUT_DIR exists. Removing its contents."
    rm -rf $PDAF_OUT_DIR/*
fi


rm $PDAF_OUT_DIR/*
export MPMD_OUTDIR=$RECOM_MAIN_DIR/mpmd_outputs

if [ ! -d "$MPMD_OUTDIR" ]; then
    echo "Directory $MPMD_OUTDIR does not exist. Creating now."
    mkdir -p $MPMD_OUTDIR
else
    echo "Directory $MPMD_OUTDIR exists. Removing its contents."
    rm -rf $MPMD_OUTDIR/*
fi

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

np=1    # process per model
NENS=512  # Ensemble size

cd $RECOM_MAIN_DIR

echo ' '
echo 'creating directories...'
for((i=1;i<=$NENS;i++))
  do
    ENSstr=`printf %03d $i`
    if [ ! -d bin_${ENSstr} ]; then
      mkdir -p bin_${ENSstr}
    else
      rm -rf bin_${ENSstr}/*
    fi 
    bindir=`pwd`/bin_${ENSstr}
    export bindir
    cp $RECOM_MAIN_DIR/results/test.clock $bindir/
    cp $RECOM_MAIN_DIR/bin/REcoM1D.x $bindir/
    echo 'Run directory: ' $wdir
done
echo ' '
#start model run
# Create MPMD configuration file
if [ -f mpmd.conf ]; then
  rm mpmd.conf
fi
#
for((i=1;i<=$NENS;i++))
  do
    ENSstr=`printf %03d $i`
    echo 'preparing mpmd.conf...'${ENSstr}
    echo '#!/bin/sh' > recom1d_${ENSstr}.sh    
    bindir=`pwd`/bin_${ENSstr}
    export bindir
    echo "cd "$bindir >> recom1d_${ENSstr}.sh

    echo "export RECOM_NAMELIST_PATH="${RECOM_MAIN_DIR}/config/ >> recom1d_${ENSstr}.sh
    echo "export RECOM_GRID_PATH="${RECOM_MAIN_DIR}/grid/ >> recom1d_${ENSstr}.sh
    echo "export RECOM_DATA_PATH="${RECOM_MAIN_DIR}/data/ >> recom1d_${ENSstr}.sh
    echo "export RECOM_FORCING_PATH="${RECOM_MAIN_DIR}/forcing/ >> recom1d_${ENSstr}.sh
    echo "export RECOM_RESULT_PATH="$bindir/ >> recom1d_${ENSstr}.sh
    echo "./REcoM1D.x -dim_ens $NENS" >> recom1d_${ENSstr}.sh
    chmod +x recom1d_${ENSstr}.sh
    echo $((i-1))' ./recom1d_'${ENSstr}.sh >> mpmd.conf
done

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# run 
date
srun --multi-prog ./mpmd.conf
date

cd $RECOM_MAIN_DIR

for((i=1;i<=$NENS;i++))
  do
    ENSstr=`printf %03d $i`
    echo 'moving mpmd outputs...'
    mv bin_${ENSstr}/REcoM1d_diagnostics.nc $MPMD_OUTDIR/REcoM1d_diagnostics_${ENSstr}.nc
    mv bin_${ENSstr}/REcoM1d_outputs.nc $MPMD_OUTDIR/REcoM1d_outputs_${ENSstr}.nc
done

rm mpmd.conf
rm -r bin_*
rm recom1d_*

