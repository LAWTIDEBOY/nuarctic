#!/bin/bash 

#----------------------------------------------------------------------------------------------
#
# script to run REcoM1D along a selected track/trajectory
# 
#----------------------------------------------------------------------------------------------
cd bin
# compile path
. set_path_REcoM.sh

# run REcoM1D
./REcoM1D.x

cd ..
