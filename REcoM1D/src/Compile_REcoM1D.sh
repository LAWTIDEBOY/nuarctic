#!/bin/bash 

#
# script that compiles REcom1D using recom utilities (gasx,...) and clean the directory

# clean directory
make clean

# link to the auxiliary modules
path_lib=/home/fbirrien/NuArctic/nuarctic/REcoM1D/src/recom/

ln -s $path_lib/gasx.mod .

# compile the code
make

# clean the directory
#rm *.mod

