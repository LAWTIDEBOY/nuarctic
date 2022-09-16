#!/usr/bin/env bash

set -e

# source environment fike
source env.sh 

# create build  and bin directories if they do not exist
#build
if [ ! -d "build" ]
then 
	mkdir build
fi
#bin
if [ ! -d "bin" ]
then 
	mkdir bin
	cp set_path_REcoM.sh bin/.
fi

# clean build directory, clear cache
cd build
if [ -f "CMakeCache.txt" ]
then 
	rm -r *
fi
cmake ..
# start with cleared compiled modules
cmake --build . --target clean
# perform compiling
cmake --build .
# create executable and move it to bin (create bin directory or overwrite files)
cmake --install .

cd ..
