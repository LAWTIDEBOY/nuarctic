#!/usr/bin/env bash

set -e

# source environment file
source env.sh 

# create build  and bin directories from scratch
if [ ! -d "build" ]
then 
	mkdir build
else
	# clear build directory
	rm -r build
	mkdir build
fi
#bin
if [ ! -d "bin" ]
then 
	mkdir bin
	cd bin
	ln -s ../set_path_REcoM.sh .
	cd ..
else 
	# clear directory
	rm -r bin
	mkdir bin
	cd bin
	ln -s ../set_path_REcoM.sh .
	cd ..
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
