#!/bin/bash -f

# run python scripts that read the time configuration file 'time.recom' and create temporary file
./write_time_configuration_REcoM.py

# change the information in the time configuration REcoM namelist 
./parse_namelist.sh ../namelist.config time_configuration_file.txt
