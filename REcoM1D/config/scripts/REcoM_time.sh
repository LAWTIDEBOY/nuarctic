#!/bin/bash -f

rm time_configuration_file.txt test.clock
# run python scripts that read the time configuration file 'time.recom' and create temporary time file and time.clock
./write_time_configuration_REcoM.py

# change the information in the time configuration REcoM namelist 
./parse_namelist.sh ../namelist.config time_configuration_file.txt

# update time.clock in results directory
rm ../../results/test.clock
cp test.clock ../../results/.
