#!/usr/bin/python3
# -*- coding: utf-8 -*
#
# Load modules
#
from glob import glob as gg
#
import numpy as np 

#
from datetime import datetime,date
#
# path and filename of the main time configuration file 'time.recom'
#
path = '../'
filename_time = path + 'time.recom'

#
# load time properties from time.recom
#
class REcoM_time:
    def __init__(self, filename):
        # load time properties
        self.estimate_time_axis(filename)
        # estimate REcoM time initialization variables
        self.estimate_initialization_variables()
        
    def estimate_time_axis(self, filename):
        # read data
        with open(filename) as f:
            lines = f.readlines()
        for l in lines:
            tmp = l.replace(' ','').replace('\t','').split('=')[-1].strip()
            if 'start_date' in l:
                date_start = datetime.strptime(tmp,'%d-%m-%Y').toordinal()
            elif 'end_date' in l:
                date_end = datetime.strptime(tmp,'%d-%m-%Y').toordinal()
            elif 'dt' in l:
                dt = float(tmp)/24.
            elif 'spinup' in l:
                spinup = int(tmp)
        
        self.date_start, self.date_end = np.copy(date_start), np.copy(date_end)
        
        # account for spinup
        date_start = date_start - spinup
        
        # create time axis
        npt = int(abs(date_end - date_start)/dt)
        self.dates = np.linspace(date_start, date_end, npt+1)
        
    def estimate_initialization_variables(self):
        
        # preliminary computation
        dt = np.sum(np.diff(self.dates))/(len(self.dates)-1)
        date_start, date_end = np.min(self.dates), np.max(self.dates) 
        ndays = np.floor(date_end) - np.floor(date_start)
        
        time_start = date_start - np.floor(date_start)
        year_start = datetime.fromordinal(int(date_start)).year
        date_ref = datetime.strptime('01-01-' + str(year_start),'%d-%m-%Y').toordinal()
        day_start = int(date_start) - int(date_ref) + 1
       
        # timestep related variables
        self.step_per_day = int(1/dt)
        self.run_length = int(ndays)
        self.run_length_unit = 'd'
        

        # initial time 
        self.timenew = time_start
        self.daynew = day_start
        self.yearnew = year_start
        
#
recom_time = REcoM_time(filename_time)


#
# write time configuration information to temporary info file to be parsed in namelist.config
#
var=recom_time

fname = 'time_configuration_file.txt'
f = open(fname, 'w')

# info about time step
f.write('step_per_day=' + str(var.step_per_day) + '\n')
f.write('run_length=' + str(var.run_length) + '\n')
f.write('run_length_unit=' + '\'' + var.run_length_unit + '\'' +'\n')

# info about starting date
f.write('timenew=' + str(var.timenew) + '\n')
f.write('daynew=' + str(var.daynew) + '\n')
f.write('yearnew=' + str(var.yearnew) + '\n')
f.close()

#
# write time configuration in test.clock
#
fname = 'test.clock'
f = open(fname, 'w')
f.write(str(int(var.timenew)) + '\t' + str(int(var.daynew)) + '\t' + str(int(var.yearnew)) + '\n')
f.write(str(int(var.timenew)) + '\t' + str(int(var.daynew)) + '\t' + str(int(var.yearnew)) + '\n')
f.close()

