#!/home/ollie/fbirrien/miniconda3/envs/pyfesom2/bin/python3
# -*- coding: utf-8 -*-

"""
scripts to merge REcoM-related variables from FESOM runs

"""
#
from glob import glob as gg
#
import xarray as xr
#
from netCDF4 import Dataset 
#
#-------------------------------------------------------
# path to transect data
path_FESOM='/work/ollie/fbirrien/NuArctic/FESOM_Outputs/'
path_output = path_FESOM + 'Transects/'

# variable of interest
#var_name = ['a_ice', 'Alk', 'benC', 'benN', 'benSi', 'pCO2','DetCalc', 'DetC', 'DetN', 'DetSi', 'DFe' , 'DiaC', 'DiaChl', 'DiaN', 'DiaSi', 'DIC', 'DIN', 'DOC' , 'DON' , 'HetC', 'HetN', 'idetz2calc', 'idetz2c', 'idetz2n', 'idetz2si', 'Kv', 'O2', 'PAR', 'pCO2s', 'PhyCalc', 'PhyC', 'PhyChl', 'PhyN', 'runoff', 'salt', 'temp', 'u', 'uice', 'v','vice','w', 'Zoo2C', 'Zoo2N']

var_name = ['a_ice', 'temp']

for vname in var_name:
  vlist = sorted(gg(path_output + vname + '*'))
  # load files corresponding to a selected variables for the different output years
  for i,filename in enumerate(vlist):
    dataset = xr.open_mfdataset(filename, combine="by_coords")
    # merge data sets
    array=dataset if i<1 else xr.merge([array, dataset])

  # write merged (year) output to file
  filename_output=vlist[0].replace('.' + vlist[0].split('.')[-2], '')
  array.to_netcdf(filename_output)
