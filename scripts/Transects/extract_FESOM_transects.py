#!/home/ollie/fbirrien/miniconda3/envs/pyfesom2/bin/python3
# -*- coding: utf-8 -*-

"""
scripts to extract REcoM-related variables from FESOM runs

"""
#
import os
import argparse
from glob import glob as gg
#
import numpy as np
import xarray as xr
from datetime import datetime
#
import pyfesom2 as pf
#
from netCDF4 import Dataset 
#
import pickle
#
#--------------------------------------------------------------------------------------
class Trajectory:
  def __init__(self, filename):
    # read MOSAiC trajectory data
    self.read_trajectory(filename)
    # month & year information
    self.dates_information()
  
  def read_trajectory(self,filename):
    ncid = Dataset(filename, "r", format="NETCDF4")
    # read dates/work/ollie/fbirrien/trajectories
    self.dates = ncid.variables['dates'][:]
    # read coordinates
    self.longitude, self.latitude = ncid.variables['longitude'][:], ncid.variables['latitude'][:]
    ncid.close()

  def dates_information(self):
    # store month & year information along the track
    month, year = [], []
    for dt in self.dates:
      tmp=datetime.fromordinal(int(dt)) 
      month.append(tmp.month), year.append(tmp.year)
    self.month, self.year = np.asarray(month), np.asarray(year)
   
#--------------------------------------------------------------------------------------

def load_nodes(path, data):
  filename_nodes = path_mesh + 'pickle_MOSAiC_transect'
  # look if nodes corresponding to MOSAiC trajectory has been store in a pickle file
  if os.path.exists(filename_nodes):
    print('The pickle file with MOSAiC mesh nodes exists, nodes with be loaded from this file')
    # load pickle file
    input_file = open(filename_nodes, "rb")
    nodes = pickle.load(input_file)
    input_file.close()
    
  else:
    print('The pickle file with MOSAiC mesh nodes does not exist but will be created')
    # extract transect and store data as pickle file
    lonlat = np.vstack((data.longitude, data.latitude))
    nodes = pf.transect_get_nodes(lonlat,mesh)
    output_file = open(filename_nodes, "wb")
    pickle.dump(nodes, output_file)
    output_file.close()

  return nodes


#--------------------------------------------------------------------------------------
def get_time(result_path, variable, years, runid="fesom", records=-1,ncfile=None, naming_convention="fesom", naming_template=None):
    """
    Get  output time from files.
    Parameters
    ----------
    result_path : string
        path to the data folder.
    variable : string
        variable name
    years : int, list
        year or list of years to open
    records: int, slice, list
        number of time steps to be considered for aggregation.
        If -1 (default), all timesteps will be taken in to account.
        If 0, only the first record will be taken
        If [0,5,7], only time steps with indexes 0,5 and 7 will be taken
        If slice(2,120,12), every 12th time step starting from the third one will be selected.
    depth: float
        The model depth closest to provided depth will be taken.
        If None, 3d field will be returned. Default = None.
    ncfile: str
        if provided, the netCDF file will be opened directly.
        Some dummy data have to be provided for result_path and years
    runid: str
        For esm-tools naming convention, use the experiment id name (e.g. test, PI-CTRL, LGM, ...)
    naming_convention : str
        The naming convention to be used. Can either be "fesom" for classic
        infrastructure, "esm-tools" for esm-tools infrastructure, or "custom",
        in which case a template string must be provided.
     naming_template : None or str
        Required if a customized naming convention is to be used. Replaced variables will be (in order) variable, runid, year.
     

    Returns
    -------
    data: xarray
       time.
    """
    paths = []
    if ncfile:
        paths = ncfile
    elif isinstance(years, int):
        if naming_convention == "fesom":
            fname = "{}.{}.{}.nc".format(variable, runid, years)
        elif naming_convention == "esm_tools":
            fname = f"{runid}.{years}.{variable}01.01.nc"
        elif naming_convention == "custom":
            fname = naming_template.format(variable, runid, years)
        else:
            raise ValueError(
                "You must have fesom, esm_tools, or custom as naming_convention!"
            )
        paths = os.path.join(result_path, fname)
    else:
        raise ValueError("year can be integer, list or one dimentional numpy array")

    dataset = xr.open_mfdataset(paths, combine="by_coords")
    time_data = np.asarray(dataset['time'])

    # transform time to python ordinal time  
    dates, time = [], []
    for t in time:
      tmp = datetime.strptime(str(t).split('.')[0], '%Y-%m-%dT%H:%M:%S')
      dates.append(tmp), time.append(tmp.toordinal())
    return np.asarray(dates), np.asarray(time)
#
class FESOM_variables:
  def __init__(self, filename, mesh):
    path_output = '/'.join(filename.split('/')[:-1])
    self.varname, self.year = filename.split('.')[0].split('/')[-1], int(filename.split('.')[-2])
   
    # extract output variable
    self.extract_variable(path_output, mesh)  
    
    # extract time
    self.extract_time(path_output)

  def extract_variable(self, path_output, mesh):
    # extract variable on user-defined mesh
    self.var = pf.get_data(path_output, self.varname, self.year, mesh, how='ori', compute=False) 

  def extract_time(self, path_output):
    # extract associated output time (python ordinal dates)
    self.dates = get_time(path_output, self.varname, self.year)

#--------------------------------------------------------------------------------------
# define list of output to extract transect from
def define_list_of_output_files(path):
  # variable of interest
  var_name = ['a_ice', 'Alk', 'benC', 'benN', 'benSi', 'pCO2','DetCalc', 'DetC', 'DetN', 'DetSi', 'DFe' , 'DiaC', 'DiaChl', 'DiaN', 'DiaSi', 'DIC', 'DIN', 'DOC' , 'DON' , 'HetC', 'HetN', 'idetz2calc', 'idetz2c', 'idetz2n', 'idetz2si', 'Kv', 'O2', 'PAR', 'pCO2s', 'PhyCalc', 'PhyC', 'PhyChl', 'PhyN', 'runoff', 'salt', 'temp', 'u', 'uice', 'v','vice','w', 'Zoo2C', 'Zoo2N']

  var_name = ['a_ice']
  # list of output file o be processed

  for vr in var_name:
    nm = vr +'.fesom.'
    try:
      filelist=np.hstack((filelist, gg(path + nm + '*')))
    except:
      filelist = gg(path + vr + '*')
  return filelist

# monthly outputs	
def extract_daily_transect_data_from_monthly_outputs(output, month, nodes):
  # extract daily data along the MOSAiC track from monthly FESOM outputs
  dates, dim = trajectory.dates, len(output.shape)
  for i, mnth in enumerate(month):
    # extract output along the track for the selected year
    tmp = output[int(mnth)-1, nodes[i],:] if dim>2 else output[int(mnth)-1, nodes[i]]
    transect = xr.concat([transect,tmp],dim='time') if i>0 else tmp
  return transect


# daily outputs
def change_time_axis(xarray, dates):
  # convert ordinal track-related dates to string
  times=[]
  for i, dt in enumerate(dates):
    tmp = datetime.fromordinal(int(dt))
    tmp = tmp.replace(hour=int(np.round(24.*(dt-float(int(dt))))))
    #tmp = datetime.strftime(tmp, '%Y-%m-%d %H:%M:%S')
    times = np.hstack((times,tmp)) if i>0 else tmp
  # change time values in xarray
  xarray['time'] = times
  return xarray
 
# extraction main routine
def extract_transect_outputs(fesom_outputs, trajectory_data, nodes, output_frequency='daily'):

  # extract FESOM outputs along the MOSAiC transect
  
  # initialization, crop trajectory dates according to the output year (cf. FESOM outputs)   
  var, year_output = fesom_outputs.var, fesom_outputs.year
  dates, year, month = trajectory_data.dates, trajectory_data.year, trajectory_data.month
  indices = np.where(year==year_output)[0]
  dates, nodes, month = dates[indices], nodes[indices], month[indices]

  # extract transect outputs
  if  var.shape[0]>12:
    print('daily FESOM outputs to be processed')
    if 'daily' in output_frequency:
      print('not coded yet') 
      transect_output = []
    elif 'monthly' in output_frequency:
      print('not coded yet')
      transect_output = []
  else:
    print ('monthly FESOM outputs to be processed')
    if 'daily' in output_frequency:
      # extract transect data
      transect_output = extract_daily_transect_data_from_monthly_outputs(var, month, nodes)
      # change time axis
      transect_output = change_time_axis(transect_output, dates)
    elif 'monthly' in output_frequency:
      print ('not coded yet')
      transect_output = []
  return transect_output
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
#
# load MOSAiC trajectory
#
path_trajectory = '/work/ollie/fbirrien/trajectories/'
filename_trajectory = path_trajectory + 'MOSAiC_track_20191005_20200731.nc'
trajectory = Trajectory(filename_trajectory)

#
# load FESOM fArc mesh
#
path_mesh='/work/ollie/fbirrien/mesh/farc/'
mesh = pf.load_mesh(path_mesh, abg=[50,15,-90])

#
# load mesh point indices that correspond to MOSAiC trajectory
#
nodes = load_nodes(path_mesh, trajectory)

#
# extract FESOM outputs along the MOSAiC track (make transects)
#
# load output
path_FESOM='/work/ollie/fbirrien/NuArctic/FESOM_Outputs/'
path_output = path_FESOM + 'Outputs/'
filelist = define_list_of_output_files(path_output)

for filename_FESOM in filelist:
  print ('variable & year to be processed', filename_FESOM.split('/')[-1].split('.')[0],  filename_FESOM.split('.')[-2])
  fesom_outputs = FESOM_variables(filename_FESOM, mesh)

  #
  # extract output along the MOSAiC transect
  # 
  transect_outputs = extract_transect_outputs(fesom_outputs, trajectory, nodes)
  #
  # store transect data
  #
  filename_transect = filename_FESOM.replace('/Outputs/','/Transects/').replace('.fesom.', '.transect.')
  transect_outputs.to_netcdf(filename_transect)

  


