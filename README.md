# <span style="color:blue"> :large_blue_circle: nuArctic :large_blue_circle: </span>

<p align="center">
    <img src="logos/logo4_whitebackground.png" alt="Drawing" style="width: 300px;">
</p>

The nuArctic project is a BMBF project that aims at exploiting data aquired during the MOSAiC expedition to improve / implement biogeochemical models with regards to the remineralization / export of organic material (nutrients and carbon).

nuArctic develops a 1D version of the biogeochemical model REcoM, combined with the assimilation framework PDAF and field data from MOSAiC.

contact:
@LaurentOziel: laurent.oziel@awi.de
@FlorentBirrien: florent.birrien@awi.de


### Table of Contents

* [0. Where to get it from](#0.)
* [1. Required pre-processing steps](#1.)
    * [1.1 How to define the time frame of the simulation](#1.1)
    * [1.2 mesh creation and configuration](#1.2)
    * [1.3 initial tracer conditions](#1.3)
    * [1.4 Atmospheric deposition](#1.4)
    * [1.5 In-situ forcing](#1.5)
     * [1.5.1 Photo-synthetically Active Radiation (PAR)] (#1.5.1)
     * [1.5.2 Temperature and salinity] (#1.5.2)
     * [1.5.3 Vertical diffusion Kz] (#1.5.3)
     * [1.5.4 Atmospheric conditions] (#1.5.4)
     * [1.5.5 Ice cover] (#1.5.5)
     * [1.5.6 Merging of all forcings in one file] (#1.5.6)
* [2. Composition of the source code](#2.)
    * [2.1 What’s new?](#2.1)
    * [2.2 Composition of the difference source code files](#2.2)
* [3. How to compile REcoM1D?](#3.)
    * [3.1 Environment and specific modules](#3.1)
    * [3.2 Compiling with Cmake](#3.1)
    * [3.3 The compiling command](#3.1)
* [4. How to run REcoM1D?](#4.)
    * [4.1 Environment and path](#4.1)
    * [4.2 Compiling with Cmake](#4.2)
    * [4.3 Namelist files](#4.3)
    * [4.4 How to run the model?](#4.4)
    * [4.5 Outputs](#4.5)

    

# <span style="color:blue">Getting started with REcoM1D</span>

## <span style="color:blue">0. Where to get it from:</span><a class="anchor" id="0."></a>

REcoM1D can be cloned/downloaded from the Github repository:
https://github.com/LAWTIDEBOY/nuarctic

The ‘nuarctic’ repository is composed of 3 main folders: 
- data and scripts, which include data and scripts to perform the pre-processing steps
- REcoM1D, which contains all the source code and scripts to configure, compile and run the model.

## <span style="color:blue">1. Required pre-processing steps:</span><a class="anchor" id="1."></a>

All the pre-processing steps prepare and format the various respective data set to be used by REcoM1D. It is performed beforehand, i.e. not directly in the model (as in the coupled version with FESOM), to avoid overloading simulation runs 
The main steps consists in providing:
- mesh configuration (folder MESH and TRACK),
- initial tracer conditions (folder initialization),
- atmospheric deposition (Carbon C, Nitrogen N and Iron Fe,  in folder atm_deposition)
- forcing, which is here derived from in-situ observations (folders ITP_CTD_Ben, ITP_PAR, MSS, and METEO, ...).

#### <span style="color:red"> :small_red_triangle_down: remarks:</span>
- Paths do not need to be changed in the following scripts, they are defined according to the current cloned/downloaded repository architecture.
- the following extraction algorithms are based on the nearest neighbor approximation to estimate the track/mesh correspondence.

### <span style="color:blue">1.1. How to define the time frame of the simulation: the file time.recom</span><a class="anchor" id="1.1"></a>
The file time.recom (in folder REcoM1D/config/) defines the time frame of the simulation. All variables need to be tuned by the user according to simulation time requirement. The time configuration file also includes a spin-up variable, which is most of the time required in simulations.

The various variables in time.recom are:
- start_date : effective starting date of the simulation (excluding spin-up)
- end_date :  last date of simulation
- dt: simulation time step (in hours)
- spinup: spin-up duration (in days)

#### <span style="color:red">Important:</span> 
- Whenever changes are performed for any of these time variables all the pre-processing steps need to be run once more, starting with the mesh scripts , to have consistent mesh/forcing data sets for the simulation.
- Time variables from time.recom are systematically overwriting default values in the simulation namelists (namelist.config) whenever the model is run (see appropriate section below) in order to synchronize the changes.

### <span style="color:blue">1.2. Mesh creation and configuration</span><a class="anchor" id="1.2"></a>
Folder: scripts/Create_Mesh/
Jupyter python notebook: Create_REcoM1D_mesh.ipynb
Data sets required: 	- observed instrument/vessel track
			- predefined mesh (for example FESOM Farc or Core2)
Output files:  data/MESH/REcoM1D_mesh.nc and data/MESH/REcoM1D_daily_mesh.nc

The scripts is dedicated to the extraction of the vertical mesh along the observed track using a predefined mesh. Depth evolves along the track. To account for this change that affects vertical discretization in Lagrangian simulations, the model allows for changing number of effective vertical level through time. Vertical mesh discretization is subsequently derived at each track point.

The script needs to be run to twice, switching the dedicate flag ‘flag_daily’ from True or False, according to the needs:
- A daily mesh is required for further processing of the forcing variables that are first derived on a daily base.
whereas
- a discretized mesh si required to format all datasets, according to the simulation time frame, in order to create the final model forcing file.

#### <span style="color:red">:small_red_triangle_down: Remark:</span> The mesh-related script needs to be run first because all the other pre-processing steps are based on along-track extracted mesh.

### <span style="color:blue">1.3. Initial tracer conditions</span><a class="anchor" id="1.3"></a>
Folder: scripts/Create_Init/)
Jupyter python notebook: 
Create_init_climatology.ipynb
Create_init_MOSAiC.ipynb
Output files: data/initialization/tracer_initialization.nc

The model offers the choice to intialize the initial contions of REcoM tracers (DO2, DAlk, DIC, DIN, DSi, DFe) either from a climatology (World Ocean Atlas 2018) or from in situ profiles for example from the MOSAiC expedition.

### <span style="color:blue">1.4. Atmospheric deposition</span><a class="anchor" id="1.4"></a>
Folder: scripts/atm_deposition/
Jupyter python notebook: Create_atm_deposition_forcing.ipynb
Data sets required: - REcoM discretized mesh REcoM1D_mesh.nc
		       - 3D mesh (farc/CORE2) related deposition data sets
		       - Carbon: MonthlyAtmCO2_gcb2021.nc (monthly global data)
		       - Nitrogen: AeolianNitrogenDep.nc (gridded monthly data sets, data from 2009)
		       -  Iron:  DustClimMonthlyAlbani.nc (default, gridded monthly climatology, other datasets are available).
Output file:  data/atm_deposition/atm_deposition.nc


Atmospheric deposition of iron, carbon and nitrogen are used by the model whenever open water conditions appears.

For each component (C, N, Fe), the deposition is extracted along the track accounting for data discretization in space (global vs gridded) and time (mainly monthly data).

### <span style="color:blue">1.5. In-situ forcing</span><a class="anchor" id="1.5"></a>

The processing and formatting of the in-situ forcing gathers several observational data sets from different sources/instruments. Each data set required its own specific processing script. Here we classify the processing scripts available according to the observed variables that are required for the forcing.

#### <span style="color:blue">1.5.1. Photo-synthetically Active Radiation (PAR):</span><a class="anchor" id="1.5.1"></a>

The PAR can be given to the model in three different ways:
1) above water shortwave radiations : in which case the model derives to air-to-sea transmission and will eventually require snow and sea-ice concentration to derive attenuation.
2) surface water PAR (defaut) : in which case the model will derive attenuation at depth.
3) Water column profile PAR : in which case the model will not calculate anything.

The model automatically "understand" which case by reading the variable name in the forcing:
if variable's name  = 'shortwave' (case 1),
if variable's name  = 'PAR_surface' (case 2)
if variable's name  = 'PAR' (case 3),

#### <span style="color:blue">1.5.2. Temperature and salinity:</span><a class="anchor" id="1.5.2"></a>

Folder: scripts/ITP_CTD/
Jupyter python notebook: Process_CTD_data.ipynb
Data sets required: - REcoM daily mesh : data/MESH/REcoM1D_daily_mesh.nc
		      - Polarstern track (for verification): data/TRACK/Polarstern_daily_track.nc
		      - files containing the ITP measurements:
				data/ITP_Ben/mosaic.whoiitpmerged.newgrid.mat  
Output file: data/ITP_CTD_Ben/ITP_daily_data.nc

The script processes the data from ITP 111 (cf. Rabe et al., 2021, Ocean MOSAiC paper) using a 3 steps approach:
- Salinity and temperature profiles are extrapolated to the surface and high depth.
- Outliers profiles are processed using a temporal interpolation.
- A moving average filter is performed (per level) to smooth day-to_day strong variations
 
The resulting file gathers the time series of daily salinity and temperature profiles

When other data sets are available for temperature and salinity profiles time series an other script
needs to be written.

#### <span style="color:blue">1.5.3. Vertical diffusion Kz:</span><a class="anchor" id="1.5.3"></a>

Folder: scripts/MSS/
Jupyter python notebook: Process_Kz_daily.ipynb
Data sets required: - REcoM daily mesh : data/MESH/REcoM1D_daily_mesh.nc
		       - files containing the MSS measurements:
				data/MSS/MOSAiC_MSS_daily_avg_Kz.nc
Output file: data/MSS/MOSAiC_MSS_daily_Kz.nc

The script interpolate Kz on a fixed vertical axis and filled the missing values with 1e-7 at the bottom. The resulting Kz is a daily vertically discretized field.    

#### <span style="color:blue">1.5.4. Atmospheric conditions:</span><a class="anchor" id="1.5.4"></a>

atmospheric conditions consist of wind and surface pressure data and have been extracted along the track from ERA5 hourly reanalysis data. The file data/METEO/ERA5_forcing_Polarstern.nc gathers all this information.

#### <span style="color:blue">1.5.5. Ice cover:</span><a class="anchor" id="1.5.5"></a>

Ice cover (or fraction) has (for now) been set to 1, i.e, cell fully covered by ice. If in-situ ice and snow data are available, they can be easily added by writing a specific pre-preprocessing script and by adding some lines to the . When doing so, think also of changing some lines related to ice cover in the source code (cf. Shortwave transmission, etc...)

#### <span style="color:blue">1.5.6. Merging of all forcings in one file:</span><a class="anchor" id="1.5.6"></a>

The REcoM forcing specific script gathers all information about forcing and interpolate the data according to the simulation time frame (time.recom). All data are then stored in REcoM1D_forcing.nc to be directly used as forcing by the model.

REcoM forcing file:
Folder: scripts/Create_Forcing/
Jupyter python notebook:Process_Forcing_REcoM.ipynb
Data sets required: 
- REcoM discretized mesh : data/MESH/REcoM1D_mesh.nc
- PAR observations: data/ITP_PAR/PAR_forcing.nc
- Atm reanalysis data: data/METEO/ ERA5_forcing_Polarstern.nc
- Temperature and salinity profiles: data/ITP_CTD_Ben/ITP_daily_data.nc
- Vertical diffusion profiles: data/MSS/MOSAiC_MSS_daily_Kz.nc
- Ice coverage (set to 1 for now).

Output File: data/REcoM_forcing_data/REcoM1D_forcing.nc

## <span style="color:blue">2. Composition of the source code:</span><a class="anchor" id="2."></a>
All code source are contained in the folder: nuarctic/REcoM1D/src/

### <span style="color:blue">2.1. What’s new?</span><a class="anchor" id="2.1"></a>
REcoM1D has been adapted from the original REcoM model and keep most of its architecture.
Some of the modules have been rewritten accounting for Lagrangian requirements (both space and time evolution). Some modules dedicated to in-situ forcing and describing tracers evolution (inherited from FESOM) have been added to complete the stand-alone 1D model. Some variables have been rewritten according to 1D Lagrangian configuration. The different Fortran files have been reorganized as dedicated modules, except REcoM1D_main.F90 which is the main running program.

### <span style="color:blue">2.2. Composition of the difference source code files</span><a class="anchor" id="2.2"></a>

The structure of REcoM remains mainly unchanged compare to the 3D version:

https://recom.readthedocs.io/en/latest/intro.html
https://github.com/FESOM/fesom2/tree/fesom2.1_recom

However 3 files have been added:

1) recom_clock.F90 (added) managed simulation time properties
2) forcing_modules.F90 (added) gathers all modules related to model forcing (ocean, ice and atmospheric components) as well for atmospheric deposition. The different atm, ocean and ice dedicated modules deal with variables (including tracers) declaration, allocation and initialization. The forcing modules itself provides a routine to read global forcing from netcdf files (cf. Forcing preprocessing) and redistribute forcing at each time step. The atmospheric deposition dedicated module deals with deposition variables initialization, reading from pre-processed files and redistribution at each time step.
3) recom_tracers.F90 (added) gathers routines for ocean tracer update that were originally located in FESOM  oce_ale_tracer.F90 and recom_sinking.F90. The mixing scheme has been adapted from the model MOPS. Tracer evolution includes vertical mixing, sinking and remineralization of tracers in the water column.

2 files have been intensely modified:

1) recom_diagnostics.F90 (rewritten) estimates diagnostics and write 1D and 2D outputs and diagnostics in netcdf file for post-processing.
2) REcoM1D_main.F90 (program) coordinates all steps of the simulation: initialization, forcing reading, main loop and outputs writing.

The other source code files are identical with 3D REcoM:

- recom_modules.F90 gathers basis modules for constant and namelists definition as well as local and  global (cf. FESOM code) variable declarations.
- recom_init.F90 (partially rewritten) is dedicated to the declaration, allocation and initialization of original REcoM global/local variables and initialization of tracers.

#### <span style="color:red"> :small_red_triangle_down: remark</span>: most tracers are initialized with very small values except DIN(tracer 3), DIC(tracer 4), Alk(tracer 5), Dsi(tracer 20), Dfe(tracer 21), O2(tracer 24) which are initialized from climatology or in-situ initial conditions.</span>

- recom_setup.F90 performs the reading of the namelists and mesh and estimate time and vertical mesh properties at each time step.
- recom_extra.F90 gathers REcoM useful computation tool routines.
- recom_sms.F90 gathers the model equation and sources minus sink of state variables.
- recom_forcing.F90 calculates some additional diagnostics and handle CO2 fluxes
- recom_main.F90 coordinates the main loop of the computation: gathers forcing and deposition at each time step and performs subsequent computation.

## <span style="color:blue">3. How to compile REcoM1D?</span><a class="anchor" id="3."></a>

### <span style="color:blue">3.1. environment and specific modules</span><a class="anchor" id="3.1"></a>
The work environment is defined using the script env.sh. The script recognizes the machine with which the user is working (local AWI laptop, Ollie, hlrn, etc….). The script then loads the required environment-related modules using the file shell located in the folder env/ollie/ (for example if Ollie is used, other environment are available).
shell gathers command lines to load modules related to processing (python) and compiling (cmake, compiler, netcdf libraries) and defines compiler to be used.

#### <span style="color:red"> :small_red_triangle_down: Remarks:</span>
- only 2 environments are for now provided (Ollie, AWI laptop), the user might create environment files according to its needs and the available modules (command module avail to know more about modules that are available on each machine). An example, which any user can follow is available for Ollie.
- 2 or more compilers can be available for one machine (GNU, intel), the user might choose 1 and its related compiling modules and libraries, shell needs be rewritten according to the changes done by the user.
- The version of  the compilers often come with related libraries and options. Therefore, the compatibility in between modules needs to be checked.

### <span style="color:blue">3.2. Compiling with Cmake</span><a class="anchor" id="3.2"></a>
 The source is compiled using Cmake and 2 CmakeLists.txt files, which defines path, executable file name, code dependency, libraries and compiling & linking options.
The main CmakeLists.txt is located in the REcoM1D folder directly, compiling options such as verbose or debugging mode can be tuned on/off according to user needs and for verification (if bugs occur). More options can be added in the file in the same way. 
This Cmake scripts calls another CmakeLists.txt file located in the  REcoM1D/src. In the latter file, path, executable file name, code dependency, libraries and compiling & linking options are effectively defined. Compiling options are specific to each compiler.

### <span style="color:blue">3.3. The compiling command</span><a class="anchor" id="3.3"></a>
When all those optional checking steps are performed (according to user/machine requirements), the code is simply compiled using the command line: ./configure.sh, which is related to the configuration scripts.

The scripts performs the compiling and create to subfolder:
- bin where the executable is stored.
- build which gathers all details of the compiling steps

When no error message or warning, the model is ready to run. 

## <span style="color:blue">4. How to run REcoM1D?</span><a class="anchor" id="4."></a>

### <span style="color:blue">4.1. Environment and path</span><a class="anchor" id="4.1"></a>
Simulation paths are exported using the script set_path_REcoM.sh. 
The folders that contains required data sets are also created here, if they do not already exist.
Remark: 
No change is required in set_path_REcoM.sh. The different paths are related to the standard model architecture.

### <span style="color:blue">4.2. Namelist files</span><a class="anchor" id="4.2"></a>
In addition to time.recom, 2 namelist files (in REcoM1D/config/) are required to run a simulation:
namelist.config gathers all information about simulation set-up time discretization, forcing and data file names and output frequencies. Time information from the time.recom systematically overwrite original time information from the namelist, as a matter of synchronization.
namelist.recom gathers all usual REcoM1D parameters which deals with simulation options and physical processes.

### <span style="color:blue">4.3. Data sets</span><a class="anchor" id="4.3"></a>
Several observation data sets are required to run the model (cf. pre-processing step).  The data are not stored in the different folders but a symbolic link (command ln -s ) is created that points in the respective data directories:
- The discretized mesh REcoM1D_mesh.nc  (in REcoM/grid/) points to the data/MESH/ folder.
-  Atmospheric deposition information (atm_deposition.nc) and tracer initialization files (tracer_initialization.nc) in the REcoM/data/ are pointing to data/atm_deposition/ and data/initialization/ respectively.
-  The files containing all the model forcing REcoM1D_forcing.nc (in REcoM1D/forcing/) points to the folder data/REcoM_forcing_data/.

### <span style="color:blue">4.4. How to run the model?</span><a class="anchor" id="4.4"></a>
The model is run simply by launching the script ./run_REcoM1D.sh. 
Result files are stored in the  REcoM/results folders and the simulation summary REcoM1D.out is stored in RecoM/bin/.

### <span style="color:blue">4.5. Outputs</span><a class="anchor" id="4.5"></a>
REcoM1d_diagnostics.nc and REcoM1d_outputs.nc are the 2 output files which gathers simulation diagnostics and 1D & 2D fields and state variables respectively.
