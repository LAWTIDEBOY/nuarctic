Implementing Ensemble Data Assimilation
---

### 1. Obtaining and Compiling the PDAF Source Code

To begin implementing Ensemble Data Assimilation into the REcoM1D model, You first need to obtain the source code for PDAF (Parallel Data Assimilation Framework) and compiling it for use. To get the most recent version of PDAF, use the following command in your terminal:
`git clone https://github.com/PDAF/PDAF.git`

Before compiling, you need to set an environment variable (PDAF_ARCH) specific to your system architecture. On 'albedo', I use linux_mpiifort_impi. To set the environment variable, execute:
`export PDAF_ARCH=linux_mpiifort_impi`

Change your directory to the source code of PDAF and initiate the compilation process:
`cd PDAF/src`
`make`

For additional guidance on compiling PDAF, visit the PDAF First Steps webpage (http://pdaf.awi.de/trac/wiki/FirstSteps). 



### 2. Cloning the PDAF Bindings Branch and Accessing

You first need to clone the specific branch of the repository that contains the PDAF bindings (pdaf_bindings). You can do this by executing the following commands in your terminal:

`git clone -b pdaf_bindings git@github.com:LAWTIDEBOY/nuarctic.git`

Then change your current directory to `./nuarctic/REcoM1D`
`cd nuarctic//REcoM1D/`

### 3. Copying Files from the PDAF Bindings Directory to the Source Directory

The next step involves copying the contents from the `./pdaf_bindings` directory to the `./src` directory.
`cp -r ./pdaf_bindings/* ./src/`

**Important Note:** The file REcoM1D_main.F90 in the ./pdaf_bindings directory will replace the existing REcoM1D_main.F90 in the ./src directory. The rest of the files from ./pdaf_bindings are user-supplied routines specifically for PDAF<!--, except `shell_pdaf`-->.

### 4. Adding PDAF Library to Environment and Build Configuration
To integrate the PDAF library with your REcoM1D model, you need to update your environment settings and modify the build configuration. This involves editing the environment file and the CMake configuration.

Add the following lines to `./env/albedo/shell`.
`NC_LIB="-L${NETCDFFROOT}/lib -lnetcdff"`
`NC_INC="-I${NETCDFFROOT}/include"`
`LIBS="$LIBS -L/path/to/PDAF/directory/lib -lpdaf-d -mkl"`

Replace `/path/to/PDAF/directory` with the actual path where you have cloned the PDAF source code. 

Navigate to `nuarctic/REcoM1D/src` and open `CMakeLists.txt`. Modify this file to include PDAF in the build process: 

* Add `-I/path/to/PDAF/directory/include` to `target_compile_options(${PROJECT_NAME}`.
* Include `-I/path/to/PDAF/directory/include` in `target_include_directories(${PROJECT_NAME}`.
* Add `-L/path/to/PDAF/directory/lib/ -lpdaf-d -DUSE_PDAF` to `target_link_libraries(${PROJECT_NAME}`.

Again, ensure that `/path/to/PDAF/directory` is replaced with the correct path to your PDAF installation.

For an example of how these configurations should look, refer to the file `nuarctic/REcoM1D/src/CMakeLists_pdaf.txt`.

### 5. Compiling the REcoM1D Model with PDAF Bindings

Once you have set up the environment and updated the build configuration, the next step is to compile the REcoM1D model with the PDAF bindings. This is achieved by running the existing configuration script. Navigate to the directory `nuarctic/REcoM1D` and locate the script named `configure.sh`. Execute the script by running the following command:
`./configure.sh` or 
`bash configure.sh`

### 6. Running the REcoM1D Model with the Job Script

After successfully compiling the REcoM1D model with PDAF bindings, the next step is to run the model. Navigate to the directory `nuarctic/REcoM1D` and locate the job script named `run_REcoM1D_mpmd.slurm`. This script is used for submitting the job to AWI albedo mechine. Ensure the `#SBATCH --ntasks` directive is set to `dim_ens`. This setting specifies the number of tasks to be used that is 1 task per ensemble member. Similarly, set the variable `NENS` to `dim_ens` within the script.

Once the script is correctly configured, submit the job to the SLURM scheduler using the following command:
`sbatch run_REcoM1D_mpmd.slurm`

Make sure that you have Python3 and Numpy available, either through a Conda environment or directly from your Python installation. To use a Conda environment, activate it within the `nuarctic/REcoM1D/config/scripts/REcoM_time.sh` script.


### Configuring Parameter Selection 
You can configure which parameters to perturb and estimate by a namelist file. Navigate to the directory `nuarctic/REcoM1D/config` and find the file named `parameter_selection.nml`. 

The values assigned to the `is_[parameter]` in this file indicate their selection status. A value of 0 means the parameter is not selected and will not be perturbed. Any integer greater than 0 indicates the parameter is selected for perturbation and estimation.

### Configuring Data Assimilation Options
`nuarctic/REcoM1D/config/pdaf.nml` contains various options that control what observations are assimilated and other PDAF options. `assim_[obs]` TRUE/FALSE specify whether the observation is assimilated. 

The `forget` option configures the degree of ensemble inflation. By default, it is set to 1, which means no inflation. Adjusting this value below 1 increases the inflation rate. The smaller the `forget` value, the higher the inflation. However, setting it below 0.9 is generally not recommended as it can lead to overly inflated ensembles.

