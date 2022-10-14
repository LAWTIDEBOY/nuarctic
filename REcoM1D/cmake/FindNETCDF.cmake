# NETCDF_Fortran_INCLUDE_DIRECTORIES
# NETCDF_Fortran_LIBRARIES

# include file netcdf.inc and search symbol e.g. NCCRE for Fortran 77 NetCDF
# include file netcdf.mod and search symbol e.g. nf90_create for Fortran 90 NetCDF (this symbol is not found on hlrn with intel ftn)
if(CMAKE_Fortran_COMPILER_LOADED)
	include(CheckFortranFunctionExists)
	check_fortran_function_exists(NCCRE HAVE_Fortran_NETCDF)
	
	if(HAVE_Fortran_NETCDF)
		set(NETCDF_Fortran_INCLUDE_DIRECTORIES "")
		set(NETCDF_Fortran_LIBRARIES "")
	else()
		find_path(NETCDF_Fortran_INCLUDE_DIRECTORIES netcdf.inc HINTS $ENV{NETCDF_DIR}/include ENV NETCDF_Fortran_INCLUDE_DIRECTORIES)
		find_library(NETCDF_Fortran_LIBRARIES netcdff HINTS ${NETCDF_Fortran_INCLUDE_DIRECTORIES}/../lib)
	endif()	
endif()

