cmake_minimum_required(VERSION 3.10)
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")

# project name and computing language in use
project(REcoM1D Fortran)

# define work directory from which REcoM1D is compiled from (containing CMakeList.txt)
set(src_home ${CMAKE_CURRENT_LIST_DIR}) 
set(bin_dir  ${CMAKE_CURRENT_LIST_DIR}/../bin)
#message(STATUS "bin dir: ${bin_dir}")

# define path to and gather source file
file(GLOB sources_Fortran ${src_home}/*.F90  
                          ${src_home}/recom/*.F90)

add_subdirectory(${src_home}/async_threads_cpp)

# include netcdf libraries                          
include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindNETCDF.cmake)
                 
# create executable
add_executable(${PROJECT_NAME} ${sources_Fortran})

# precision definition
target_compile_definitions(${PROJECT_NAME} PRIVATE USE_PRECISION=2)

# add verbose, i.e. comments during simulation
if(${VERBOSE})
   target_compile_definitions(${PROJECT_NAME} PRIVATE VERBOSE)
endif()

# define the compiling option according to the compiler in use
# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
	target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip -init=zero -no-wrap-margin)

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
	# allow debugging
	if(${DEBUGGING})
		target_compile_options(${PROJECT_NAME} PRIVATE -O0 -g -fbacktrace -ffree-line-length-none)
	else()
	# optimized compiling
		target_compile_options(${PROJECT_NAME} PRIVATE -O3 -finit-local-zero  -finline-functions -march=native -fimplicit-none  -fdefault-real-8 -ffree-line-length-none)
	endif()
	
	if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10 )
		target_compile_options(${PROJECT_NAME} PRIVATE -fallow-argument-mismatch)
	endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray )
       target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64)
endif()

#librairies
# include path to required librairies (netcdf,...) 
target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES})
# link to librairies
target_link_libraries(${PROJECT_NAME} ${NETCDF_Fortran_LIBRARIES})
target_link_libraries(${PROJECT_NAME} async_threads_cpp)

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

# executable
# set path and
#set(RECOM1D_INSTALL_FILEPATH "${CMAKE_CURRENT_LIST_DIR}/../bin/REcoM1D.exe" CACHE FILEPATH "file path to executable required")
set(RECOM1D_INSTALL_FILEPATH "${bin_dir}/REcoM1D.exe" CACHE FILEPATH "file path to executable required")
get_filename_component(RECOM1D_INSTALL_PATH  ${RECOM1D_INSTALL_FILEPATH} DIRECTORY)
get_filename_component(RECOM1D_INSTALL_NAME ${RECOM1D_INSTALL_FILEPATH} NAME)
#create executable
install(PROGRAMS ${PROJECT_BINARY_DIR}/${PROJECT_NAME} DESTINATION ${RECOM1D_INSTALL_PATH} RENAME ${RECOM1D_INSTALL_NAME})
