cmake_minimum_required(VERSION 3.10)

# set default build type cache entry (do so before project(...) is called, which would create this cache entry on its own)
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "setting default build type: Release")
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif()

project(REcoM1D)
set(TOPLEVEL_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CRAY OFF CACHE BOOL "compile with cray ftn")
set(VERBOSE ON CACHE BOOL "toggle debug output")
set(DEBUGGING ON CACHE BOOL "debbuging compiling options")
add_subdirectory(src)
