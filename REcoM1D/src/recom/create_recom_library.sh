#!/bin/bash 
#
# script that creates a library for recom using all the functional Fortran routines 
# created using:makedepf90 -coco -o foobar *.F90 > .depend, which summarize all fortran dependencies (Makefile oriented)
#

# name of the library
library_name=recomlib.so

# list of all the functional routines
path=/home/fbirrien/NuArctic/nuarctic/REcoM1D/src/recom

#
# compile individual Fortran files based on dependencies
#
# no dependencies
gfortran -fpic -c singledouble.F90 gsw_mod_kinds.F90 
# 1 level dependency
gfortran -fpic -c f2pCO2.F90 p80.F90 DNAD.F90 p2fCO2.F90
gfortran -fpic -c gsw_add_mean.F90 gsw_mod_baltic_data.F90 gsw_mod_error_functions.F90 gsw_mod_saar_data.F90
gfortran -fpic -c gsw_mod_specvol_coefficients.F90 gsw_mod_teos10_constants.F90 gsw_mod_toolbox.F90 gsw_util_indx.F90
# 2 level dependency
gfortran -fpic -c depth2press.F90 eos.F90 phsolvers.F90 rho.F90
gfortran -fpic -c gsw_add_barrier.F90 gsw_ct_from_pt.F90 gsw_ct_from_t.F90 gsw_entropy_part.F90 
gfortran -fpic -c gsw_entropy_part_zerop.F90 gsw_gibbs.F90 gsw_gibbs_pt0_pt0.F90 gsw_pt0_from_t.F90
gfortran -fpic -c gsw_pt_from_ct.F90 gsw_pt_from_t.F90 gsw_rho.F90 gsw_saar.F90 gsw_sa_from_sp_baltic.F90
gfortran -fpic -c gsw_sa_from_sp.F90 gsw_specvol.F90 gsw_sp_from_sa_baltic.F90 gsw_sp_from_sa.F90 
gfortran -fpic -c gsw_t_from_ct.F90 gsw_util_xinterp1.F90
gfortran -fpic -c sw_adtg.F90
# 3 level dependency
gfortran -fpic -c rhoinsitu.F90 sw_ptmp.F90 
# 4 level dependency
gfortran -fpic -c sw_temp.F90 tpot.F90 varsolver.F90 
# 5 level dependency
gfortran -fpic -c constants.F90 tis.F90
# 6 level dependency
gfortran -fpic -c derivauto.F90 vars.F90 
# 7 level dependency
gfortran -fpic -c buffesm.F90 derivnum.F90 gasx.F90
# 8 (last) level dependency
gfortran -fpic -c errors.F90

#
# link objects into library
#

# list all the objects .o files created by compiling 
list=$(ls $path/*.o)
# link all the objects .o into a dedicated dynamic library
gfortran -shared -o $library_name $list


#
# clean directory
#
#rm *.o *.mod
