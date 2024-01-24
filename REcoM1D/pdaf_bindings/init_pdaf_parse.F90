!$Id: init_pdaf_parse.F90 1105 2023-02-18 18:11:14Z lnerger $
!>  Parse command line options for PDAF
!!
!! This routine calls the command line parser to initialize
!! variables for the data assimilation with PDAF.
!!
!! Using the parser is optional and shows one possibility
!! to modify the variables of the compiled program. An
!! alternative to this might be Fortran namelist files.
!!
!! __Revision history:__
!! * 2011-15 - Lars Nerger - Initial code extracted from init_pdaf
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf_parse()

  USE parser, &           ! Parser function
      ONLY: parse
  USE mod_parallel_pdaf, &
      ONLY: mype_world
  USE mod_assimilation, & ! Variables for assimilation
      ONLY: dim_ens, screen, filtertype, subtype, type_trans,       &
      type_forget, forget, type_sqrt, perturb_scale,                & 
      init_delt_obs, delt_obs, write_ens, bgc_layer,                &
      parameter_estimation

  USE obs_chla_pdafomi, &
      ONLY: assim_chla, rms_obs_chla
  USE obs_din_pdafomi, &
      only: assim_din, rms_obs_din
  USE obs_dsi_pdafomi, &
      only: assim_dsi, rms_obs_dsi



  IMPLICIT NONE

! *** Local variables ***
  CHARACTER(len=32) :: handle  ! handle for command line parser

! variable in the namelist file
  NAMELIST  /pdaf_nml/dim_ens, screen, filtertype, subtype, type_trans,       &
      type_forget, forget, type_sqrt, bgc_layer, perturb_scale,               &
      init_delt_obs, delt_obs, write_ens, assim_din, rms_obs_din, assim_dsi,  &
      rms_obs_dsi, assim_chla, rms_obs_chla, parameter_estimation
      
  ! *** Read namelist for PDAF config ***
  IF (mype_world==0) WRITE (*,'(5x,a)') 'Read namelist for PDAF configuration'
  OPEN (10,file='../config/pdaf.nml')
  READ (10,NML=pdaf_nml)
  CLOSE (10)

      
! **********************************
! *** Parse command line options ***
! **********************************
  ! General settings for PDAF
  handle = 'dim_ens'                 ! set ensemble size/rank of covar matrix
  CALL parse(handle, dim_ens)
  handle = 'screen'                  ! set verbosity of PDAF
  CALL parse(handle, screen)
  handle = 'filtertype'              ! Choose filter algorithm
  CALL parse(handle, filtertype)
  handle = 'subtype'                 ! Set subtype of filter
  CALL parse(handle, subtype)

  ! Filter-specific settings
  handle = 'type_trans'              ! Type of ensemble transformation in SEIK/ETKF/LSEIK/LETKF
  CALL parse(handle, type_trans)
  handle = 'type_forget'             ! Set type of forgetting factor
  CALL parse(handle, type_forget)
  handle = 'forget'                  ! Set forgetting factor
  CALL parse(handle,forget)
  handle = 'type_sqrt'               ! Set type of transformation square-root (SEIK-sub4, ESTKF)
  CALL parse(handle, type_sqrt)



END SUBROUTINE init_pdaf_parse
