!$Id: mod_obs_pdaf.F90 nmamnun $
!> callback_obs_pdafomi
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specific routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!! When adding an observation type, one has to add one module
!! obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
!! In addition one has to add a call to the different routines include
!! in this file. It is recommended to keep the order of the calls
!! consistent over all files.
!!
!! __Revision history:__
!! * Later revisions - see repository log
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

  ! ! Include functions for different observations
  USE obs_chla_pdafomi, ONLY: assim_chla, init_dim_obs_chla
  USE obs_din_pdafomi, ONLY: assim_din, init_dim_obs_din
  USE obs_dsi_pdafomi, ONLY: assim_dsi, init_dim_obs_dsi

  USE PDAFomi,         ONLY: PDAFomi_set_debug_flag

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_chla ! Observation dimensions
  INTEGER :: dim_obs_din ! Observation dimensions
  INTEGER :: dim_obs_dsi ! Observation dimensions

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************
  ! CALL PDAFomi_set_debug_flag(1)
  ! Initialize number of observations
  dim_obs_chla = 0
  dim_obs_din = 0
  dim_obs_dsi = 0
  ! dim_obs_B = 0
  ! dim_obs_C = 0

  ! ! Call observation-specific routines
  ! ! The routines are independent, so it is not relevant
  ! ! in which order they are called
  IF (assim_chla) CALL init_dim_obs_chla(step, dim_obs_chla)
  IF (assim_din) CALL init_dim_obs_din(step, dim_obs_din)
  IF (assim_dsi) CALL init_dim_obs_dsi(step, dim_obs_dsi)

  dim_obs = dim_obs_chla + dim_obs_din + dim_obs_dsi ! + dim_obs_B + dim_obs_C

END SUBROUTINE init_dim_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_TYPE.
!!
SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! ! Include functions for different observations
  USE obs_chla_pdafomi, ONLY: obs_op_chla
  USE obs_din_pdafomi, ONLY: obs_op_din
  USE obs_dsi_pdafomi, ONLY: obs_op_dsi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
  REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL, INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi

  CALL obs_op_chla(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_din(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_dsi(dim_p, dim_obs, state_p, ostate)

END SUBROUTINE obs_op_pdafomi



