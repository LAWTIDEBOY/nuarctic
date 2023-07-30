!$Id: next_observation_pdaf.F90 332 2019-12-30 09:37:03Z lnerger $
!>  Initialize information on next observation
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters
!!
!! The subroutine is called before each forecast phase
!! by PDAF_get_state. It has to initialize the number
!! of time steps until the next available observation
!! (nsteps) and the current model time (time). In
!! addition the exit flag (exit) has to be initialized.
!! It indicates if the data assimilation process is
!! completed such that the ensemble loop in the model
!! routine can be exited.
!!
!! The routine is called from PDAF_get_state by all processes
!!
!! Version for the 2D tutorial model.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

  USE mod_assimilation, &     ! Assimilation variables
       ONLY: step_null, step_assim, init_delt_obs, delt_obs
  USE mod_parallel_pdaf, &    ! Parallelization variables
       ONLY: mype_world
  USE obs_chla_pdafomi, &
      ONLY: chla_steps, dim_step


  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: stepnow  !< Number of the current time step
  INTEGER, INTENT(out) :: nsteps   !< Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   !< Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     !< Current model (physical) time
  ! Local variables 
  INTEGER   :: i, id 

! *******************************************************
! *** Set number of time steps until next observation ***
! *******************************************************
! C
!   IF (stepnow - step_null==0) THEN
!     step_assim = 0
!     nsteps = init_delt_obs  
!   ELSE
!     step_assim = step_assim + 1
!     nsteps = delt_obs
!   END IF




  id = 0
  Do i = 1, dim_step
    IF(chla_steps(i)==stepnow) id = i 
  END DO 

  IF (stepnow - step_null==0) THEN
    step_assim = 0
    nsteps = init_delt_obs        
    IF (mype_world == 0) WRITE (*, '(i7, 3x, a, i7)') & 
      stepnow, 'First observation at time step=', init_delt_obs
  ELSE
    step_assim = step_assim + 1
    nsteps = chla_steps(id+1) - chla_steps(id) 
    IF (mype_world == 0) WRITE (*, '(i7, 3x, a, i7)') &
      stepnow, 'Next observation at time step', stepnow + nsteps
  END IF

  IF (mype_world == 0) WRITE(*,*) 'step_assim= ', step_assim

  doexit = 0 ! Do not exit assimilation
  



END SUBROUTINE next_observation_pdaf
