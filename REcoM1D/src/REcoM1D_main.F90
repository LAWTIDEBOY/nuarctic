!=============================================================================!
!
!                 Lagrangian SCM version of REcoM2
!
!=============================================================================!
!                      Main driving routine
!=============================================================================!    
program main

! modules
use general_config
use mod_mesh
use REcoM_config
use REcoM_clock
use ocean_module
use ice_module
use atm_module
use forcing_module
use atm_deposition_module
use REcoM_init
use REcoM_setup
use REcoM_main

!=============================================================================!    
! local variables
type(t_mesh),             target, save :: mesh
integer				       :: istep
!=============================================================================! 
!
! initialization step
! 
! Read namelist
call read_namelist 

! create time axis
call clock_init

call get_run_steps(nsteps)


! Read mesh data (time and vertical grid)
call read_mesh(mesh)


! declare, allocate and initialize ocean arrays
call ocean_setup(mesh)

! declare, allocate and initialize ice arrays
call ice_setup

! declare, allocate and initialize atm arrays
call atm_setup

! initialize REcoM arrays
call recom_initialization(mesh)

! initial conditions

! Read (ice, turbulence, PAR, wind, ...) forcing
call read_forcing(mesh)

! read atm (CO2, Fe, N) deposition
call read_deposition(mesh)

! check if new year
call clock_newyear

!=============================================================================!
! 
! main loop of the model 
! 
do istep=1,nsteps
	! update time
	call clock
	! run main 
	call recom(istep,mesh)
	! update tracers
	! 1) perform mixing (due to turbulence)
	call recom_mixing
	! 2) vertical diffusion and sinking

	
enddo


!=============================================================================!
! 
! postprocessing diagnostics and deallocation steps
! 


end program
