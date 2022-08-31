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
use REcoM_forcing
use REcoM_tracers
use REcoM_diagnostics
use REcoM_main

!=============================================================================!    
! local variables
type(t_mesh),             target, save :: mesh
integer				       :: istep
!=============================================================================! 
print*, 'start initialization'
!
! initialization step
! 
! Read namelist
call read_namelist 

! create time axis
call clock_init

call get_run_steps

! Read mesh data (time and vertical grid)
call read_mesh(mesh)

! 
! declare, allocate and initialize arrays 
! 
! ocean
call ocean_setup(mesh, .True.)

! ice 
call ice_setup

! atm
call atm_setup

! diagnostics
call setup_diagnostics(.True.)

! initialize REcoM arrays
call recom_initialization(mesh)

! initial conditions
! Read (ice, turbulence, PAR, wind, ...) forcing
call read_forcing(mesh)

! read atm (CO2, Fe, N) deposition
call read_deposition(mesh)

! check if new year
call clock_newyear

print*, 'initialization done'
!=============================================================================!
! 
! main loop of the model 
! 
print*, 'start computation loop'
do istep=1,nsteps
        !print*, istep, nsteps
	! update time
	call clock

	! run main 
	!!!! double check is here
	call recom(istep,mesh)

	! update tracers
	! 1) perform mixing (due to turbulence)
	call recom_mixing
	! 2) remineralisation and sinking
	call update_tracers
	
	! perform diagnostics
	if (mask_diagnostic(istep)) call store_diagnostics(index_diagnostic(istep))
	
enddo

print*, 'computation done'
!=============================================================================!
print*, 'start diagnostics'
! 
! write diagnostics to netcdf and deallocate arrays
! 
! 
! write diagnostics
! 
call write_diagnostics
print*, 'diagnostics done'
! 
! deallocate arrays
! 
! diagnostic
call setup_diagnostics(.False.)
! fluxes
call deallocate_flux
! forcing
!call deallocate_ocean
call forcing_setup(.False.)
call deposition_setup(.False.)
call ocean_setup(mesh, .False.)
end program
