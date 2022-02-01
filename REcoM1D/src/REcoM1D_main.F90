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
use ocean_module
use ice_module
use forcing_module
use REcoM_setup
use REcoM_clock

!=============================================================================!    
! local variables
type(t_mesh),             target, save :: mesh
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
call ocean_array_setup(mesh)

! declare, allocate and initialize ice arrays
call ice_array_setup(mesh)

! initialize REcoM arrays
call recom_init(mesh)

! initial conditions

! Read forcing
call read_forcing(mesh)

!=============================================================================!
!
! main loop of the model 
!
istep=7000
call get_spatial_configuration_step(istep,mesh)
print*, 'nb of levels', nlevels
call get_forcing(istep)




!=============================================================================!
!
! postprocessing diagnostics and deallocation steps
!


end program
