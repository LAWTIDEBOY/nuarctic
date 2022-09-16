! ==============================================================
module REcoM_setup

  use general_config
  use REcoM_config 
  use ocean_module
  use REcoM_clock
  
  implicit none
  
  integer		:: nsteps
  
  contains
  
! ==============================================================
subroutine read_namelist

  !use REcoM_diag, only:ldiag_carbon, ldiag_silicate, recom_diag_freq, &
  !                      recom_diag_freq_unit, recom_logfile_outfreq, precom_diag_list
                
  implicit none

  character(len=4096)   :: nml_file, nml_path
  
  
  namelist /clockinit/ timenew, daynew, yearnew
  
  !----------------------------------------
  call get_environment_variable("RECOM_NAMELIST_PATH", nml_path)
  ! configuration namelist
  nml_file =trim(nml_path)//trim('namelist.config')    ! name of general configuration namelist file
  open (20,file=nml_file)
  read (20,NML=simulationname)
  read (20,NML=timestep)
  read (20,NML=clockinit)
  read (20,NML=meshproperties)
  read (20,NML=forcingproperties)  
  read (20,NML=diagnostics)
  read (20,NML=calendar)  
  !read (20,NML=paths)
  close (20)
  dt=86400./real(step_per_day)
  
  ! ocean namelist (tracers related)
  !nml_file =trim(nml_path)//trim('namelist.ocean')
  !open (20,file=nml_file)
  !read (20,NML=ocean_tracers)
  !sclose(20)  
  
  ! recom specific namelist
  nml_file =trim(nml_path)//trim('namelist.recom') ! name of recom namelist file   !read (20,NML=precom_diag_list)
  open (20,file=nml_file)
  read (20,NML=pavariables)
  read (20,NML=pasinking)
  read (20,NML=painitialization_N)
  read (20,NML=paArrhenius)
  read (20,NML=palimiter_function)
  read (20,NML=palight_calculations)
  read (20,NML=paphotosynthesis)
  read (20,NML=paassimilation)
  read (20,NML=pairon_chem)
  read (20,NML=pazooplankton)
  read (20,NML=pasecondzooplankton)
  read (20,NML=pagrazingdetritus)
  read (20,NML=paaggregation)
  read (20,NML=padin_rho_N)
  read (20,NML=padic_rho_C1)
  read (20,NML=paphytoplankton_N)
  read (20,NML=paphytoplankton_C)
  read (20,NML=paphytoplankton_ChlA)
  read (20,NML=padetritus_N)
  read (20,NML=padetritus_C)
  read (20,NML=paheterotrophs)
  read (20,NML=paseczooloss)
  read (20,NML=pairon)
  read (20,NML=pacalc)
  read (20,NML=pabenthos_decay_rate)
  read (20,NML=paco2_flux_param)
  read (20,NML=paalkalinity_restoring)
  close (20)

write(*,*) '#--------------------------------------------------------------------#'
write(*,*) 'Namelist files read'
write(*,*) 'Check model setup'
write(*,*) 'simulation name: ', trim(runid)
write(*,*) 'time (initial and step): ', timenew, daynew, yearnew, dt
write(*,*) 'mesh filename: ', trim(meshname)
write(*,*) '#--------------------------------------------------------------------#'
end subroutine

! ==============================================================
!
subroutine read_mesh(mesh)

  use mod_mesh
  !use read_mesh_interface
  
  implicit none

#include "netcdf.inc" 

  type(t_mesh), intent(inout), target :: mesh
  character(len=4096)	:: filename
  Real(kind=8),dimension(:),allocatable  :: tmp
  integer		:: ncid, status
  integer		:: lon_varid, lat_varid, dpth_varid, dim_id
  integer		:: zbar_varid, Z_varid, nlevels_varid 
  !real (kind=8), dimension(:), allocatable :: lon, lat, dpth
  character(len=4096) :: grid_path
  character(len=4096) :: dname='time', lname='nl'
   !---------------------------------------------
  call get_environment_variable("RECOM_GRID_PATH", grid_path)
  filename = trim(grid_path) // trim(meshname)
  
  ! general 1D mesh characteristics
  !mesh%nb_of_nodes = 1
  mesh%ocean_area = 0
  mesh%ulevels=1 
  ! open file
  status=nf_open(trim(filename), nf_nowrite, ncid)
  
  !
  ! horizontal mesh (representing time)
  !
  ! inquire dimensions id and dimensions
  status=nf_inq_dimid(ncid,trim(dname),dim_id)
  status=nf_inq_dim(ncid,dim_id,trim(dname),mesh%npt)

  ! allocate coordinates arrays
  allocate(mesh%geo_coords(2,mesh%npt), mesh%depth(mesh%npt), tmp(mesh%npt))

  ! inquire variable ids
  status=nf_inq_varid(ncid, 'longitude', lon_varid)
  status=nf_inq_varid(ncid, 'latitude', lat_varid)
  status=nf_inq_varid(ncid, 'depth', dpth_varid)
  
  ! read variables
  status=nf_get_vara_double(ncid,lon_varid,1,mesh%npt,tmp)
  mesh%geo_coords(1,:) = tmp
  status=nf_get_vara_double(ncid,lat_varid,1,mesh%npt,tmp)
  mesh%geo_coords(2,:) = tmp
  status=nf_get_vara_double(ncid,dpth_varid,1,mesh%npt,mesh%depth)

  !
  ! vertical mesh (representing time)
  !
  ! inquire dimensions id and dimensions
  status=nf_inq_dimid(ncid,trim(lname),dim_id)
  status=nf_inq_dim(ncid,dim_id,trim(lname),mesh%nl)
  
  ! allocate vertical mesh arrays
  allocate(mesh%zbar(mesh%nl), mesh%Z(mesh%nl-1), mesh%nlevels(mesh%npt))
  
  ! inquire variable ids
  status=nf_inq_varid(ncid, 'zbar', zbar_varid)
  status=nf_inq_varid(ncid, 'Z', Z_varid)
  status=nf_inq_varid(ncid, 'nlevels', nlevels_varid)

  ! read variables
  status=nf_get_vara_double(ncid,zbar_varid,1,mesh%nl,mesh%zbar)
  status=nf_get_vara_double(ncid,Z_varid,1,mesh%nl-1,mesh%Z)
  status=nf_get_vara_int(ncid,nlevels_varid,1,mesh%npt, mesh%nlevels) 
  
  status=nf_close(ncid)
  deallocate (tmp)
end subroutine

! ==============================================================
subroutine get_run_steps
  ! Coded by Qiang Wang
  !--------------------------------------------------------------
  use recom_clock
  
  implicit none
  integer      :: i, temp_year, temp_mon, temp_fleapyear

  ! clock should have been inialized before calling this routine

  if(run_length_unit=='s') then
     nsteps=run_length
  elseif(run_length_unit=='d') then
     nsteps=step_per_day*run_length
  elseif(run_length_unit=='m') then
     nsteps=0
     temp_mon=month-1
     temp_year=yearnew
     temp_fleapyear=fleapyear
     do i=1,run_length
        temp_mon=temp_mon+1
        if(temp_mon>12) then
           temp_year=temp_year+1
           temp_mon=1
           call check_fleapyr(temp_year, temp_fleapyear)
        end if
        nsteps=nsteps+step_per_day*num_day_in_month(temp_fleapyear,temp_mon)
     end do
  elseif(run_length_unit=='y') then
     nsteps=0
     do i=1,run_length
        temp_year=yearnew+i-1
        call check_fleapyr(temp_year, temp_fleapyear)
        nsteps=nsteps+step_per_day*(365+temp_fleapyear)
     end do
  end if
end subroutine
! ==============================================================
subroutine get_spatial_configuration_step(istep,mesh)
  use mod_mesh
  use REcoM_Glovar
  use ocean_module


  !use get_spatial_configuration_step_interface
  implicit none 
  
  integer, intent(in) :: istep
  type(t_mesh), intent(in), target :: mesh

  ! set current geo-coordinates
  geo_coords(:) = mesh%geo_coords(:,istep)
  ! set number of vertical levels concerned for vertical column computation (varies along the track)
  ulevel=1
  nlevel=mesh%nlevels(istep)
  ! set the depth of the different vertical nodes
  hnode=0.d0
  hnode(1:nlevel-1) = zbar(1:nlevel-1) - zbar(2:nlevel)

end subroutine

! ==============================================================
end module
  
