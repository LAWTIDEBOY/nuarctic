module ocean_array_setup_interface
  interface
    subroutine ocean_array_setup(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module ice_array_setup_interface
  interface
    subroutine ice_array_setup(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module read_forcing_interface
  interface
    subroutine read_forcing(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
! ==============================================================
module ocean_module
  ! Ocean module that gather ocean variables and related routines that are required  REcOM1D (derived from FESOM2)
  implicit none
  save
  
  integer :: num_tracers=30
  integer, dimension(30) :: tracer_ID 
  real(kind=8), allocatable :: tr_arr(:,:)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! REcoM/ocean related tracers list:
  ! 1)  Temperature 									(ID:   0)
  ! 2)  Salinity    									(ID:   1)

  !!!! Dissolved  matter
  ! inorganic
  ! 3)  Dissolved inorganic Nitrogen 					(DIN)		(ID:1001)
  ! 4)  Dissolved inorganic Carbon   					(DIC)		(ID:1002)
  ! 5)  Total Alkalinity 	    					(Alk)		(ID:1003)
  ! organic
  ! 13) Dissolved organic Nitrogen 					(DON)		(ID:1011)
  ! 14) Dissolved organic Carbon 					(DOC)		(ID:1012)

  !!! Phytoplankton
  ! Phytoplankton (1st group, small (nano) phytoplankton)
  ! 6)  concentration (intracellular) of Nitrogen in small phytoplankton(PhyN)		(ID:1004)
  ! 7)  concentration (intracellular) of Carbon in small phytoplankton	(PhyC)		(ID:1005)
  ! 8)  Current concentration of Chlorophyl A in small phytoplankton	(PhyChl)	(ID:1006)
  ! Diatoms (2nd group of phytoplankton, larger phytoplankton)
  ! 15) concentration of Nitrogen in diatoms				(DiaN)		(ID:1013)
  ! 16) concentration of Carbon in diatoms				(DiaC)		(ID:1014)
  ! 17) concentration of Chlorophyl A in diatoms			(DiaChl)	(ID:1015)
  ! 18) concentration of Silicate in diatoms				(DiaSi)		(ID:1016)

  !!! Zooplanton
  ! heterotrophs (1st group of phytoplankton)
  ! 11) concentration of Nitrogen in heterotrophs			(HetN)		(ID:1009)
  ! 12) concentration of Carbon in heterotrophs				(HetC)		(ID:1010)
  ! 2nd group of phytoplankton
  ! 25) concentration of Nitrogen in the 2nd group of zooplankton	(Zoo2N)		(ID:1023)
  ! 26) concentration of Carbon in the 2nd group of zooplankton		(Zoo2C)		(ID:1024)
    
  !!! Detritus
  ! from Heterotroph (1st group of zooplankton)
  ! 9)  Concentration of Nitrogen in detritus				(DetN)		(ID:1007)
  ! 10) Concentration of Carbon in detritus				(DetC)		(ID:1008)
  ! from 2nd group of zooplankton
  ! 27) Concentration of Nitrogen in detritus				(idetz2n)	(ID:1025)
  ! 28) Concentration of Carbon in detritus				(idetz2c)	(ID:1026)
  ! 29) Concentration of Silicate in detritus				(idetz2Si)	(ID:1027)
  
  !!! 
  ! 22) Calcification of Phytoplankton 					(PhyCalc)	(ID:1020)
  ! 23) Calcification of Detritus (zooplankton 1) 			(DetCalc)	(ID:1021)
  ! 30) Calcification of Detritus (zooplankton 2)			(idetz2calc)	(ID:1028)
  
  !!! miscealleneous
  ! 20) concentration of Dissolved Silicate				(DSi) 		(ID:1018)
  ! 21)	concentration of Dissolved Iron					(DFe)		(ID:1019)
  ! 24) concentration in dioxygen					(O2)		(ID:1022)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! depth related ALE arrays
  integer :: nlevels
  real(kind=8), allocatable,dimension(:)   :: hnode, zbar_3d_n, Z_3d_n

  ! geo coordinates along the simulated track
  real(kind=8), allocatable,dimension(:) :: geo_coords	
  
  ! forcing variables
  real(kind=8), allocatable, dimension(:) :: temperature, salinity, Kz, PAR
  real(kind=8), allocatable  :: shortwave
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ocean namelist
  namelist /ocean_tracers/ num_tracers, tracer_ID 
  contains
  !
  !--------------------------------------------------------------------------------
  !

subroutine ocean_array_setup(mesh)
  ! 
  ! allocate and initialize ocean related arrays
  ! 
  use mod_mesh
  use REcoM_config
  use REcoM_GloVar
  
  implicit none
  
  type(t_mesh), intent(in), target :: mesh
  
  integer :: i
  !
  !
  nl=mesh%nl
  !nb_of_nodes=mesh%nb_of_nodes
  
  !!!!!!!!!!!!!!
  ! tracers
  !!!!!!!!!!!!!!
  ! allocate tracers arrays
  allocate(tr_arr(nl-1, num_tracers))
  
  ! initialize tracer arrays
  tr_arr=0.d0
  write(*,*) 'initializing tracers'
  
  !!!!!!!!!!!!!!
  ! alkalinity
  !!!!!!!!!!!!!!
  if (restore_alkalinity) then
    ! initialization
    Alk_surf=0.d0
    relax_alk=0.d0
    virtual_alk=0.d0
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! geo coordinates (lon/lat) of the track point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(geo_coords(2))
  geo_coords = 0.d0     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! depth-related ALE arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) 'ALE arrays'
  ! allocation
  allocate(hnode(nl-1), zbar_3d_n(nl), Z_3d_n(nl-1))
  ! initialization
  zbar_3d_n=mesh%zbar
  Z_3d_n=mesh%Z
  hnode = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! forcing ocean data: temperature, salinity, PAR, Kz
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! allocation
  allocate(temperature(nl), salinity(nl), PAR(nl), Kz(nl))
  ! initialization
  temperature 	= 0.d0
  salinity 	= 0.d0
  PAR		= 0.d0
  Kz		= 0.d0
  
end subroutine  

!

end module


module ice_module
  ! Ocean module that gather ocean variables and related routines that are required  REcOM1D (derived from FESOM2)
  implicit none
  save
  
   real(kind=8)		   :: aice

  contains

subroutine ice_array_setup(mesh)
  ! 
  ! allocate and initialize ocean related arrays
  ! 
  use mod_mesh
  use REcoM_config
  use REcoM_GloVar
  
  type(t_mesh), intent(in), target :: mesh
  
  nl=mesh%nl
  !nb_of_nodes=mesh%nb_of_nodes
  
  aice=1.d0
end subroutine

end module

module forcing_module
  ! Ocean module that gather ocean variables and related routines that are required  REcOM1D (derived from FESOM2)
  implicit none
  save
   
   ! temperature and salinity
   real(kind=8), allocatable, dimension(:,:)  :: temperature_forcing, salinity_forcing
   ! ice cover
   real(kind=8), allocatable, dimension(:)  :: aice_forcing
   ! turbulent forcing and Photosynthetically Active Radiation
   real(kind=8), allocatable, dimension(:,:)  :: Kz_forcing, PAR_forcing
   
   ! shortwave radiation at surface or ice/ocean interface
   real(kind=8), allocatable, dimension(:) :: shortwave_forcing
   !  
   contains
   
subroutine read_forcing(mesh)
   ! 
   ! allocate and initialize ocean related arrays
   ! 
   use mod_mesh
   use general_config
   use REcoM_config
 
#include "netcdf.inc" 

   type(t_mesh), intent(in), target :: mesh
   real(kind=8), allocatable, dimension(:,:) :: tmp
   integer		:: start1(1), count1(1), start2(2), count2(2)
   character(len=4096)	:: filename
   integer		:: i, ncid, status, nptf, nlf
   integer		:: ndims_in, nvars_in, ngatts_in, unlimdimid_in
   integer		:: aice_varid, T_varid, S_varid, dim_id
   integer		:: Kz_varid, PAR_varid, sw_varid 
   integer		:: ndims, nvars
   character(len=4096) :: variable_name
   character(len=4096) :: dname='time', lname='level'
   character(len=4096) :: iname = 'aice', Tname='temperature', Sname = 'salinity'
   character(len=4096) :: Pname = 'PAR', Swname='shortwave', Kname='Kz'
   !---------------------------------------------
   filename = trim(data_path) // trim(forcingname)
   

   ! initialization 
   ! dimensions
   nl=mesh%nl
   npt=mesh%npt
   print*, 'nl', nl
   ! netcdf counters
   start1=1
   count1=npt
   start2=(/1, 1/)
   count2=(/npt, nl/)
   allocate(tmp(npt,nl))
   ! initialize/allocate arrays
   call forcing_array_setup
   
   ! read forcing from netcdf dedicated forcing file
   status=nf_open(trim(filename), nf_nowrite, ncid)
   !
   ! check forcing array dimension
   !
   ! inquire time dimension
   status=nf_inq_dimid(ncid,trim(dname),dim_id)
   status=nf_inq_dim(ncid,dim_id,trim(dname),nptf)
   print*, nptf
   if (nptf .ne. npt) then
   	write(*,*) 'forcing and mesh time frame are different, please check', nptf, npt
   	return
   endif
   ! inquire level dimension
   status=nf_inq_dimid(ncid,trim(lname),dim_id)
   status=nf_inq_dim(ncid,dim_id,trim(lname),nlf)
   if (nlf .ne. nl) then
   	write(*,*) 'forcing and mesh nb of vertical levels are different, please check', nlf, nl
   	return
   endif  
   
   ! read ice cover  
   status=nf_inq_varid(ncid, trim(iname), aice_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,aice_varid, start1,count1,aice_forcing)
   else
   	!call handle_err(status)
   	write(*,*) 'no ice cover specified as forcing'
   endif
   
   ! read water column temperature
   status=nf_inq_varid(ncid, trim(Tname), T_varid)
   print*, trim(Tname), T_varid, status
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,T_varid, start2, count2, tmp)
   	temperature_forcing=transpose(tmp)
   	!print*, 'test test', temperature_forcing(:,1)
   else
   	!call handle_err(status)
   	write(*,*) 'no temperature profile specified as forcing'
   endif
   
   ! read salinity in the water column
   status=nf_inq_varid(ncid, Sname, S_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,S_varid,start2,count2,tmp)
	salinity_forcing=transpose(tmp)
   else
   	!call handle_err(status)
   	write(*,*) 'no salinity profile specified as forcing'
   endif   
   
   
   ! read turbulence in the water column
   status=nf_inq_varid(ncid, Kname, Kz_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,Kz_varid,start2,count2,tmp)
	Kz_forcing = transpose(tmp)
   else
   	!call handle_err(status)
   	write(*,*) 'no turbulence profile specified as forcing'
   endif    
   
   ! read PAR in the water column
   status=nf_inq_varid(ncid, Pname, PAR_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,PAR_varid,start2,count2,tmp)
   	PAR_forcing = transpose(tmp)
   else
   	!call handle_err(status)
   	write(*,*) 'no PAR profile specified as forcing'
   endif    
   
   ! read shortwave radiation at the surface
   status=nf_inq_varid(ncid, Swname, sw_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,sw_varid,start1,count1,shortwave_forcing) 
   	! estimate PAR in the water column from surface shortwave radiation
   else
   	!call handle_err(status)
   	write(*,*) 'no surface shortwave radiation specified as forcing'
   endif    
     
   status=nf_close(ncid)
   deallocate(tmp)
end subroutine

subroutine get_forcing(istep)
    use ice_module
    use ocean_module
    
    implicit none
    
    integer, intent(in) :: istep
    
    
    ! get forcing at selected time step
    aice = aice_forcing(istep)
    print*, size(temperature), size(temperature_forcing(:,istep))
    temperature(:) = temperature_forcing(:, istep)
    salinity(:) = salinity_forcing(:, istep)
    Kz(:) = Kz_forcing(:, istep)
    PAR(:) = PAR_forcing(:, istep)
    !shortwave = shortwave_forcing(istep)
    print*, 'temperature check', temperature
end subroutine

subroutine forcing_array_setup
  ! 
  ! allocate and initialize ocean related arrays
  ! 
  use mod_mesh
  use REcoM_config
  use REcoM_GloVar
  
  
  ! allocation of the forcing variables
  
  allocate(aice_forcing(npt))
  allocate(temperature_forcing(nl, npt), salinity_forcing(nl, npt))
  allocate(Kz_forcing(nl, npt), PAR_forcing(nl,npt))
  allocate(shortwave_forcing(npt))
 
  ! array initialization
  aice_forcing = 0.d0
  temperature_forcing = 0.d0
  salinity_forcing = 0.d0
  Kz_forcing = 0.d0
  PAR_forcing = 0.d0
  shortwave_forcing = 0.d0 

end subroutine

end module

