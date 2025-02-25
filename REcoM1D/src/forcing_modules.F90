! ==============================================================
module ocean_module
  ! Ocean module that gather ocean variables and related routines that are required  REcOM1D (derived from FESOM2)

  use mod_mesh
  use general_config
  use REcoM_config
  use REcoM_GloVar

  implicit none

  
  integer, parameter :: num_tracers=30
  integer, dimension(num_tracers) :: tracer_ID 
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
  ! 19) Concentration of Silicate in detritus                           (DetSi)         (ID:1017)
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
  integer :: nlevel, ulevel
  real(kind=8), allocatable,dimension(:)   :: hnode, zbar, Z

  ! geo coordinates along the simulated track
  real(kind=8), allocatable,dimension(:) :: geo_coords	
  
  ! forcing variables
  real(kind=8), allocatable, dimension(:) :: temperature, salinity, Kz, PAR
  real(kind=8)				  :: PAR_surface, shortwave
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! initialization variables for tracers
  real(kind=8), allocatable, dimension(:) :: DIN_init, DIC_init, Alk_init, DSi_init, DFe_init, DO2_init 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ocean namelist
  !namelist /ocean_tracers/ num_tracers, tracer_ID 
  data tracer_ID / 0, 1, 1001, 1002, 1003, 1004, 1005, 1006, &
                 1007, 1008, 1009, 1010, 1011, 1012, 1013, &
                 1014, 1015, 1016, 1017, 1018, 1019, 1020, &
                 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028 /
  contains
  !
  !--------------------------------------------------------------------------------
  !

subroutine ocean_setup(mesh, boolean)
  ! 
  ! allocate and initialize ocean related arrays
  ! 
 
  implicit none
  
  type(t_mesh), intent(in), target :: mesh
  logical, intent(in)		   :: boolean
  integer :: i
  !
  !
  if (boolean) then
  nl=mesh%nl
  !nb_of_nodes=mesh%nb_of_nodes
  
  !!!!!!!!!!!!!!
  ! tracers
  !!!!!!!!!!!!!!
  ! allocate tracers arrays
  allocate(tr_arr(nl-1, num_tracers))
  ! allocate initialization of specific tracers
  allocate (DIN_init(nl-1), DIC_init(nl-1), Alk_init(nl-1))
  allocate (DSi_init(nl-1), DFe_init(nl-1), DO2_init(nl-1))
  tr_arr=0.d0
  DIN_init=0.d0
  DIC_init=0.d0
  Alk_init=0.d0
  DSi_init=0.d0
  DFe_init=0.d0
  DO2_init=0.d0

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
  allocate(hnode(nl-1), zbar(nl), Z(nl-1))
  ! initialization
  zbar=mesh%zbar
  Z=mesh%Z
  hnode = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! forcing ocean data: temperature, salinity, PAR, Kz
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! allocation
  allocate(temperature(nl-1), salinity(nl-1), PAR(nl-1), Kz(nl))
  ! initialization
  temperature 	= 0.d0
  salinity 	= 0.d0
  PAR		= 0.d0
  Kz		= 0.d0
  
  else
  	deallocate(tr_arr, geo_coords)
  	deallocate(DIN_init, DIC_init, Alk_init)
  	deallocate(DSi_init, DFe_init, DO2_init)
	deallocate(hnode, zbar, Z)
	deallocate(temperature)
	deallocate(salinity)
	deallocate(PAR)
	deallocate(Kz)
 endif
end subroutine  
! 
subroutine read_tracer_initialization(mesh)
 
#include "netcdf.inc" 


   type(t_mesh), intent(in), target :: mesh
   integer		:: start1(1), count1(1)
   character(len=4096)	:: filename
   integer		:: i, ncid, status, nlf, dim_id
   integer		:: DIN_varid, DIC_varid, Alk_varid
   integer		:: DSi_varid, DFe_varid, DO2_varid
   character(len=4096) :: data_path
   character(len=4096) :: lname='index'
   character(len=4096) :: DO2_name='DO2', DIN_name='DIN', DIC_name='DIC'
   character(len=4096) :: DSi_name='DSi', DFe_name='DFe', Alk_name='Alk'
 
   call get_environment_variable("RECOM_DATA_PATH", data_path)
   filename = trim(data_path) // trim(tracerinitname)
   print*, 'initialization tracer file: ', trim(filename)
   ! initialization 
   ! dimensions
   nl=mesh%nl
   ! netcdf counters
   start1=1
   count1=nl-1

   ! read initialization from netcdf dedicated tracer initialization file (pre-processed before simulation)
   ! 
   ! open file 
   status=nf_open(trim(filename), nf_nowrite, ncid)
   ! check forcing array dimension
   ! inquire level dimension
   status=nf_inq_dimid(ncid,trim(lname),dim_id)
   status=nf_inq_dim(ncid,dim_id,trim(lname),nlf)
   if (nlf .ne. nl-1) then
   	write(*,*) 'forcing and mesh nb of vertical levels are different, please check', nlf, nl-1
   	return
   endif  
   !----------------------------------------------------------------------------
   ! read DIN (tracer:3, ID:1001)
   status=nf_inq_varid(ncid, trim(DIN_name), DIN_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid, DIN_varid, start1, count1, DIN_init)
   else
   	!call handle_err(status)
   	write(*,*) 'no DIN specified for initialization'
   endif

   ! read DIC (tracer:4, ID:1002)
   status=nf_inq_varid(ncid, trim(DIC_name), DIC_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid, DIC_varid, start1, count1, DIC_init)
   else
   	!call handle_err(status)
   	write(*,*) 'no DIC specified for initialization'
   endif 

   ! read Alk (tracer:5, ID:1003)
   status=nf_inq_varid(ncid, trim(Alk_name), Alk_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid, Alk_varid, start1, count1, Alk_init)
   else
   	!call handle_err(status)
   	write(*,*) 'no Alk specified for initialization'
   endif   
   
   ! read DSi (tracer:20, ID:1018)
   status=nf_inq_varid(ncid, trim(DSi_name), DSi_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid, DSi_varid, start1, count1, DSi_init)
   else
   	!call handle_err(status)
   	write(*,*) 'no DSi specified for initialization'
   endif 
   
   ! read DFe (tracer:21, ID:1019)
   status=nf_inq_varid(ncid, trim(DFe_name), DFe_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid, DFe_varid, start1, count1, DFe_init)
   else
   	!call handle_err(status)
   	write(*,*) 'no DFe specified for initialization'
   endif 

   ! read DO2 (tracer:24, ID:1022)
   status=nf_inq_varid(ncid, trim(DO2_name), DO2_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid, DO2_varid, start1, count1, DO2_init)
   else
   	!call handle_err(status)
   	write(*,*) 'no O2 specified for initialization'
   endif    
   
   status=nf_close(ncid)

end subroutine


end module

! ==============================================================
module ice_module
  ! ice module that gather ice variables
  implicit none
  save
  
   real(kind=8)		   :: aice

  contains

subroutine ice_setup
  ! 
  ! initialize ice related arrays
  ! 
  aice=1.d0
end subroutine

end module

! ==============================================================
module atm_module
  ! atm-related module
  implicit none
  save
  
   real(kind=8) :: u_wind, v_wind, press_air	! wind components and mean sea level pressure

  contains

subroutine atm_setup
  ! 
  ! allocate and initialize ocean related arrays
  ! 
  ! wind
  u_wind=0.d0
  v_wind=0.d0
  ! mean sea level pressure
  press_air=0.d0
end subroutine

end module

! ==============================================================
module forcing_module
  ! Ocean module that gather ocean variables and related routines that are required  REcOM1D (derived from FESOM2)
  
  use mod_mesh
  use general_config
  use REcoM_config
  use ice_module
  use ocean_module
  use atm_module
  
  implicit none
   
   ! temperature and salinity
   real(kind=8), allocatable, dimension(:,:)  :: temperature_forcing, salinity_forcing
   ! ice cover
   real(kind=8), allocatable, dimension(:)  :: aice_forcing
   ! turbulent forcing and Photosynthetically Active Radiation
   real(kind=8), allocatable, dimension(:,:)  :: Kz_forcing, PAR_forcing
   
   ! shortwave radiation at ice surface (shortwave) or ice/ocean interface (PAR_surface)
   real(kind=8), allocatable, dimension(:) :: shortwave_forcing, PAR_surface_forcing
 
   ! atmospheric forcing
   real(kind=8), allocatable, dimension(:) :: pressure_forcing, uwind_forcing, vwind_forcing 
   
   ! flag whether PAR is provided as forcing or not
   logical				   :: flag_PAR, flag_PAR_surface
   
   contains
   
subroutine read_forcing(mesh)
   ! 
   ! allocate and initialize ocean related arrays
   ! 
 
#include "netcdf.inc" 

   type(t_mesh), intent(in), target :: mesh
   real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2
   integer		:: start1(1), count1(1), start2(2), count2(2), start3(2), count3(2)
   character(len=4096)	:: filename
   integer		:: i, ncid, status, nptf, nlf
   integer		:: aice_varid, T_varid, S_varid, dim_id
   integer		:: Kz_varid, PAR_varid, PARSF_varid, sw_varid 
   integer		:: slp_varid, uw_varid, vw_varid
   character(len=4096) :: forcing_path
   character(len=4096) :: dname='time', lname='level_node'
   character(len=4096) :: iname = 'aice', Tname='temperature', Sname = 'salinity'
   character(len=4096) :: Pname = 'PAR', Swname='shortwave', PSname='PAR_surface', Kname='Kz'
   character(len=4096) :: uwname= 'uwind', vwname='vwind', slpname='sea_level_pressure'
   !---------------------------------------------
   call get_environment_variable("RECOM_FORCING_PATH", forcing_path)
   filename = trim(forcing_path) // trim(forcingname)
   print*, 'filename forcing: ', trim(filename)
   ! initialization 
   ! dimensions
   nl=mesh%nl
   npt=mesh%npt
   ! netcdf counters
   start1=1
   count1=npt
   start2=(/1, 1/)
   count2=(/npt, nl-1/)
   start3=(/1, 1/)
   count3=(/npt, nl/)  
   allocate(tmp(npt,nl-1), tmp2(npt,nl))
   ! initialize/allocate arrays
   call forcing_setup(.True.)
   
   ! read forcing from netcdf dedicated forcing file (pre-processed before simulation)
   !
   ! open file 
   status=nf_open(trim(filename), nf_nowrite, ncid)
   ! check forcing array dimension
   ! 
   ! inquire time dimension
   status=nf_inq_dimid(ncid,trim(dname),dim_id)
   status=nf_inq_dim(ncid,dim_id,trim(dname),nptf)
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
   !----------------------------------------------------------------------------
   ! read ice cover  
   status=nf_inq_varid(ncid, trim(iname), aice_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,aice_varid, start1,count1,aice_forcing)
   else
   	!call handle_err(status)
   	write(*,*) 'no ice cover specified as forcing'
   endif
   !----------------------------------------------------------------------------
   ! read wind
   ! u-wind    
   status=nf_inq_varid(ncid, trim(uwname), uw_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,uw_varid, start1, count1, uwind_forcing)
   else
   	!call handle_err(status)
   	write(*,*) 'no wind specified as forcing'
   endif
   ! v-wind
   status=nf_inq_varid(ncid, trim(vwname), vw_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,vw_varid, start1, count1, vwind_forcing)
   else
   	!call handle_err(status)
   	write(*,*) 'no wind specified as forcing'
   endif   
   
   ! read mean sea level pressure
   status=nf_inq_varid(ncid, trim(slpname), slp_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid, slp_varid, start1, count1, pressure_forcing)
   else
   	!call handle_err(status)
   	write(*,*) 'mean sea level pressure specified as forcing'
   endif   
   !----------------------------------------------------------------------------   
   ! read water column temperature
   status=nf_inq_varid(ncid, trim(Tname), T_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,T_varid, start2, count2, tmp)
   	temperature_forcing=transpose(tmp)
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
   	status=nf_get_vara_double(ncid,Kz_varid,start3,count3,tmp2)
	Kz_forcing = transpose(tmp2)
   else
   	!call handle_err(status)
   	write(*,*) 'no turbulence profile specified as forcing'
   endif    
   
   ! read PAR in the water column
   status=nf_inq_varid(ncid, Pname, PAR_varid)
   if (status .eq. NF_NOERR) then
   	flag_PAR=.True.
   	status=nf_get_vara_double(ncid,PAR_varid,start2,count2,tmp)
   	PAR_forcing = transpose(tmp)
   	write(*,*) 'PAR profiles are specified as forcing, surface PAR or above-ice shortwave is to be provided'
   else
        flag_PAR=.False.
        
        ! read shortwave radiation at the atm/snow(ice) surface or under ice PAR (ice-ocean interface)
        status=nf_inq_varid(ncid, PSname, PARSF_varid)
        ! PAR at the ocean surface (under-ice)
	if (status .eq. NF_NOERR) then
		flag_PAR_surface=.True.
   		status=nf_get_vara_double(ncid, PARSF_varid, start1, count1, PAR_surface_forcing)
   		write(*,*) 'PAR (ocean surface) time series is specified as forcing'
   	else
   		flag_PAR_surface=.False.
   		!call handle_err(status)
   		
   	endif
        ! shortwave
        status=nf_inq_varid(ncid, Swname, sw_varid)
	if (status .eq. NF_NOERR) then
   		status=nf_get_vara_double(ncid, sw_varid, start1, count1, shortwave_forcing)
   		write(*,*) 'shortwave radiation (above ice) time series is specified as forcing'
   	endif
   	
   endif    
     
   status=nf_close(ncid)
   deallocate(tmp)
   
end subroutine

subroutine get_forcing(istep)
    
    implicit none
    
    integer, intent(in) :: istep
    
    
    ! get forcing at selected time step
    ! ice related forcing
    aice = aice_forcing(istep)
    
    ! ocean related forcing
    temperature(:) = temperature_forcing(:, istep)
    salinity(:) = salinity_forcing(:, istep)
    Kz(:) = Kz_forcing(:, istep)
    
    ! just for output of salinity and temperature
    tr_arr(:,1) = temperature_forcing(:, istep) ! temperature
    tr_arr(:,2) = salinity_forcing(:, istep) ! salinity
    
    
    
    if (flag_PAR) then
    	PAR(:) = PAR_forcing(:, istep)
    elseif(.not. flag_PAR .and. flag_PAR_surface) then
    	PAR_surface = PAR_surface_forcing(istep) * 4
    else
    	shortwave = shortwave_forcing(istep)
    endif
    
    ! atmospheric forcing
    u_wind = uwind_forcing(istep)
    v_wind = vwind_forcing(istep)
    press_air = pressure_forcing(istep)
    
end subroutine

subroutine forcing_setup(boolean)
  ! 
  ! allocate and initialize ocean related arrays
  ! 
  use mod_mesh
  use REcoM_config
  use REcoM_GloVar
  
  implicit none
  
  logical, intent(in)	:: boolean
  
  if (boolean) then
  	! allocation of the forcing variables
  	allocate(aice_forcing(npt))
  	allocate(temperature_forcing(nl-1, npt), salinity_forcing(nl-1, npt))
  	allocate(Kz_forcing(nl, npt), PAR_forcing(nl-1,npt))
  	allocate(shortwave_forcing(npt), PAR_surface_forcing(npt))
  	allocate(uwind_forcing(npt), vwind_forcing(npt), pressure_forcing(npt))
  	! array initialization
  	aice_forcing = 0.d0
  	temperature_forcing = 0.d0
  	salinity_forcing = 0.d0
  	Kz_forcing = 0.d0
  	PAR_forcing = 0.d0
  	PAR_surface_forcing = 0.d0
  	shortwave_forcing = 0.d0 
  	uwind_forcing = 0.d0
  	vwind_forcing = 0.d0  
  	pressure_forcing = 0.d0
  else
  	deallocate(aice_forcing, temperature_forcing, salinity_forcing)
  	deallocate(Kz_forcing, PAR_forcing, PAR_surface_forcing, shortwave_forcing)
  	deallocate(uwind_forcing, vwind_forcing, pressure_forcing)
  endif
end subroutine

end module
! ==============================================================
module atm_deposition_module
  use mod_mesh
  use general_config
  use REcoM_config
  use REcoM_GloVar
  
  
  implicit none
  save
   
   ! C02, Iron (Fe) and Nitrogen (N) atm deposition forcing 
   real(kind=8), allocatable, dimension(:) :: CO2_forcing, Fe_forcing, N_forcing

contains
subroutine read_deposition(mesh)
   ! 
   ! allocate and initialize ocean related arrays
   ! 
 
#include "netcdf.inc" 

   type(t_mesh), intent(in), target :: mesh
   real(kind=8), allocatable, dimension(:) :: tmp
   integer		:: start1(1), count1(1)
   character(len=4096)	:: filename
   integer		:: i, ncid, status, nptf
   integer		:: CO2_varid, Fe_varid, N_varid, dim_id
   character(len=4096) :: data_path
   character(len=4096) :: dname='time'
   character(len=4096) :: CO2name = 'CO2_deposition', Fename='Fe_deposition', Nname = 'N_deposition'

   !---------------------------------------------
   call get_environment_variable("RECOM_DATA_PATH", data_path)
   filename = trim(data_path) // trim(atmdepositionname)

   ! initialization 
   ! dimensions
   npt=mesh%npt
   ! netcdf counters
   start1=1
   count1=npt
   allocate(tmp(npt))
   ! initialize/allocate arrays
   call deposition_setup(.True.)
   
   ! read atm deposition from netcdf dedicated deposition file (pre-processed before simulation)
   status=nf_open(trim(filename), nf_nowrite, ncid)
   ! 
   ! check forcing array dimension
   ! 
   ! inquire time dimension
   status=nf_inq_dimid(ncid,trim(dname),dim_id)
   status=nf_inq_dim(ncid,dim_id,trim(dname),nptf)
   if (nptf .ne. npt) then
   	write(*,*) 'deposition and mesh time frame are different, please check', nptf, npt
   	return
   endif
   !----------------------------------------------------------------------------
   ! read CO2 deposition
   status=nf_inq_varid(ncid, trim(CO2name), CO2_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,CO2_varid, start1,count1,CO2_forcing)
   else
   	!call handle_err(status)
   	write(*,*) 'no CO2 deposition specified as forcing'
   endif
   !----------------------------------------------------------------------------
   ! read Iron (Fe) deposition   
   status=nf_inq_varid(ncid, trim(Fename), Fe_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,Fe_varid, start1, count1, Fe_forcing)
   else
   	!call handle_err(status)
   	write(*,*) 'no Iron (Fe) deposition specified as forcing'
   endif
   !----------------------------------------------------------------------------
   ! read Nitrogen (N) deposition   
   status=nf_inq_varid(ncid, trim(Nname), N_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,N_varid, start1, count1, N_forcing)
   else
   	!call handle_err(status)
   	write(*,*) 'no wind specified as forcing'
   endif    
   status=nf_close(ncid)
   deallocate(tmp)
end subroutine 


subroutine get_atm_deposition(istep)
    
    implicit none
    
    integer, intent(in) :: istep
    
    
    ! get deposition at selected time step
    ! CO2 deposition
    ! AtmCO2 is not needed, only to scale isotope 13 and 14
    if (constant_CO2) then
    	AtmCO2 = CO2_for_spinup
    else
    	AtmCO2 = CO2_forcing(istep)
    endif
    
    ! Iron (Fe) deposition 
    if (UseFeDust) then
    	GloFeDust = Fe_forcing(istep)
    else
    	GloFeDust = 0.d0
    endif
    
    ! Nitrogen (N) deposition
     if (useAeolianN) then
    	GloNDust = N_forcing(istep)
    else
    	GloNDust = 0.d0
    endif   

end subroutine

subroutine deposition_setup(boolean)

  use REcoM_config
  
  implicit none
  logical, intent(in) :: boolean
  ! 
  ! allocate and initialize arrays related to atm deposition
  ! 
  
  if (boolean) then
  	! allocation of the forcing variables
  	allocate(CO2_forcing(npt), Fe_forcing(npt), N_forcing(npt))

  	! array initialization
  	CO2_forcing = 0.d0
  	Fe_forcing = 0.d0
 	N_forcing = 0.d0
  else
	deallocate(CO2_forcing, Fe_forcing, N_forcing)
  endif  
  	
end subroutine

end module

