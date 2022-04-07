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
  real(kind=8), allocatable,dimension(:)   :: hnode, zbar, Z

  ! geo coordinates along the simulated track
  real(kind=8), allocatable,dimension(:) :: geo_coords	
  
  ! forcing variables
  real(kind=8), allocatable, dimension(:) :: temperature, salinity, Kz, PAR
  real(kind=8)				  :: shortwave
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ocean namelist
  namelist /ocean_tracers/ num_tracers, tracer_ID 
  contains
  !
  !--------------------------------------------------------------------------------
  !

subroutine ocean_setup(mesh)
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
  allocate(hnode(nl-1), zbar(nl), Z(nl-1))
  ! initialization
  zbar=mesh%zbar
  Z=mesh%Z
  hnode = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! forcing ocean data: temperature, salinity, PAR, Kz
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! allocation
  allocate(temperature(nl-1), salinity(nl-1), PAR(nl-1), Kz(nl-1))
  ! initialization
  temperature 	= 0.d0
  salinity 	= 0.d0
  PAR		= 0.d0
  Kz		= 0.d0
  
end subroutine  

end module


module ice_module
  ! ice module that gather ice variables
  implicit none
  save
  
   real(kind=8)		   :: aice

  contains

subroutine ice_setup
  ! 
  ! allocate and initialize ocean related arrays
  ! 
  aice=1.d0
end subroutine

end module

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
 
   ! atmospheric forcing
   real(kind=8), allocatable, dimension(:) :: pressure_forcing, uatm_forcing, vatm_forcing 
   
   ! flag whether PAR is provided as forcing or not
   logical				   :: flag_PAR
   
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
   integer		:: slp_varid, uw_varid, vw_varid
   integer		:: ndims, nvars
   character(len=4096) :: variable_name
   character(len=4096) :: dname='time', lname='level'
   character(len=4096) :: iname = 'aice', Tname='temperature', Sname = 'salinity'
   character(len=4096) :: Pname = 'PAR', Swname='shortwave', Kname='Kz'
   character(len=4096) :: uwname= 'uatm', vwname='vatm', slpname='sea_level_pressure'
   !---------------------------------------------
   filename = trim(forcing_path) // trim(forcingname)
   

   ! initialization 
   ! dimensions
   nl=mesh%nl
   npt=mesh%npt
   ! netcdf counters
   start1=1
   count1=npt
   start2=(/1, 1/)
   count2=(/npt, nl/)
   allocate(tmp(npt,nl))
   ! initialize/allocate arrays
   call forcing_setup
   
   ! read forcing from netcdf dedicated forcing file (pre-processed before simulation)
   ! 
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
   	status=nf_get_vara_double(ncid,uw_varid, start1, count1, uatm_forcing)
   else
   	!call handle_err(status)
   	write(*,*) 'no wind specified as forcing'
   endif
   ! v-wind
   status=nf_inq_varid(ncid, trim(vwname), vw_varid)
   if (status .eq. NF_NOERR) then
   	status=nf_get_vara_double(ncid,vw_varid, start1, count1, vatm_forcing)
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
   	status=nf_get_vara_double(ncid,Kz_varid,start2,count2,tmp)
	Kz_forcing = transpose(tmp)
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
   else
   	write(*,*) 'no PAR profile specified as forcing, surface shortwave is to be provided'
        flag_PAR=.False.
        ! read shortwave radiation at the surface
        status=nf_inq_varid(ncid, Swname, sw_varid)
	if (status .eq. NF_NOERR) then
   		status=nf_get_vara_double(ncid,sw_varid,start1,count1,shortwave_forcing)
   	else
   		!call handle_err(status)
   		write(*,*) 'no surface shortwave radiation and/or PAR are specified as forcing'
   	endif
   endif    
     
     
     
   status=nf_close(ncid)
   deallocate(tmp)
end subroutine

subroutine get_forcing(istep)
    use ice_module
    use ocean_module
    use atm_module
    
    implicit none
    
    integer, intent(in) :: istep
    
    
    ! get forcing at selected time step
    ! ice related forcing
    aice = aice_forcing(istep)

    ! ocean related forcing
    temperature(:) = temperature_forcing(:, istep)
    salinity(:) = salinity_forcing(:, istep)
    Kz(:) = Kz_forcing(:, istep)
    
    if (flag_PAR) then
    	PAR(:) = PAR_forcing(:, istep)
    else
    	shortwave = shortwave_forcing(istep)
    endif
    
    ! atmospheric forcing
    u_wind = uatm_forcing(istep)
    v_wind = vatm_forcing(istep)
    press_air = pressure_forcing(istep)
    
end subroutine

subroutine forcing_setup
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

module atm_deposition_module
  implicit none
  save
   
   ! C02, Iron (Fe) and Nitrogen (N) atm deposition forcing 
   real(kind=8), allocatable, dimension(:) :: CO2_forcing, Fe_forcing, N_forcing

contains
subroutine read_deposition(mesh)
   ! 
   ! allocate and initialize ocean related arrays
   ! 
   use mod_mesh
   use general_config
   use REcoM_config
 
#include "netcdf.inc" 

   type(t_mesh), intent(in), target :: mesh
   real(kind=8), allocatable, dimension(:) :: tmp
   integer		:: start1(1), count1(1)
   character(len=4096)	:: filename
   integer		:: i, ncid, status, nptf, nlf
   integer		:: ndims_in, nvars_in, ngatts_in, unlimdimid_in
   integer		:: CO2_varid, Fe_varid, N_varid, dim_id
   integer		:: ndims, nvars
   character(len=4096) :: variable_name
   character(len=4096) :: dname='time'
   character(len=4096) :: CO2name = 'CO2_deposition', Fename='Fe_deposition', Nname = 'N_deposition'

   !---------------------------------------------
   filename = trim(data_path) // trim(atmdepositionname)
   

   ! initialization 
   ! dimensions
   npt=mesh%npt
   ! netcdf counters
   start1=1
   count1=npt
   allocate(tmp(npt))
   ! initialize/allocate arrays
   call deposition_setup
   
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
   ! read CO2 depsosition
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
    use REcoM_config
    use REcoM_GloVar
    implicit none
    
    integer, intent(in) :: istep
    
    
    ! get deposition at selected time step
    ! CO2 deposition
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

subroutine deposition_setup
  ! 
  ! allocate and initialize arrays related to atm deposition
  ! 
  use REcoM_config
  
  ! allocation of the forcing variables
  allocate(CO2_forcing(npt), Fe_forcing(npt), N_forcing(npt))

  ! array initialization
  CO2_forcing = 0.d0
  Fe_forcing = 0.d0
  N_forcing = 0.d0

end subroutine

end module

