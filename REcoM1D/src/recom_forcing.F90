module REcoM_forcing


contains
!===============================================================================
! REcoM_Forcing
!===============================================================================
subroutine REcoM_computation(Nn, state, SurfSW, Loc_slp, Temp, Sali, PARc)


  use REcoM_declarations
  use REcoM_LocVar
  use REcoM_GloVar
  use general_config
  use recom_config
  use recom_clock
  use ocean_module
  use ice_module
  use atm_module
  use forcing_module
  use recom_extra
  use gasx, only : pistonvel, flxco2, o2flux
  use recom_sms


  Implicit none
  
  !external pistonvel, flxco2, o2flux
  
  ! Arguments
  integer, intent (in)			     :: Nn            ! number of vertical computation nodes 
  real(kind=8),dimension(:,:), intent(inout) :: state 
  real(kind=8), intent(in)                   :: SurfSW        ! [W/m2] ShortWave radiation at surface
  real(kind=8), intent(in)                   :: Loc_slp       ! [Pa] sea-level pressure 
  Real(kind=8), dimension(:), intent(in)     :: Temp          ! [degrees C] Ocean temperature
  real(kind=8), dimension(:),intent(inout)      :: PARc 
  real(kind=8), intent(in)                   :: Sali                 ! Salinity of current surface layer
  
  ! local variables
  Real(kind=8)                              :: Latr          
 
! Subroutine Depth

  Real(kind=8),dimension(:),allocatable         :: zF                   ! [m] Depth of fluxes
  Real(kind=8),dimension(:,:),allocatable       :: SinkVel              ! [m/day]
  Real(kind=8),dimension(:), allocatable        :: thick                ! [m] Vertical distance between two nodes = Thickness 
  Real(kind=8),dimension(:), allocatable        :: recipthick           ! [1/m] reciprocal of thick

! Subroutine CO2Flux /mocsy
  Real(kind=8)                         :: REcoM_DIC            ! [mmol/m3] Conc of DIC in the surface water, used to calculate CO2 flux
  Real(kind=8)                         :: REcoM_Alk            ! [mmol/m3] Conc of Alk in the surface water, used to calculate CO2 flux
  Real(kind=8)                         :: REcoM_Si             ! [mol/m3] Conc of Si in the surface water, used to calculate CO2 flux
  Real(kind=8)                         :: REcoM_Phos           ! [mol/m3] Conc of Phos in the surface water, used to calculate the CO2 flux
  Real(kind=8)                         :: Latd                 ! latitude in degree
  Real(kind=8)                         :: Lond                 ! longitude in degree
  Real(kind=8)                         :: REcoM_T              ! temperature again, for mocsy minimum defined as -2
  Real(kind=8)                         :: REcoM_S              ! temperature again, for mocsy minimum defined as 21
! atm pressure, now read in as forcing!!
  Real(kind=8)                         :: Patm                 ! atmospheric pressure [atm]

! Subroutine o2flux /mocsy 
  Real(kind=8)                          :: ppo                 ! atmospheric pressure, divided by 1 atm 
  Real(kind=8)                          :: REcoM_O2            ! [mmol/m3] Conc of O2 in the surface water, used to calculate O2 flux

! misceallaneous
  Real(kind=8)				:: tmp, diags1d_tmp			

! Subroutine REcoM_sms_computation
  Real(kind=8),dimension(:,:), allocatable :: sms, aux                ! matrix that entail changes in tracer concentrations

! units for gas exchanges (! with character length, ! interface not possible with library)
  character(2), parameter :: optP='m ', optKf='dg'
  character(3), parameter :: optB='u74', optK1K2='l  '
  character(4), parameter :: optS='Sprc'
  character(6), parameter :: optCON='mol/m3'
  character(7), parameter :: optT='Tpot   ', optGAS='Pinsitu'

!Diagnostics
  integer                              :: idiags, j
!===============================================================================
  ! allocate and initialize arrays 
  allocate(zF(nl), SinkVel(nl,4), thick(nl-1), recipthick(nl-1), sms(nl-1,bgc_num), aux(nl-1,bgc_num))
  zF=0.d0
  SinkVel=0.d0
  thick=0.d0
  recipthick=0.d0
  sms=0.d0
  aux=0.d0
  
  ! variable initialization
  tiny_N   = tiny_chl/chl2N_max      ! 0.00001/ 3.15d0   Chl2N_max [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
  tiny_N_d = tiny_chl/chl2N_max_d    ! 0.00001/ 4.2d0
  tiny_C   = tiny_N  /NCmax          ! NCmax   = 0.2d0   [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
  tiny_C_d = tiny_N_d/NCmax_d        ! NCmax_d = 0.2d0 
  tiny_Si  = tiny_C_d/SiCmax         ! SiCmax = 0.8d0

 
  !!!!!!!!!!!!!!!!! change here -> look ine recom_extra
  !!!!!!!!!!!!!!!!!! sinking velocity and tracer evolution needs to be included at the early stage in the loop
  ! compute 
  call Cobeta      
  call Depth_calculations(Nn, SinkVel, zF, thick, recipthick)

!! ----- mocsy -------! 

!! convert from mmol/m3 to mol/m3
  REcoM_DIC  = max(tiny*1e-3, state(one,idic)*1e-3)
  REcoM_Alk  = max(tiny*1e-3, state(one,ialk)*1e-3)
  REcoM_Si   = max(tiny*1e-3, state(one,isi) *1e-3)

!! convert N to P with Redfield ratio
  REcoM_Phos = max(tiny*1e-3, state(one,idin)*1e-3) /16
!! minimum set to 2 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
  REcoM_T    = max(2.d0, Temp(1))
!! maximum set to 40 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
  REcoM_T    = min(REcoM_T, 40.d0) 
!! minimum set to 21: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble in regions with S between 19 and 21 and ice conc above 97%
  REcoM_S    = max(21.d0, Sali) 
!! maximum set to 43: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble   REcoM_S    = min(REcoM_S, 43.d0)

!! convert from Pa to atm.
  Patm = Loc_slp/Pa2atm  

!! point coordinates
!! longitude
  Lond=geo_coords(1)/rad !! convert from rad to degree
!! latitutde
  Latr=geo_coords(2)
  Latd=geo_coords(2)/rad !! convert from rad to degree



  call pistonvel(ULoc, Loc_ice_conc, kw660)
  call flxco2(co2flux, co2ex, dpco2surf,                                                    &
                  ph, pco2surf, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
                  REcoM_T, REcoM_S, REcoM_Alk, REcoM_DIC, REcoM_Si, REcoM_Phos, kw660, LocAtmCO2, Patm, thick(One), Lond,Latd,     &
                  optCON, optT, optP, optB, optK1K2, optKf, optGAS, optS)

! changed optK1K2='l  ' to 'm10'
  if((co2flux>1.e10) .or. (co2flux<-1.e10)) then
!     co2flux=0.0  
      print*,'ERROR: co2 flux !'
      print*,'pco2surf: ',pco2surf
      print*,'co2: ',co2
      print*,'rhoSW: ', rhoSW
      print*,'temp: ',REcoM_T
      print*,'tempis: ',tempis
      print*,'REcoM_S: ', REcoM_S
      print*, 'REcoM_Alk: ', REcom_Alk
      print*, 'REcoM_DIC: ', REcoM_DIC
      print*, 'REcoM_Si: ', REcoM_Si
      print*, 'REcoM_Phos: ', REcoM_Phos
      print*, 'kw660: ',kw660
      print*, 'LocAtmCO2: ', LocAtmCO2
      print*, 'Patm: ', Patm
      print*, 'thick(One): ',thick(One) 
      print*, 'Nmocsy: ', Nmocsy
      print*, 'Lond: ', Lond
      print*, 'Latd: ', Latd   
      print*, 'ULoc: ', ULoc
      print*, 'Loc_ice_conc: ', Loc_ice_conc
      stop
    endif

! ice-fraction is already considered in piston-velocity, so don't apply it here
   dflux     = co2flux *1.e3 *SecondsPerDay  !* (1.d0 - Loc_ice_conc)
   co2flux_seaicemask = co2flux * 1.e3 !  [mmol/m2/s]  * (1.d0 - Loc_ice_conc)

! then oxygen
   ppo = Loc_slp/Pa2atm !1 !slp divided by 1 atm
   REcoM_O2 = max(tiny*1e-3,state(one,ioxy)*1e-3) ! convert from mmol/m3 to mol/m3 for mocsy

   call  o2flux(REcoM_T, REcoM_S, kw660, ppo, REcoM_O2, o2ex)
     
   o2flux_seaicemask = o2ex * 1.e3 ! back to mmol here [mmol/m2/s] 

! Source-Minus-Sinks


!  addtiny(1:nn,1) = state(1:nn,isi)
!  addtiny(1:nn,2) = state(1:nn,idetsi) 
!  addtiny(1:nn,3) = state(1:nn,idiasi)
!  addtiny(1:nn,4) = state(1:nn,idetz2si)


!!!!!!!! keep on changing from here on
!!!!!!!!
!!!!!! recom sms needs to be integrated correcly
   call REcoM_sms_computation(Nn, state, thick, recipthick, SurfSW, sms, Temp, zF, PARc)

!  addtiny(1:nn,1) = (state(1:nn,isi)           - aux(1:nn,isi))
!  addtiny(1:nn,2) = (state(1:nn,idetsi)        - aux(1:nn,idetsi))
!  addtiny(1:nn,3) = (state(1:nn,idiasi)        - aux(1:nn,idiasi)) 
!  addtiny(1:nn,4) = (state(1:nn,idetz2si)      - aux(1:nn,idetz2si))

!  aux=0.0d0
!  aux(1:nn,:)        = state(1:nn,:) + sms(1:nn,:)

  state(1:nn,:)      = max(tiny,state(1:nn,:) + sms(1:nn,:))
  state(1:nn,ipchl)  = max(tiny_chl,state(1:nn,ipchl))
  state(1:nn,iphyn)  = max(tiny_N,  state(1:nn,iphyn))
  state(1:nn,iphyc)  = max(tiny_C,  state(1:nn,iphyc))
  state(1:nn,idchl)  = max(tiny_chl,state(1:nn,idchl))
  state(1:nn,idian)  = max(tiny_N_d,state(1:nn,idian))
  state(1:nn,idiac)  = max(tiny_C_d,state(1:nn,idiac))
  state(1:nn,idiasi) = max(tiny_Si, state(1:nn,idiasi))

!  addtiny(1:nn,5) = (state(1:nn,isi)           - aux(1:nn,isi))
!  addtiny(1:nn,6) = (state(1:nn,idetsi)        - aux(1:nn,idetsi))
!  addtiny(1:nn,7) = (state(1:nn,idiasi)        - aux(1:nn,idiasi)) 
!  addtiny(1:nn,8) = (state(1:nn,idetz2si)      - aux(1:nn,idetz2si))

!  addtiny(1:nn,5) = state(1:nn,isi)
!  addtiny(1:nn,6) = state(1:nn,idetsi) 
!  addtiny(1:nn,7) = state(1:nn,idiasi)
!  addtiny(1:nn,8) = state(1:nn,idetz2si)


!-------------------------------------------------------------------------------
! Diagnostics
  if (Diags) then
	do idiags = one,8
	  diags1d_tmp = 0
	  do j=1,nn
	  	tmp = diags2Dloc(j,idiags) * thick(j)
	  	if (.not. isnan(tmp)) diags1d_tmp = diags1d_tmp + tmp
	  enddo
	  LocDiags1D(idiags) = diags1d_tmp
	end do
  end if
! array deallocation
deallocate(zF, SinkVel, thick, recipthick, sms, aux)

end subroutine REcoM_computation

end module REcoM_forcing
