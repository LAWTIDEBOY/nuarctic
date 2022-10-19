module recom_main

  use REcoM_declarations
  use REcoM_LocVar
  use REcoM_GloVar
  use general_config
  use recom_config
  use mod_mesh
  use recom_clock
  use ocean_module
  use ice_module
  use atm_module
  use atm_deposition_module
  use forcing_module
  use recom_setup
  use recom_extra
  use recom_sms
  use recom_forcing



contains
! CONTENT:
! ------------
!    subroutine recom_init
!    subroutine recom_forcing
!    subroutine recom_sms
!
! initially written by REcoM group, adapted by  O. Gurses 02.03.2020
! rewritten for REcoM1D by F. Birrien

subroutine recom(istep,mesh)
  
  implicit none
  type(t_mesh), intent(in) , target :: mesh
  integer, intent(in)		    :: istep
! ======================================================================================
!! Depth information

!! zbar (depth of layers) and Z (mid depths of layers)
!! zbar is negative 
!! zbar(nl) allocate the array for storing the standard depths
!! Z(nl-1)  mid-depths of cells

!! max. number of levels at node n
!! nzmax = nlevels_nod2D(n)
!! u_ice and v_ice are at nodes
!! u_w, v_w are at nodes (interpolated from elements)
!! u_wind and v_wind are always at nodes
! ======================================================================================

  real(kind=8)               :: SW, Loc_slp
  integer                    :: nzmax  
  integer                    :: idiags

  !real(kind=8)               :: Sali
  !real (kind=8), allocatable :: Temp(:),  zr(:), PAR(:)
  real(kind=8),  allocatable :: C(:,:)
! ======================================================================================
  ! allocate local tracer array
  allocate(C(nl-1,bgc_num))

! ======================================================================================
!************************* SURFACE BOUNDARY  *********************************			
  ! atmospheric deposition (CO2, Fe, N)
  call get_atm_deposition(istep) 

  ! restore alkalinity
  relax_alk = surf_relax_Alk*(Alk_surf - tr_arr(1,2+ialk)) ! 1 temp, 2 salt

! ======================================================================================
!************************* get forcing and update steps *********************************
  !!!!!!!!!!!!!!!!!!!!!!!!!!! test continue here
  call get_spatial_configuration_step(istep,mesh)
  call get_forcing(istep)	

! ======================================================================================
!**************************gather variables for recom_forcing *****************************************				
  ! number of computation level
  nzmax = nlevel-1
  
  ! local ice concentration and mean sea level pressure
  Loc_ice_conc = aice 
  Loc_slp = press_air

  !!---- Benthic layers
  LocBenthos(1:benthos_num) = Benthos(1:benthos_num)

  !!---- Local conc of [H+]-ions from last time time step. Stored in LocVar
  !!---- used as first guess for H+ conc.in subroutine CO2flux (provided by recom_init)
  Hplus = GloHplus                                  

  !!---- Interpolated wind from atmospheric forcing 
  !!---- temporarily stored in module LocVar
  ULoc = sqrt(u_wind**2+v_wind**2)

  !!---- Atmospheric CO2 in LocVar                                                                        
  LocAtmCO2  = AtmCO2   


  !!---- Shortwave penetration and PAR initialization
  if (flag_PAR) then
     ! if PAR profile is available 
     SW = 0.d0
  else
     if (flag_PAR_surface) then
     	! PAR is available at the ocean-ice interface
     	SW = PAR_surface
     else
     	! shortwave is available at the atm-snow/ice interface and require attenuation
     	SW = parFrac * shortwave
     	SW = SW * (1.d0 - aice) !! this is a simplification, how to account for snow and ice thickness
     endif
     PAR(1:nzmax) = 0.d0 ! PAR initialization (computation in sms)
  endif
  
  
  !!!! not needed as well as PAR3D
  !!---- Temperature in water column and surface salinity
  !!Temp(1:nzmax) = tr_arr(1:nzmax, 1)
  !!Sali = tr_arr(1, 2)
  !!---- Depth of the nodes in the water column 
  !!zr(1:nzmax) = Z_3d_n(1:nzmax, n)   
     
  !!---- Biogeochemical tracers
  C(1:nzmax,1:bgc_num) = tr_arr(1:nzmax, 3:num_tracers)             

                       
  !!---- a_ice(row): Ice concentration in the local node
  FeDust = GloFeDust * (1 - aice) * dust_sol    
  NDust = GloNDust  * (1 - aice)

  allocate(Diags3Dloc(nzmax,8))
  Diags3Dloc(:,:) = 0.d0



! ======================================================================================
!******************************** RECOM FORCING ****************************************

  call REcoM_computation(nzmax, C, SW, Loc_slp, temperature(1:nzmax), salinity(1), PAR(1:nzmax))
     
  ! update tracers subsequently
  tr_arr(1:nzmax, 3:num_tracers)  = C(1:nzmax, 1:bgc_num)

! =====================================================================================================
!******************************** saved local variables****************************************
  Benthos(1:benthos_num)         = LocBenthos(1:benthos_num)         ! Updating Benthos values
  Diags2D(1:8)                   = LocDiags2D(1:8)                   ! Updating diagnostics
  GloPCO2surf                    = pco2surf
  GlodPCO2surf                   = dpco2surf

  GloCO2flux                     = dflux
  GloCO2flux_seaicemask          = co2flux_seaicemask                 !  [mmol/m2/s]
  GloO2flux_seaicemask           = o2flux_seaicemask                  !  [mmol/m2/s]

  GloHplus                       = hplus
  AtmFeInput                     = FeDust
  AtmNInput                      = NDust 
  DenitBen                       = LocDenit

  GlodecayBenthos(1:benthos_num) = decayBenthos(1:benthos_num)/SecondsPerDay ! convert from [mmol/m2/d] to [mmol/m2/s]  

  PAR3D(1:nzmax)             = PAR(1:nzmax)
   
  do idiags = 1,diags3d_num
    Diags3D(1:nzmax,idiags)  = Diags3Dloc(1:nzmax,idiags) ! 1=NPPnano, 2=NPPdia
  end do

  deallocate(Diags3Dloc,C)

end subroutine recom

end module recom_main
