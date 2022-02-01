!>
!! @par Copyright
!! This code is subject to the FESOM-REcoM - License - Agreement in it's most recent form.
!! Please see URL xxx
!! 
!! @brief Module for defining variables used in REcoM, ex constant sinking velocity and 
!! local time step dt
!!
!! @remarks This module contains namelist for recom
!! @author xxx, FESOM-REcoM, Bremerhaven (2019)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   2014: original implementation (V. Schourup-Kristensen)
!!   2019: ongoing
!
!==========================================================================================================
module general_config
  ! general configuration from FESOM
  implicit none
  save
  ! general variables and parameters
  integer		        :: mstep
  real(kind=8), parameter      :: pi=3.14159265358979
  real(kind=8), parameter      :: rad=pi/180.0
  real(kind=8), parameter      :: g=9.81
  !_____________________________________________________________________________
  ! *** Simulation name ***
  character(len=4096)           :: runid='test1'       ! a model/setup name
  namelist /simulationname/ runid
  
  !_____________________________________________________________________________
  ! *** time ***
  !real(kind=8)           :: dt=3600.     ! number of steps per day
  !integer                :: nt=24*365	! run length
  !character              :: time_unit='s' ! unit: y, d, s
  !namelist /timestep/ dt, nt, time_unit

  ! *** time step ***
  integer                :: step_per_day=24     ! number of steps per day
  integer                :: run_length=1	! run length
  character              :: run_length_unit='y' ! unit: y, d, s
  namelist /timestep/ step_per_day, run_length, run_length_unit
  !_____________________________________________________________________________
  ! *** track ***
  integer		:: nx=1
  integer		:: nt_track=365
  real(kind=8)		:: dt_track=1.
  character		:: time_track_unit='d' 
  character(len=4096)   :: meshname='MOSAiC_mesh'       ! a model/setup name
  namelist /meshproperties/ meshname
  !_____________________________________________________________________________
  ! *** forcing ***  
  character(len=4096)   :: forcingname='MOSAiC_forcing'  ! file name which contains REcoM forcing
  namelist /forcingproperties/ forcingname
  !_____________________________________________________________________________
  ! *** Paths for all in and out ***
  character(len=4096)        :: grid_path='./grid/'
  character(len=4096)        :: data_path='./data/'
  character(len=4096)        :: result_path='./result/'
  namelist /paths/  grid_path, data_path, result_path
  
  !_____________________________________________________________________________
  ! *** fleap_year ***
  logical                       :: include_fleapyear=.false.
  namelist /calendar/ include_fleapyear
  
       
  !_____________________________________________________________________________
  ! *** restart_log ***
  !integer                :: logfile_outfreq=1      ! logfile info. outp. freq., # steps
  !integer                :: restart_length=1
  !character              :: restart_length_unit='m'
  
  !namelist /restart_log/   restart_length, restart_length_unit, logfile_outfreq

  !_____________________________________________________________________________
  ! *** others ***
  logical                       :: r_restart
  real(kind=8)			:: dt, nt
end module

!==========================================================================================================
module mod_mesh

type t_mesh

! horizontal mesh characteristics
!integer                                     :: nb_of_nodes          ! only 1 node at a time

integer					    :: npt		  ! number of nodes along the track (= number of time steps)

! vertical mesh characteritics
integer                                     :: nl		  ! total number of layers

real(kind=8), allocatable, dimension(:)     :: zbar, Z 		  ! initial layer, mid-depth layer and element depth

! number of levels at elem and vertices considering bottom topography
integer, allocatable, dimension(:) :: nlevels

! nodes coordinates and depth
real(kind=8), allocatable, dimension(:,:) :: geo_coords
real(kind=8), allocatable, dimension(:)   :: depth

! miscealleaneous (to be potentially removed)
! upper boudnary index of all vertical vertice/element loops, default==1 but when 
integer				      :: ulevels
real(kind=8)                          :: ocean_area


end type t_mesh
end module

!==========================================================================================================
module recom_config
  implicit none
  save

!! *** General constants ***

  Integer :: idin = 1, idic = 2, ialk = 3, iphyn = 4, iphyc = 5, &
             ipchl = 6, idetn = 7, idetc = 8, ihetn = 9,         &
             ihetc = 10, idon = 11, idoc = 12, idian = 13,       &
             idiac = 14, idchl = 15, idiasi = 16, idetsi = 17,   &
             isi = 18, ife = 19, iphycal = 20, idetcal = 21,     &
             ioxy= 22, izoo2n = 23, izoo2c =24, idetz2n = 25,   &
             idetz2c = 26, idetz2si = 27, idetz2calc = 28    

  Integer :: ivphy = 1, ivdia = 2, ivdet = 3, ivdetsc = 4
  
  Integer :: nl, npt 

  Real(kind=8)                 :: zero           = 0.d0
  Integer                      :: one            = 1
  Real(kind=8)                 :: tiny           = 2.23D-16
  Real(kind=8)                 :: tiny_chl       = 0.00001 
  Real(kind=8)                 :: SecondsPerDay  = 86400.d0     ! [s/day]
  Real(kind=8)                 :: Pa2atm         = 101325.d0    ! [Pa/atm] 
  Real(kind=8)                 :: redO2C         = 1.453        ! O2:C ratio Anderson and Sarmiento, 1994

!! *** General configuration ***
  logical                :: REcoM_restart        = .false.
  logical                :: recom_binary_write   = .false.  ! Determines if tracervalue snapshots are saved. For fine grids it may crash the model to set this to true

  logical                :: recom_binary_init    = .false.  ! Restart from binary
  Integer                :: bgc_num               = 22
  Integer                :: diags3d_num           = 2       ! Number of diagnostic 3d tracers to be saved
  Real(kind=8)           :: VDet                  = 20.d0   ! Sinking velocity, constant through the water column and positive downwards
  Real(kind=8)           :: VDet_zoo2             = 200.d0   ! Sinking velocity, constant through the water column 
  Real(kind=8)           :: VPhy                  = 0.d0    !!! If the number of sinking velocities are different from 3, code needs to be changed !!!
  Real(kind=8)           :: VDia                  = 0.d0 
  Logical                :: allow_var_sinking     = .true.   
  Integer                :: biostep               = 1                    ! Number of times biology should be stepped forward for each time step		 
  Logical                :: REcoM_Geider_limiter  = .false.              ! Decides what routine should be used to calculate limiters in sms
  Logical                :: REcoM_Grazing_Variable_Preference = .true.  ! Decides if grazing should have preference for phyN or DiaN
  Logical                :: REcoM_Second_Zoo      = .false.    ! Decides whether having macrozooplankton and second detritus or not
  Logical                :: Grazing_detritus     = .false.    ! Decides grazing on detritus                            
  Logical                :: zoo2_fecal_loss     = .false.    ! Decides fecalloss for the second zooplankton            
  Logical                :: zoo2_initial_field     = .false.    ! Decides initialization of secondzoo        ! NOT CODED YET OG
  Logical                :: het_resp_noredfield    = .true.    ! Decides respiratation of copepods              
  Logical                :: diatom_mucus    = .true.    ! Effect of nutrient limitation on the aggregation
  Logical                :: Graz_pref_new    = .true.    ! If it is true Fasham et 1990, otherwise original recom variable preference
  Logical                :: Diags                 = .true.    !!!!!!!!!!!!!!!!!!!!!!Change in recom.F90 Diagnostics -> Diags
  Logical                :: constant_CO2          = .true.
  Logical                :: use_Fe2N              = .true.               ! use Fe2N instead of Fe2C, as in MITgcm version
  Logical                :: use_photodamage       = .false.             ! use Alvarez et al (2018) for chlorophyll degradation
  logical                :: HetRespFlux_plus      = .true.     !MB More stable computation of zooplankton respiration fluxes adding a small number to HetN
  character(100)         :: REcoMDataPath         = ''
  logical                :: restore_alkalinity    = .true.
  logical                :: NitrogenSS            = .false.   ! This one only activates rivers! And in principle denitrification, but denitrification is commented out. When set to true, external sources and sinks of nitrogen are activated (Riverine, aeolian and denitrification)
  logical                :: useAeolianN           = .false.   ! When set to true, aeolian nitrogen deposition is activated
  integer                :: firstyearoffesomcycle = 1948      ! The first year of the actual physical forcing (e.g. JRA-55) used
  integer                :: lastyearoffesomcycle  = 2009      ! Last year of the actual physical forcing used
  integer                :: numofCO2cycles        = 1         ! Number of cycles of the forcing planned 
  integer                :: currentCO2cycle       = 1         ! Which CO2 cycle we are currently running
  Logical                :: DIC_PI                = .true.
  integer                :: Nmocsy                = 1         ! Length of the vector that is passed to mocsy (always one for recom)
  logical                :: recom_debug           =.true.
  integer                :: benthos_num           = 4

  namelist /pavariables/ REcoM_restart,         recom_binary_write,      &
                       recom_binary_init,                 bgc_num,               diags3d_num,             &
                       VDet,          VDet_zoo2,     &
                       VPhy,                              VDia,                    &
                       allow_var_sinking,                 biostep,               REcoM_Geider_limiter,    &
                       REcoM_Grazing_Variable_Preference, REcoM_Second_Zoo,      Grazing_detritus,        &
                       zoo2_fecal_loss,                   zoo2_initial_field,    het_resp_noredfield,     &
                       diatom_mucus,                      Graz_pref_new,                &
                       Diags      ,                       constant_CO2,            &
                       use_Fe2N,                          use_photodamage,       HetRespFlux_plus,        &
                       REcoMDataPath,                     restore_alkalinity,    			  &
                       NitrogenSS,                        useAeolianN,           firstyearoffesomcycle,   &
                       lastyearoffesomcycle,              numofCO2cycles,        currentCO2cycle,         &
                       DIC_PI,                            Nmocsy,                recom_debug,             &
                       benthos_num
!!------------------------------------------------------------------------------
!! *** Sinking ***
  Real(kind=8)                 :: Vdet_a         = 0.0288       ! [1/day]
                                                                                                                           
  namelist /pasinking/ Vdet_a
!!------------------------------------------------------------------------------
!! *** Initialization ***
  Real(kind=8)                 :: cPhyN          = 0.2d0
  Real(kind=8)                 :: cHetN          = 0.2d0
  Real(kind=8)                 :: cZoo2N          = 0.2d0

  namelist /painitialization_N/ cPhyN, cHetN, cZoo2N
!!------------------------------------------------------------------------------
!! *** Arrhenius function ***
  Real(kind=8)                 :: recom_Tref     = 288.15d0       ! [K]
  Real(kind=8)                 :: C2K            = 273.15d0       !     Conversion from degrees C to K
  Real(kind=8)                 :: Ae             = 4500.d0        ! [K] Slope of the linear part of the Arrhenius function
  Real(kind=8)                 :: reminSi        = 0.02d0
  namelist /paArrhenius/ recom_Tref, C2K, Ae, reminSi  
!!------------------------------------------------------------------------------
!! *** For limiter function ***
  Real(kind=8)                 :: NMinSlope      = 50.d0 
  Real(kind=8)                 :: SiMinSlope     = 1000.d0
  Real(kind=8)                 :: NCmin          = 0.04d0
  Real(kind=8)                 :: NCmin_d        = 0.04d0
  Real(kind=8)                 :: SiCmin         = 0.04d0
  Real(kind=8)                 :: k_Fe           = 0.04d0
  Real(kind=8)                 :: k_Fe_d         = 0.12d0
  Real(kind=8)                 :: k_si           = 4.d0
  Real(kind=8)                 :: P_cm           = 3.0d0          ! [1/day]   the rate of C-specific photosynthesis
  Real(kind=8)                 :: P_cm_d         = 3.5d0
  namelist /palimiter_function/ NMinSlope, SiMinSlope, NCmin, NCmin_d, SiCmin, k_Fe, k_Fe_d, k_si, P_cm, P_cm_d   
!!------------------------------------------------------------------------------
!! *** For light calculations ***
  Real(kind=8)                 :: k_w            = 0.04d0         ! [1/m]              Light attenuation coefficient
  Real(kind=8)                 :: a_chl          = 0.03d0         ! [1/m * 1/(mg Chl)] Chlorophyll specific attenuation coefficients
  namelist /palight_calculations/ k_w, a_chl   
!!------------------------------------------------------------------------------
!! *** Photosynthesis ***
  Real(kind=8)                 :: alfa           = 0.14d0	  ! [(mmol C*m2)/(mg Chl*W*day)] 
  Real(kind=8)                 :: alfa_d         = 0.19d0         ! An initial slope of the P-I curve
  Real(kind=8)                 :: parFrac        = 0.43d0
  namelist /paphotosynthesis/ alfa, alfa_d, parFrac   
!!------------------------------------------------------------------------------
!! *** Assimilation ***
  Real(kind=8)                 :: V_cm_fact     = 0.7d0          ! scaling factor for temperature dependent maximum of C-specific N-uptake
  Real(kind=8)                 :: V_cm_fact_d   = 0.7d0  
  Real(kind=8)                 :: NMaxSlope     = 1000.d0       ! Max slope for limiting function
  Real(kind=8)                 :: SiMaxSlope    = 1000.d0
  Real(kind=8)                 :: NCmax         = 0.2d0         ! [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
  Real(kind=8)                 :: NCmax_d       = 0.2d0
  Real(kind=8)                 :: SiCmax        = 0.8d0
  Real(kind=8)                 :: NCuptakeRatio = 0.2d0         ! [mmol N/mmol C] Maximum uptake ratio of N:C
  Real(kind=8)                 :: NCUptakeRatio_d = 0.2d0
  Real(kind=8)                 :: SiCUptakeRatio= 0.2d0
  Real(kind=8)                 :: k_din         = 0.55d0          ! [mmol N/m3] Half-saturation constant for nitrate uptake
  Real(kind=8)                 :: k_din_d       = 1.0d0
  Real(kind=8)                 :: Chl2N_max     = 3.15d0           ! [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
  Real(kind=8)                 :: Chl2N_max_d   = 4.2d0
  Real(kind=8)                 :: res_phy       = 0.01d0          ! [1/day] Maintenance respiration rate constant
  Real(kind=8)                 :: res_phy_d     = 0.01d0
  Real(kind=8)                 :: biosynth      = 2.33d0          ! [mmol C/mmol N] Cost of biosynthesis
  Real(kind=8)                 :: biosynthSi    = 0.d0
  namelist /paassimilation/ V_cm_fact, V_cm_fact_d, NMaxSlope, SiMaxSlope, NCmax, NCmax_d, SiCmax,        &
                       NCuptakeRatio, NCUptakeRatio_d, SiCUptakeRatio, k_din, k_din_d,                    &
                       Chl2N_max, Chl2N_max_d, res_phy, res_phy_d, biosynth, biosynthSi
!!------------------------------------------------------------------------------
!! *** Iron chemistry ***
  Real(kind=8)                 :: totalligand     = 1.d0        ! [mumol/m3] order 1. Total free ligand
  Real(kind=8)                 :: ligandStabConst = 100.d0      ! [m3/mumol] order 100. Ligand-free iron stability constant
  namelist /pairon_chem/ totalligand, ligandStabConst
!!------------------------------------------------------------------------------
!! *** Zooplankton ***
  Real(kind=8)                 :: graz_max      = 2.4d0           ! [mmol N/(m3 * day)] Maximum grazing loss parameter 
  Real(kind=8)                 :: epsilonr       = 0.35d0          ! [(mmol N)2 /m6] Half saturation constant for grazing loss 
  Real(kind=8)                 :: res_het       = 0.01d0          ! [1/day] Respiration by heterotrophs and mortality (loss to detritus)
  Real(kind=8)                 :: Redfield      = 6.625           ! [mmol C/mmol N] Redfield ratio of C:N = 106:16
  Real(kind=8)                 :: loss_het      = 0.05d0          ! [1/day] Temperature dependent N degradation of extracellular organic N (EON)
  Real(kind=8)                 :: pzDia         = 0.5d0           ! Maximum diatom preference
  Real(kind=8)                 :: sDiaNsq       = 0.d0
  Real(kind=8)                 :: pzPhy         = 1.0d0           ! Maximum nano-phytoplankton preference                                                                                           
  Real(kind=8)                 :: sPhyNsq       = 0.d0
  real(kind=8)                 :: tiny_het       = 1.d-5          ! for more stable computation of HetRespFlux (_plus). Value can be > tiny because HetRespFlux ~ hetC**2.
  namelist /pazooplankton/ graz_max, epsilonr, res_het, Redfield, loss_het, pzDia, sDiaNsq, pzPhy, sPhyNsq, tiny_het
!-------------------------------------------------------------------------------                                   
                                                                                                                         
! Second Zooplankton                                                                                                                                                                                    
                             
  Real(kind=8)                 :: graz_max2      = 0.1d0          ! [mmol N/(m3 * day)] Maximum grazing loss parameter                                                                                        
  Real(kind=8)                 :: epsilon2       = 0.0144d0          ! [(mmol N)2 /m6] Half saturation constant for grazing loss                                                                              
  Real(kind=8)                 :: res_zoo2       = 0.0107d0          ! [1/day] Respiration by heterotrophs and mortality (loss to detritus)                                                            
! Real(kind=8)                 :: recip_res_zoo2 = 1./res_zoo2                                                                                                                                               
  Real(kind=8)                 :: loss_zoo2      = 0.003d0          ! [1/day] Temperature dependent N degradation of extracellular organic N
  Real(kind=8)                 :: fecal_rate_n      = 0.13d0                                                         
  Real(kind=8)                 :: fecal_rate_c      = 0.295d0                                                          
  Real(kind=8)                 :: pzDia2         = 1.d0           ! Maximum diatom preference                                                                                                                 
  Real(kind=8)                 :: sDiaNsq2       = 0.d0
  Real(kind=8)                 :: pzPhy2         = 0.5d0           ! Maximum diatom preference                                                                                                                
  Real(kind=8)                 :: sPhyNsq2       = 0.d0
  Real(kind=8)                 :: pzHet          = 0.8d0           ! Maximum diatom preference                                                                                                               
  Real(kind=8)                 :: sHetNsq        = 0.d0
  Real(kind=8)                 :: t1_zoo2        = 28145.d0           ! Krill temp. function constant1                                                                                                       
  Real(kind=8)                 :: t2_zoo2        = 272.5d0          ! Krill temp. function constant2                                                                                                         
  Real(kind=8)                 :: t3_zoo2        = 105234.d0          ! Krill temp. function constant3                                                                                                      
  Real(kind=8)                 :: t4_zoo2        = 274.15d0          ! Krill temp. function constant3                                                                                              
  namelist /pasecondzooplankton/ graz_max2, epsilon2, res_zoo2, & !recip_res_zoo2 &                                                                                                                         
                                 loss_zoo2, fecal_rate_n, fecal_rate_c, pzDia2, pzPhy2, pzHet,        &
                                 sDiaNsq2, sPhyNsq2, sHetNsq,            &
                                 t1_zoo2, t2_zoo2, t3_zoo2, t4_zoo2
!-------------------------------------------------------------------------------                                                                                                                          
!! *** Detritus Grazing Params ***                                                                                                                                                                        
  Real(kind=8)                 :: pzDet         = 1.d0           ! Maximum small detritus prefence by first zooplankton                                                                                   
  Real(kind=8)                 :: sDetNsq       = 0.d0  
  Real(kind=8)                 :: pzDetZ2       = 1.d0           ! Maximum large detritus preference by first zooplankton                                                                                 
  Real(kind=8)                 :: sDetZ2Nsq     = 0.d0
  Real(kind=8)                 :: pzDet2         = 1.d0           ! Maximum small detritus prefence by second zooplankton                                                                                
  Real(kind=8)                 :: sDetNsq2       = 0.d0  
  Real(kind=8)                 :: pzDetZ22       = 1.d0           ! Maximum large detritus preference by second zooplankton                                                                              
  Real(kind=8)                 :: sDetZ2Nsq2     = 0.d0
  namelist /pagrazingdetritus/ pzDet, sDetNsq, pzDetZ2, sDetZ2Nsq, &
                                 pzDet2, sDetNsq2, pzDetZ22, sDetZ2Nsq2
!!------------------------------------------------------------------------------
!! *** Aggregation ***
  Real(kind=8)                 :: agg_PD        = 0.165d0          ! [m3/(mmol N * day)] Maximum aggregation loss parameter for DetN
  Real(kind=8)                 :: agg_PP        = 0.015d0          ! [m3/(mmol N * day)] Maximum aggregation loss parameter for PhyN and DiaN (plankton)
  namelist /paaggregation/ agg_PD, agg_PP
!!------------------------------------------------------------------------------
!! *** DIN ***
  Real(kind=8)                 :: rho_N         = 0.11d0           ! [1/day] Temperature dependent N degradation of extracellular organic N (EON) (Remineralization of DON)
  namelist /padin_rho_N/ rho_N
!!------------------------------------------------------------------------------
!! *** DIC ***
  Real(kind=8)                 :: rho_C1         = 0.1d0           ! [1/day] Temperature dependent C degradation of extracellular organic C (EOC)
  namelist /padic_rho_C1/ rho_C1
!!------------------------------------------------------------------------------
!! *** Phytoplankton N ***
  Real(kind=8)                 :: lossN         = 0.05d0          ! [1/day] Phytoplankton loss of organic N compounds
  Real(kind=8)                 :: lossN_d       = 0.05d0
  namelist /paphytoplankton_N/ lossN, lossN_d
!!------------------------------------------------------------------------------
!! *** Phytoplankton C ***
  Real(kind=8)                 :: lossC         = 0.10d0          ! [1/day] Phytoplankton loss of carbon 
  Real(kind=8)                 :: lossC_d       = 0.10d0
  namelist /paphytoplankton_C/ lossC, lossC_d
!!------------------------------------------------------------------------------
!! *** Phytoplankton ChlA ***
  Real(8)                      :: deg_Chl       = 0.25d0        ! [1/day]
  Real(kind=8)                 :: deg_Chl_d     = 0.25d0
  namelist /paphytoplankton_ChlA/ deg_Chl, deg_Chl_d
!!------------------------------------------------------------------------------
!! *** Detritus N ***
  Real(kind=8)                 :: grazEff       = 0.4d0         ! [] Grazing efficiency (fraction of grazing flux into zooplankton pool) 
  Real(kind=8)                 :: grazEff2      = 0.8d0         ! [] Grazing efficiency (fraction of grazing flux into second zooplankton pool)
  Real(kind=8)                 :: reminN        = 0.165d0        ! [1/day] Temperature dependent remineralisation rate of detritus	
  namelist /padetritus_N/ grazEff, grazEff2, reminN
!!------------------------------------------------------------------------------
!! *** Detritus C ***
  Real(kind=8)                 :: reminC        = 0.15d0        ! [1/day] Temperature dependent remineralisation rate of detritus
  Real(kind=8)                 :: rho_c2        = 0.1d0        ! [1/day] Temperature dependent C degradation of TEP-C
  namelist /padetritus_C/ reminC, rho_c2
!!------------------------------------------------------------------------------
!! *** Heterotrophs ***
  Real(kind=8)                 :: lossN_z       = 0.15d0
  Real(kind=8)                 :: lossC_z       = 0.15d0
  namelist /paheterotrophs/ lossN_z, lossC_z
!!------------------------------------------------------------------------------
! Second Zooplankton                                                                                                                                                                                        
                          
  Real(kind=8)                 :: lossN_z2       = 0.02d0
  Real(kind=8)                 :: lossC_z2       = 0.02d0
  namelist /paseczooloss/ lossN_z2, lossC_z2
!!------------------------------------------------------------------------------
!! *** Iron ***
!! only Fe2C or Fe2N is used, but I am not allowed to introduce an if-statement here
  Real(kind=8)                 :: Fe2N          = 0.033d0     ! Fe2C * 6.625
  Real(kind=8)                 :: Fe2N_benthos  = 0.15d0        ! test, default was 0.14 Fe2C_benthos * 6.625 - will have to be tuned. [umol/m2/day]
  Real(kind=8)                 :: Fe2C          = 0.005d0
  Real(kind=8)                 :: Fe2C_benthos  = 0.02125       !0.68d0/32.d0       ! [umol/m2/day]
  Real(kind=8)                 :: kScavFe       = 0.07d0
  namelist /pairon/ Fe2N, Fe2N_benthos, Fe2C, Fe2C_benthos, kScavFe
!!------------------------------------------------------------------------------
!! *** Calcification ***
  Real(kind=8)                 :: calc_prod_ratio = 0.02d0
  Real(kind=8)                 :: calc_diss_guts  = 0.0d0
  Real(kind=8)                 :: calc_diss_rate  = 0.005714d0   !20.d0/3500.d0
  Real(kind=8)                 :: calc_diss_rate2 = 0.005714d0
  namelist /pacalc/ calc_prod_ratio, calc_diss_guts, calc_diss_rate, calc_diss_rate2
!!------------------------------------------------------------------------------
!! *** Benthos ***
  Real(kind=8)                 :: decayRateBenN   = 0.005d0
  Real(kind=8)                 :: decayRateBenC   = 0.005d0
  Real(kind=8)                 :: decayRateBenSi  = 0.005d0
  Real(kind=8)                 :: q_NC_Denit      = 0.86d0         ! N:C quota of the denitrification process
  namelist /pabenthos_decay_rate/ decayRateBenN, decayRateBenC, decayRateBenSi, q_NC_Denit
!!------------------------------------------------------------------------------
!! *** CO2-flux ***
  Real(kind=8)                 :: permil          = 0.000000976    ! 1.e-3/1024.5d0 ! Converting DIC from [mmol/m3] to [mol/kg]
  Real(kind=8)                 :: permeg          = 1.e-6          ! [atm/uatm] Changes units from uatm to atm
  Real(kind=8)                 :: Xacc            = 1.e-12         ! Accuracy for ph-iteration (phacc)
!  Real(kind=8)                :: pCO2a           = 380.d0         ! [uatm] Atmospheric partial pressure of CO2
  Real(kind=8)                 :: CO2_for_spinup  = 278.d0         !  
  namelist /paco2_flux_param/ permil, permeg, Xacc, CO2_for_spinup
!!------------------------------------------------------------------------------
!! *** Alkalinity restoring ***
  Real(kind=8)                 :: surf_relax_Alk = 3.2e-07 !10.d0/31536000.d0
  namelist /paalkalinity_restoring/ surf_relax_Alk

end module recom_config
!
!==========================================================================================================
!
Module REcoM_declarations
  implicit none
  save

  Integer       :: save_count_recom
  Real(kind=8)  :: tiny_N                 ! Min PhyN
  Real(kind=8)  :: tiny_N_d               ! Min DiaN
  Real(kind=8)  :: tiny_C                 ! Min PhyC
  Real(kind=8)  :: tiny_C_d               ! Min DiaC
  Real(kind=8)  :: tiny_Si                ! Min DiaSi
!-------------------------------------------------------------------------------
! Temperature dependence of rates
  Real(kind=8)  :: rTref                  ! [1/K] Reciproque value of reference temp for Arrhenius function
  Real(kind=8)  :: rTloc                  ! [1/K] Reciproque of local ocean temp
  Real(kind=8)  :: arrFunc                ! []    Temp dependence of rates 
!  if (REcoM_Second_Zoo) then
  Real(kind=8)  :: arrFuncZoo2           ! []    Temperature function for krill     
!  endif
  Real(kind=8)  :: reminSiT
!-------------------------------------------------------------------------------
! Quotas
  Real(kind=8)  :: quota, quota_dia           ! [mmol N/mmol C]  Quota between phytoplankton N and C
  Real(kind=8)  :: recipQuota, recipQuota_dia ! [mmol C/mmol N]  Reciproque of 'quota'
  Real(kind=8)  :: Chl2C, Chl2C_dia           ! [mg ChlA/mmol C] Quota between phytoplankton ChlA and C
  Real(kind=8)  :: Chl2C_plast, Chl2C_plast_dia   ! [mg ChlA/mmol C] needed for photodamage
  Real(kind=8)  :: Chl2N, Chl2N_dia           ! [mg ChlA/mmol N] Quota between phytoplankton ChlA and N
  Real(kind=8)  :: qSiC
  Real(kind=8)  :: qSiN
  Real(kind=8)  :: recipQZoo                  ! [mmol C/mmol N]  Quota between heterotrophic C and N 
!  if (REcoM_Second_Zoo) then
  Real(kind=8)  :: recipQZoo2                  ! [mmol C/mmol N]  Quota between second zoo  C and N                                                                                          
!  endif
!!! Grazing detritus Quotas for converting                                                                                                                                                                
  Real(kind=8)  :: recipDet                  ! [mmol C/mmol N]  Quota between second zoo  C and N                                                                                                         
  Real(kind=8)  :: recipDet2                  ! [mmol C/mmol N]  Quota between second zoo  C and N 

!-------------------------------------------------------------------------------
! For limiter function
  Real(kind=8)          :: qlimitFac, qlimitFacTmp !            Factor that regulates photosynthesis
  Real(kind=8),external :: recom_limiter           !            Function calculating qlimitFac
  Real(kind=8)          :: FeLimitFac              ! [Mumol/m3] Half sat constant for iron  
  Real(kind=8)          :: pMax, pMax_dia          ! [1/day]    Maximum rate of C-specific photosynthesis 
!-------------------------------------------------------------------------------
! Light
  Real(kind=8)  :: kappar                  ! [1/m]  Light attenuation coefficient modified by chla
  Real(kind=8)  :: kappastar              ! []
  Real(kind=8)  :: kdzUpper,kdzLower      ! []     light attenuation * deltaZ at lower and upper control volume border
  Real(kind=8)  :: chl_upper,chl_lower    ! [mg/m3]     chl  at lower and upper control volume border
  Real(kind=8)  :: Chlave                 ! [mg/m3]     vertical average chl between two nodes
  Real(kind=8)  :: Upperlight, Lowerlight ! [?]    light at upper and lower border of control volume
  Real(kind=8)  :: PARave                 ! [?]    Average light in the control volumes
!-------------------------------------------------------------------------------
! Photosynthesis
  Real(kind=8)  :: Cphot, Cphot_dia       ! [1/day] C-specific rate of photosynthesis
!-------------------------------------------------------------------------------
! Assimilation
  Real(kind=8)  :: V_cm                   ! scaling factor for temperature dependent maximum of C-specific N-uptake
  Real(kind=8)  :: limitFacN,limitFacN_dia! Factor that regulates N-assimilation. Calc from function recom_limiter
  Real(kind=8)  :: limitFacSi
  Real(kind=8)  :: N_assim, N_assim_dia   ! [mmol N/(mmol C * day)] C specific N utilization rate
  Real(kind=8)  :: Si_assim
!-------------------------------------------------------------------------------
! Chlorophyll
  Real(kind=8)  :: ChlSynth, ChlSynth_dia           ! [mg CHL/ mmol N] CHL a synthesis regulation term 
  Real(kind=8)  :: phyRespRate, phyRespRate_dia     ! [1/day] Phytoplankton respiration rate
  Real(kind=8)  :: KOchl, KOchl_dia                 ! coefficient for damage to the photosynthetic apparatus 
!-------------------------------------------------------------------------------
! Iron chemistry
  Real(kind=8),external :: iron_chemistry 
!-------------------------------------------------------------------------------
! Zooplankton
  Real(kind=8)  :: DiaNsq  
  Real(kind=8)  :: varpzdia, fDiaN                  ! Part of Diatoms available for food
  Real(kind=8)  :: PhyNsq  
  Real(kind=8)  :: varpzPhy, fPhyN                  ! Part of Nano available for food
  Real(kind=8)  :: food, foodsq                     ! [(mmol N)2/m6]
  Real(kind=8)  :: grazingFlux_phy, grazingFlux_Dia ! [mmol N / (m3 * day)]
  Real(kind=8)  :: grazingFlux
  Real(kind=8)  :: HetRespFlux                      ! Zooplankton respiration
  Real(kind=8)  :: HetLossFlux                      ! [(mmol N)2/(m6 * day)] Zooplankton mortality (quadratic loss)
!-------------------------------------------------------------------------------
!  if (REcoM_Second_Zoo) then
! Second Zooplankton                                                                                           
       
  Real(kind=8)  :: DiaNsq2, PhyNsq2, HetNsq   
  Real(kind=8)  :: varpzDia2, fDiaN2, varpzPhy2, fPhyN2, varpzHet, fHetN ! Part of Diatoms available for food   
  Real(kind=8)  :: food2, foodsq2                     ! [(mmol N)2/m6]                                   
  Real(kind=8)  :: grazingFlux_phy2, grazingFlux_Dia2, grazingFlux_het2 ! [mmol N / (m3 * day)]        
  Real(kind=8)  :: grazingFlux2
  Real(kind=8)  :: Zoo2RespFlux     ! Zooplankton respiration                   
  Real(kind=8)  :: Zoo2LossFlux    ! [(mmol N)2/(m6 * day)] Zooplankton mortality (quadratic loss)  
  Real(kind=8)  :: Zoo2fecalloss_n    ! [(mmol N)/(m3*day)] Second zoo fecal pellet                        
  Real(kind=8)  :: Zoo2fecalloss_c    ! [(mmol N)/(m3*day)] Second zoo fecal pellet                               
  Real(kind=8)  :: recip_res_zoo22 
!  endif
! Grazing Detritus                                                                                                                                                                                       \
                                                                                                                                                                                                          
  Real(kind=8)  :: DetNsq, DetZ2Nsq, DetNsq2, DetZ2Nsq2  
  Real(kind=8)  :: varpzDet, varpzDetZ2, varpzDet2, varpzDetZ22         ! Part of Diatoms available for food                                                                                              
  Real(kind=8)  :: fDetN, fDetZ2N, fDetN2, fDetZ2N2                                                                                                                         
  Real(kind=8)  :: grazingFlux_Det, grazingFlux_DetZ2 ! [mmol N / (m3 * day)]                                                                                                                            
  Real(kind=8)  :: grazingFlux_Det2, grazingFlux_DetZ22 ! [mmol N / (m3 * day)]                                                                                                                          \
                                                                                   
!---------------------------------------------------------------------------------
! Aggregation
  Real(kind=8)  :: AggregationRate                  ! [1/day] AggregationRate (of nitrogen)
!-------------------------------------------------------------------------------
! Calcification
  Real(kind=8)  :: calcification
  Real(kind=8)  :: calc_loss_agg
  Real(kind=8)  :: calc_loss_gra
  Real(kind=8)  :: calc_diss
  Real(kind=8)  :: calc_loss_gra2 !zoo2 detritus
  Real(kind=8)  :: calc_diss2     !zoo2 detritus
!-------------------------------------------------------------------------------
! Diagnostics
  Real(kind=8)  :: recipbiostep                         ! 1/number of steps per recom cycle
  Real(kind=8),allocatable,dimension(:,:) :: Diags3Dloc
!-------------------------------------------------------------------------------
! Benthos
  Real(kind=8),allocatable,dimension(:) :: decayBenthos ! [1/day] Decay rate of detritus in the benthic layer
  Real(kind=8),allocatable,dimension(:) :: wFluxDet     ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
  Real(kind=8),allocatable,dimension(:) :: wFluxPhy     ! [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of phytoplankton
  Real(kind=8),allocatable,dimension(:) :: wFluxDia     ! [mmol/(m2 * day)] Flux of N,C, Si and chl through sinking of diatoms 	
  Real(kind=8)              :: Vben_det     ! [m/day] speed of sinking into benthos from water column
  Real(kind=8)              :: Vben_det_seczoo !second zooplankton sinking benthos  
  Real(kind=8)              :: Vben_phy
  Real(kind=8)              :: Vben_dia 
  Real(kind=8)              :: Ironflux     ! [umol Fe/(m2*day)] Flux of Fe from sediment to water
!_______________________________________________________________________________
! Arrays added for RECOM implementation:
!!---- PAR
!real(kind=8),allocatable,dimension(:)     :: PAR

end module REcoM_declarations

!==========================================================================================================
! For arrays needed for the whole 2D or 3D domain, but only needed in REcoM
!-------------------------------------------------------------------------------
Module REcoM_GloVar
  implicit none
  save
	
  Real(kind=8),allocatable,dimension(:)   :: Benthos          ! 4 types of benthos-tracers with size [4 n2d]
  Real(kind=8),dimension(12)              :: AtmCO2           ! [uatm] Atmospheric CO2 partial pressure. One value for the whole planet for each month

  Real(kind=8) 				  :: GloPCO2surf      ! [uatm] Surface ocean CO2 partial pressure
  Real(kind=8)				  :: GloCO2flux       ! [mmol/m2/day] Positive downwards
  Real(kind=8)				  :: GloCO2flux_seaicemask       ! [mmol/m2/day] Positive downwards
  Real(kind=8)				  :: GloO2flux_seaicemask       ! [mmol/m2/day] Positive downwards
  Real(kind=8)				  :: GloHplus         ! [mol/kg] Concentrations of H-plus ions in the surface ocean
  Real(kind=8)				  :: GlodPCO2surf       ! [mmol/m2/day] ocean-atmosphere  
  Real(kind=8),allocatable,dimension(:)   :: GlodecayBenthos  ! [1/day] Decay rate of detritus in the benthic layer saved for oce_ale_tracer.F90

  Real(kind=8),allocatable,dimension(:)	  :: GlowFluxDet    ! 
  Real(kind=8),allocatable,dimension(:)	  :: GlowFluxPhy    ! 
  Real(kind=8),allocatable,dimension(:)	  :: GlowFluxDia    ! 

  Real(kind=8),allocatable,dimension(:,:) :: addtiny
  Real(kind=8),allocatable,dimension(:,:) :: Gloaddtiny
  Real(kind=8),allocatable,dimension(:,:) :: auxy 

!  Real(kind=8),allocatable,dimension(:,:)   :: GlowFlux         ! 
  Real(kind=8),allocatable,dimension(:)   :: diags2D          ! Diagnostics in 2D [8 n2d]
  Real(kind=8),allocatable,dimension(:,:) :: diags3D          ! Diagnostics in 3D [2 nl-1 n2d]
  Real(kind=8)				  :: DenitBen         ! Benthic denitrification Field in 2D [n2d 1]

  Real(kind=8)				  :: Alk_surf         ! Surface alkalinity field used for restoring
  Real(kind=8)				  :: relax_alk
  Real(kind=8)                            :: virtual_alk
  real(kind=8), allocatable,dimension(:)  :: PAR3D           ! Light in the water column [nl-1 n2d]

!! Cobeta, Cos(Angle of incidence)
  Real(kind=8) ::  cosAI
   
end module REcoM_GloVar

!==========================================================================================================
! For variables saved locally for each column and then used in REcoM
!-------------------------------------------------------------------------------
Module REcoM_locVar

  Real(kind=8),allocatable,dimension(:) :: LocBenthos ! Storing the values for benthos in current watercolumn: N,C,Si and Calc
  Real(kind=8) :: Hplus                     ! [mol/kg] Concentrations of H-plus ions in the surface node
  Real(kind=8) :: pCO2surf(1)                  ! [uatm] Partial pressure of CO2 in surface layer at current 2D node	
  Real(kind=8) :: dflux(1)                     ! [mmol/m2/day] Flux of CO2 into the ocean
  Real(kind=8) :: o2ex(1)                     ! [mmol/m2/s] Flux of O2 into the ocean
  Real(kind=8) :: ULoc(1)                      ! Wind strength above current 2D node, change array size if used with mocsy input vector longer than one
  Real(kind=8) :: dpCO2surf(1)              ! [uatm] difference of oceanic pCO2 minus atmospheric pCO2

! mocsy output -----------------------------------------------------------------------------------------------------------------------------
  Real(kind=8) :: co2flux(1)                   ! air-to-sea flux of CO2 [mol/(m^2 * s)]
  Real(kind=8) :: co2ex(1)                     ! time rate of change of surface CO2 due to gas exchange [mol/(m^3 * s)]
  Real(kind=8) :: dpco2(1)                     ! difference of oceanic pCO2 minus atmospheric pCO2 [uatm]
  Real(kind=8) :: ph(1)                        ! pH on total scale
  Real(kind=8) :: pco2(1)                      ! oceanic partial pressure of CO2 (uatm)
  Real(kind=8) :: fco2(1)                      ! oceanic fugacity of CO2 (uatm)
  Real(kind=8) :: co2(1)                       ! aqueous CO2 concentration [mol/m^3]
  Real(kind=8) :: hco3(1)                      ! bicarbonate (HCO3-) concentration [mol/m^3]
  Real(kind=8) :: co3(1)                       ! carbonate (CO3--) concentration [mol/m^3]
  Real(kind=8) :: OmegaA(1)                    ! Omega for aragonite, i.e., the aragonite saturation state
  Real(kind=8) :: OmegaC(1)                    ! Omega for calcite, i.e., the   calcite saturation state
  Real(kind=8) :: BetaD(1)                     ! BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
  Real(kind=8) :: rhoSW(1)                     ! rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
  Real(kind=8) :: p(1)                         ! pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
  Real(kind=8) :: tempis(1)                    ! in-situ temperature [degrees C]
  Real(kind=8) :: kw660(1)                     ! gas transfer velocity (piston velocity) for CO2 [m/s] 
  Real(kind=8) :: co2flux_seaicemask(1)        ! air-to-sea flux of CO2 [mmol/m2/s]
  Real(kind=8) :: o2flux_seaicemask(1)        ! air-to-sea flux of CO2 [mmol/m2/s]
!-------------------------------------------------------------------------------

  Real(kind=8) :: bt, dic_molal, talk_molal ! Common block: Species
  Real(kind=8) :: k1, k2, kw, kb, ff        ! Common block: Equilibrium_constants
  Real(kind=8) :: Loc_ice_conc(1)           ! Used to calculate flux of DIC in REcoM 0 -> 1
  Real(kind=8) :: LocAtmCO2(1)              ! [uatm]
  Real(kind=8) :: LocDiags2D(8)
  Real(kind=8) :: LocDenit
!  if (REcoM_Second_Zoo) then
  Real(kind=8) :: res_zoo2_a, res_zoo2_f
  Real(kind=8) :: grazingFluxcarbonzoo2                      ! grazingfluxcarbon 
!  endif
  Integer      :: currentCO2year
	
end module REcoM_LocVar
