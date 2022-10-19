module recom_init

  use mod_mesh
  use recom_clock
  use recom_declarations
  use recom_GloVar
  use recom_locVar
  use recom_config
  use ocean_module
! CONTENT:
! ------------
!    subroutine recom_init
!
! written by V. Schourup-Kristensen, 
! adapted to fesom2.0 by ozgur gurses, 22.05.2020
! adapted to REcoM1D by F.Birrien 21.12.2021   

contains
!===============================================================================
! allocate & initialise arrays for REcoM
subroutine recom_initialization(mesh)


    ! fesom modules
    !use o_ARRAYS
    !use o_MESH
    implicit none
!#include "netcdf.inc"

    type(t_mesh), intent(in) , target :: mesh
!==============================================================================================================================================
   ! nb_of_nodes=mesh%nb_of_nodes
    nl=mesh%nl
      
    !___allocate________________________________________________________________
    allocate(GlodecayBenthos(benthos_num))
    allocate(PAR3D(nl-1)) 
    allocate(Benthos(benthos_num))
    allocate(LocBenthos(benthos_num))
    allocate(decayBenthos(benthos_num)) ! [1/day] Decay rate of detritus in the benthic layer

!    allocate(wFluxDet(benthos_num))     ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
    allocate(wFluxPhy(benthos_num))     ! [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of phytoplankton
    allocate(wFluxDia(benthos_num))     ! [mmol/(m2 * day)] Flux of N,C, Si and chl through sinking of diatoms 	

!    allocate(GlowFlux(2))     ! 
!    GlowFlux=0.0d0

if (REcoM_Second_Zoo) then
    allocate(GlowFluxDet(benthos_num*2))
    allocate(wFluxDet(benthos_num*2))     ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
else
    allocate(GlowFluxDet(benthos_num))
    allocate(wFluxDet(benthos_num))     ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
end if 
    allocate(GlowFluxPhy(benthos_num))
    allocate(GlowFluxDia(benthos_num))

    GlowFluxDet=0.0d0
    GlowFluxPhy=0.0d0
    GlowFluxDia=0.0d0    

    allocate(addtiny(nl-1,8))
    addtiny=0.0d0
    allocate(Gloaddtiny(nl-1,8))
    allocate(auxy(nl-1,bgc_num))
    Gloaddtiny=0.0d0
    auxy=0.0d0


    !___initialize______________________________________________________________
    ! Dust surface deposition
    GloFeDust = 0.d0
    AtmFeInput = 0.d0

    GloNDust = 0.d0
    AtmNInput = 0.d0
    
    ! 
    cosAI = 0.d0

    ! initialise 2d field of CO2 related diagnostics
    GloPCO2surf = 0.d0
    GloCO2flux = 0.d0
    GloCO2flux_seaicemask = 0.0d0 
    GloO2flux_seaicemask = 0.0d0 
    GlodPCO2surf = 0.d0
    GlodecayBenthos = 0.0d0

    ! initialise 2d field of CO2 related variables 
    ! atmospheric CO2 is read from a file if constant_CO2 is false 
    AtmCO2 = 0.d0 
    Hplus = 0.d0
    pco2surf = 0.d0
    dflux = 0.d0
    co2flux_seaicemask = 0.d0
    o2flux_seaicemask = 0.d0
    dpco2surf= 0.d0

    if (Diags) then
        allocate(diags2D(8))
        diags2D(:)      = 0.d0
        allocate(diags3D(nl-1,diags3d_num))
        diags3D(:,:)    = 0.d0
    end if  

    PAR3D(:) = 0.d0
    DenitBen = 0.d0
    !__________________________________________________________________________
    ! Initialization of benthos
    ! Benthic layer consists of Benthos(1) = N, Benthos(2)=C, Benthos(3)=Si, Benthos(4)=calc
  
    Benthos(:) = 0.d0 !tiny
    GloHplus = exp(-8.d0 * log(10.d0)) ! = 10**(-8)
    !___________________________________________________________________________
    write(*,*) 'Benthic layers are set'

! initialize tracer arrays (previously define in ocean_module)
! DIN(3), DIC(4), Alk(5), DSi(20), DFe, O2(24) are initialize from climatology or MOSAiC initial conditions
! if no specific values are defined then it stays equal to 0
    call read_tracer_initialization(mesh)
    if (sum(abs(DIN_init))>0) tr_arr(:,3)=DIN_init
    if (sum(abs(DIC_init))>0) tr_arr(:,4)=DIC_init
    if (sum(abs(Alk_init))>0) tr_arr(:,5)=Alk_init
    if (sum(abs(DSi_init))>0) tr_arr(:,20)=DSi_init
    if (sum(abs(DFe_init))>0) tr_arr(:,21)=DFe_init
    if (sum(abs(DO2_init))>0) tr_arr(:,24)=DO2_init      
! tracer 6  = PhyN   -> Intracellular conc of Nitrogen in small phytoplankton
    tr_arr(:,6)  = tiny_chl/chl2N_max   !tiny                  
    
! tracer 7  = PhyC   -> Intracellular conc of Carbon in small phytoplankton
    tr_arr(:,7)  = tiny_chl/chl2N_max/NCmax !tiny * Redfield       
    
! tracer 8  = PhyChl -> Current intracellular ChlA conc
    tr_arr(:,8)  = tiny_chl  !tiny * 1.56d0         

! Detritus: Nitrogen, Carbon and Silicate concentration (DetN(9), DetC(10), DetSi(19))
    tr_arr(:,9)  = tiny                  			
    tr_arr(:,10) = tiny               
    tr_arr(:,19) = tiny                     			
! Heterotroph: Nitrogen and Carbon concentration (HetN(11), HetC(12))
    tr_arr(:,11) = tiny                  			
    tr_arr(:,12) = tiny * Redfield       			
! Dissolved inorganic nitrogen and carbon concentration (DON(13), DOC(14))
    tr_arr(:,13) = tiny                  			
    tr_arr(:,14) = tiny                  			
! Diatoms nitrogen, carbon, Chlorophyll and silicate concentration (DiaN(15), DiaC(16), DiaChl(17), DiaSi(18))
    tr_arr(:,15) = tiny_chl/chl2N_max !tiny                  	
    tr_arr(:,16) = tiny_chl/chl2N_max/NCmax !tiny * Redfield 	
    tr_arr(:,17) = tiny_chl !tiny * 1.56d0        		 
    tr_arr(:,18) = tiny_chl/chl2N_max_d/NCmax_d/SiCmax !tiny          

! concentration in Dissolved Iron (DFe) 
    !tr_arr(:,21) = tr_arr(:,21) * 1.e9 ! Fe [mol/L] => [umol/m3] Check the units again!

! calcification: concentration from phytoplankton and Detritus (PhyCalc(22) and DetCalc(23))
    tr_arr(:,22) = tiny !cPhyN * 0.25d0    
    tr_arr(:,23) = tiny


if (REcoM_Second_Zoo) then
    ! 2nd zooplankton group: concentration in Nitrogen and Carbon (Zoo2N(25), Zoo2N(26))
    tr_arr(:,25) = tiny                   
    tr_arr(:,26) = tiny * Redfield
    ! associated Detritus: concentration in Nitrogen, Carbon and Silicate (DetZ2N(27), DetZ2C(28), DetZ2Si(29))  
    tr_arr(:,27) = tiny                               
    tr_arr(:,28) = tiny                                               
    tr_arr(:,29) = tiny
    ! calcification of Detritus (DetZ2Calc(30))                  
    tr_arr(:,30) = tiny   
endif

!---------------------------------------------------------------------------------------------------------
write(*,*) 'Tracers have been initialized as spinup from WOA/glodap netcdf files'

end subroutine recom_initialization

subroutine deallocate_flux

    deallocate(GlodecayBenthos, PAR3D, Benthos, LocBenthos, decayBenthos)   
    deallocate(wFluxPhy, wFluxDia, GlowFluxDet, wFluxDet)
    deallocate(GlowFluxPhy, GlowFluxDia)
    deallocate(addtiny, Gloaddtiny, auxy)

end subroutine

end module recom_init

