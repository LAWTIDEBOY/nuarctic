module REcoM_diagnostics

 use general_config
 use recom_config
 use ocean_module
 
 implicit none
 ! unit definition
 character(len=*), parameter :: m_unit='m', d_unit='days' 
 character(len=*), parameter :: P_unit='uatm', f_unit='mmolC/m2/d'
 character(len=*), parameter :: aFe_unit='umolFe/m2/s', aN_unit='mmolN/m2/s'
 character(len=*), parameter :: qm_unit='mol/kg', qa_unit='mmol/m2'
 character(len=*), parameter :: T_unit= 'C' , S_unit ='psu', W_unit = 'W/m2'
 character(len=*), parameter :: Conc_unit = 'mmol/m3'
 
 ! storage array definition
 integer, dimension(:), allocatable 	  :: day_diag, year_diag, nlevel_diag
 real(kind=8), dimension(:) , allocatable :: zcell_diag, znode_diag
 
 ! C02 related outputs
 real(kind=8), dimension(:) , allocatable :: dpCO2s_diag, pCO2s_diag, CO2f_diag
 ! H+
 real(kind=8), dimension(:) , allocatable :: Hp_diag
 ! atm input 
 real(kind=8), dimension(:) , allocatable :: aFe_diag, aN_diag
 
 ! benthos-related outputs
 real(kind=8), dimension(:) , allocatable :: denb_diag, benN_diag, benC_diag, &
 						benSi_diag, benCalc_diag
 ! NPP and GPP diatoms and nano-PhytoPlankton
 real(kind=8), dimension(:) , allocatable :: NPPn_diag, NPPd_diag, 		&
 						GPPn_diag, GPPd_diag
 ! N-assimilation diatoms and nano-Phytoplankton
 real(kind=8), dimension(:) , allocatable :: NNAn_diag, NNAd_diag,		&
  						GNAn_diag, GNAd_diag
  						
 !tracers output
 ! water column temperature and salinity, PAR 
 real(kind=8), dimension(:,:) , allocatable :: temp_diag, sali_diag, PAR_diag
 ! DIN and DIC
 real(kind=8), dimension(:,:) , allocatable :: DIN_diag, DIC_diag
 ! alkalinity
 real(kind=8), dimension(:,:) , allocatable :: Alk_diag
 ! phytoplankton
 real(kind=8), dimension(:,:) , allocatable :: PhyN_diag, PhyC_diag,		&
 						  PhyChl_diag, PhyCalc_diag
 ! Detritus
 real(kind=8), dimension(:,:) , allocatable :: DetN_diag, DetC_diag,		&
 						  DetSi_diag, DetCalc_diag
 ! Heterotrophs
 real(kind=8), dimension(:,:) , allocatable :: HetN_diag, HetC_diag
 ! DOC and DON
 real(kind=8), dimension(:,:) , allocatable :: DON_diag, DOC_diag
 ! other dissolved matters
 real(kind=8), dimension(:,:) , allocatable :: DSi_diag, DFe_diag
 ! Diatoms
 real(kind=8), dimension(:,:) , allocatable :: DiaC_diag, DiaN_diag, 		&
 						  DiaChl_diag, DiaSi_diag
 ! oxygen
 real(kind=8), dimension(:,:) , allocatable :: O2_diag
 ! zooplankton
 real(kind=8), dimension(:,:) , allocatable :: Zoo2N_diag, Zoo2C_diag
 ! 2nd class of zooplankton
 real(kind=8), dimension(:,:) , allocatable :: idetz2n_diag, idetz2c_diag,	&
 						  idetz2si_diag, idetz2calc_diag
 						  
 ! 
 integer				    :: ndiag
 logical, dimension(:), allocatable         :: mask_diagnostic
 integer, dimension(:), allocatable	    :: index_diagnostic
 						  
contains

!----------------------------------------------------------------------------------
subroutine diagnostics_properties
! 
! compute number of diagnostics steps and mask
! 
! estimate computational time in seconds
  use recom_setup
  
  implicit none
  integer     :: fac, it, is, fac_day
  
  allocate(mask_diagnostic(nsteps), index_diagnostic(nsteps))
  mask_diagnostic = .False.
  index_diagnostic = 0
  fac_day=int(step_per_day/24.)
  if (diag_freq_unit=='h') then
	fac=1*fac_day
  elseif (diag_freq_unit=='d') then
  	fac=24*fac_day
  else
  	print*, 'not coded yet only daily or sub-daily diagnostics available'
  	fac=24
  endif 
  write(*,*) 'simulation time step (h):', 24./step_per_day, 'diagnostic frequency (h)', fac/fac_day
  ! if comptutational timestep is in seconds 
  if(run_length_unit=='s') fac = fac * 3600
  it = 0
  do is = 1,nsteps
  	if (mod(is, fac)==0) then
  		it=it+1
  		mask_diagnostic(is)=.True.
  		index_diagnostic(is)=it
  	endif
  enddo
  ndiag = it
end subroutine
					  					  
subroutine setup_diagnostics(boolean)
! 
! (de)allocate and initialize diagnostics
! 

 logical, intent(in) :: boolean
 
 if (boolean) then
	call diagnostics_properties
 	! allocate 
 	allocate(day_diag(ndiag), year_diag(ndiag), nlevel_diag(ndiag))
 	allocate(zcell_diag(nl-1), znode_diag(nl))
 	allocate(dpCO2s_diag(ndiag), pCO2s_diag(ndiag), CO2f_diag(ndiag))
 	allocate(Hp_diag(ndiag), aFe_diag(ndiag), aN_diag(ndiag))
 	allocate(denb_diag(ndiag), benN_diag(ndiag), benC_diag(ndiag),		&
 		benSi_diag(ndiag), benCalc_diag(ndiag))
 	allocate(NPPn_diag(ndiag), NPPd_diag(ndiag), GPPn_diag(ndiag), 		&
 		GPPd_diag(ndiag))
 	allocate(NNAn_diag(ndiag), NNAd_diag(ndiag), GNAn_diag(ndiag), 		&
 		GNAd_diag(ndiag))
 	! 
 	allocate(temp_diag(ndiag, nl-1), sali_diag(ndiag, nl-1), 			&
 		PAR_diag(ndiag, nl-1), O2_diag(ndiag, nl-1))
 	allocate(DIN_diag(ndiag, nl-1), DIC_diag(ndiag, nl-1), Alk_diag(ndiag, nl-1))
 	allocate(PhyN_diag(ndiag, nl-1), PhyC_diag(ndiag, nl-1), 			&
 		PhyChl_diag(ndiag, nl-1), PhyCalc_diag(ndiag, nl-1))
 	allocate(DetN_diag(ndiag, nl-1), DetC_diag(ndiag, nl-1), 			&
 		DetSi_diag(ndiag, nl-1), DetCalc_diag(ndiag, nl-1))
 	allocate(HetN_diag(ndiag, nl-1), HetC_diag(ndiag, nl-1))
 	allocate(DON_diag(ndiag, nl-1), DOC_diag(ndiag, nl-1))
  	allocate(DSi_diag(ndiag, nl-1), DFe_diag(ndiag, nl-1))
 	allocate(DiaC_diag(ndiag, nl-1), DiaN_diag(ndiag, nl), 			&
 		DiaChl_diag(ndiag, nl-1), DiaSi_diag(ndiag, nl-1))
 	if (REcoM_Second_Zoo) allocate(Zoo2N_diag(ndiag, nl-1), Zoo2C_diag(ndiag, nl-1))	
  	if (REcoM_Second_Zoo) allocate(idetz2n_diag(ndiag, nl), idetz2c_diag(ndiag, nl-1), 	&
  			idetz2si_diag(ndiag, nl-1), idetz2calc_diag(ndiag, nl-1))	 	
 	
 	
 	! initialize to zero
 	day_diag        = 0
 	year_diag       = 0 
 	nlevel_diag     = 0
 	
 	znode_diag      = 0.d0
 	zcell_diag      = 0.d0
 	
 	dpCO2s_diag	= 0.d0
 	pCO2s_diag	= 0.d0
 	CO2f_diag	= 0.d0
 	Hp_diag		= 0.d0
 	aFe_diag	= 0.d0
 	aN_diag		= 0.d0
 	denb_diag	= 0.d0
 	benN_diag	= 0.d0
 	benC_diag	= 0.d0
 	benSi_diag	= 0.d0
 	benCalc_diag	= 0.d0
	NPPn_diag	= 0.d0
	NPPd_diag	= 0.d0
	GPPn_diag	= 0.d0
	GPPd_diag	= 0.d0
	NNAn_diag	= 0.d0
	NNAd_diag	= 0.d0
	GNAn_diag	= 0.d0
	GNAd_diag	= 0.d0
	
	temp_diag	= 0.d0
	sali_diag	= 0.d0
	PAR_diag	= 0.d0
	O2_diag		= 0.d0
	DIN_diag	= 0.d0
	DIC_diag	= 0.d0
	Alk_diag	= 0.d0
	PhyN_diag	= 0.d0
	PhyC_diag	= 0.d0
	PhyChl_diag	= 0.d0
	PhyCalc_diag	= 0.d0
	DetN_diag	= 0.d0
	DetC_diag	= 0.d0
	DetSi_diag	= 0.d0
	DetCalc_diag	= 0.d0
	HetN_diag	= 0.d0
	HetC_diag	= 0.d0
	DON_diag	= 0.d0
	DOC_diag	= 0.d0
	DSi_diag	= 0.d0
	DFe_diag	= 0.d0
	DiaC_diag	= 0.d0
	DiaN_diag	= 0.d0
	DiaChl_diag	= 0.d0
	DiaSi_diag	= 0.d0

	if (REcoM_Second_Zoo) then
		Zoo2N_diag	= 0.d0
		Zoo2C_diag	= 0.d0
		idetz2n_diag	= 0.d0
		idetz2c_diag	= 0.d0
		idetz2si_diag	= 0.d0
		idetz2calc_diag	= 0.d0
	endif
  else
  	deallocate(day_diag, year_diag, nlevel_diag, zcell_diag, znode_diag)
  	deallocate(dpCO2s_diag, pCO2s_diag, CO2f_diag)
  	deallocate(Hp_diag, aFe_diag, aN_diag)
  	deallocate(denb_diag, benN_diag, benC_diag, benSi_diag, benCalc_diag)
  	deallocate(NPPn_diag, NPPd_diag, GPPn_diag, GPPd_diag)
  	deallocate(NNAn_diag, NNAd_diag, GNAn_diag, GNAd_diag)
  	deallocate(temp_diag, sali_diag, PAR_diag, O2_diag)
  	deallocate(DIN_diag, DIC_diag, Alk_diag)
  	deallocate(PhyN_diag, PhyC_diag, PhyChl_diag, PhyCalc_diag)
  	deallocate(DetN_diag, DetC_diag, DetSi_diag, DetCalc_diag)
  	deallocate(HetN_diag, HetC_diag)
  	deallocate(DON_diag, DOC_diag)
  	deallocate(DSi_diag, DFe_diag) 
  	deallocate(DiaC_diag, DiaN_diag, DiaChl_diag, DiaSi_diag)
  	if (REcoM_Second_Zoo) deallocate(Zoo2N_diag, Zoo2C_diag, idetz2n_diag, 		&
  					idetz2c_diag, idetz2si_diag, idetz2calc_diag)
  	deallocate(mask_diagnostic, index_diagnostic)
 endif   	

end subroutine
!----------------------------------------------------------------------------------

subroutine store_diagnostics(idiag)
   use REcoM_GloVar
   use recom_clock
   
   implicit none
   integer, intent(in) :: idiag
   integer 	       :: nlvl, flaglyr
   	
   nlvl = nlevel-1
   nlevel_diag(idiag)		= nlvl
   year_diag(idiag)		= yearnew
   day_diag(idiag)		= daynew + int(timenew/86400.)
   	
   ! deal with leap year
   call check_fleapyr(year_diag(idiag), flaglyr)
   if (day_diag(idiag)>365 .and. flaglyr<1) then
   	day_diag(idiag) = 1
   	year_diag(idiag) = year_diag(idiag) + 1
   endif
   	
   ! deal with depth and level
   if (idiag.eq.1) then
   	zcell_diag		= Z
   	znode_diag		= zbar
   endif
   !
   !zcell_diag, znode_diag
   	
   dpCO2s_diag(idiag)	= GlodPCO2surf
   pCO2s_diag(idiag)	= GloPCO2surf
   CO2f_diag(idiag)	= GloCO2flux
   Hp_diag(idiag)		= GloHplus
   aFe_diag(idiag)		= AtmFeInput
   aN_diag(idiag)		= AtmNInput
   denb_diag(idiag)	= DenitBen
   benN_diag(idiag)	= Benthos(1)
   benC_diag(idiag)	= Benthos(2)
   benSi_diag(idiag)	= Benthos(3)
   benCalc_diag(idiag)	= Benthos(4)
   NPPn_diag(idiag)	= diags2D(1)
   NPPd_diag(idiag)	= diags2D(2)
   GPPn_diag(idiag)	= diags2D(3)
   GPPd_diag(idiag)	= diags2D(4)
   NNAn_diag(idiag)	= diags2D(5)
   NNAd_diag(idiag)	= diags2D(6)
   GNAn_diag(idiag)	= diags2D(7)
   GNAd_diag(idiag)	= diags2D(8)
	
   temp_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,1)
   sali_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,2)
   PAR_diag(idiag,1:nlvl)		= PAR(1:nlvl)
   O2_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,24)
   DIN_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,3)
   DIC_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,4)
   Alk_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,5)
   PhyN_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,6)
   PhyC_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,7)
   PhyChl_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,8)
   PhyCalc_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,22)
   DetN_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,9)
   DetC_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,10)
   DetSi_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,19)
   DetCalc_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,23)
   HetN_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,11)
   HetC_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,12)
   DON_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,13)
   DOC_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,14)
   DSi_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,20)
   DFe_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,21)
   DiaC_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,16)
   DiaN_diag(idiag,1:nlvl)		= tr_arr(1:nlvl,15)
   DiaChl_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,17)
   DiaSi_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,18)

   if (REcoM_Second_Zoo) then
	Zoo2N_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,25)
	Zoo2C_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,26)
	idetz2n_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,27)
	idetz2c_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,28)
	idetz2si_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,29)
	idetz2calc_diag(idiag,1:nlvl)	= tr_arr(1:nlvl,30)
    endif

end subroutine
!----------------------------------------------------------------------------------
subroutine write_diagnostics
	
  implicit none
  
  character(len=4096),parameter	:: diag_name='REcoM1d_outputs.nc', units="units", lname="long_name"
  character(len=4096)		:: result_path, filename
  integer, dimension(4)		:: dd
  integer			:: status, fileid
  integer			:: dim_time, dim_cell, dim_node
  ! Id of output variables
  integer			:: year_varid, day_varid, zcell_varid, znode_varid, nlvl_varid  
  integer			:: dpCO2_varid, pCO2_varid, CO2f_varid
  integer			:: Hp_varid, aFe_varid, aN_varid
  integer			:: denb_varid, benN_varid, benC_varid, benSi_varid, benCalc_varid
  integer			:: NPPn_varid, NPPd_varid, GPPn_varid, GPPd_varid
  integer			:: NNAn_varid, NNAd_varid, GNAn_varid, GNAd_varid
  integer			:: T_varid, S_varid, Alk_varid, PAR_varid
  integer			:: O2_varid, DIN_varid, DIC_varid, DOC_Varid, DON_varid
  integer			:: PhyN_varid, PhyC_varid, PhyChl_varid, PhyCalc_varid
  integer			:: DetN_varid, DetC_varid, DetSi_varid, DetCalc_varid
  integer			:: HetN_varid, HetC_varid, DSi_Varid, DFe_varid
  integer			:: DiaN_varid, DiaC_varid, DiaChl_varid, DiaSi_varid
  integer			:: Zoo2N_varid, Zoo2C_varid
  integer			:: detz2n_varid,detz2c_varid, detz2si_varid, detz2calc_varid
  ! name of output variables	
  character(len=20), parameter	:: time_name='time', cell_name='level_cell', node_name='level_node'
  character(len=20), parameter	:: day_name='day', year_name='year', nlvl_name='level'
  character(len=20), parameter	:: znode_name='z_node', zcell_name='z_cell'
  character(len=20), parameter	:: dpCO2s_name = 'dpCO2s', pCO2s_name = 'pCO2s', CO2f_name='CO2f'
  character(len=20), parameter	:: Hp_name='Hp', aFe_name='aFe', aN_name='aN'
  character(len=20), parameter	:: denb_name='denb', benN_name='benN', benC_name = 'benC'
  character(len=20), parameter	:: benSi_name='benSi', benCalc_name='benCalc'
  character(len=20), parameter	:: NPPn_name='NPPn', NPPd_name='NPPd', GPPn_name='GPPn', GPPd_name='GPPd'
  character(len=20), parameter	:: NNAn_name='NNAn', NNAd_name='NNAd', GNAd_name='GNAd', GNAn_name='GNAn'
  character(len=20), parameter	:: temp_name='temp', salt_name='salt', PAR_name='PAR'
  character(len=20), parameter	:: DIN_name='DIN', DIC_name='DIC', ALk_name='Alk', O2_name='O2'
  character(len=20), parameter	:: PhyN_name='PhyN', PhyC_name='PhyC', PhyChl_name='PhyChl', PhyCalc_name='PhyCalc'
  character(len=20), parameter	:: DetN_name='DetN', DetC_name='DetC', DetSi_name='DetSi', DetCalc_name='DetCalc'
  character(len=20), parameter	:: HetN_name='HetN', HetC_name='HetC', DSi_name='DSi', DFe_name='DFe'
  character(len=20), parameter	:: DON_name='DON', DOC_name='DOC'
  character(len=20), parameter	:: DiaN_name='DiaN', DiaC_name='DiaC', DiaChl_name='DiaChl', DiaSi_name='DiaSi'
  character(len=20), parameter	:: Zoo2N_name='Zoo2N', Zoo2C_name='Zoo2C'
  character(len=20), parameter	:: detz2n_name='idetz2n', detz2c_name='idetz2c'
  character(len=20), parameter	:: detz2si_name='idetz2si', detz2calc_name='idetz2calc'
  ! variable long names
  character(len=12), parameter	:: day_lname="calendar day", year_lname="year"
  character(len=12), parameter	:: nlvl_lname="series of nb of active depth level in the computation "
  character(len=50), parameter	:: znode_lname="z at vertical nodes", zcell_lname="z at center of vertical cells"
  character(len=50), parameter	:: dpCO2s_lname = "Difference of oceanic pCO2 minus atmospheric pCO2"
  character(len=50), parameter	:: pCO2s_lname = "Partial pressure of oceanic CO2"
  character(len=50), parameter	:: CO2f_lname="CO2-flux into the surface water"
  character(len=50), parameter	:: Hp_lname="Mean of H-plus ions in the surface water"
  character(len=50), parameter	:: aFe_lname="Atmospheric iron input", aN_lname="Atmospheric DIN input"
  character(len=50), parameter	:: denb_lname="Benthic denitrification rate"
  character(len=50), parameter	:: benN_lname="Benthos Nitrogen", benC_lname = "Benthos Carbon"
  character(len=50), parameter	:: benSi_lname="Benthos silicon", benCalc_lname="Benthos calcite"
  character(len=50), parameter	:: NPPn_lname="Mean NPP nanophytoplankton", NPPd_lname="Mean NPP diatoms"
  character(len=50), parameter	:: GPPn_lname="Mean GPP nanophytoplankton", GPPd_lname="Mean GPP diatoms"
  character(len=50), parameter	:: NNAn_lname="Net N-assimilation nanophytoplankton", NNAd_lname="Net N-assimilation diatoms"
  character(len=50), parameter	:: GNAd_lname="Gross N-assimilation diatoms", GNAn_lname="Gross N-assimilation nanophytoplankton"
  character(len=50), parameter	:: temp_lname="Temperature profile", salt_lname="Salinity profile"
  character(len=50), parameter	:: PAR_lname="Photosynthecally Active Radiation"
  character(len=50), parameter	:: DIN_lname="Dissolved Inorganic Nitrogen", DIC_lname="Dissolved Inorganic Carbon"
  character(len=50), parameter	:: Alk_lname="Total Alkalinity", O2_lname="concentration of dioxygen in water"
  character(len=50), parameter	:: PhyN_lname="Intracellular concentration of Nitrogen in small phytoplankton"
  character(len=50), parameter	:: PhyC_lname="Intracellular concentration of Carbon in small phytoplankton" 
  character(len=50), parameter	:: PhyChl_lname="Current intracellular Chlorophyl A concentration"
  character(len=50), parameter	:: PhyCalc_lname="Current intracellular Calcite concentration"
  character(len=50), parameter	:: DetN_lname="Concentration of Nitrogen in Detritus", DetC_lname="Concentration of Carbon in Detritus"
  character(len=50), parameter	:: DetSi_lname="Concentration of Silicon in Detritus", DetCalc_lname="Concentration of Calcite in Detritus"
  character(len=50), parameter	:: HetN_lname="Concentration of Nitrogen in heterotrophs", HetC_lname="Concentration of Carbon in heterotrophs"
  character(len=50), parameter	:: DSi_lname="concentration of Dissolved Silicate", DFe_lname="concentration of Dissolved Iron"
  character(len=50), parameter	:: DON_lname="Dissolved organic Nitrogen in water", DOC_lname="Dissolved organic Carbon in water"
  character(len=50), parameter	:: DiaN_lname="concentration of Nitrogen in diatoms", DiaC_lname="concentration of Carbon in diatoms"
  character(len=50), parameter	:: DiaChl_lname="concentration of Chlorophyll A in diatoms", DiaSi_lname="concentration of Silicon in diatoms"
  character(len=50), parameter	:: Zoo2N_lname="Intracellular concentration of Nitrogen in second zooplankton"
  character(len=50), parameter	:: Zoo2C_lname="Intracellular concentration of Carbon in second zooplankton"
  character(len=50), parameter	:: detz2n_lname="Concentration of Nitrogen in detritus from second zooplankton"
  character(len=50), parameter	:: detz2c_lname="Concentration of Carbon in detritus from second zooplankton"
  character(len=50), parameter	:: detz2si_lname="Concentration of Silicon in detritus from second zooplankton"
  character(len=50), parameter	:: detz2calc_lname="Concentration of Calcite in detritus from second zooplankton"


#include "netcdf.inc" 
  
  call get_environment_variable("RECOM_RESULT_PATH", result_path)
  filename = trim(result_path) // trim(diag_name)
  dd=0
  
  print*,'number of diagnostic steps',ndiag
  ! open netCDF file
  status=nf_create(filename, NF_CLOBBER, fileid)

  ! define output array dimensions
  status = nf_def_dim(fileid,trim(time_name),ndiag, dim_time)
  status = nf_def_dim(fileid, trim(cell_name),nl-1,dim_cell)
  status = nf_def_dim(fileid, trim(node_name),nl,dim_node)


  ! define variable
  !dates
  dd(1)= dim_time
  status = nf_def_var(fileid,trim(day_name), nf_int,1, dd, day_varid)
  status = nf_def_var(fileid,trim(year_name), nf_int,1, dd, year_varid) 
  status = nf_def_var(fileid,trim(nlvl_name), nf_int,1, dd, nlvl_varid) 
  ! level and depth information: depth proxies z at center of cell and at cell nodes
  dd(1) = dim_node
  status = nf_def_var(fileid,trim(znode_name), nf_real,1, dd, znode_varid)	
  dd(1) = dim_cell
  status = nf_def_var(fileid,trim(zcell_name), nf_real,1, dd, zcell_varid)
  !print*, 'test', ndiag, nl, nl-1, dim_time, dim_cell, dim_node
  
  ! time dependent variables
  dd(1)=dim_time
  ! time series
  status=nf_def_var(fileid, trim(dpCO2s_name), nf_real, 1, dd, dpCO2_varid)
  status=nf_def_var(fileid, trim(pCO2s_name), nf_real, 1, dd, pCO2_varid)
  status=nf_def_var(fileid, trim(CO2f_name), nf_real, 1, dd, CO2f_varid)
  status=nf_def_var(fileid, trim(Hp_name), nf_real, 1, dd, Hp_varid)
  status=nf_def_var(fileid, trim(aFe_name), nf_real, 1, dd, aFe_varid)
  status=nf_def_var(fileid, trim(aN_name), nf_real, 1, dd, aN_varid)
  
  status=nf_def_var(fileid, trim(denb_name), nf_real, 1, dd, denb_varid)
  status=nf_def_var(fileid, trim(benN_name), nf_real, 1, dd, benN_varid)
  status=nf_def_var(fileid, trim(benC_name), nf_real, 1, dd, benC_varid)
  status=nf_def_var(fileid, trim(benSi_name), nf_real, 1, dd, benSi_varid)
  status=nf_def_var(fileid, trim(benCalc_name), nf_real, 1, dd, benCalc_varid)
  
  ! diagnostics
  status=nf_def_var(fileid, trim(NPPn_name), nf_real, 1, dd, NPPn_varid)
  status=nf_def_var(fileid, trim(NPPd_name), nf_real, 1, dd, NPPd_varid)
  status=nf_def_var(fileid, trim(GPPn_name), nf_real, 1, dd, GPPn_varid)
  status=nf_def_var(fileid, trim(GPPd_name), nf_real, 1, dd, GPPd_varid)
  
  status=nf_def_var(fileid, trim(NNAn_name), nf_real, 1, dd, NNAn_varid)
  status=nf_def_var(fileid, trim(NNAd_name), nf_real, 1, dd, NNAd_varid)
  status=nf_def_var(fileid, trim(GNAn_name), nf_real, 1, dd, GNAn_varid)
  status=nf_def_var(fileid, trim(GNAd_name), nf_real, 1, dd, GNAd_varid)
  
  ! depth at middle of cell
  dd(1) = dim_time
  dd(2) = dim_cell
  
  ! state variables(time and depth dependent, mainly tracers)
  status=nf_def_var(fileid, trim(PAR_name), nf_real, 2, dd, PAR_varid)
  
  status=nf_def_var(fileid, trim(temp_name), nf_real, 2, dd, T_varid)
  status=nf_def_var(fileid, trim(salt_name), nf_real, 2, dd, S_varid)  
 
  status=nf_def_var(fileid, trim(DIN_name), nf_real, 2, dd, DIN_varid) 
  status=nf_def_var(fileid, trim(DIC_name), nf_real, 2, dd, DIC_varid)
  status=nf_def_var(fileid, trim(Alk_name), nf_real, 2, dd, Alk_varid) 
  
  status=nf_def_var(fileid, trim(PhyN_name), nf_real, 2, dd, PhyN_varid)  
  status=nf_def_var(fileid, trim(PhyC_name), nf_real, 2, dd, PhyC_varid)
  status=nf_def_var(fileid, trim(PhyChl_name), nf_real, 2, dd, PhyChl_varid)
  
  status=nf_def_var(fileid, trim(DetN_name), nf_real, 2, dd, DetN_varid)
  status=nf_def_var(fileid, trim(DetC_name), nf_real, 2, dd, DetC_varid)
  
  status=nf_def_var(fileid, trim(HetN_name), nf_real, 2, dd, HetN_varid)
  status=nf_def_var(fileid, trim(HetC_name), nf_real, 2, dd, HetC_varid)  
  
  status=nf_def_var(fileid, trim(DON_name), nf_real, 2, dd, DON_varid)
  status=nf_def_var(fileid, trim(DOC_name), nf_real, 2, dd, DOC_varid)
  
  status=nf_def_var(fileid, trim(DiaN_name), nf_real, 2, dd, DiaN_varid)
  status=nf_def_var(fileid, trim(DiaC_name), nf_real, 2, dd, DiaC_varid) 
  status=nf_def_var(fileid, trim(DiaChl_name), nf_real, 2, dd, DiaChl_varid)
  status=nf_def_var(fileid, trim(DiaSi_name), nf_real, 2, dd, DiaSi_varid) 
  
  status=nf_def_var(fileid, trim(DetSi_name), nf_real, 2, dd, DetSi_varid)
  
  status=nf_def_var(fileid, trim(DSi_name), nf_real, 2, dd, DSi_varid)
  status=nf_def_var(fileid, trim(DFe_name), nf_real, 2, dd, DFe_varid)
  
  status=nf_def_var(fileid, trim(PhyCalc_name), nf_real, 2, dd, Phycalc_varid) 
  status=nf_def_var(fileid, trim(DetCalc_name), nf_real, 2, dd, DetCalc_varid)
  
  status=nf_def_var(fileid, trim(O2_name), nf_real, 2, dd, O2_varid)
 
  if (REcoM_Second_Zoo) then
    	status=nf_def_var(fileid, trim(Zoo2N_name), nf_real, 2, dd, Zoo2N_varid)
  	status=nf_def_var(fileid, trim(Zoo2C_name), nf_real, 2, dd, Zoo2C_varid) 
  	status=nf_def_var(fileid, trim(detz2n_name), nf_real, 2, dd, detz2n_varid)
  	status=nf_def_var(fileid, trim(detz2c_name), nf_real, 2, dd, detz2c_varid) 
  	status=nf_def_var(fileid, trim(detz2si_name), nf_real, 2, dd, detz2si_varid)
  	status=nf_def_var(fileid, trim(detz2calc_name), nf_real, 2, dd, detz2calc_varid)   
  endif   

  ! define attributes (units, variable names) long_name, unit
  ! dates, depth proxies day and year
  ! dates

  status = nf_put_att_text(fileid, day_varid, trim(lname), len(trim(day_lname)), trim(day_lname))
  status = nf_put_att_text(fileid, year_varid, trim(lname), len(trim(year_lname)), trim(year_lname))
  status = nf_put_att_text(fileid, nlvl_varid, trim(lname), len(trim(nlvl_lname)), trim(nlvl_lname))  
  	 
  ! depth proxy at middle of cells
  status = nf_put_att_text(fileid, zcell_varid, trim(units), len(trim(m_unit)), m_unit)
  status = nf_put_att_text(fileid, zcell_varid, trim(lname), len(trim(zcell_lname)), trim(zcell_lname)) 
  
  ! depth at vertical nodes
  status = nf_put_att_text(fileid, znode_varid, trim(units), len(trim(m_unit)), m_unit)
  status = nf_put_att_text(fileid, znode_varid, trim(lname), len(trim(znode_lname)), trim(znode_lname))
  	  
  ! dpCO2
  status = nf_put_att_text(fileid, dpCO2_varid, trim(units), len(trim(P_unit)), P_unit)
  status = nf_put_att_text(fileid, dpCO2_varid, trim(lname), len(trim(dpCO2s_lname)), trim(dpCO2s_lname)) 

  ! pCO2
  status = nf_put_att_text(fileid, pCO2_varid, trim(units), len(trim(P_unit)), P_unit)
  status = nf_put_att_text(fileid, pCO2_varid, trim(lname), len(trim(pCO2s_lname)), trim(pCO2s_lname)) 

  ! CO2f
  status = nf_put_att_text(fileid, CO2f_varid, trim(units), len(trim(f_unit)), f_unit)
  status = nf_put_att_text(fileid, CO2f_varid, trim(lname), len(trim(CO2f_lname)), trim(CO2f_lname)) 
  ! Hp
  status = nf_put_att_text(fileid, Hp_varid, trim(units), len(trim(qm_unit)), qm_unit)
  status = nf_put_att_text(fileid, Hp_varid, trim(lname), len(trim(Hp_lname)), trim(Hp_lname))
  ! aFe
  status = nf_put_att_text(fileid, aFe_varid, trim(units), len(trim(aFe_unit)), aFe_unit)
  status = nf_put_att_text(fileid, aFe_varid, trim(lname), len(trim(aFe_lname)), trim(aFe_lname))   
  
  ! aN
  status = nf_put_att_text(fileid, aN_varid, trim(units), len(trim(aN_unit)), aN_unit)
  status = nf_put_att_text(fileid, aN_varid, trim(lname), len(trim(aN_lname)), trim(aN_lname)) 
  
  ! denb
  status = nf_put_att_text(fileid, denb_varid, trim(units), len(trim(qa_unit)), qa_unit)
  status = nf_put_att_text(fileid, denb_varid, trim(lname), len(trim(denb_lname)), trim(denb_lname)) 
  
  ! benN
  status = nf_put_att_text(fileid, benN_varid, trim(units), len(trim(qa_unit)), qa_unit)
  status = nf_put_att_text(fileid, benN_varid, trim(lname), len(trim(benN_lname)), trim(benN_lname)) 
    
  ! benC
  status = nf_put_att_text(fileid, benC_varid, trim(units), len(trim(qa_unit)), qa_unit)
  status = nf_put_att_text(fileid, benC_varid, trim(lname), len(trim(benC_lname)), trim(benC_lname)) 
    
  ! benSi
  status = nf_put_att_text(fileid, benSi_varid, trim(units), len(trim(qa_unit)), qa_unit)
  status = nf_put_att_text(fileid, benSi_varid, trim(lname), len(trim(benSi_lname)), trim(benSi_lname)) 
    
  ! benCalc
  status = nf_put_att_text(fileid, benCalc_varid, trim(units), len(trim(qa_unit)), qa_unit)
  status = nf_put_att_text(fileid, benCalc_varid, trim(lname), len(trim(benCalc_lname)), trim(benCalc_lname)) 
  
  ! NPPn
  status = nf_put_att_text(fileid, NPPn_varid, trim(units), len(trim(f_unit)), f_unit)
  status = nf_put_att_text(fileid, NPPn_varid, trim(lname), len(trim(NPPn_lname)), trim(NPPn_lname))
  
  ! NPPd
  status = nf_put_att_text(fileid, NPPd_varid, trim(units), len(trim(f_unit)), f_unit)
  status = nf_put_att_text(fileid, NPPd_varid, trim(lname), len(trim(NPPd_lname)), trim(NPPd_lname))
    
  ! GPPn
  status = nf_put_att_text(fileid, GPPn_varid, trim(units), len(trim(f_unit)), f_unit)
  status = nf_put_att_text(fileid, GPPn_varid, trim(lname), len(trim(GPPn_lname)), trim(GPPn_lname))
    
  ! GPPd
  status = nf_put_att_text(fileid, GPPd_varid, trim(units), len(trim(f_unit)), f_unit)
  status = nf_put_att_text(fileid, GPPd_varid, trim(lname), len(trim(GPPd_lname)), trim(GPPd_lname))
    
  ! NNAn
  status = nf_put_att_text(fileid, NNAn_varid, trim(units), len(trim(f_unit)), f_unit)
  status = nf_put_att_text(fileid, NNAn_varid, trim(lname), len(trim(NNAn_lname)), trim(NNAn_lname))
    
  ! NNAd
  status = nf_put_att_text(fileid, NNAd_varid, trim(units), len(trim(f_unit)), f_unit)
  status = nf_put_att_text(fileid, NNAd_varid, trim(lname), len(trim(NNAd_lname)), trim(NNAd_lname))  

  ! GNAn
  status = nf_put_att_text(fileid, GNAn_varid, trim(units), len(trim(f_unit)), f_unit)
  status = nf_put_att_text(fileid, GNAn_varid, trim(lname), len(trim(GNAn_lname)), trim(GNAn_lname))
    
  ! GNAd
  status = nf_put_att_text(fileid, GNAd_varid, trim(units), len(trim(f_unit)), f_unit)
  status = nf_put_att_text(fileid, GNAd_varid, trim(lname), len(trim(GNAd_lname)), trim(GNAd_lname))  
  
  ! temp
  status = nf_put_att_text(fileid, T_varid, trim(units), len(trim(T_unit)), T_unit)
  status = nf_put_att_text(fileid, T_varid, trim(lname), len(trim(temp_lname)), trim(temp_lname))   
  
  ! salt
  status = nf_put_att_text(fileid, S_varid, trim(units), len(trim(S_unit)), S_unit)
  status = nf_put_att_text(fileid, S_varid, trim(lname), len(trim(salt_lname)), trim(salt_lname)) 
    
  ! PAR
  status = nf_put_att_text(fileid, PAR_varid, trim(units), len(trim(W_unit)), W_unit)
  status = nf_put_att_text(fileid, PAR_varid, trim(lname), len(trim(PAR_lname)), trim(PAR_lname))   
  
  ! DIN
  status = nf_put_att_text(fileid, DIN_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DIN_varid, trim(lname), len(trim(DIN_lname)),trim(DIN_lname))  
   
  ! DIC
  status = nf_put_att_text(fileid, DIC_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DIC_varid, trim(lname), len(trim(DIC_lname)), trim(DIC_lname))  
  
  ! Alk
  status = nf_put_att_text(fileid, Alk_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, Alk_varid, trim(lname), len(trim(Alk_lname)), trim(Alk_lname))   

  ! PhyN
  status = nf_put_att_text(fileid, PhyN_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, PhyN_varid, trim(lname), len(trim(PhyN_lname)), trim(PhyN_lname))   
  
  ! PhyC
  status = nf_put_att_text(fileid, PhyC_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, PhyC_varid, trim(lname), len(trim(PhyC_lname)), trim(PhyC_lname))   
  
  ! PhyChl
  status = nf_put_att_text(fileid, PhyChl_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, PhyChl_varid, trim(lname), len(trim(PhyChl_lname)), trim(PhyChl_lname))   

  ! PhyCalc
  status = nf_put_att_text(fileid, PhyCalc_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, PhyCalc_varid, trim(lname), len(trim(PhyCalc_lname)), trim(PhyCalc_lname))  
    
  ! DetN
  status = nf_put_att_text(fileid, DetN_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DetN_varid, trim(lname), len(trim(DetN_lname)), trim(DetN_lname))  
  
  ! DetC
  status = nf_put_att_text(fileid, DetC_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DetC_varid, trim(lname), len(trim(DetC_lname)), trim(DetC_lname))  

  ! DetSi
  status = nf_put_att_text(fileid, DetSi_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DetSi_varid, trim(lname), len(trim(DetSi_lname)), trim(DetSi_lname))  
  
  ! DetCalc
  status = nf_put_att_text(fileid, DetCalc_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DetCalc_varid, trim(lname), len(trim(DetCalc_lname)), trim(DetCalc_lname))  
 
  ! HetN
  status = nf_put_att_text(fileid, HetN_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, HetN_varid, trim(lname), len(trim(HetN_lname)), trim(HetN_lname))
  
  ! HetC
  status = nf_put_att_text(fileid, HetC_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, HetC_varid, trim(lname), len(trim(HetC_lname)), trim(HetC_lname))  
   
  ! DON
  status = nf_put_att_text(fileid, DON_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DON_varid, trim(lname), len(trim(DON_lname)), trim(DON_lname)) 
    
  ! DOC
  status = nf_put_att_text(fileid, DOC_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DOC_varid, trim(lname), len(trim(DOC_lname)), trim(DOC_lname)) 
  
  ! O2
  status = nf_put_att_text(fileid, O2_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, O2_varid, trim(lname), len(trim(O2_lname)), trim(O2_lname)) 

  ! DiaN
  status = nf_put_att_text(fileid, DiaN_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DiaN_varid, trim(lname), len(trim(DiaN_lname)), trim(DiaN_lname))   
  
  ! DiaC
  status = nf_put_att_text(fileid, DiaC_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DiaC_varid, trim(lname), len(trim(DiaC_lname)), trim(DiaC_lname))   
  
  ! DiaChl
  status = nf_put_att_text(fileid, DiaChl_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DiaChl_varid, trim(lname), len(trim(DiaChl_lname)), trim(DiaChl_lname))   
  
  ! DiaSi
  status = nf_put_att_text(fileid, DiaSi_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DiaSi_varid, trim(lname), len(trim(DiaSi_lname)), trim(DiaSi_lname))   
  
  ! DSi
  status = nf_put_att_text(fileid, DSi_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DSi_varid, trim(lname), len(trim(DSi_lname)), trim(DSi_lname))   
    
  ! DFe
  status = nf_put_att_text(fileid, DFe_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  status = nf_put_att_text(fileid, DFe_varid, trim(lname), len(trim(DFe_lname)), trim(DFe_lname))   


  
  if (REcoM_Second_Zoo) then
  
    	! Zoo2N
  	status = nf_put_att_text(fileid, Zoo2N_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  	status = nf_put_att_text(fileid, Zoo2N_varid, trim(lname), len(trim(Zoo2N_lname)), trim(Zoo2N_lname))   
    
  	! Zoo2C
  	status = nf_put_att_text(fileid, Zoo2C_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  	status = nf_put_att_text(fileid, Zoo2C_varid, trim(lname), len(trim(Zoo2C_lname)), trim(Zoo2C_lname))   
  	
  	! idetz2n
  	status = nf_put_att_text(fileid, detz2n_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  	status = nf_put_att_text(fileid, detz2n_varid, trim(lname), len(trim(detz2n_lname)), trim(detz2n_lname))   
    
  	! idetz2c
  	status = nf_put_att_text(fileid, detz2c_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  	status = nf_put_att_text(fileid, detz2c_varid, trim(lname), len(trim(detz2c_lname)), trim(detz2c_lname))   
    
  	! idetz2si
  	status = nf_put_att_text(fileid, detz2si_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  	status = nf_put_att_text(fileid, detz2si_varid, trim(lname), len(trim(detz2si_lname)), trim(detz2si_lname))     

  	! idetz2calc
  	status = nf_put_att_text(fileid, detz2calc_varid, trim(units), len(trim(Conc_unit)), Conc_unit)
  	status = nf_put_att_text(fileid, detz2calc_varid, trim(lname), len(trim(detz2calc_lname)), trim(detz2calc_lname))   
  endif
 
  ! end definition mode
   status = nf_enddef(fileid)
   
  ! fill and store variable in the netcdf file
  status = nf_put_var_int(fileid, day_varid, day_diag)
  status = nf_put_var_int(fileid, year_varid, year_diag)
  status = nf_put_var_int(fileid, nlvl_varid, nlevel_diag)
  status = nf_put_var_double(fileid, zcell_varid, zcell_diag)
  status = nf_put_var_double(fileid, znode_varid, znode_diag)
  
  status = nf_put_var_double(fileid, dpCO2_varid, dpCO2s_diag)
  status = nf_put_var_double(fileid, pCO2_varid, pCO2s_diag)
  status = nf_put_var_double(fileid, CO2f_varid, CO2f_diag)
  status = nf_put_var_double(fileid, Hp_varid, Hp_diag)
  status = nf_put_var_double(fileid, aFe_varid, aFe_diag)
  status = nf_put_var_double(fileid, aN_varid, aN_diag)    

  status = nf_put_var_double(fileid, denb_varid, denb_diag)
  status = nf_put_var_double(fileid, benN_varid, benN_diag)  
  status = nf_put_var_double(fileid, benC_varid, benC_diag)
  status = nf_put_var_double(fileid, benSi_varid, benSi_diag)    
  status = nf_put_var_double(fileid, benCalc_varid, benCalc_diag)
  
  status = nf_put_var_double(fileid, NPPn_varid, NPPn_diag)
  status = nf_put_var_double(fileid, NPPd_varid, NPPd_diag)  
  status = nf_put_var_double(fileid, GPPn_varid, GPPn_diag)
  status = nf_put_var_double(fileid, GPPd_varid, GPPd_diag) 
  
  status = nf_put_var_double(fileid, NNAn_varid, NNAn_diag)
  status = nf_put_var_double(fileid, NNAd_varid, NNAd_diag)  
  status = nf_put_var_double(fileid, GNAn_varid, GNAn_diag)
  status = nf_put_var_double(fileid, GNAd_varid, GNAd_diag) 
  
  status = nf_put_var_double(fileid, T_varid, temp_diag)  
  status = nf_put_var_double(fileid, S_varid, sali_diag)
  status = nf_put_var_double(fileid, PAR_varid, PAR_diag) 

  status = nf_put_var_double(fileid, DIN_varid, DIN_diag)
  status = nf_put_var_double(fileid, DIC_varid, DIC_diag)  
  status = nf_put_var_double(fileid, Alk_varid, Alk_diag)
  status = nf_put_var_double(fileid, O2_varid, O2_diag) 
  print*, 'O2', O2_diag(10,:)
  
  status = nf_put_var_double(fileid, PhyN_varid, PhyN_diag)
  status = nf_put_var_double(fileid, PhyC_varid, PhyC_diag)  
  status = nf_put_var_double(fileid, PhyChl_varid, PhyChl_diag)
  status = nf_put_var_double(fileid, Phycalc_varid, PhyCalc_diag) 

  status = nf_put_var_double(fileid, DetN_varid, DetN_diag)
  status = nf_put_var_double(fileid, DetC_varid, DetC_diag)  
  status = nf_put_var_double(fileid, DetSi_varid, DetSi_diag)
  status = nf_put_var_double(fileid, DetCalc_varid, DetCalc_diag) 

  status = nf_put_var_double(fileid, HetN_varid, HetN_diag)
  status = nf_put_var_double(fileid, HetC_varid, HetC_diag) 
   
  status = nf_put_var_double(fileid, DSi_varid, DSi_diag)
  status = nf_put_var_double(fileid, DFe_varid, DFe_diag) 
  
  status = nf_put_var_double(fileid, DON_varid, DON_diag)
  status = nf_put_var_double(fileid, DOC_varid, DOC_diag) 
  
  status = nf_put_var_double(fileid, DiaN_varid, DiaN_diag)
  status = nf_put_var_double(fileid, DiaC_varid, DiaC_diag)  
  status = nf_put_var_double(fileid, DiaChl_varid, DiaChl_diag)
  status = nf_put_var_double(fileid, DiaSi_varid, DiaSi_diag)
   
  
  if (REcoM_Second_Zoo) then
    	status = nf_put_var_double(fileid, Zoo2N_varid, Zoo2N_diag)
  	status = nf_put_var_double(fileid, Zoo2C_varid, Zoo2C_diag)  
  	status = nf_put_var_double(fileid, detz2n_varid, idetz2n_diag)
  	status = nf_put_var_double(fileid, detz2c_varid, idetz2c_diag) 
  	status = nf_put_var_double(fileid, detz2si_varid, idetz2si_diag)
  	status = nf_put_var_double(fileid, detz2calc_varid, idetz2calc_diag) 
  endif   

  ! close file
  status=nf_close(fileid)

end subroutine

!----------------------------------------------------------------------------------

end module

