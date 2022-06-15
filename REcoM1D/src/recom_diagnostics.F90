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
 real(kind=8), dimension(:) , allocatable :: dates_diag
 real(kind=8), dimension(:,:) , allocatable :: depth_diag, bound_diag
 
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
  integer     :: fac, it, is
  
  allocate(mask_diagnostic(nsteps), index_diagnostic(nsteps))
  mask_diagnostic = .False.
  index_diagnostic = 0
  
  if (diag_freq_unit=='h') then
	fac=1
  elseif (diag_freq_unit=='d') then
  	fac=24
  else
  	print*, 'not coded yet only daily or sub-daily diagnostics available'
  	fac=24
  endif 
   
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
 	allocate(dates_diag(ndiag), depth_diag(ndiag, nl-1), bound_diag(ndiag,nl))
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
  	deallocate(dates_diag, depth_diag, bound_diag)
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
   !!!!!!!!
   !! continue her to propose the daily/monthly mean
   implicit none
   integer, intent(in) :: idiag
   integer 	       :: nlvl
   
   	nlvl = nlevel-1
   	
   	dates_diag(idiag)		= yearnew*1e8 + daynew*1e5 + timenew
   	depth_diag(idiag,1:nlvl)	= hnode(1:nlvl)
   	bound_diag(idiag,1:nlvl+1)	= zbar(1:nlvl+1)
   	
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
  
  character(len=4096),parameter	:: diag_name='diagnostics.nc', units="units", lname="long_name"
  character(len=4096)		:: result_path, filename
  integer, dimension(4)		:: dd
  integer, parameter		:: namelength = 50, unitlength=6
  integer			:: status, fileid
  integer			:: dim_level, dim_time, dim_bound
  integer			:: dates_varid, depth_varid, bound_varid  
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
  logical 			:: test
#include "netcdf.inc" 
  
  call get_environment_variable("RECOM_RESULT_PATH", result_path)
  filename = trim(result_path) // trim(diag_name)
  dd=0
  test=.False.
  print*,'number of diagnostic steps',ndiag
  ! open netCDF file
  status=nf_create(filename, NF_CLOBBER, fileid)

  ! define output array dimensions
  status = nf_def_dim(fileid,'time',ndiag,dim_time)
  status = nf_def_dim(fileid,'level',nl-1,dim_level)
  status = nf_def_dim(fileid,'level_bounds',nl,dim_bound)
  ! define variable
  !dates
  dd(1)= dim_time
  status = nf_def_var(fileid,'dates', nf_real,1, dd, dates_varid)
  
  ! time dependent variables
  status=nf_def_var(fileid, 'dpCO2s', nf_real, 1, dd, dpCO2_varid)
  status=nf_def_var(fileid, 'pCO2s', nf_real, 1, dd, pCO2_varid)
  status=nf_def_var(fileid, 'CO2f', nf_real, 1, dd, CO2f_varid)
  status=nf_def_var(fileid, 'Hp', nf_real, 1, dd, Hp_varid)
  status=nf_def_var(fileid, 'aFe', nf_real, 1, dd, aFe_varid)
  status=nf_def_var(fileid, 'aN', nf_real, 1, dd, aN_varid)
  
  status=nf_def_var(fileid, 'denb', nf_real, 1, dd, denb_varid)
  status=nf_def_var(fileid, 'benN', nf_real, 1, dd, benN_varid)
  status=nf_def_var(fileid, 'benC', nf_real, 1, dd, benC_varid)
  status=nf_def_var(fileid, 'benSi', nf_real, 1, dd, benSi_varid)
  status=nf_def_var(fileid, 'benCalc', nf_real, 1, dd, benCalc_varid)
  
  status=nf_def_var(fileid, 'NPPn', nf_real, 1, dd, NPPn_varid)
  status=nf_def_var(fileid, 'NPPd', nf_real, 1, dd, NPPd_varid)
  status=nf_def_var(fileid, 'GPPn', nf_real, 1, dd, GPPn_varid)
  status=nf_def_var(fileid, 'GPPd', nf_real, 1, dd, GPPd_varid)
  
  status=nf_def_var(fileid, 'NNAn', nf_real, 1, dd, NNAn_varid)
  status=nf_def_var(fileid, 'NNAd', nf_real, 1, dd, NNAd_varid)
  status=nf_def_var(fileid, 'GNAn', nf_real, 1, dd, GNAn_varid)
  status=nf_def_var(fileid, 'GNAd', nf_real, 1, dd, GNAd_varid)
  
  ! level and depth information: depth at center of cell and depth of cell nodes
  dd(1) = dim_time
  dd(2) = dim_level
  status = nf_def_var(fileid,'depth', nf_real,2, dd, depth_varid)	
  dd(1) = dim_time
  dd(2) = dim_bound	
  status = nf_def_var(fileid,'depth_bounds', nf_real,2, dd, bound_varid)
  
  ! time and depth dependent variables (mainly tracers)
  dd(1) = dim_time
  dd(2) = dim_level
  status=nf_def_var(fileid, 'temp', nf_real, 2, dd, T_varid)
  status=nf_def_var(fileid, 'salt', nf_real, 2, dd, S_varid)  
  status=nf_def_var(fileid, 'PAR', nf_real, 2, dd, PAR_varid)
 
  status=nf_def_var(fileid, 'DIN', nf_real, 2, dd, DIN_varid) 
  status=nf_def_var(fileid, 'DIC', nf_real, 2, dd, DIC_varid)
  status=nf_def_var(fileid, 'Alk', nf_real, 2, dd, Alk_varid) 
  status=nf_def_var(fileid, 'O2', nf_real, 2, dd, O2_varid)
  status=nf_def_var(fileid, 'PhyN', nf_real, 2, dd, PhyN_varid)  
  status=nf_def_var(fileid, 'PhyC', nf_real, 2, dd, PhyC_varid)
  status=nf_def_var(fileid, 'PhyChl', nf_real, 2, dd, PhyChl_varid)
  status=nf_def_var(fileid, 'PhyCalc', nf_real, 2, dd, Phycalc_varid) 
  
  status=nf_def_var(fileid, 'DetN', nf_real, 2, dd, DetN_varid)
  status=nf_def_var(fileid, 'DetC', nf_real, 2, dd, DetC_varid) 
  status=nf_def_var(fileid, 'DetSi', nf_real, 2, dd, DetSi_varid)
  status=nf_def_var(fileid, 'DetCalc', nf_real, 2, dd, DetCalc_varid)
    
  status=nf_def_var(fileid, 'HetN', nf_real, 2, dd, HetN_varid)
  status=nf_def_var(fileid, 'HetC', nf_real, 2, dd, HetC_varid) 
  
  status=nf_def_var(fileid, 'DSi', nf_real, 2, dd, DSi_varid)
  status=nf_def_var(fileid, 'DFe', nf_real, 2, dd, DFe_varid) 
  
  status=nf_def_var(fileid, 'DON', nf_real, 2, dd, DON_varid)
  status=nf_def_var(fileid, 'DOC', nf_real, 2, dd, DOC_varid)
    
  status=nf_def_var(fileid, 'DiaN', nf_real, 2, dd, DiaN_varid)
  status=nf_def_var(fileid, 'DiaC', nf_real, 2, dd, DiaC_varid) 
  status=nf_def_var(fileid, 'DiaChl', nf_real, 2, dd, DiaChl_varid)
  status=nf_def_var(fileid, 'DiaSi', nf_real, 2, dd, DiaSi_varid) 
  
 
  if (REcoM_Second_Zoo) then
    	status=nf_def_var(fileid, 'Zoo2N', nf_real, 2, dd, Zoo2N_varid)
  	status=nf_def_var(fileid, 'Zoo2C', nf_real, 2, dd, Zoo2C_varid) 
  	status=nf_def_var(fileid, 'idetz2n', nf_real, 2, dd, detz2n_varid)
  	status=nf_def_var(fileid, 'idetz2c', nf_real, 2, dd, detz2c_varid) 
  	status=nf_def_var(fileid, 'idetz2si', nf_real, 2, dd, detz2si_varid)
  	status=nf_def_var(fileid, 'idetz2calc', nf_real, 2, dd, detz2calc_varid)   
  endif   

  ! define attributes (units, variable names) long_name, unit
  ! dates, depth and depth bounds
  ! dates

  print*, 'test', fileid, dates_varid, depth_varid, trim(units), d_unit
  status = nf_put_att_text(fileid, dates_varid, trim(units), unitlength, d_unit)
  status = nf_put_att_text(fileid, dates_varid, trim(lname), namelength,		&
  	"days of simulation")
  test=.False.
  if (test) then 
  ! depth
  status = nf_put_att(fileid, depth_varid, trim(units), unitlength, m_unit)
  status = nf_put_att(fileid, depth_varid, trim(lname), namelength, 		&
  	"depth (center of vertical cell)") 
  
  ! depth cell boundaries 
   status = nf_put_att(fileid, bound_varid, trim(units), unitlength, m_unit)
  status = nf_put_att(fileid, bound_varid, trim(lname), namelength, 		&
  	"depth (boundaries of vertical cell)")
  	  
  ! dpCO2
  status = nf_put_att(fileid, dpCO2_varid, trim(units), unitlength, P_unit)
  status = nf_put_att(fileid, dpCO2_varid, trim(lname), namelength, 		&
  	"Difference of oceanic pCO2 minus atmospheric pCO2") 
  ! pCO2
  status = nf_put_att(fileid, pCO2_varid, trim(units), unitlength, P_unit)
  status = nf_put_att(fileid, pCO2_varid, trim(lname), namelength, 		&
  	"Partial pressure of oceanic CO2") 
  ! CO2f
  status = nf_put_att(fileid, CO2f_varid, trim(units), unitlength, f_unit)
  status = nf_put_att(fileid, CO2f_varid, trim(lname), namelength, 		&
  	"CO2-flux into the surface water") 
  ! Hp
  status = nf_put_att(fileid, Hp_varid, trim(units), unitlength, qm_unit)
  status = nf_put_att(fileid, Hp_varid, trim(lname), namelength, 			&
  	"Mean of H-plus ions in the surface water") 
  ! aFe
  status = nf_put_att(fileid, aFe_varid, trim(units), unitlength, aFe_unit)
  status = nf_put_att(fileid, aFe_varid, trim(lname), namelength, "Atmospheric iron input")   
  
  ! aN
  status = nf_put_att(fileid, aN_varid, trim(units), unitlength, aN_unit)
  status = nf_put_att(fileid, aN_varid, trim(lname), namelength, "Atmospheric DIN input") 
  
  ! denb
  status = nf_put_att(fileid, denb_varid, trim(units), unitlength, qa_unit)
  status = nf_put_att(fileid, denb_varid, trim(lname), namelength, "Benthic denitrification rate") 
  
  ! benN
  status = nf_put_att(fileid, benN_varid, trim(units), unitlength, qa_unit)
  status = nf_put_att(fileid, benN_varid, trim(lname), namelength, "Benthos Nitrogen") 
    
  ! benC
  status = nf_put_att(fileid, benC_varid, trim(units), unitlength, qa_unit)
  status = nf_put_att(fileid, benC_varid, trim(lname), namelength, "Benthos Carbon") 
    
  ! benSi
  status = nf_put_att(fileid, benSi_varid, trim(units), unitlength, qa_unit)
  status = nf_put_att(fileid, benSi_varid, trim(lname), namelength, "Benthos silicon") 
    
  ! benCalc
  status = nf_put_att(fileid, benCalc_varid, trim(units), unitlength, qa_unit)
  status = nf_put_att(fileid, benCalc_varid, trim(lname), namelength, "Benthos calcite") 
  
  ! NPPn
  status = nf_put_att(fileid, NPPn_varid, trim(units), unitlength, f_unit)
  status = nf_put_att(fileid, NPPn_varid, trim(lname), namelength, "Mean NPP nanophytoplankton")
  
  ! NPPd
  status = nf_put_att(fileid, NPPd_varid, trim(units), unitlength, f_unit)
  status = nf_put_att(fileid, NPPd_varid, trim(lname), namelength, "Mean NPP diatoms")
    
  ! GPPn
  status = nf_put_att(fileid, GPPn_varid, trim(units), unitlength, f_unit)
  status = nf_put_att(fileid, GPPn_varid, trim(lname), namelength, "Mean GPP nanophytoplankton")
    
  ! GPPd
  status = nf_put_att(fileid, GPPd_varid, trim(units), unitlength, f_unit)
  status = nf_put_att(fileid, GPPd_varid, trim(lname), namelength, "Mean GPP diatoms")
    
  ! NNAn
  status = nf_put_att(fileid, NNAn_varid, trim(units), unitlength, f_unit)
  status = nf_put_att(fileid, NNAn_varid, trim(lname), namelength, "Net N-assimilation nanophytoplankton")
    
  ! NNAd
  status = nf_put_att(fileid, NNAd_varid, trim(units), unitlength, f_unit)
  status = nf_put_att(fileid, NNAd_varid, trim(lname), namelength, "Net N-assimilation diatoms")  

  ! GNAn
  status = nf_put_att(fileid, GNAn_varid, trim(units), unitlength, f_unit)
  status = nf_put_att(fileid, GNAn_varid, trim(lname), namelength, "Gross N-assimilation nanophytoplankton")
    
  ! GNAd
  status = nf_put_att(fileid, GNAd_varid, trim(units), unitlength, f_unit)
  status = nf_put_att(fileid, GNAd_varid, trim(lname), namelength, "Gross N-assimilation diatoms")  
  
  ! temp
  status = nf_put_att(fileid, T_varid, trim(units), unitlength, T_unit)
  status = nf_put_att(fileid, T_varid, trim(lname), namelength, "Temperature profile")   
  
  ! salt
  status = nf_put_att(fileid, S_varid, trim(units), unitlength, S_unit)
  status = nf_put_att(fileid, S_varid, trim(lname), namelength, "Salinity profile") 
    
  ! PAR
  status = nf_put_att(fileid, PAR_varid, trim(units), unitlength, W_unit)
  status = nf_put_att(fileid, PAR_varid, trim(lname), namelength, "PAR")   
  
  ! DIN
  status = nf_put_att(fileid, DIN_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DIN_varid, trim(lname), namelength, "Dissolved Inorganic Nitrogen")  
   
  ! DIC
  status = nf_put_att(fileid, DIC_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DIC_varid, trim(lname), namelength, "Dissolved Inorganic Carbon")  
  
  ! Alk
  status = nf_put_att(fileid, Alk_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, Alk_varid, trim(lname), namelength, "Total Alkalinity")   

  ! PhyN
  status = nf_put_att(fileid, PhyN_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, PhyN_varid, trim(lname), namelength, 			&
  "Intracellular concentration of Nitrogen in small phytoplankton")   
  
  ! PhyC
  status = nf_put_att(fileid, PhyC_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, PhyC_varid, trim(lname), namelength, 			&
  "Intracellular concentration of Carbon in small phytoplankton")   
  
  ! PhyChl
  status = nf_put_att(fileid, PhyChl_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, PhyChl_varid, trim(lname), namelength, "Current intracellular Chlorophyl A concentration")   

  ! PhyCalc
  status = nf_put_att(fileid, PhyCalc_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, PhyCalc_varid, trim(lname), namelength, "Current intracellular Calcite concentration")  
    
  ! DetN
  status = nf_put_att(fileid, DetN_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DetN_varid, trim(lname), namelength, "Concentration of Nitrogen in Detritus")  
  
  ! DetC
  status = nf_put_att(fileid, DetC_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DetC_varid, trim(lname), namelength, "Concentration of Carbon in Detritus")  

  ! DetSi
  status = nf_put_att(fileid, DetSi_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DetSi_varid, trim(lname), namelength, "Concentration of Silicon in Detritus")  
  
  ! DetCalc
  status = nf_put_att(fileid, DetCalc_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DetCalc_varid, trim(lname), namelength, "Concentration of Calcite in Detritus")  
 
  ! HetN
  status = nf_put_att(fileid, HetN_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, HetN_varid, trim(lname), namelength, "Concentration of Nitrogen in heterotrophs")   
  
  ! HetC
  status = nf_put_att(fileid, HetC_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, HetC_varid, trim(lname), namelength, "Concentration of Carbon in heterotrophs")  
   
  ! DON
  status = nf_put_att(fileid, DON_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DON_varid, trim(lname), namelength, "Dissolved organic Nitrogen in water") 
    
  ! DOC
  status = nf_put_att(fileid, DOC_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DOC_varid, trim(lname), namelength, "Dissolved organic Carbon in water") 
  
  ! O2
  status = nf_put_att(fileid, O2_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, O2_varid, trim(lname), namelength, "concentration of dioxygen in water") 

  ! DiaN
  status = nf_put_att(fileid, DiaN_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DiaN_varid, trim(lname), namelength, "concentration of Nitrogen in diatoms")   
  
  ! DiaC
  status = nf_put_att(fileid, DiaC_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DiaC_varid, trim(lname), namelength, "concentration of Carbon in diatoms")   
  
  ! DiaChl
  status = nf_put_att(fileid, DiaChl_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DiaChl_varid, trim(lname), namelength, "concentration of Chlorophyll A in diatoms")   
  
  ! DiaSi
  status = nf_put_att(fileid, DiaSi_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DiaSi_varid, trim(lname), namelength, "concentration of Silicon in diatoms")   
  
  ! DSi
  status = nf_put_att(fileid, DSi_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DSi_varid, trim(lname), namelength, "concentration of Dissolved Silicate")   
    
  ! DFe
  status = nf_put_att(fileid, DFe_varid, trim(units), unitlength, Conc_unit)
  status = nf_put_att(fileid, DFe_varid, trim(lname), namelength, "concentration of Dissolved Iron")   


  
  if (REcoM_Second_Zoo) then
  
    	! Zoo2N
  	status = nf_put_att(fileid, Zoo2N_varid, trim(units), unitlength, Conc_unit)
  	status = nf_put_att(fileid, Zoo2N_varid, trim(lname), namelength, 				&
  	"Intracellular concentration of Nitrogen in second zooplankton")   
    
  	! Zoo2C
  	status = nf_put_att(fileid, Zoo2C_varid, trim(units), unitlength, Conc_unit)
  	status = nf_put_att(fileid, Zoo2C_varid, trim(lname), namelength, 				&
  	"Intracellular concentration of Carbon in second zooplankton")   
  	
  	! idetz2n
  	status = nf_put_att(fileid, detz2n_varid, trim(units), unitlength, Conc_unit)
  	status = nf_put_att(fileid, detz2n_varid, trim(lname), namelength, 				&
  	"Concentration of Nitrogen in detritus from second zooplankton")   
    
  	! idetz2c
  	status = nf_put_att(fileid, detz2c_varid, trim(units), unitlength, Conc_unit)
  	status = nf_put_att(fileid, detz2c_varid, trim(lname), namelength,  				&
  	"Concentration of Carbon in detritus from second zooplankton")   
    
  	! idetz2si
  	status = nf_put_att(fileid, detz2si_varid, trim(units), unitlength, Conc_unit)
  	status = nf_put_att(fileid, detz2si_varid, trim(lname), namelength,  				&
  	"Concentration of Silicon in detritus from second zooplankton")     

  	! idetz2calc
  	status = nf_put_att(fileid, detz2calc_varid, trim(units), unitlength, Conc_unit)
  	status = nf_put_att(fileid, detz2calc_varid, trim(lname), namelength,  				&
  	"Concentration of Calcite in detritus from second zooplankton")   
  endif
  endif
  ! end definition mode
   status = nf_enddef(fileid)
   
  ! fill and store variable in the netcdf file
  status = nf_put_var(fileid, dates_varid, dates_diag)
  status = nf_put_var(fileid, depth_varid, Z)
  status = nf_put_var(fileid, bound_varid, zbar)
  
  status = nf_put_var(fileid, dpCO2_varid, dpCO2s_diag)
  status = nf_put_var(fileid, pCO2_varid, pCO2s_diag)
  status = nf_put_var(fileid, CO2f_varid, CO2f_diag)
  status = nf_put_var(fileid, Hp_varid, Hp_diag)
  status = nf_put_var(fileid, aFe_varid, aFe_diag)
  status = nf_put_var(fileid, aN_varid, aN_diag)    

  status = nf_put_var(fileid, denb_varid, denb_diag)
  status = nf_put_var(fileid, benN_varid, benN_diag)  
  status = nf_put_var(fileid, benC_varid, benC_diag)
  status = nf_put_var(fileid, benSi_varid, benSi_diag)    
  status = nf_put_var(fileid, benCalc_varid, benCalc_diag)
  
  status = nf_put_var(fileid, NPPn_varid, NPPn_diag)
  status = nf_put_var(fileid, NPPd_varid, NPPd_diag)  
  status = nf_put_var(fileid, GPPn_varid, GPPn_diag)
  status = nf_put_var(fileid, GPPd_varid, GPPd_diag) 
  
  status = nf_put_var(fileid, NNAn_varid, NNAn_diag)
  status = nf_put_var(fileid, NNAd_varid, NNAd_diag)  
  status = nf_put_var(fileid, GNAn_varid, GNAn_diag)
  status = nf_put_var(fileid, GNAd_varid, GNAd_diag) 
  
  status = nf_put_var(fileid, T_varid, temp_diag)  
  status = nf_put_var(fileid, S_varid, sali_diag)
  status = nf_put_var(fileid, PAR_varid, PAR_diag) 

  status = nf_put_var(fileid, DIN_varid, DIN_diag)
  status = nf_put_var(fileid, DIC_varid, DIC_diag)  
  status = nf_put_var(fileid, Alk_varid, Alk_diag)
  status = nf_put_var(fileid, O2_varid, O2_diag) 
  
  status = nf_put_var(fileid, PhyN_varid, PhyN_diag)
  status = nf_put_var(fileid, PhyC_varid, PhyC_diag)  
  status = nf_put_var(fileid, PhyChl_varid, PhyChl_diag)
  status = nf_put_var(fileid, Phycalc_varid, PhyCalc_diag) 

  status = nf_put_var(fileid, DetN_varid, DetN_diag)
  status = nf_put_var(fileid, DetC_varid, DetC_diag)  
  status = nf_put_var(fileid, DetSi_varid, DetSi_diag)
  status = nf_put_var(fileid, DetCalc_varid, DetCalc_diag) 

  status = nf_put_var(fileid, HetN_varid, HetN_diag)
  status = nf_put_var(fileid, HetC_varid, HetC_diag) 
   
  status = nf_put_var(fileid, DSi_varid, DSi_diag)
  status = nf_put_var(fileid, DFe_varid, DFe_diag) 
  
  status = nf_put_var(fileid, DON_varid, DON_diag)
  status = nf_put_var(fileid, DOC_varid, DOC_diag) 
  
  status = nf_put_var(fileid, DiaN_varid, DiaN_diag)
  status = nf_put_var(fileid, DiaC_varid, DiaC_diag)  
  status = nf_put_var(fileid, DiaChl_varid, DiaChl_diag)
  status = nf_put_var(fileid, DiaSi_varid, DiaSi_diag)
   
  
  if (REcoM_Second_Zoo) then
    	status = nf_put_var(fileid, Zoo2N_varid, Zoo2N_diag)
  	status = nf_put_var(fileid, Zoo2C_varid, Zoo2C_diag)  
  	status = nf_put_var(fileid, detz2n_varid, idetz2n_diag)
  	status = nf_put_var(fileid, detz2c_varid, idetz2c_diag) 
  	status = nf_put_var(fileid, detz2si_varid, idetz2si_diag)
  	status = nf_put_var(fileid, detz2calc_varid, idetz2calc_diag) 
  endif   

  ! close file
  status=nf_close(fileid)

end subroutine

!----------------------------------------------------------------------------------

end module

