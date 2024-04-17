!$Id: init_pdaf.F90 1093 2023-02-17 13:32:58Z nmamnun $
!>  Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! The initialization routine PDAF_init is called
!! to perform the internal initialization of PDAF.
!!
!! This subroutine call read_param_selection where the selected parameters
!! are read from a 
!!
!! In addition, the observations are read from the NetCDF files. 
!!
!! __Revision history:__
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf()

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAF_init, PDAF_get_state

  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel

  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, incremental,  &
       type_forget, forget, rank_analysis_enkf, locweight, cradius, sradius,  &
       filename, type_trans, type_sqrt, delt_obs, init_delt_obs, type_winf,   &
       limit_winf, pf_res_type, pf_noise_type, pf_noise_amp, type_hyb,        &
       hyb_gamma, hyb_kappa, n_fields_1d, n_fields_0d, n_params, off_fields,  & 
       dim_fields, dim_field_1d, f_id, tr_id, step_null, perturb_scale,       &
       write_ens, parameter_estimation, n_fields, bgc_layer, step_assim,      &
       obs_dir

  USE mod_utils, &
      ONLY: file_exist, error_handler, get_unit,  e_warn, e_err

  USE obs_chla_pdafomi, &
      ONLY: assim_chla, rms_obs_chla, chla_steps, chla_depths, dim_step_chla, dim_depth_chla, &
      MOSAiC_Chla



  USE obs_din_pdafomi, &
      only: assim_din, rms_obs_din, din_steps, din_depths, dim_step_din, dim_depth_din, &
      MOSAiC_DIN


  USE obs_dsi_pdafomi, &
      only: assim_dsi, rms_obs_dsi, dsi_steps, dsi_depths, dim_step_dsi, dim_depth_dsi, &
      MOSAiC_DSi



  USE ocean_module, ONLY: tr_arr   ! Array containing all the tracers

  USE recom_clock, ONLY: timeold, timenew



  USE netcdf


  IMPLICIT NONE

! *** Local variables ***
  INTEGER         :: i, j, k      ! counters  
  INTEGER         :: filter_param_i(7) ! Integer parameter array for filter
  REAL(kind=8)    :: filter_param_r(3) ! Real parameter array for filter
  INTEGER         :: status_pdaf       ! PDAF status flag
  INTEGER         :: doexit, steps     ! Not used in this implementation
  REAL(kind=8)    :: timenow           ! Not used in this implementation

  INTEGER         :: dbg_id  ! Debugging flag: >0 for debug output; =0 for no debug output
  
  CHARACTER(len=110)          :: chla_obs_file, din_obs_file, dsi_obs_file 
  INTEGER                     :: ncid, dimid_step, dimid_depth, varid_chla, varid_din, varid_dsi
  INTEGER                     :: ierr, nc_status
   
  INTEGER                     :: startv(2), cntv(2)
  CHARACTER(len=256)          :: msgstring      ! String for error handling message 

! *** External subroutines ***
  EXTERNAL :: init_ens_pdaf            ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time, and dimension of next observation
              distribute_state_pdaf, & ! Routine to distribute a state vector to model fields
              prepoststep_ens_pdaf     ! User supplied pre/poststep routine


! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0)  WRITE (*,'(/1x,a)') '--- INITIALIZE PDAF ---'

! set debug flag 
  dbg_id = 0 !0 for no debug output
  CALL PDAF_set_debug_flag(dbg_id) 


  parameter_estimation = .TRUE.
  step_null = 0
  step_assim = step_null

  assim_chla = .TRUE.
  rms_obs_chla = 3.0D-01

  assim_din = .TRUE.
  rms_obs_din = 3.0D-01



  assim_dsi = .TRUE.
  rms_obs_dsi = 3.0D-01



! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************
! *** IO options ***
  screen      = 2  ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 6    ! Type of filter
                    !   (1) SEIK
                    !   (2) EnKF
                    !   (3) LSEIK
                    !   (4) ETKF
                    !   (5) LETKF
                    !   (6) ESTKF
                    !   (7) LESTKF
                    !   (8) localized EnKF
                    !   (9) NETF
                    !  (10) LNETF
                    !  (12) PF
                    !  (100) GENOBS
  dim_ens = 9       ! Size of ensemble for all ensemble filters
  subtype = 0       ! subtype of filter:
                    !   SEIK:
                    !     (0) mean forecast; new formulation
                    !     (1) mean forecast; old formulation
                    !     (2) fixed error space basis
                    !     (3) fixed state covariance matrix
                    !     (4) SEIK with ensemble transformation
                    !   EnKF:
                    !     (0) analysis for large observation dimension
                    !     (1) analysis for small observation dimension
                    !   LSEIK:
                    !     (0) mean forecast;
                    !     (2) fixed error space basis
                    !     (3) fixed state covariance matrix
                    !     (4) LSEIK with ensemble transformation
                    !   ETKF:
                    !     (0) ETKF using T-matrix like SEIK
                    !     (1) ETKF following Hunt et al. (2007)
                    !       There are no fixed basis/covariance cases, as
                    !       these are equivalent to SEIK subtypes 2/3
                    !   LETKF:
                    !     (0) LETKF using T-matrix like SEIK
                    !     (1) LETKF following Hunt et al. (2007)
                    !       There are no fixed basis/covariance cases, as
                    !       these are equivalent to LSEIK subtypes 2/3
                    !   ESTKF:
                    !     (0) Standard form of ESTKF
                    !     (2) fixed ensemble perturbations
                    !     (3) fixed state covariance matrix
                    !   LESTKF:
                    !     (0) Standard form of LESTKF
                    !     (2) fixed ensemble perturbations
                    !     (3) fixed state covariance matrix
                    !   NETF:
                    !     (0) Standard form of NETF
                    !   LNETF:
                    !     (0) Standard form of LNETF
                    !   PF:
                    !     (0) Standard form of PF
  type_trans = 0    ! Type of ensemble transformation
                    !   SEIK/LSEIK and ESTKF/LESTKF:
                    !     (0) use deterministic omega
                    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   ETKF/LETKF:
                    !     (0) use deterministic symmetric transformation
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   NETF/LNETF:
                    !     (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                    !     (1) use identity transformation
  forget  = 1.0     ! Forgetting factor
  type_forget = 0   ! Type of forgetting factor
                    ! SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
                    ! NETF/LNETF/PF
                    !   (0) apply inflation on forecast ensemble
                    !   (2) apply inflation on analysis ensemble
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
  rank_analysis_enkf = 0   ! rank to be considered for inversion of HPH
                    ! in analysis of EnKF; (0) for analysis w/o eigendecomposition
  type_winf = 0     ! NETF/LNETF: Type of weights inflation: (1) use N_eff/N>limit_winf
  limit_winf = 0.0  ! Limit for weights inflation
  type_hyb = 0      ! LKNETF: Type of hybrid weight:
                    !   (0) use fixed hybrid weight hyb_gamma
                    !   (1) use gamma_lin: (1 - N_eff/N_e)*hyb_gamma
                    !   (2) use gamma_alfa: hybrid weight from N_eff/N>=hyb_gamma
                    !   (3) use gamma_ska: 1 - min(s,k)/sqrt(hyb_kappa) with N_eff/N>=hyb_gamma
                    !   (4) use gamma_sklin: 1 - min(s,k)/sqrt(hyb_kappa) >= 1-N_eff/N>=hyb_gamma
  hyb_gamma =  1.0  ! Hybrid filter weight for state (1.0: LETKF, 0.0: LNETF)
  hyb_kappa = 30.0  ! Hybrid norm for using skewness and kurtosis (type_hyb 3 or 4)
  pf_res_type = 1   ! Resampling type for particle filter
                    !   (1) probabilistic resampling
                    !   (2) stochastic universal resampling
                    !   (3) residual resampling
  pf_noise_type = 0 ! Type of pertubing noise in PF: (0) no perturbations
                    ! (1) constant stddev, (2) amplitude of stddev relative of ensemble variance
  pf_noise_amp = 0.0 ! Noise amplitude for particle filter


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs = 24      ! This should be set according to the data availability

! *** Which observation type to assimilate
  ! assim_OBSTYPE = .true.

! *** specifications for observations ***
  ! rms_obs_OBSTYPE = 0.5    ! Observation error standard deviation

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  cradius = 2.0     ! Cut-off radius in grid points for observation domain in local filters
  sradius = cradius ! Support radius for 5th-order polynomial
                    ! or radius for 1/e for exponential weighting

  perturb_scale = 0.5D+00
  bgc_layer     = 47


! *** Parse command line options   ***
! ***              OR              ***
! ***   Read from namelist file    ***

  call init_pdaf_parse()



! ***         set dimessions         ***
! assign tracer ids
  f_id%DIN        = 1
  f_id%DIC        = 2
  f_id%DSi        = 3
  f_id%NanoN      = 4
  f_id%NanoC      = 5
  f_id%NanoChl    = 6
  f_id%DiaN       = 7
  f_id%DiaC       = 8
  f_id%DiaChl     = 9
  f_id%DiaSi      = 10
  f_id%NanoCaCO3  = 11
  f_id%DON        = 0
  f_id%DOC        = 0
  f_id%DetN       = 0
  f_id%DetC       = 0
  f_id%TotChl     = 12
  f_id%NPP        = 0


  n_fields_1d = 12
  n_fields_0d = 0
  n_fields = n_fields_1d + n_fields_0d



  CALL read_param_selection() 


! assign tracer ids
  tr_id%Temp      = 1
  tr_id%Salt      = 2
  tr_id%DIN       = 3
  tr_id%DIC       = 4
  tr_id%ALK       = 5
  tr_id%NanoN     = 6
  tr_id%NanoC     = 7
  tr_id%NanoChl   = 8
  tr_id%DetN      = 9
  tr_id%DetC      = 10
  tr_id%HetN      = 11
  tr_id%HetC      = 12
  tr_id%DON       = 13
  tr_id%DOC       = 14
  tr_id%DiaN      = 15
  tr_id%DiaC      = 16
  tr_id%DiaChl    = 17
  tr_id%DiaSi     = 18
  tr_id%DetSi     = 19
  tr_id%DSi       = 20
  tr_id%DFe       = 21
  tr_id%NanoCaCO3 = 22
  tr_id%DetCaCO3  = 23
  tr_id%DO2       = 24
  tr_id%ZooN      = 24
  tr_id%ZooC      = 26
  tr_id%DetZooN   = 27
  tr_id%DetZooC   = 28
  tr_id%DetZooSi  = 29
  tr_id%DetZooCalCO3 = 30

! ***         Define state dimension         ***  
  dim_field_1d = bgc_layer

  ALLOCATE(dim_fields(n_fields + n_params))
  IF (n_fields_1d /= 0) THEN 
    DO i = 1, n_fields_1d
      dim_fields(i) =  dim_field_1d
    END DO
  END IF 

  IF (n_fields_0d /= 0) THEN 
    DO i = n_fields_1d + 1, n_fields_1d + n_fields_0d
      dim_fields(i) =  1
    END DO
  END IF 

  DO i = n_fields + 1, n_fields + n_params
    dim_fields(i) = 1
  END DO 



  ALLOCATE(off_fields(n_fields + n_params))

  IF (n_fields_1d /= 0) THEN 
    DO i = 1, n_fields_1d
      off_fields(i) =  dim_field_1d * (i - 1)
    END DO
  END IF 

  IF (n_fields_0d /= 0) THEN 
    DO i = 1, n_fields_0d
      off_fields(n_fields_1d + i) =  dim_field_1d*n_fields_1d + i - 1
    END DO
  END IF 

  DO i = 1, n_params
    off_fields(n_fields + i) = dim_field_1d*n_fields_1d + n_fields_0d + i -1
  END DO 



! state dimenssion 
  dim_state_p = n_fields_1d * dim_field_1d + n_fields_0d + n_params
  IF (mype_world==0) WRITE(*,*) 'dim_state_p = ', dim_state_p



! *** Initial Screen output ***
! *** This is optional      ***
  IF (mype_world==0) THEN
    WRITE (*,*) '-- Overview of PDAF configuration ------'
    WRITE (*,*) 'PDAF [../config/pdaf.nml]:'
    WRITE (*,*) 'dim_ens: ',        dim_ens
    WRITE (*,*) 'screen: ',         screen
    WRITE (*,*) 'filtertype: ',     filtertype
    WRITE (*,*) 'subtype: ',        subtype
    WRITE (*,*) 'type_trans: ',     type_trans
    WRITE (*,*) 'type_forget: ',    type_forget
    WRITE (*,*) 'forget: ',         forget
    WRITE (*,*) 'perturb_scale: ',  perturb_scale
    WRITE (*,*) 'init_delt_obs: ',  init_delt_obs
    WRITE (*,*) 'delt_obs: ',       delt_obs
    WRITE (*,*) 'assim_din: ',      assim_din
    WRITE (*,*) 'rms_obs_din: ',    rms_obs_din
    WRITE (*,*) 'assim_dsi: ',      assim_dsi
    WRITE (*,*) 'rms_obs_dsi: ',    rms_obs_dsi
    WRITE (*,*) 'write_ens: ',      write_ens
    WRITE (*,*) '-- End of PDAF configuration -----------'
  ENDIF




  ! *************************************
  ! ***   Read observations here ***
  ! *************************************
  obs_dir = '../../data'
  chla_obs_file  = TRIM(obs_dir)//TRIM('/MOSAiC_Chla_forLaurent_20220905.nc')

  IF ( file_exist(chla_obs_file) ) THEN
    IF(mype_world==0) WRITE(*,*) '--- read observation from file ', chla_obs_file
    ! Open the MOSAiC_Chla NetCDF file
    nc_status = NF90_OPEN(chla_obs_file, nf90_nowrite, ncid)
    IF (nc_status /= NF90_NOERR) THEN
      WRITE(msgstring, * ) "Error in opening ", TRIM(chla_obs_file)
      CALL error_handler( e_err, "init_pdaf", msgstring )      
    ENDIF

  ELSE 

    WRITE(msgstring, * ) TRIM(chla_obs_file), " does not exist"
    CALL error_handler( e_err, "init_pdaf", msgstring )

  ENDIF

    dim_step_chla = 49
  
  ! Allocate the data_array with the appropriate size
  ALLOCATE( chla_steps(dim_step_chla) )

  ! Get dimension ID 
  nc_status = NF90_INQ_DIMID(ncid, "step", dimid_step)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error getting step dimension ID", TRIM(chla_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )    
  ENDIF


  ! Get the dimension values
  nc_status = nf90_get_var(ncid, dimid_step, chla_steps)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error getting dimension values for step", TRIM(chla_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )      
  ENDIF


  dim_depth_chla = 20
  ALLOCATE( chla_depths(dim_depth_chla) )



  nc_status = NF90_INQ_DIMID(ncid, "depth", dimid_depth)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error getting depth dimension ID", TRIM(chla_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )    
  ENDIF

  nc_status = NF90_GET_VAR(ncid, dimid_depth, chla_depths)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error getting dimension values for depth", TRIM(chla_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )      
  ENDIF





  
  IF(mype_world==0) WRITE(*,*) 'chla_steps:', chla_steps
  IF(mype_world==0) WRITE(*,*) 'chla_depths:', chla_depths


  ! read chl-a data 
  ALLOCATE(MOSAiC_Chla( dim_step_chla, dim_depth_chla ))

  ! Get variable ID 
  nc_status = NF90_INQ_VARID(ncid, "Chl_a", varid_chla)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error in getting variable ID for Chl_a", TRIM(chla_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )    
  ENDIF
  
  ! Read the variable
  startv(1) = 1 ! depth 
  startv(2) = 1 ! step 
  cntv(1) = dim_depth_chla
  cntv(2) = dim_step_chla
       
  nc_status = NF90_GET_VAR(ncid, varid_chla, MOSAiC_Chla, start=startv, count=cntv)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = nf90_close(ncid) 
    WRITE(msgstring, * ) "Error in getting variable ", TRIM(chla_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )
  endif

  nc_status = nf90_close(ncid)
  IF (nc_status /= NF90_NOERR) THEN
  nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error in getting variable ", TRIM(chla_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )
  ENDIF


  IF(mype_world==0) WRITE(*,*) 'MOSAiC_Chla', MOSAiC_Chla 


  din_obs_file = TRIM(obs_dir)//TRIM('/PS122_NUTRIENTS.nc')
  IF ( file_exist(din_obs_file) ) THEN
    IF(mype_world==0) WRITE(*,*) '--- read observation from file ', din_obs_file
    ! Open the PS122_NUTRIENTS.nc  file
    nc_status = NF90_OPEN(din_obs_file, nf90_nowrite, ncid)
    IF (nc_status /= NF90_NOERR) THEN
      WRITE(msgstring, * ) "Error in opening ", TRIM(din_obs_file)
      CALL error_handler( e_err, "init_pdaf", msgstring )      
    ENDIF

  ELSE 

    WRITE(msgstring, * ) TRIM(din_obs_file), " does not exist"
    CALL error_handler( e_err, "init_pdaf", msgstring )

  ENDIF

  dim_step_din = 63
  dim_step_dsi = 63



  ! Allocate the data_array with the appropriate size
  ALLOCATE( din_steps(dim_step_din) )
  ALLOCATE( dsi_steps(dim_step_dsi) )


  ! Get dimension ID 
  nc_status = NF90_INQ_DIMID(ncid, "step", dimid_step)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error getting step dimension ID", TRIM(din_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )    
  ENDIF


  ! Get the dimension values
  nc_status = NF90_GET_VAR(ncid, dimid_step, din_steps)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error getting dimension values for step", TRIM(din_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )      
  ENDIF



  dim_depth_din = 47
  ALLOCATE( din_depths(dim_depth_din) )
  dim_depth_dsi = 47
  ALLOCATE( dsi_depths(dim_depth_dsi) )



  nc_status = NF90_INQ_DIMID(ncid, "depth", dimid_depth)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error getting depth dimension ID", TRIM(din_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )    
  ENDIF

  nc_status = nf90_get_var(ncid, dimid_depth, din_depths)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error getting dimension values for depth ", TRIM(din_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )      
  ENDIF


  
  IF(mype_world==0) WRITE(*,*) 'chla_steps:', chla_steps
  IF(mype_world==0) WRITE(*,*) 'chla_depths:', chla_depths
  
  IF(mype_world==0) WRITE(*,*) 'din_steps:', din_steps
  IF(mype_world==0) WRITE(*,*) 'din_depths:', din_depths


  ! read DIN data 
  ALLOCATE(MOSAiC_DIN( dim_step_din, dim_depth_din ))
  ! read DSi data 
  ALLOCATE(MOSAiC_dsi( dim_step_dsi, dim_depth_dsi ))

  ! Get variable ID 
  nc_status = NF90_INQ_VARID(ncid, "DIN", varid_DIN)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error in getting variable ID for DIN ", TRIM(din_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )    
  ENDIF
  
  ! Read the variable
  startv(1) = 1 ! depth 
  startv(2) = 1 ! step 
  cntv(1) = dim_depth_din
  cntv(2) = dim_step_din
       
  nc_status = NF90_GET_VAR(ncid, varid_DIN, MOSAiC_DIN, start=startv, count=cntv)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = nf90_close(ncid) 
    WRITE(msgstring, * ) "Error in getting variable ", TRIM(din_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )
  endif




  ! Get variable ID 
  nc_status = NF90_INQ_VARID(ncid, "DSi", varid_dsi)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error in getting variable ID for DSi ", TRIM(din_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )    
  ENDIF
  
  ! Read the variable
  startv(1) = 1 ! depth 
  startv(2) = 1 ! step 
  cntv(1) = dim_depth_dsi
  cntv(2) = dim_step_dsi
       
  nc_status = NF90_GET_VAR(ncid, varid_dsi, MOSAiC_dsi, start=startv, count=cntv)
  IF (nc_status /= NF90_NOERR) THEN
    nc_status = nf90_close(ncid) 
    WRITE(msgstring, * ) "Error in getting variable ", TRIM(din_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )
  endif






  nc_status = nf90_close(ncid)
  IF (nc_status /= NF90_NOERR) THEN
  nc_status = NF90_CLOSE(ncid) 
    WRITE(msgstring, * ) "Error in getting variable ", TRIM(din_obs_file)
    CALL error_handler( e_err, "init_pdaf", msgstring )
  ENDIF


  IF(mype_world==0) WRITE(*,*) 'MOSAiC_DIN', MOSAiC_DIN 
  IF(mype_world==0) WRITE(*,*) 'MOSAiC_dsi', MOSAiC_dsi 









! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** Here, the full selection of filters is        ***
! *** implemented. In a real implementation, one    ***
! *** reduce this to selected filters.              ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

  whichinit: IF (filtertype == 2) THEN
     ! *** EnKF with Monte Carlo init ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = 0           ! Smoother lag (not implemented here)
     filter_param_r(1) = forget      ! Forgetting factor

     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 6,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSEIF (filtertype == 9) THEN
     ! *** NETF ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = 0           ! Size of lag in smoother
     filter_param_i(4) = 0           ! Not used for NETF (Whether to perform incremental analysis)
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_winf   ! Type of weights inflation
     filter_param_r(1) = forget      ! Forgetting factor
     filter_param_r(2) = limit_winf  ! Limit for weights inflation

     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSEIF (filtertype == 10) THEN
     ! *** LNETF ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = 0           ! Size of lag in smoother
     filter_param_i(4) = 0           ! Not used for NETF (Whether to perform incremental analysis)
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_winf   ! Type of weights inflation
     filter_param_r(1) = forget      ! Forgetting factor
     filter_param_r(2) = limit_winf  ! Limit for weights inflation

     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSEIF (filtertype == 11) THEN
     ! *** Hybrid filter LKNETF                    ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = 0           ! Smoother lag (not implemented here)
     filter_param_i(4) = 0           ! Whether to perform incremental analysis (not implemented for LKNETF)
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_hyb    ! Type of hybrid weight
     filter_param_r(1) = forget      ! Forgetting factor
     filter_param_r(2) = hyb_gamma   ! Hybrid filter weight for state
     filter_param_r(3) = hyb_kappa   ! Normalization factor for hybrid weight

     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7,&
          filter_param_r, 3, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSEIF (filtertype == 12) THEN
     ! *** Particle Filter ***
     filter_param_i(1) = dim_state_p     ! State dimension
     filter_param_i(2) = dim_ens       ! Size of ensemble
     filter_param_r(1) = pf_noise_amp  ! Noise amplitude
     ! Optional parameters
     filter_param_i(3) = pf_res_type   ! Resampling type
     filter_param_i(4) = pf_noise_type ! Perturbation type
     filter_param_i(5) = type_forget   ! Type of forgetting factor
     filter_param_i(6) = type_winf     ! Type of weights inflation
     filter_param_r(2) = forget        ! Forgetting factor
     filter_param_r(3) = limit_winf    ! Limit for weights inflation

     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 6, &
          filter_param_r, 3, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSE
     ! *** All other filters                       ***
     ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = 0           ! Smoother lag (not implemented here)
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
     filter_param_r(1) = forget      ! Forgetting factor

     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  END IF whichinit


! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF


! **********************************
! *** Prepare ensemble forecasts ***
! **********************************

  CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
       distribute_state_pdaf, prepoststep_ens_pdaf, status_pdaf)

END SUBROUTINE init_pdaf
