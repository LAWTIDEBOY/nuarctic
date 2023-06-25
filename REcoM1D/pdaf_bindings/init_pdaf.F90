!$Id: init_pdaf.F90 1093 2023-02-17 13:32:58Z lnerger $
!>  Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! to perform the internal initialization of PDAF.
!!
!! This variant is for the online mode of PDAF.
!!
!! This routine is generic. However, it assumes a constant observation
!! error (rms_obs). Further, with parallelization the local state
!! dimension dim_state_p is used.
!!
!! __Revision history:__
!! * 2008-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf()

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAF_init, PDAF_get_state

  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel

  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
       incremental, type_forget, forget, &
       rank_analysis_enkf, locweight, cradius, sradius, &
       filename, type_trans, type_sqrt, delt_obs, &
       type_winf, limit_winf, pf_res_type, pf_noise_type, pf_noise_amp, &
       type_hyb, hyb_gamma, hyb_kappa, n_fields_1d, n_fields_0d, &
       n_params, off_fields, dim_fields, dim_field_1d, f_id, tr_id, &
       step_null, perturb_scale

  USE recom_config, &             ! REcoM parameters
      ONLY: NCuptakeRatio, NCUptakeRatio_d, k_din, k_din_d, &
      Chl2N_max, Chl2N_max_d, deg_Chl, deg_Chl_d, &
      graz_max, graz_max2, grazEff, grazEff2, lossN, lossN_d, &
      lossN_z, lossC_z, lossN_z2, lossC_z2, reminN, reminC

  USE ocean_module, ONLY: tr_arr   ! Array containing all the tracers

  USE recom_clock, ONLY: timeold, timenew

  USE mod_perturbation_pdaf, ONLY:  perturb_lognorm

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: i, j, k      ! counters
  INTEGER :: iseed(4)               ! Seed for random number
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(3) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps     ! Not used in this implementation
  REAL    :: timenow           ! Not used in this implementation

! *** External subroutines ***
  EXTERNAL :: init_ens_pdaf            ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time,
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_ens_pdaf                ! User supplied pre/poststep routine


! ***************************
! ***   Initialize PDAF   ***
! ***************************

!  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - ONLINE MODE'
!  END IF

  perturb_scale = 0.25

  dim_ens = 9

WRITE(*,*) 'dim_ens = ', dim_ens

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
  f_id%TotChl     = 12
  f_id%DON        = 0
  f_id%DOC        = 0
  f_id%DetN       = 0
  f_id%DetC       = 0
  f_id%NPP        = 0



  WRITE(*,*) f_id


! assign tracer ids
  tr_id%Temp          = 1
  tr_id%Salt          = 2
  tr_id%DIN           = 3
  tr_id%DIC           = 4
  tr_id%ALK           = 5
  tr_id%NanoN         = 6
  tr_id%NanoC         = 7
  tr_id%NanoChl       = 8
  tr_id%DetN          = 9
  tr_id%DetC          = 10
  tr_id%HetN          = 11
  tr_id%HetC          = 12
  tr_id%DON           = 13
  tr_id%DOC           = 14
  tr_id%DiaN          = 15
  tr_id%DiaC          = 16
  tr_id%DiaChl        = 17
  tr_id%DiaSi         = 18
  tr_id%DetSi         = 19
  tr_id%DSi           = 20
  tr_id%DFe           = 21
  tr_id%NanoCaCO3      = 22
  tr_id%DetCaCO3      = 23
  tr_id%DO2           = 24
  tr_id%ZooN          = 24
  tr_id%ZooC          = 26
  tr_id%DetZooN       = 27
  tr_id%DetZooC       = 28
  tr_id%DetZooSi      = 29
  tr_id%DetZooCalCO3  = 30

  WRITE(*,*) 'tr_id = ', tr_id




! *** Define state dimension ***


  n_fields_1d = 14
  n_fields_0d = 0
  n_params = 0

  dim_field_1d = 25


  WRITE(*,*) 'dim_field_1d = ', dim_field_1d

  ALLOCATE(dim_fields(n_fields_1d))

  DO i = 1, n_fields_1d
    dim_fields(i) =  dim_field_1d
  END DO

  WRITE(*,*) 'dim_fields = ', dim_fields

  ! DO i = 1, n_fields_0d
  !   dim_fields(n_fields_1d + i) =  1
  ! END DO

  ! DO i = 1, n_params
  !   dim_fields(n_fields_1d + n_fields_0d + i) =  1
  ! END DO

  ALLOCATE(off_fields(n_fields_1d))
  DO i = 1, n_fields_1d
    off_fields(i) =  dim_field_1d * (i - 1)
  END DO

  WRITE(*,*) 'off_fields = ', off_fields

  ! DO i =  1, n_fields_0d
  !   off_fields(n_fields_1d + i) =  dim_field_1d*n_fields_1d + i - 1
  ! END DO

  ! DO i =  1, n_params
  !   off_fields(n_fields_1d + n_fields_0d + i) =  dim_field_1d*n_fields_1d + i - 1
  ! END DO



  dim_state_p = n_fields_1d * dim_field_1d + n_fields_0d + n_params

  WRITE(*,*) 'dim_state_p = ', dim_state_p


!  NCuptakeRatio
  iseed(1)=1
  iseed(2)=3
  iseed(3)=20*task_id + 1
  iseed(4)=41
  NCuptakeRatio = perturb_lognorm(NCuptakeRatio, perturb_scale, iseed)

!  NCUptakeRatio_d
  iseed(1)=2
  iseed(2)=5
  iseed(3)=19*task_id + 2
  iseed(4)=39
  NCuptakeRatio = perturb_lognorm(NCUptakeRatio_d, perturb_scale, iseed)

!  k_din
  iseed(1)=3
  iseed(2)=7
  iseed(3)=18*task_id + 3
  iseed(4)=37
  NCuptakeRatio = perturb_lognorm(k_din, perturb_scale, iseed)

!  k_din_d
  iseed(1)=4
  iseed(2)=11
  iseed(3)=17*task_id + 4
  iseed(4)=35
  NCuptakeRatio = perturb_lognorm(k_din_d, perturb_scale, iseed)

!  Chl2N_max
  iseed(1)=5
  iseed(2)=13
  iseed(3)=16*task_id + 5
  iseed(4)=33
  NCuptakeRatio = perturb_lognorm(Chl2N_max, perturb_scale, iseed)

!  Chl2N_max_d
  iseed(1)=6
  iseed(2)=17
  iseed(3)=15*task_id + 6
  iseed(4)=31
  NCuptakeRatio = perturb_lognorm(Chl2N_max_d, perturb_scale, iseed)

!  deg_Chl
  iseed(1)=7
  iseed(2)=19
  iseed(3)=14*task_id + 7
  iseed(4)=29
  NCuptakeRatio = perturb_lognorm(deg_Chl, perturb_scale, iseed)

!  deg_Chl_d
  iseed(1)=8
  iseed(2)=23
  iseed(3)=13*task_id + 8
  iseed(4)=27
  NCuptakeRatio = perturb_lognorm(deg_Chl_d, perturb_scale, iseed)

!  graz_max
  iseed(1)=9
  iseed(2)=29
  iseed(3)=12*task_id + 9
  iseed(4)=25
  NCuptakeRatio = perturb_lognorm(graz_max, perturb_scale, iseed)

!  graz_max2
  iseed(1)=10
  iseed(2)=1
  iseed(3)=11*task_id + 10
  iseed(4)=23
  NCuptakeRatio = perturb_lognorm(graz_max2, perturb_scale, iseed)

!  grazEff
  iseed(1)=11
  iseed(2)=31
  iseed(3)=10*task_id + 11
  iseed(4)=21
  NCuptakeRatio = perturb_lognorm(grazEff, perturb_scale, iseed)

!  grazEff2
  iseed(1)=12
  iseed(2)=37
  iseed(3)=9*task_id + 12
  iseed(4)=19
  NCuptakeRatio = perturb_lognorm(grazEff2, perturb_scale, iseed)

!  lossN
  iseed(1)=13
  iseed(2)=41
  iseed(3)=8*task_id + 13
  iseed(4)=17
  NCuptakeRatio = perturb_lognorm(lossN, perturb_scale, iseed)

!  lossN_d
  iseed(1)=14
  iseed(2)=43
  iseed(3)=7*task_id + 14
  iseed(4)=15
  NCuptakeRatio = perturb_lognorm(lossN_d, perturb_scale, iseed)

!  lossN_z
  iseed(1)=15
  iseed(2)=47
  iseed(3)=6*task_id + 15
  iseed(4)=13
  NCuptakeRatio = perturb_lognorm(lossN_z, perturb_scale, iseed)

!  lossC_z
  iseed(1)=16
  iseed(2)=53
  iseed(3)=5*task_id + 16
  iseed(4)=11
  NCuptakeRatio = perturb_lognorm(lossC_z, perturb_scale, iseed)

!  lossN_z2
  iseed(1)=17
  iseed(2)=59
  iseed(3)=4*task_id + 17
  iseed(4)=9
  NCuptakeRatio = perturb_lognorm(lossN_z2, perturb_scale, iseed)

!  lossC_z2
  iseed(1)=18
  iseed(2)=61
  iseed(3)=3*task_id + 18
  iseed(4)=7
  NCuptakeRatio = perturb_lognorm(lossC_z2, perturb_scale, iseed)

!  reminN
  iseed(1)=19
  iseed(2)=67
  iseed(3)=2*task_id + 19
  iseed(4)=5
  NCuptakeRatio = perturb_lognorm(reminN, perturb_scale, iseed)

!  reminC
  iseed(1)=20
  iseed(2)=71
  iseed(3)=1*task_id + 20
  iseed(4)=3
  NCuptakeRatio = perturb_lognorm(reminC, perturb_scale, iseed)

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
                    !   (2) use gamma_alpha: hybrid weight from N_eff/N>=hyb_gamma
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

! *** File names
  filename = 'output.dat'


! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

!  call init_pdaf_parse()


! *** Initial Screen output ***
! *** This is optional      ***

! IF (mype_world == 0) call init_pdaf_info()






  WRITE(*,*) "dim_ens = ", dim_ens
  WRITE(*,*) "filtertype = ", filtertype
  WRITE(*,*) "subtype  = ", subtype 








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
