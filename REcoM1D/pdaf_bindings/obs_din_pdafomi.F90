!$Id: obs_din_pdafomi.F90 1012 2022-03-31 12:31:09Z lnerger $
!> PDAF-OMI observation module for type A observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! TYPE = A
!!
!! __Observation type A:__
!! The observation type A in this tutorial are 28 observations at specified 
!! model grid points.
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs)
!! and for the observation operator (obs_op).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! In addition, there are two optional routine, which are required if filters 
!! with localization are used:
!! * init_dim_obs_l_TYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!! * localize_covar_TYPE \n
!!           Only required if the localized EnKF is used:
!!           Apply covariance localization in the LEnKF.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_din_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_world, mype_filter    ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l   ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_din        !< Whether to assimilate this data type
  REAL(kind=8)    :: rms_obs_din  !< Observation error standard deviation (for constant errors)

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.
  INTEGER                     :: dim_step_din, dim_depth_din
  INTEGER, ALLOCATABLE        :: din_steps(:), din_depths(:)
  REAL(kind=8), ALLOCATABLE   :: MOSAiC_DIN(:, :)

! ***********************************************************************
! *** The following two data types are used in PDAFomi                ***
! *** They are declared in PDAFomi and only listed here for reference ***
! ***********************************************************************

! Data type to define the full observations by internally shared variables of the module
!   TYPE obs_f
!           Mandatory variables to be set in INIT_DIM_OBS
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!                                           ! (0) Cartesian, (1) Cartesian periodic
!                                           ! (2) simplified geographic, (3) geographic haversine function
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)
!           
!           Optional variables - they can be set in INIT_DIM_OBS
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!
!           Variables with predefined values - they can be changed in INIT_DIM_OBS
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!      INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs. 
!                                           ! or (0) obs. restricted to those relevant for a process domain
!
!           The following variables are set in the routine PDAFomi_gather_obs
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      INTEGER :: dim_obs_g                 ! global number of observations
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!      INTEGER :: off_obs_g                 ! Offset of this observation in overall global obs. vector
!      INTEGER :: obsid                     ! Index of observation over all assimilated observations
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
!                                           ! (only if full obs. are restricted to process domain))
!   END TYPE obs_f

! Data type to define the local observations by internally shared variables of the module
!   TYPE obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!      INTEGER :: locweight                 ! Specify localization function
!      REAL :: lradius                      ! localization radius
!      REAL :: sradius                      ! support radius for localization function
!   END TYPE obs_l
! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
  TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

CONTAINS

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
  SUBROUTINE init_dim_obs_din(step, dim_obs)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs
    USE mod_assimilation, &
         ONLY: filtertype, cradius, off_fields, f_id
         

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER                   :: i, j, index                      ! Counters
    REAL(kind=8)              :: din_file 
    REAL(kind=8)              :: rmsd, sigma2
    INTEGER                   :: dim_obs_p                 ! Number of process-local observations
    REAL(kind=8), ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL(kind=8), ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL(kind=8), ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates 
    INTEGER                   :: cnt, cnt0                 ! Counters
    REAL(kind=8), ALLOCATABLE :: obs_field(:)  ! Observation field read from file
    CHARACTER(len=2)          :: stepstr          ! String for time step

    INTEGER                   :: obs_loc(1)
    LOGICAL                   :: do_assim_now
! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_world==0) &
         WRITE (*,'(8x,a)') 'Assimilate MOSAiC_DIN'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_din) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2


! **********************************
! *** Read  observations ***
! **********************************

  do_assim_now = .FALSE.
  Do i = 1, dim_step_din
    IF(din_steps(i)==step) THEN 
      do_assim_now = .TRUE. 
      index = i
    END IF 
  END DO 

  WRITE(*,*) 'step: ', step, 'DIN assim ', do_assim_now

 
  ! dim_obs_p = 0
  ! dim_obs = 0   

  IF (do_assim_now) THEN
    ALLOCATE(obs_field(dim_depth_din) )
    obs_field(:) = MOSAiC_DIN(index, :)
    IF (mype_world==0) WRITE(*,*) "DIN obs at step ", step, ':', obs_field
! *** Count valid observations that lie at thi step ***
    cnt = 0
    DO i= 1, dim_depth_din
      IF (obs_field(i) > 0.0D+00) cnt = cnt + 1
    END DO
    dim_obs_p = cnt
    dim_obs = cnt

  ELSE 
    dim_obs_p = 0
    dim_obs = 0
  END IF   


  IF(mype_filter==0) &
    WRITE (*,'(8x, a, i6)') 'DIN number of full observations', dim_obs



  IF (dim_obs_p > 0) THEN
    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations on the process sub-domain ***
    ! Allocate process-local observation arrays
    ALLOCATE(obs_p(dim_obs_p))
    ALLOCATE(ivar_obs_p(dim_obs_p))
    ALLOCATE(ocoord_p(1, dim_obs_p))
    ! Allocate process-local index array
    ! This array has a many rows as required for the observation operator
    ! 1 if observations are at grid points; >1 if interpolation is required
    ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))
    cnt = 0
    cnt0 = 0
    DO i= 1, dim_depth_din
      cnt0 = cnt0 + 1
      IF (obs_field(i) > 0.0D+00) THEN
        cnt = cnt + 1
        thisobs%id_obs_p(1, cnt) = cnt0
        din_file = obs_field(i)
        obs_p(cnt) = din_file 
        ! Define observation errors for process-local observations
        ! We assume 30% RMSD of actual concentration (rms_obs_din = 0.3)
        rmsd = din_file * rms_obs_din         
        ! *** Set inverse observation error variances ***
        ivar_obs_p(cnt) = 1.0D+00 / rmsd
        ocoord_p(1, cnt) = REAL(i, 8)
      END IF
    END DO
   ! add total DIN off-set to the state vector 
    thisobs%id_obs_p = thisobs%id_obs_p + off_fields(f_id%DIN)
    
  ELSE 
    ALLOCATE(obs_p(1))
    ALLOCATE(ivar_obs_p(1))
    ALLOCATE(ocoord_p(thisobs%ncoord, 1))
    ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))

  END IF 

! ****************************************
! *** Gather global observation arrays ***
! ****************************************
  CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
       thisobs%ncoord, cradius, dim_obs)



! ********************
! *** Finishing up ***
! ********************


  ! Deallocate all local arrays
  IF(do_assim_now) DEALLOCATE(obs_field)
  DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)





    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_din



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_din(dim_p, dim_obs, state_p, ostate)

    USE PDAFomi, &
         ONLY: PDAFomi_obs_op_gridpoint

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL(kind=8), INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL(kind=8), INTENT(inout) :: ostate(dim_obs)       !< Full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    ! observation operator for observed grid point values
    CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)

  END SUBROUTINE obs_op_din



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
! C 
!   SUBROUTINE init_dim_obs_l_din(domain_p, step, dim_obs, dim_obs_l)

!     ! Include PDAFomi function
!     USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l

!     ! Include localization radius and local coordinates
!     USE mod_assimilation, &   
!          ONLY: coords_l, local_range, locweight, srange

!     IMPLICIT NONE

! ! *** Arguments ***
!     INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
!     INTEGER, INTENT(in)  :: step         !< Current time step
!     INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
!     INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector


! ! **********************************************
! ! *** Initialize local observation dimension ***
! ! **********************************************

!     CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
!          locweight, local_range, srange, dim_obs_l)

!   END SUBROUTINE init_dim_obs_l_din



! !-------------------------------------------------------------------------------
! !> Perform covariance localization for local EnKF on the module-type observation
! !!
! !! The routine is called in the analysis step of the localized
! !! EnKF. It has to apply localization to the two matrices
! !! HP and HPH of the analysis step for the module-type
! !! observation.
! !!
! !! This routine calls the routine PDAFomi_localize_covar
! !! for each observation type. The call allows to specify a
! !! different localization radius and localization functions
! !! for each observation type.
! !!
!   SUBROUTINE localize_covar_din(dim_p, dim_obs, HP_p, HPH, coords_p)

!     ! Include PDAFomi function
!     USE PDAFomi, ONLY: PDAFomi_localize_covar

!     ! Include localization radius and local coordinates
!     USE mod_assimilation, &   
!          ONLY: local_range, locweight, srange

!     IMPLICIT NONE

! ! *** Arguments ***
!     INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
!     INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
!     REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
!     REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
!     REAL, INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


! ! *************************************
! ! *** Apply covariance localization ***
! ! *************************************

!     CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, local_range, srange, &
!          coords_p, HP_p, HPH)

!   END SUBROUTINE localize_covar_din

END MODULE obs_din_pdafomi
