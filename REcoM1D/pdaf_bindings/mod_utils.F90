! $Id$
MODULE mod_utils
  ! USES
  USE mod_constants, ONLY: rdk, pi

  IMPLICIT NONE

  SAVE

  ! Module variables 
  INTEGER :: e_warn, e_err


  CONTAINS


  SUBROUTINE error_handler(level, routine, text )
    !-----------------------------------------------------------------------
    !> write warnings and errors to standard output
    !> If error terminarte the program
    IMPLICIT NONE
      INTEGER, INTENT(in)           :: level
      CHARACTER(len=*), INTENT(in)  :: routine, text
      CHARACTER(len=16)             :: msgtype
      CHARACTER(len=256)            :: wherefrom
      e_warn = 0
      e_err = 1
      wherefrom = TRIM(routine)
      IF (level == e_warn) THEN
        msgtype = 'WARNING FROM:'
        WRITE(*,*) TRIM(TRIM(msgtype)//' '//TRIM(wherefrom)//' '//TRIM(text))
      ELSE IF (level == e_err)  THEN
        msgtype = 'ERROR FROM:'
        WRITE(*,*) TRIM(TRIM(msgtype)//' '//TRIM(wherefrom)//' '//TRIM(text))
        STOP
      END IF
  END SUBROUTINE error_handler


    
  FUNCTION file_exist (file_name)
    ! inquire whether a file exist
    ! if exist return .TRUE.
    ! if not return .FALSE.
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: file_name
    LOGICAL :: file_exist
    INQUIRE (file=file_name, exist=file_exist)
  END FUNCTION file_exist


  FUNCTION get_unit()
    ! get available file unit number
    IMPLICIT NONE
    INTEGER :: get_unit
    INTEGER :: i, iunit
    LOGICAL :: open
    CHARACTER(len=256) :: msgstring
    iunit = -1
    DO i = 11, 99
      INQUIRE (i, OPENED=open) ! if a file is opened open=true
      IF (.not. open) THEN
        get_unit = i
        RETURN
      ENDIF
    ENDDO
    ! IF you get here it is an error
    WRITE(msgstring, *) "Unable to find an available unit number between 11 and 99"
    CALL error_handler(e_err,"get_unit", msgstring)
  END FUNCTION get_unit


  FUNCTION is_file_open(iunit)
  !> Function that returns .true. if this unit number refers to an open file.
  IMPLICIT NONE
    INTEGER, INTENT(in) :: iunit
    LOGICAL :: is_file_open
    INTEGER :: ios
    LOGICAL :: open
    CHARACTER(len=256) :: msgstring
    INQUIRE (UNIT=iunit, OPENED=open, IOSTAT=ios)
    IF ( ios /= 0 ) THEN
      WRITE(msgstring,*)'Unable to determine status of file unit ', iunit
      CALL error_handler(e_warn, 'is_file_open: ', msgstring )
    ENDIF
    is_file_open = open
  END FUNCTION is_file_open



  SUBROUTINE calculate_ens_mean(dim_state, dim_ens, ens_state, mean_state) 
    IMPLICIT NONE
    INTEGER, INTENT(in)         :: dim_state !< Size of state vector 
    INTEGER, INTENT(in)         :: dim_ens !< Size of ensemble
    REAL(kind=rdk), INTENT(in)  :: ens_state(dim_state, dim_ens)
    REAL(kind=rdk), INTENT(out) :: mean_state(dim_state)
    ! local variables 
    INTEGER       :: i, member 
    REAL(kind=rdk)  :: invdim_ens

    invdim_ens    = 1.0 / REAL(dim_ens, rdk)

    mean_state = 0.0
    DO member = 1, dim_ens
      DO i = 1, dim_state
        mean_state(i) = mean_state(i) + ens_state(i, member)
      END DO
    END DO
    mean_state(:) = invdim_ens * mean_state(:)


  END SUBROUTINE calculate_ens_mean







  SUBROUTINE calculate_ens_variance(dim_state, dim_ens, ens_state, variance) 
    IMPLICIT NONE
    INTEGER, INTENT(in)         :: dim_state !< Size of state vector 
    INTEGER, INTENT(in)         :: dim_ens !< Size of ensemble
    REAL(kind=rdk), INTENT(in)  :: ens_state(dim_state, dim_ens)
    REAL(kind=rdk), INTENT(out) :: variance(dim_state)
    ! local variables 
    INTEGER       :: i, member 
    REAL(kind=rdk)  :: mean_state(dim_state)
    REAL(kind=rdk)  :: invdim_ens    
    REAL(kind=rdk)  :: invdim_ensm1

    invdim_ens    = 1.0 / REAL(dim_ens, rdk)
    invdim_ensm1    = 1.0 / REAL(dim_ens-1, rdk)

    
    mean_state = 0.0
    DO member = 1, dim_ens
      DO i = 1, dim_state
        mean_state(i) = mean_state(i) + ens_state(i, member)
      END DO
    END DO
    mean_state(:) = invdim_ens * mean_state(:)


    variance = 0.0
    DO member = 1, dim_ens
      DO i = 1, dim_state
        variance(i) = variance(i) + (mean_state(i) - ens_state(i, member))**2
      END DO
    END DO
    variance(:) = invdim_ensm1 * variance(:)


  END SUBROUTINE calculate_ens_variance





  REAL FUNCTION log_transformed_mean(mean, stddev) 
    ! performs log transformation on actual concentration data, primarily applied in the context 
    ! of chlorophyll-a observation. Given that chlorophyll-a concentrations are distributed 
    ! according to a lognormal distribution, log transformation proves essential to align with 
    ! the assumptions of the Kalman filter—a necessary tool for data assimilation and estimation.
    ! However, the log transformation process is not a simple, straightforward procedure. 
    ! Despite its complexity, employing log transformation allows for a more suitable 
    ! representation of the chlorophyll-a data. It enables better compatibility with the ensemble 
    ! Kalman filter's Gaussian assumptions. 
    IMPLICIT NONE     
    REAL(kind=rdk), INTENT(in) :: mean    ! mean of the original (non-log-transformed) data
    REAL(kind=rdk), INTENT(in) :: stddev ! standard deviation (σ) of the original (non-log-transformed) data
    ! Convert mean and standard deviation to parameters of the log-normal distribution
    REAL(kind=rdk) :: sigma_lognormal_squared ! variance of the log-normal distribution
    REAL(kind=rdk) :: mu_lognormal            ! log-transformed mean of the log-normal distribution
    ! Calculate the variance of the log-normal distribution (sigma_lognormal_squared):
    ! σ^2 = ln(1.0 + (stddev / μ^2))
    sigma_lognormal_squared = log(1.0D+00 + (stddev**2 / mean**2)) 
    ! Calculate the log-transformed mean of the log-normal distribution
    ! μ  = ln(mean) - (sigma_lognormal_squared / 2.0)
    log_transformed_mean = log(mean) - (sigma_lognormal_squared / 2.0D+00)    

  END FUNCTION log_transformed_mean

      
  REAL FUNCTION exponential_transform_mean(mu, sigma) 
  ! Exponential transform mean = exp(μ+σ^2/2)
  ! μ is the mean of the corresponding normal distribution (the mean of the logarithm samples) 
  ! σ is the standard deviation of the corresponding normal distribution.
    IMPLICIT NONE
    REAL(kind=rdk), INTENT(IN) :: mu, sigma   ! Parameters of the lognormal distribution
    ! Compute the expected value of the lognormal distribution
    exponential_transform_mean = EXP(mu + (sigma**2) / 2.0D+00)
  END FUNCTION exponential_transform_mean




  REAL FUNCTION lognorm_expected_value(mean, stddev)
    ! calculate the expected value lognormal distribution 
    ! from sample mean and standar deviation 
    IMPLICIT NONE     
    REAL(kind=rdk), INTENT(in) :: mean    ! mean of the original (non-log-transformed) data
    REAL(kind=rdk), INTENT(in) :: stddev ! standard deviation (σ) of the original (non-log-transformed) data
    ! Convert mean and standard deviation to parameters of the log-normal distribution
    REAL(kind=rdk) :: sigma2 ! variance of the log-normal distribution
    REAL(kind=rdk) :: mu            ! log-transformed mean of the log-normal distribution
    ! Calculate the variance of the log-normal distribution (sigma_lognormal_squared):
    ! σ^2 = ln(1.0 + (stddev / μ^2))
    sigma2 = log(1.0D+00 + (stddev**2 / mean**2)) 
    ! Calculate the log-transformed mean of the log-normal distribution
    ! μ  = ln(mean) - (sigma_lognormal_squared / 2.0)
    mu = log(mean) - (sigma2 / 2.0D+00)   
    ! Compute the expected value of the lognormal distribution
    lognorm_expected_value = EXP(mu + sigma2 / 2.0D+00)
    
  END FUNCTION lognorm_expected_value


      
  REAL FUNCTION beta_expected(samplemean, samplestd )   
    ! calculate the expected value beta distribution 
    ! from sample mean and standar deviation 
    IMPLICIT NONE
    REAL(kind=rdk), INTENT(IN) :: samplemean   ! sample mean 
    REAL(kind=rdk), INTENT(IN) :: samplestd   ! sample standard deviation
    ! local variable 
    REAL(kind=rdk) :: alpha, beta             ! Parameters of the beta distribution
    ! Estimate parameters alpha and beta using method of moments
    alpha = (samplemean * (1.0D+00 - samplemean) - samplestd**2 * samplemean) / samplestd**2
    beta  = alpha * (1.0D+00 - samplemean) / samplemean
    ! Calculate the expected value of the beta distribution
    beta_expected = alpha / (alpha + beta)
  END FUNCTION beta_expected










END MODULE mod_utils
