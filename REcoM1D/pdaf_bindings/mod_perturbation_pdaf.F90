! $Id: mod_perturbation_pdaf.F 2021-06-15 $
! Purpose:
!   holding procedures to generate random perturbation of inputs
! Initial code Nabir Mamnun (nabir.mamnun@awi.de) June 2021
! BOP
MODULE mod_perturbation_pdaf
  USE mod_constants, ONLY: rdk, pi

  IMPLICIT NONE
  SAVE
  CONTAINS


! ---------------------------------
    SUBROUTINE set_random_seed(seed)
      INTEGER, INTENT(IN) :: seed
      INTEGER :: iseed(4)     
      ! Convert the scalar seed to a seed array with four elements
      iseed = [seed, seed, seed, seed]
      ! Set the seed for the random number generator
      CALL RANDOM_SEED(PUT=iseed)
    END SUBROUTINE set_random_seed



! ---------------------------------
    !   procedure for log-normal parameters
    REAL FUNCTION perturb_lognorm (val_in, relvar, rnd_seed)
      !   FUNCTION - generate random perturbation of parameters in a lognormal distribution
      !   the function is written particularly for REcoM2 model parameters in mind.
      IMPLICIT NONE
      !   Data dictionary
      REAL(kind=rdk), INTENT(IN) :: val_in ! param value
      REAL(kind=rdk), INTENT(IN) :: relvar ! perturbation scale
      INTEGER, INTENT(IN), OPTIONAL :: rnd_seed(4) ! seed for random number generation
      !   Local variable
      REAL(kind=rdk)  :: rnd_num  ! output of the random number generation
      REAL(kind=rdk)  :: logval !
      REAL(kind=rdk)  ::sigma2 ! for perturbation
      INTEGER       :: f_seed(4)
      !   define seed for random number genaration
      IF (PRESENT(rnd_seed)) THEN
        f_seed = rnd_seed
      ELSE
        f_seed(1)=19
        f_seed(2)=23
        f_seed(3)=143
        f_seed(4)=17
      END IF
      !   generate random number
      !   the seed of the random number generator; the array elements must be between 0 and 4095,
      !   and ISEED(4) must be odd.
      CALL DLARNV(3, f_seed, 1, rnd_num)
      !   in terms of μ and σ
      !   E(X)    = e^(μ + σ^2/2)
      !   Var(X)  = e^(2μ + σ^2) . (e^(σ^2) - 1)
      !   If the mean E and variance V for the lognormal distribution are given,
      !   then the corresponding μ and σ^2 for the normal distribution are given by
      !   σ^2 = log(1+V/E^2)  ! for relative varience E = 1
      !   μ = log(E) - σ^2/2
      sigma2 = LOG(1.0D+00 + relvar*relvar)
      logval = LOG(val_in) - 5.0D-01 * sigma2
      !   Given a zero mean,
      !   unit deviation normally distributed random variable x,
      !   then a log-normally distributed is simply
      !   y = exp(E+sigma*x)
      !   evaluate the perturbed parameters
      perturb_lognorm = EXP(logval + SQRT(sigma2) * rnd_num)
    END FUNCTION perturb_lognorm

! ---------------------------------
  FUNCTION perturb_lognorm_ens(val_in, relvar, N, rnd_seed) result(perturbed_array)
      ! This function generates an ensemble of perturbed values based on a lognormal distribution.

      IMPLICIT NONE

      ! Input parameters
      REAL(kind=rdk), INTENT(IN) :: val_in      ! Original parameter value
      REAL(kind=rdk), INTENT(IN) :: relvar      ! Perturbation scale (relative variance)
      INTEGER, INTENT(IN)       :: N            ! Ensemble size
      INTEGER, INTENT(IN), OPTIONAL :: rnd_seed(4)  ! Seed for random number generation

      ! Output parameter
      REAL(kind=rdk) :: perturbed_array(N)      ! Output array with lognormal-distributed perturbed values

      ! Local variables
      REAL(kind=rdk) :: rnd_num(N)              ! Array to store generated random numbers
      REAL(kind=rdk) :: logval                 ! Log-transformed input value
      REAL(kind=rdk) :: sigma2                 ! Variance used for perturbation
      INTEGER       :: f_seed(4)               ! Final seed for random number generation
      INTEGER       :: i                       ! Loop index

      ! Set the seed for random number generation
      IF (PRESENT(rnd_seed)) THEN
          f_seed = rnd_seed
      ELSE
          f_seed = (/ 19, 23, 143, 17 /)       ! Default seed values
      END IF

      ! Generate random numbers
      CALL DLARNV(3, f_seed, N, rnd_num)

      ! Calculate variance and log-transformed value for perturbation
      sigma2 = LOG(1.0D+00 + relvar*relvar)
      logval = LOG(val_in) - 0.5D0 * sigma2

      ! Generate the perturbed ensemble
      DO i = 1, N
          perturbed_array(i) = EXP(logval + SQRT(sigma2) * rnd_num(i))
      END DO

  END FUNCTION perturb_lognorm_ens


















! ---------------------------------
    REAL FUNCTION perturbbeta (mean, stddev) ! PerturbBeta(mean, stddev) , rnd_seed
      IMPLICIT NONE
      REAL(kind=rdk), INTENT(IN) :: mean, stddev   ! Mean and standard deviation of the beta distribution

      REAL(kind=rdk) :: alpha, beta               ! Parameters of the beta distribution
      REAL(kind=rdk) :: u

      ! Estimate the alpha and beta parameters
      alpha = mean * (mean * (1.0D+00 - mean) - stddev**2) / stddev**2
      beta = (1.0D+00 - mean) * (mean * (1.0D+00 - mean) - stddev**2) / stddev**2
    
      ! Generate a uniform random number between 0 and 1
      CALL RANDOM_NUMBER(u)
    
      ! Transform the uniform random number to beta-distributed random variable
      perturbbeta = beta_transform(u, alpha, beta)
    
  END FUNCTION perturbbeta 

! ---------------------------------
  REAL FUNCTION perturb_beta(mean, stddev, rnd_seed) 
    IMPLICIT NONE
    REAL(kind=rdk), INTENT(IN) :: mean, stddev   ! Mean and standard deviation of the beta distribution
    INTEGER, INTENT(IN), OPTIONAL :: rnd_seed(4) ! seed for random number generation
    ! local variables 
    REAL(kind=rdk) :: alpha, beta               ! Parameters of the beta distribution
    REAL(kind=rdk) :: value                     ! Result of perturbing a value from the beta distribution
    REAL(kind=rdk)  :: rnd_num
    INTEGER :: idist, n 
    INTEGER :: iseed(4)

    ! Set the seed for the random number generator

    IF (PRESENT(rnd_seed)) THEN
        iseed = rnd_seed
    ELSE
      iseed(1)=19
      iseed(2)=23
      iseed(3)=143
      iseed(4)=17
    END IF

    ! Estimate the alpha and beta parameters
    alpha = mean * (mean * (1.0D+00 - mean) - stddev**2) / stddev**2
    beta = (1.0D+00 - mean) * (mean * (1.0D+00 - mean) - stddev**2) / stddev**2
    
    ! Use DLARNV to generate uniform random numbers between 0 and 1
    idist = 1 ! uniform random 
    n = 1     ! number of randon value 
    CALL DLARNV(1, iseed, 1, rnd_num)

    ! Transform the uniform random numbers to beta-distributed random variable
    perturb_beta = beta_transform(rnd_num, alpha, beta)
    
  END FUNCTION perturb_beta


! ---------------------------------
  FUNCTION perturb_beta_ens (mean, stddev, N, rnd_seed) RESULT(beta_perturbed_array)
    IMPLICIT NONE
    REAL(kind=rdk), INTENT(IN) :: mean, stddev   ! Mean and standard deviation of the beta distribution
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN), OPTIONAL :: rnd_seed(4) ! seed for random number generation
    ! local variables 
    REAL(kind=rdk) :: beta_perturbed_array(N)
    REAL(kind=rdk) :: alpha, beta               ! Parameters of the beta distribution
    REAL(kind=rdk) :: value                     ! Result of perturbing a value from the beta distribution
    REAL(kind=rdk)  :: rnd_num(N)
    INTEGER :: idist, i 
    INTEGER :: iseed(4)

    ! Set the seed for the random number generator
    IF (PRESENT(rnd_seed)) THEN
        iseed = rnd_seed
    ELSE
      iseed(1)=19
      iseed(2)=23
      iseed(3)=143
      iseed(4)=17
    END IF

    ! Estimate the alpha and beta parameters
    alpha = mean * (mean * (1.0D+00 - mean) - stddev**2) / stddev**2
    beta = (1.0D+00 - mean) * (mean * (1.0D+00 - mean) - stddev**2) / stddev**2
    
    ! Use DLARNV to generate uniform random numbers between 0 and 1
    idist = 1 ! uniform random 
    CALL DLARNV(idist, iseed, N, rnd_num)

    ! Transform the uniform random numbers to beta-distributed random variable
    DO i = 1, N
      beta_perturbed_array(i) = beta_transform(rnd_num(i), alpha, beta)
    END DO 
    
  END FUNCTION perturb_beta_ens

! ---------------------------------
  REAL FUNCTION beta_transform(invalue, alpha, beta) 
    IMPLICIT NONE
    REAL(kind=rdk), INTENT(IN) :: invalue             ! Uniform random number between 0 and 1
    REAL(kind=rdk), INTENT(IN) :: alpha, beta   ! Parameters of the beta distribution

    ! Inverse transform sampling to get beta-distributed random variable
    beta_transform = EXP((LOG(invalue) / alpha) / (1.0D+00 / beta))

  END FUNCTION beta_transform

! ---------------------------------
 REAL FUNCTION normal_perturbation(mean, stddev) 
  IMPLICIT NONE
    REAL(kind=rdk), INTENT(IN) :: mean, stddev   ! Mean and standard deviation of the normal distribution
    ! Local variables 
    REAL(kind=rdk) :: u1, u2, z
    
    ! Generate two uniform random numbers between 0 and 1
    CALL RANDOM_NUMBER(u1)
    CALL RANDOM_NUMBER(u2)
    
    ! Transform the uniform random numbers to standard normal random variables
    z = SQRT(-2.0D+00 * LOG(u1)) * COS (2.0D+00 * pi * u2)
    
    ! Transform the standard normal random variable to a normal random variable with mean and stddev
    normal_perturbation = mean + stddev * z
    
  END FUNCTION normal_perturbation
  
! ---------------------------------
  REAL FUNCTION norm_perturb(in_value, stddev, rnd_seed)
    IMPLICIT NONE
    !   Data dictionary
      REAL(kind=rdk), INTENT(IN) :: in_value 
      REAL(kind=rdk), INTENT(IN) :: stddev 
      INTEGER, INTENT(IN), OPTIONAL :: rnd_seed(4) ! seed for random number generation
      !   Local variable
      REAL(kind=rdk)  :: rnd_num  ! output of the random number generation
      INTEGER       :: f_seed(4)
      !   define seed for random number genaration
      IF (PRESENT(rnd_seed)) THEN
        f_seed = rnd_seed
      ELSE
        f_seed(1)=19
        f_seed(2)=23
        f_seed(3)=143
        f_seed(4)=17
      END IF

      CALL DLARNV(3, f_seed, 1, rnd_num)

      norm_perturb = in_value * rnd_num

  END FUNCTION norm_perturb

! ---------------------------------  
  FUNCTION norm_perturb_ens(in_value, stddev, N, rnd_seed) RESULT(perturbed_array)

    IMPLICIT NONE
    !   Data dictionary
      REAL(kind=rdk), INTENT(IN) :: in_value 
      REAL(kind=rdk), INTENT(IN) :: stddev 
      INTEGER, INTENT(IN) :: N     ! Ensemble size 
      INTEGER, INTENT(IN), OPTIONAL :: rnd_seed(4) ! seed for random number generation
      ! Output
      REAL(kind=rdk) :: perturbed_array(N)  ! Output array containing lognormal-distributed perturbed values
       !   Local variable
      REAL(kind=rdk)  :: rnd_num(N)  ! output of the random number generation
      INTEGER       :: f_seed(4)
      INTEGER       :: i
      !   define seed for random number genaration
      IF (PRESENT(rnd_seed)) THEN
        f_seed = rnd_seed
      ELSE
        f_seed(1)=19
        f_seed(2)=23
        f_seed(3)=143
        f_seed(4)=17
      END IF

      CALL DLARNV(3, f_seed, N, rnd_num)
      DO i = 1, N      
      perturbed_array(i) = in_value * rnd_num(i) 
      END DO 


  END FUNCTION norm_perturb_ens









END MODULE mod_perturbation_pdaf
