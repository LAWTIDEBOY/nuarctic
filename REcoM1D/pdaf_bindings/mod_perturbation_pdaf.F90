! $Id: mod_perturbation_pdaf.F 2021-06-15 $
! Purpose:
!   holding procedures to generate random perturbation of inputs
! Initial code Nabir Mamnun (nabir.mamnun@awi.de) June 2021
! BOP
MODULE mod_perturbation_pdaf
  !   USE mod_assimilation, ONLY: rdk, pi
  IMPLICIT NONE
  SAVE
  CONTAINS
    !   procedure for log-normal parameters
    REAL FUNCTION perturb_lognorm (val_in, relvar, rnd_seed)
      !   FUNCTION - generate random perturbation of parameters in a lognormal distribution
      !   the function is written particularly for REcoM2 model parameters in mind.
      IMPLICIT NONE
      !   Data dictionary
      REAL(kind=8), INTENT(IN) :: val_in ! param value
      REAL(kind=8), INTENT(IN) :: relvar ! perturbation scale
      INTEGER, INTENT(IN), OPTIONAL :: rnd_seed(4) ! seed for random number generation
      !   Local variable
      REAL(kind=8)  :: rnd_num  ! output of the random number generation
      REAL(kind=8)  :: logval !
      REAL(kind=8)  ::sigma2 ! for perturbation
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

END MODULE mod_perturbation_pdaf
