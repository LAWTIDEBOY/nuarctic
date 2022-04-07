!> \file tpot.f90
!! \BRIEF 
!>    Module with tpot subroutine - compute potential T from in situ T,S,P
MODULE mtpot
CONTAINS
!>    Compute potential temperature from arrays of in situ temp, salinity, and pressure.
!!    This subroutine is needed because sw_ptmp is a function (using scalars not arrays)
SUBROUTINE tpot(salt, tempis, press, pressref, tempot)
  !    Purpose:
  !    Compute potential temperature from arrays of in situ temp, salinity, and pressure.
  !    Needed because sw_ptmp is a function

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

  USE msingledouble
  USE msw_ptmp
  IMPLICIT NONE


! INPUT variables
  !> salinity [psu]
  REAL(kind=rx), INTENT(in) :: salt
  !> in situ temperature [C]
  REAL(kind=rx), INTENT(in):: tempis
  !> pressure [db]
  REAL(kind=rx), INTENT(in) :: press
  !> pressure reference level [db]
  REAL(kind=rx), INTENT(in) :: pressref

! OUTPUT variables:
  !> potential temperature [C] for pressref
  REAL(kind=rx), INTENT(out) :: tempot

  REAL(kind=r8) :: dsalt, dtempis, dpress, dpressref
  REAL(kind=r8) :: dtempot

  dsalt     = DBLE(salt)
  dtempis   = DBLE(tempis)
  dpress    = DBLE(press)
  dpressref = DBLE(pressref)

  dtempot   = sw_ptmp(dsalt, dtempis, dpress, dpressref)

  tempot = SGLE(dtempot)

  RETURN
END SUBROUTINE tpot
END MODULE mtpot
