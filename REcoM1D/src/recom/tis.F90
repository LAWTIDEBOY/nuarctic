!> \file tis.f90
!! \BRIEF 
!>    Module with tis subroutine - compute in situ T from S,T,P
MODULE mtis
CONTAINS
!>    Compute in situ temperature from arrays of potential temp, salinity, and pressure.
!!    This subroutine is needed because sw_temp is a function (using scalars not arrays)
SUBROUTINE tis(salt, tempot, press, pressref, tempis)
  !    Purpose:
  !    Compute in situ temperature from arrays of in situ temp, salinity, and pressure.
  !    Needed because sw_temp is a function

  USE msingledouble
  USE msw_temp
  IMPLICIT NONE


! INPUT variables
  !> salinity [psu]
  REAL(kind=rx), INTENT(in) :: salt
  !> potential temperature [C]
  REAL(kind=rx), INTENT(in) :: tempot
  !> pressure [db]
  REAL(kind=rx), INTENT(in) :: press
  !> pressure reference level [db]
  REAL(kind=rx), INTENT(in) :: pressref

! OUTPUT variables:
  !> in situ temperature [C] 
  REAL(kind=rx), INTENT(out) :: tempis

! REAL(kind=r8) :: dsalt, dtempis, dpress, dpressref
! REAL(kind=r8) :: dtempot


! REAL(kind=r8) :: sw_temp
! REAL(kind=rx) :: sw_temp
! EXTERNAL sw_temp

    !dsalt     = DBLE(salt(i))
    !dtempot   = DBLE(tempot(i))
    !dpress    = DBLE(press(i))
    !dpressref = DBLE(pressref)
    !dtempis   = sw_temp(dsalt, dtempot, dpress, dpressref)
    !tempis(i) = REAL(dtempis)

  tempis = sw_temp(salt, tempot, press, pressref)

  RETURN
END SUBROUTINE tis
END MODULE mtis
