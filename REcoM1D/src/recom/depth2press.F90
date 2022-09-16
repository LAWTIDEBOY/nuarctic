!> \file depth2press.f90
!! \BRIEF 
!> Module with depth2press subroutine - converts depth to pressure
!! with Saunders (1981) formula
MODULE mdepth2press
CONTAINS
!>     Compute pressure [db] from depth [m] & latitude [degrees north].
!!     This subroutine is needed because p80 is a function (using scalars not arrays)
SUBROUTINE depth2press(depth, lat, pdbar)

  !     Purpose:
  !     Compute pressure [db] from depth [m] & latitude [degrees north].
  !     Needed because p80 is a function 

  USE msingledouble
  USE mp80
  IMPLICIT NONE

! INPUT variables
  !> depth [m]
  REAL(kind=rx), INTENT(in) :: depth
  !> latitude [degrees]
  REAL(kind=rx), INTENT(in) :: lat

! OUTPUT variables:
  !> pressure [db]
  REAL(kind=rx), INTENT(out) :: pdbar


! REAL(kind=rx) ::  p80
! EXTERNAL p80

  pdbar = p80(depth, lat)

  RETURN
END SUBROUTINE depth2press
END MODULE mdepth2press
