module recom_extra

contains
! ==============================================================
subroutine integrate_nod(datas, int2D)
  use ocean_module
  IMPLICIT NONE
  real(kind=8), intent(in)       :: datas(:)
  real(kind=8), intent(inout)    :: int2D

  int2D=0.d0
  int2D=dot_product(datas,hnode) ! integrate on the water column


end subroutine integrate_nod


!===============================================================================
! Subroutine for calculating flux-depth and thickness of control volumes
!===============================================================================

subroutine Depth_calculations(Nn,wF,zF,thick,recipthick)
  use recom_config
  use ocean_module
 
  implicit none
  Integer,                  intent(in)     :: Nn	     ! Total number of nodes
! Output
  real(kind=8),dimension(:),intent(out)    :: zF             ! [m] Depth of vertical fluxes

  real(kind=8),dimension(:),intent(out)    :: thick          ! [m] Distance between two nodes = thickness
  real(kind=8),dimension(:),intent(out)    :: recipthick     ! [1/m] reciprocal of thickness

  real(kind=8),dimension(:,:),intent(out)    :: wF             ! [m/day] Velocities of fluxes at the border of the control volumes
  Integer                                  :: k             ! Index for depth      

! ======================================================================================
!! zbar (depth of layers) and Z (mid depths of layers)
!! zbar is negative 
!! zbar(nl) allocate the array for storing the standard depths
!! Z(nl-1)  mid-depths of cells

    ! zbar_n: depth of layers due to ale thinkness variations at every node n 
!    allocate(zbar_n(nl))
!    allocate(zbar_3d_n(nl,myDim_nod2D+eDim_nod2D))
    
    ! Z_n: mid depth of layers due to ale thinkness variations at every node n 
!    allocate(Z_n(nl-1))
!    allocate(Z_3d_n(nl-1,myDim_nod2D+eDim_nod2D)) 
! ============================================================================== modular

  wF(2:Nn,ivphy) = VPhy  
  wF(2:Nn,ivdia) = VDia
  wF(2:Nn,ivdet) = VDet
  wF(2:Nn,ivdetsc) = VDet_zoo2

  wF(1,:)         = 0.d0
  wF(Nn+1,:)      = 0.d0

!----------------------------------------------------
! calculate thickness of vertical layers

    thick   =0.0
    recipthick=0.0
    zF=0.0

   do k=1, Nn !nlevels_nod2D(n)-1
      thick(k)=hnode(k)
      if (hnode(k) > 0.) then
         recipthick(k)=1./hnode(k)
      else
         recipthick(k)=0.
      end if
   end do
!end do

! layer depth (negative)
   do k=1, Nn+1 !nlevels_nod2D(n)
      zF(k)=zbar(k)
   end do

  
end subroutine Depth_calculations

!===============================================================================
! Subroutine for calculating cos(AngleOfIncidence)
!===============================================================================
subroutine Cobeta
  use general_config
  use REcoM_GloVar
  use recom_clock
  use ocean_module
  
  Implicit none
  	
! Declarations
  Real(kind=8)                     :: yearfrac              ! The fraction of the year that has passed [0 1]
  Real(kind=8)                     :: yDay                  ! yearfrac in radians [0 2*pi]
  Real(kind=8)                     :: declination   = 0.d0  ! Declination of the sun at present lat and time
  Real(kind=8)                     :: CosAngleNoon  = 0.d0  ! Cos(Angle of Incidence) at Noon ?

! Constants
  Real(kind=8)		           :: nWater        = 1.33 
!
! find day (****NOTE for year starting in winter*****)  
! Paltridge, G. W. and C. M. R. Platt, Radiative Processes in Meteorology and Climatology, Developments in Atmospheric Sciences, vol. 5, Elsevier Scientific Publishing Company, Amsterdam, Oxford, New York, 1976, ISBN 0-444-41444-4.
  yearfrac    = mod(real(daynew),real(ndpyr))/real(ndpyr)
  yDay        = 2 * pi * yearfrac
  declination = 0.006918                   &
         	- 0.399912 * cos(    yDay) &
     	        + 0.070257 * sin(    yDay) &
          	- 0.006758 * cos(2 * yDay) &
        	+ 0.000907 * sin(2 * yDay) &
        	- 0.002697 * cos(3 * yDay) &
        	+ 0.001480 * sin(3 * yDay) 
   
   cosAngleNoon   =   sin(geo_coords(2)) * sin(declination) &
     		         + cos(geo_coords(2)) * cos(declination)
   cosAI       = sqrt(1-(1-cosAngleNoon**2)/nWater**2)

end subroutine Cobeta

!================================================================================
! Calculating second zooplankton respiration rates
!================================================================================
 subroutine krill_resp
   use REcoM_LocVar
   use recom_clock
   use ocean_module
   !use recom_forcing

   implicit none
                                                                                            
   if (geo_coords(2)<0.0) then  !SH
      if ((daynew .le. 105)) then
       res_zoo2_a = 0.d0
      else if((105 .le. daynew).and.(daynew .le. 150)) then
       res_zoo2_a = (-1./90.)*daynew +7./6.
      else if((150 .lt. daynew).and.(daynew .lt. 250)) then
       res_zoo2_a = -0.5
      else if((250 .le. daynew).and.(daynew .le. 295)) then
       res_zoo2_a = (1/90.)*daynew - 59./18.
      else if((daynew .gt. 295)) then
       res_zoo2_a = 0.d0
      end if
   else
      if ((daynew .le. 65)) then
       res_zoo2_a = -0.5
      else if((285 .le. daynew).and.(daynew .le. 330)) then
       res_zoo2_a = (-1./90.)*daynew +57./18.
      else if((330 .lt. daynew)) then
       res_zoo2_a = -0.5
      else if((65 .le. daynew).and.(daynew .le. 110)) then
       res_zoo2_a = (1/90.)*daynew - 22./18.
      else if((110 .lt. daynew).and.(daynew .lt. 285)) then
       res_zoo2_a = 0.d0
      end if
   endif



!  if ((Latr .lt. 0).and.(daynew .le. 105)) then
!       res_zoo2_a = 0.d0
!   else if((Latr .lt. 0).and.(105 .le. daynew).and.(daynew .le. 150)) then
!       res_zoo2_a = (-1./90.)*daynew +7./6.
!   else if((Latr .lt. 0).and.(150 .lt. daynew).and.(daynew .lt. 250)) then
!       res_zoo2_a = -0.5
!   else if((Latr .lt. 0).and.(250 .le. daynew).and.(daynew .le. 295)) then
!       res_zoo2_a = (1/90.)*daynew - 59./18.
!   else if((Latr .lt. 0).and.(daynew .gt. 295)) then
!       res_zoo2_a = 0.d0
!  end if

!!For Northern Hemisphere


!  if ((Latr .ge. 0).and.(daynew .le. 65)) then
!       res_zoo2_a = -0.5
!   else if((Latr .ge. 0).and.(285 .le. daynew).and.(daynew .le. 330)) then
!       res_zoo2_a = (-1./90.)*daynew +57./18.
!   else if((Latr .ge. 0).and.(330 .lt. daynew)) then
!       res_zoo2_a = -0.5
!   else if((Latr .ge. 0).and.(65 .le. daynew).and.(daynew .le. 110)) then
!       res_zoo2_a = (1/90.)*daynew - 22./18.
!   else if((Latr .ge. 0).and.(110 .lt. daynew).and.(daynew .lt. 285)) then
!       res_zoo2_a = 0.d0
!  end if

 
 end subroutine krill_resp

end module
