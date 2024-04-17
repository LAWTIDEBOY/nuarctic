
!$Id: read_param_selection.F90 13-01-2024 nmamnun $
!>  Parse selected parameters from the namelist file 
!!
!! __Revision history:__
!! * 2011-15 - Lars Nerger - Initial code extracted from init_pdaf
!! * Later revisions - see repository log
!!
SUBROUTINE read_param_selection()

  USE mod_parallel_pdaf, &
      ONLY: mype_world

  USE mod_assimilation, &
       ONLY: f_id, n_fields, n_params

  IMPLICIT NONE

! *** Local variables ***
  INTEGER ::  index 
  INTEGER ::  is_NCuptakeRatio, is_NCUptakeRatio_d, is_SiCUptakeRatio,        &
              is_k_din, is_k_din_d, is_alfa, is_alfa_d, is_P_cm, is_P_cm_d,   &
              is_Chl2N_max, is_Chl2N_max_d, is_deg_Chl, is_deg_Chl_d,         &
              is_graz_max, is_graz_max2, is_grazEff, is_grazEff2, is_lossN,   &
              is_lossN_d, is_lossN_z, is_lossN_z2, is_lossC_z, is_lossC_z2,   &
              is_reminN, is_reminC, is_VDet, is_VDet_zoo2


! assign 0 to all index  
  f_id%NCuptakeRatio    = 0
  f_id%NCUptakeRatio_d  = 0
  f_id%SiCUptakeRatio   = 0
  f_id%k_din            = 0
  f_id%k_din_d          = 0
  f_id%alfa             = 0
  f_id%alfa_d           = 0
  f_id%P_cm             = 0
  f_id%P_cm_d           = 0
  f_id%Chl2N_max        = 0
  f_id%Chl2N_max_d      = 0
  f_id%deg_Chl          = 0
  f_id%deg_Chl_d        = 0
  f_id%graz_max         = 0
  f_id%graz_max2        = 0
  f_id%grazEff          = 0
  f_id%grazEff2         = 0
  f_id%lossN            = 0
  f_id%lossN_d          = 0
  f_id%lossN_z          = 0
  f_id%lossN_z2         = 0
  f_id%lossC_z          = 0
  f_id%lossC_z2         = 0
  f_id%reminN           = 0
  f_id%reminC           = 0
  f_id%VDet             = 0
  f_id%VDet_zoo2        = 0


! read namelist file
  NAMELIST /param_nml/is_NCuptakeRatio, is_NCUptakeRatio_d, is_SiCUptakeRatio,      &
              is_k_din, is_k_din_d, is_alfa, is_alfa_d, is_P_cm, is_P_cm_d,         &
              is_Chl2N_max, is_Chl2N_max_d, is_deg_Chl, is_deg_Chl_d, is_graz_max,  &
              is_graz_max2, is_grazEff, is_grazEff2, is_lossN, is_lossN_d,          &
              is_lossN_z, is_lossN_z2, is_lossC_z, is_lossC_z2, is_reminN,          &
              is_reminC, is_VDet, is_VDet_zoo2
      
  ! *** Read namelist for PDAF config ***
  IF (mype_world==0) WRITE (*,'(5x,a)') 'Read namelist for parameter selection'
  OPEN (10,file='../config/parameter_selection.nml')
  READ (10,NML=param_nml)
  CLOSE (10)

  n_params = 0
  index = n_fields

  IF (is_NCuptakeRatio /= 0) THEN
    index = index + 1
    f_id%NCuptakeRatio = index
    n_params = n_params + 1
  END IF 

  IF (is_NCUptakeRatio_d /= 0) THEN
    index = index + 1
    f_id%NCUptakeRatio_d = index
    n_params = n_params + 1
  END IF 

  IF (is_SiCUptakeRatio /= 0) THEN
    index = index + 1
    f_id%SiCUptakeRatio = index
    n_params = n_params + 1
  END IF 

  IF (is_k_din /= 0) THEN
    index = index + 1
    f_id%k_din = index
    n_params = n_params + 1
  END IF 

  IF (is_k_din_d /= 0) THEN
    index = index + 1
    f_id%k_din_d = index
    n_params = n_params + 1
  END IF 

  IF (is_alfa  /= 0) THEN
    index = index + 1
    f_id%alfa = index
    n_params = n_params + 1
  END IF 

  IF (is_alfa_d /= 0) THEN
    index = index + 1
    f_id%alfa_d = index
    n_params = n_params + 1
  END IF 

  IF (is_P_cm  /= 0) THEN
    index = index + 1
    f_id%P_cm = index
    n_params = n_params + 1
  END IF 

  IF (is_P_cm_d /= 0) THEN
    index = index + 1
    f_id%P_cm_d = index
    n_params = n_params + 1
  END IF 

  IF (is_Chl2N_max /= 0) THEN
    index = index + 1
    f_id%Chl2N_max = index
    n_params = n_params + 1
  END IF 

  IF (is_Chl2N_max_d /= 0) THEN
    index = index + 1
    f_id%Chl2N_max_d = index
    n_params = n_params + 1
  END IF 

  IF (is_deg_Chl /= 0) THEN
    index = index + 1
    f_id%deg_Chl = index
    n_params = n_params + 1
  END IF 

  IF (is_deg_Chl_d /= 0) THEN
    index = index + 1
    f_id%deg_Chl_d = index
    n_params = n_params + 1
  END IF 

  IF (is_graz_max /= 0) THEN
    index = index + 1
    f_id%graz_max = index
    n_params = n_params + 1
  END IF 

  IF (is_graz_max2 /= 0) THEN
    index = index + 1
    f_id%graz_max2 = index
    n_params = n_params + 1
  END IF 

  IF (is_grazEff /= 0) THEN
    index = index + 1
    f_id%grazEff = index
    n_params = n_params + 1
  END IF 

  IF (is_grazEff2 /= 0) THEN
    index = index + 1
    f_id%grazEff2 = index
    n_params = n_params + 1
  END IF 

  IF (is_lossN /= 0) THEN
    index = index + 1
    f_id%lossN = index
    n_params = n_params + 1
  END IF 

  IF (is_lossN_d /= 0) THEN
    index = index + 1
    f_id%lossN_d = index
    n_params = n_params + 1
  END IF 

  IF (is_lossN_z /= 0) THEN
    index = index + 1
    f_id%lossN_z = index
    n_params = n_params + 1
  END IF 

  IF (is_lossN_z2 /= 0) THEN
    index = index + 1
    f_id%lossN_z2 = index
    n_params = n_params + 1
  END IF 

  IF (is_lossC_z /= 0) THEN
    index = index + 1
    f_id%lossC_z = index
    n_params = n_params + 1
  END IF 

  IF (is_lossC_z2 /= 0) THEN
    index = index + 1
    f_id%lossC_z2 = index
    n_params = n_params + 1
  END IF 

  IF (is_reminN /= 0) THEN
    index = index + 1
    f_id%reminN = index
    n_params = n_params + 1
  END IF 

  IF (is_reminC /= 0) THEN
    index = index + 1
    f_id%reminC = index
    n_params = n_params + 1
  END IF 

  IF (is_VDet /= 0) THEN
    index = index + 1
    f_id%VDet = index
    n_params = n_params + 1
  END IF 

  IF (is_VDet_zoo2 /= 0) THEN
    index = index + 1
    f_id%VDet_zoo2 = index
    n_params = n_params + 1
  END IF 


  IF (mype_world==0) THEN

    WRITE(*,*) "!! * * * * * * * * * * * * * * * * * * !!"
    WRITE(*,*) "f_id%NCuptakeRatio    =", f_id%NCuptakeRatio
    WRITE(*,*) "f_id%NCUptakeRatio_d  =", f_id%NCUptakeRatio_d
    WRITE(*,*) "f_id%SiCUptakeRatio   =", f_id%SiCUptakeRatio
    WRITE(*,*) "f_id%k_din            =", f_id%k_din
    WRITE(*,*) "f_id%k_din_d          =", f_id%k_din_d
    WRITE(*,*) "f_id%alfa             =", f_id%alfa
    WRITE(*,*) "f_id%alfa_d           =", f_id%alfa_d 
    WRITE(*,*) "f_id%P_cm             =", f_id%P_cm
    WRITE(*,*) "f_id%P_cm_d           =", f_id%P_cm_d
    WRITE(*,*) "f_id%Chl2N_max        =", f_id%Chl2N_max
    WRITE(*,*) "f_id%Chl2N_max_d      =", f_id%Chl2N_max_d
    WRITE(*,*) "f_id%deg_Chl          =", f_id%deg_Chl
    WRITE(*,*) "f_id%deg_Chl_d        =", f_id%deg_Chl_d
    WRITE(*,*) "f_id%graz_max         =", f_id%graz_max
    WRITE(*,*) "f_id%graz_max2        =", f_id%graz_max2
    WRITE(*,*) "f_id%grazEff          =", f_id%grazEff
    WRITE(*,*) "f_id%grazEff2         =", f_id%grazEff2
    WRITE(*,*) "f_id%lossN            =", f_id%lossN
    WRITE(*,*) "f_id%lossN_d          =", f_id%lossN_d
    WRITE(*,*) "f_id%lossN_z          =", f_id%lossN_z
    WRITE(*,*) "f_id%lossN_z2         =", f_id%lossN_z2
    WRITE(*,*) "f_id%lossC_z          =", f_id%lossC_z
    WRITE(*,*) "f_id%lossC_z2         =", f_id%lossC_z2
    WRITE(*,*) "f_id%reminN           =", f_id%reminN
    WRITE(*,*) "f_id%reminC           =", f_id%reminC
    WRITE(*,*) "f_id%VDet             =", f_id%VDet
    WRITE(*,*) "f_id%VDet_zoo2        =", f_id%VDet_zoo2

    WRITE(*,*) "* * * * * * * * * * * * * * * * * *"
    WRITE(*,*) "n_params = ", n_params
    WRITE(*,*) "!! * * * * * * * * * * * * * * * * * * !!"

  END IF 

END SUBROUTINE read_param_selection