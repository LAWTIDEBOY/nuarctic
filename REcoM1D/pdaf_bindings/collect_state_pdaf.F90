!$Id: collect_state_pdaf.F90 nmamnun $
!>  Initialize state vector from model fields
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters.
!!
!! This subroutine is called during the forecast
!! phase from PDAF_put_state_X or PDAF_assimilate_X
!! after the propagation of each ensemble member.
!! The supplied state vector has to be initialized
!! from the model fields (typically via a module).
!! With parallelization, MPI communication might be
!! required to initialize state vectors for all
!! subdomains on the model PEs.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * Later revisions - see repository log
!!
SUBROUTINE collect_state_pdaf(dim_p, state_p)

  USE ocean_module, &
      ONLY: tr_arr

  USE recom_config, &             ! REcoM parameters
      ONLY: NCuptakeRatio, NCUptakeRatio_d, k_din, k_din_d, alfa, alfa_d, P_cm,     & 
      P_cm_d, Chl2N_max, Chl2N_max_d, deg_Chl, deg_Chl_d, graz_max, graz_max2,      &
      grazEff, grazEff2, lossN, lossN_d, lossN_z, lossC_z, lossN_z2, lossC_z2,      &
      reminN, reminC,  VDet, VDet_zoo2

  USE mod_assimilation, &
      ONLY: off_fields, dim_fields, dim_field_1d, f_id, tr_id

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector

! *** local variables ***
! *************************************************


! *** Initialize state vector from model fields ***
  IF (f_id%DIN /= 0)  &
    state_p(off_fields(f_id%DIN)+1:off_fields(f_id%DIN)+dim_field_1d)             & 
    = tr_arr(1:dim_field_1d, tr_id%DIN)

  IF (f_id%DIC /= 0)  &
    state_p(off_fields(f_id%DIC)+1:off_fields(f_id%DIC)+dim_field_1d)             &
    = tr_arr(1:dim_field_1d, tr_id%DIC)

  IF (f_id%DSi /= 0)  &
    state_p(off_fields(f_id%DSi)+1:off_fields(f_id%DSi)+dim_field_1d)             &
    = tr_arr(1:dim_field_1d, tr_id%DSi)

  IF (f_id%NanoN /= 0)  &
    state_p(off_fields(f_id%NanoN)+1:off_fields(f_id%NanoN)+dim_field_1d)         &
    = LOG(tr_arr(1:dim_field_1d, tr_id%NanoN))

  IF (f_id%NanoC /= 0)  &
    state_p(off_fields(f_id%NanoC)+1:off_fields(f_id%NanoC)+dim_field_1d)         &
    = LOG(tr_arr(1:dim_field_1d, tr_id%NanoC))

  IF (f_id%NanoChl /= 0)  &
    state_p(off_fields(f_id%NanoChl)+1:off_fields(f_id%NanoChl)+dim_field_1d)     &
    = LOG(tr_arr(1:dim_field_1d, tr_id%NanoChl))

  IF (f_id%DiaN /= 0)  &
    state_p(off_fields(f_id%DiaN)+1:off_fields(f_id%DiaN)+dim_field_1d)           &
    = LOG(tr_arr(1:dim_field_1d, tr_id%DiaN))

  IF (f_id%DiaC /= 0)  &
    state_p(off_fields(f_id%DiaC)+1:off_fields(f_id%DiaC)+dim_field_1d)           &
    = LOG(tr_arr(1:dim_field_1d, tr_id%DiaC))

  IF (f_id%DiaChl /= 0)  &
    state_p(off_fields(f_id%DiaChl)+1:off_fields(f_id%DiaChl)+dim_field_1d)       &
    = LOG(tr_arr(1:dim_field_1d, tr_id%DiaChl))

  IF (f_id%DiaSi /= 0)  &
    state_p(off_fields(f_id%DiaSi)+1:off_fields(f_id%DiaSi)+dim_field_1d)         &
    = LOG(tr_arr(1:dim_field_1d, tr_id%DiaSi))

  IF (f_id%NanoCaCO3 /= 0)  &
    state_p(off_fields(f_id%NanoCaCO3)+1:off_fields(f_id%NanoCaCO3)+dim_field_1d) &
    = LOG(tr_arr(1:dim_field_1d, tr_id%NanoCaCO3))

  IF (f_id%DON /= 0)  &
    state_p(off_fields(f_id%DON)+1:off_fields(f_id%DON)+dim_field_1d)             &
    = LOG(tr_arr(1:dim_field_1d, tr_id%DON))

  IF (f_id%DOC /= 0)  &
    state_p(off_fields(f_id%DOC)+1:off_fields(f_id%DOC)+dim_field_1d)             &
    = LOG(tr_arr(1:dim_field_1d, tr_id%DOC))

  IF (f_id%DetN /= 0)  &
    state_p(off_fields(f_id%DetN)+1:off_fields(f_id%DetN)+dim_field_1d)  &
    = LOG(tr_arr(1:dim_field_1d, tr_id%DetN))

  IF (f_id%DetC /= 0)  &
    state_p(off_fields(f_id%DetC)+1:off_fields(f_id%DetC)+dim_field_1d)           &
    = LOG(tr_arr(1:dim_field_1d, tr_id%DetC))

  IF (f_id%TotChl /= 0)  &
    state_p(off_fields(f_id%TotChl)+1:off_fields(f_id%TotChl)+dim_field_1d)       &
    = LOG(  tr_arr(1:dim_field_1d, tr_id%NanoChl)                                 &
          + tr_arr(1:dim_field_1d, tr_id%DiaChl))

  ! IF (f_id%NPP /= 0) state_p(off_fields(f_id%NPP)+1) = NPP? 
  ! NPP is not a state variable, perhaps we can compute the vertically average NPP
  ! daily mean and access the value 

  IF (f_id%NCuptakeRatio /= 0) state_p(off_fields(f_id%NCuptakeRatio)+1) = NCuptakeRatio
  
  IF (f_id%NCUptakeRatio_d /= 0) state_p(off_fields(f_id%NCUptakeRatio_d)+1) = NCUptakeRatio_d

  IF (f_id%k_din /= 0) state_p(off_fields(f_id%k_din)+1) = k_din

  IF (f_id%k_din_d /= 0) state_p(off_fields(f_id%k_din_d)+1) = k_din_d

  IF (f_id%alfa /= 0) state_p(off_fields(f_id%alfa)+1) = alfa

  IF (f_id%alfa_d /= 0) state_p(off_fields(f_id%alfa_d)+1) = alfa_d

  IF (f_id%P_cm /= 0) state_p(off_fields(f_id%P_cm)+1) = P_cm

  IF (f_id%P_cm_d /= 0) state_p(off_fields(f_id%P_cm_d)+1) = P_cm_d

  IF (f_id%Chl2N_max /= 0) state_p(off_fields(f_id%Chl2N_max)+1) = Chl2N_max

  IF (f_id%Chl2N_max_d /= 0) state_p(off_fields(f_id%Chl2N_max_d)+1) = Chl2N_max_d

  IF (f_id%deg_Chl /= 0) state_p(off_fields(f_id%deg_Chl)+1) = deg_Chl

  IF (f_id%deg_Chl_d /= 0) state_p(off_fields(f_id%deg_Chl_d)+1) = deg_Chl_d

  IF (f_id%graz_max /= 0) state_p(off_fields(f_id%graz_max)+1) = graz_max

  IF (f_id%graz_max2 /= 0) state_p(off_fields(f_id%graz_max2)+1) = graz_max2

  IF (f_id%grazEff /= 0) state_p(off_fields(f_id%grazEff)+1) = grazEff

  IF (f_id%grazEff2 /= 0) state_p(off_fields(f_id%grazEff2)+1) = grazEff2

  IF (f_id%lossN /= 0) state_p(off_fields(f_id%lossN)+1) = lossN

  IF (f_id%lossN_d /= 0) state_p(off_fields(f_id%lossN_d)+1) = lossN_d

  IF (f_id%lossN_z /= 0) state_p(off_fields(f_id%lossN_z)+1) = lossN_z

  IF (f_id%lossN_z2 /= 0) state_p(off_fields(f_id%lossN_z2)+1) = lossN_z2

  IF (f_id%lossC_z /= 0) state_p(off_fields(f_id%lossC_z)+1) = lossC_z

  IF (f_id%lossC_z2 /= 0) state_p(off_fields(f_id%lossC_z2)+1) = lossC_z2

  IF (f_id%reminN /= 0) state_p(off_fields(f_id%reminN)+1) = reminN

  IF (f_id%reminC /= 0) state_p(off_fields(f_id%reminC)+1) = reminC

  IF (f_id%VDet /= 0) state_p(off_fields(f_id%VDet)+1) = VDet

  IF (f_id%VDet_zoo2 /= 0) state_p(off_fields(f_id%VDet_zoo2)+1) = VDet_zoo2



END SUBROUTINE collect_state_pdaf
