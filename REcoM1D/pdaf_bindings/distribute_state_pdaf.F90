!$Id: distribute_state_pdaf.F90 871 2021-11-22 16:44:34Z lnerger $
!>  Initialize model fields from state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters.
!!
!! During the forecast phase of the filter this
!! subroutine is called from PDAF_get_state
!! supplying a model state which has to be evolved.
!! The routine has to initialize the fields of the
!! model (typically available through a module) from
!! the state vector of PDAF. With parallelization,
!! MPI communication might be required to
!! initialize all subdomains on the model PEs.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

  USE ocean_module, &
      ONLY: tr_arr

  USE recom_config, &             ! REcoM parameters
      ONLY: NCuptakeRatio, NCUptakeRatio_d, k_din, k_din_d, alfa, alfa_d, P_cm,     & 
      P_cm_d, Chl2N_max, Chl2N_max_d, deg_Chl, deg_Chl_d, graz_max, graz_max2,      &
      grazEff, grazEff2, lossN, lossN_d, lossN_z, lossC_z, lossN_z2, lossC_z2,      &
      reminN, reminC,  VDet, VDet_zoo2

  USE mod_assimilation, &
      ONLY: parameter_estimation, off_fields, dim_fields, dim_field_1d, &
      f_id, tr_id, step_assim



  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** local variables ***


! *************************************************
! *** Initialize model fields from state vector ***
! *** for process-local model domain            ***
!**************************************************

  IF (f_id%DIN /= 0) tr_arr(1:dim_field_1d, tr_id%DIN)   &
    = state_p(off_fields(f_id%DIN)+1:off_fields(f_id%DIN)+dim_field_1d)

  IF (f_id%DIC /= 0) tr_arr(1:dim_field_1d, tr_id%DIC)   &
    = state_p(off_fields(f_id%DIC)+1:off_fields(f_id%DIC)+dim_field_1d)

  IF (f_id%DSi /= 0) tr_arr(1:dim_field_1d, tr_id%DSi)   &
    = state_p(off_fields(f_id%DSi)+1:off_fields(f_id%DSi)+dim_field_1d)

  IF (f_id%NanoN /= 0) tr_arr(1:dim_field_1d, tr_id%NanoN) &
    = EXP(state_p(off_fields(f_id%NanoN)+1:off_fields(f_id%NanoN)+dim_field_1d))

  IF (f_id%NanoC /= 0) tr_arr(1:dim_field_1d, tr_id%NanoC) &
    = EXP(state_p(off_fields(f_id%NanoC)+1:off_fields(f_id%NanoC)+dim_field_1d))

  IF (f_id%NanoChl /= 0) tr_arr(1:dim_field_1d, tr_id%NanoChl) &
    = EXP(state_p(off_fields(f_id%NanoChl)+1:off_fields(f_id%NanoChl)+dim_field_1d))

  IF (f_id%DiaN /= 0) tr_arr(1:dim_field_1d, tr_id%DiaN) &
    = EXP(state_p(off_fields(f_id%DiaN)+1:off_fields(f_id%DiaN)+dim_field_1d))

  IF (f_id%DiaC /= 0) tr_arr(1:dim_field_1d, tr_id%DiaC) &
    = EXP(state_p(off_fields(f_id%DiaC)+1:off_fields(f_id%DiaC)+dim_field_1d))

  IF (f_id%DiaChl /= 0) tr_arr(1:dim_field_1d, tr_id%DiaChl) &
    = EXP(state_p(off_fields(f_id%DiaChl)+1:off_fields(f_id%DiaChl)+dim_field_1d))

  IF (f_id%DiaSi /= 0) tr_arr(1:dim_field_1d, tr_id%DiaSi) &
    = EXP(state_p(off_fields(f_id%DiaSi)+1:off_fields(f_id%DiaSi)+dim_field_1d))

  IF (f_id%NanoCaCO3 /= 0) tr_arr(1:dim_field_1d, tr_id%NanoCaCO3) &
    = EXP(state_p(off_fields(f_id%NanoCaCO3)+1:off_fields(f_id%NanoCaCO3)+dim_field_1d))

  IF (f_id%DON /= 0) tr_arr(1:dim_field_1d, tr_id%DON)   &
    = state_p(off_fields(f_id%DON)+1:off_fields(f_id%DON)+dim_field_1d)

  IF (f_id%DOC /= 0) tr_arr(1:dim_field_1d, tr_id%DOC)   &
    = state_p(off_fields(f_id%DOC)+1:off_fields(f_id%DOC)+dim_field_1d)
    
  IF (f_id%DetN /= 0) tr_arr(1:dim_field_1d, tr_id%DetN)   &
    = state_p(off_fields(f_id%DetN)+1:off_fields(f_id%DetN)+dim_field_1d)
    
  IF (f_id%DetC /= 0) tr_arr(1:dim_field_1d, tr_id%DetC)   &
    = state_p(off_fields(f_id%DetC)+1:off_fields(f_id%DetC)+dim_field_1d)


  ! IF ( step_assim == 0 .OR. parameter_estimation) THEN 
    IF (f_id%NCuptakeRatio /= 0)    NCuptakeRatio = state_p(off_fields(f_id%NCuptakeRatio)+1)
    IF (f_id%NCUptakeRatio_d /= 0)  NCUptakeRatio_d = state_p(off_fields(f_id%NCUptakeRatio_d)+1)
    IF (f_id%k_din /= 0)            k_din = state_p(off_fields(f_id%k_din)+1)
    IF (f_id%k_din_d /= 0)          k_din_d = state_p(off_fields(f_id%k_din_d)+1)
    IF (f_id%alfa /= 0)             alfa = state_p(off_fields(f_id%alfa)+1)
    IF (f_id%alfa_d /= 0)           alfa_d = state_p(off_fields(f_id%alfa_d)+1)
    IF (f_id%P_cm /= 0)             P_cm = state_p(off_fields(f_id%P_cm)+1)
    IF (f_id%P_cm_d /= 0)           P_cm_d = state_p(off_fields(f_id%P_cm_d)+1)
    IF (f_id%Chl2N_max /= 0)        Chl2N_max = state_p(off_fields(f_id%Chl2N_max)+1)
    IF (f_id%Chl2N_max_d /= 0)      Chl2N_max_d = state_p(off_fields(f_id%Chl2N_max_d)+1)
    IF (f_id%deg_Chl /= 0)          deg_Chl = state_p(off_fields(f_id%deg_Chl)+1)
    IF (f_id%deg_Chl_d /= 0)        deg_Chl_d = state_p(off_fields(f_id%deg_Chl_d)+1)
    IF (f_id%graz_max /= 0)         graz_max = state_p(off_fields(f_id%graz_max)+1)
    IF (f_id%graz_max2 /= 0)        graz_max2 = state_p(off_fields(f_id%graz_max2)+1)
    IF (f_id%grazEff /= 0)          grazEff = state_p(off_fields(f_id%grazEff)+1)
    IF (f_id%grazEff2 /= 0)         grazEff2 = state_p(off_fields(f_id%grazEff2)+1)
    IF (f_id%lossN /= 0)            lossN = state_p(off_fields(f_id%lossN)+1)
    IF (f_id%lossN_d /= 0)          lossN_d = state_p(off_fields(f_id%lossN_d)+1)
    IF (f_id%lossN_z /= 0)          lossN_z = state_p(off_fields(f_id%lossN_z)+1)
    IF (f_id%lossN_z2 /= 0)         lossN_z2 = state_p(off_fields(f_id%lossN_z2)+1)
    IF (f_id%lossC_z /= 0)          lossC_z = state_p(off_fields(f_id%lossC_z)+1)
    IF (f_id%lossC_z2 /= 0)         lossC_z2 = state_p(off_fields(f_id%lossC_z2)+1)
    IF (f_id%reminN /= 0)           reminN = state_p(off_fields(f_id%reminN)+1)
    IF (f_id%reminC /= 0)           reminC = state_p(off_fields(f_id%reminC)+1)
    IF (f_id%VDet /= 0)             VDet = state_p(off_fields(f_id%VDet)+1)
    IF (f_id%VDet_zoo2 /= 0)           VDet_zoo2 = state_p(off_fields(f_id%VDet_zoo2)+1)
  ! END IF 

END SUBROUTINE distribute_state_pdaf
