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

  USE mod_assimilation, &
      ONLY: off_fields, dim_fields, dim_field_1d, f_id, tr_id



  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** local variables ***


! *************************************************
! *** Initialize model fields from state vector ***
! *** for process-local model domain            ***
!**************************************************

  tr_arr(1:dim_field_1d, tr_id%DIN) = EXP(state_p(off_fields(f_id%DIN)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%DIC) = EXP(state_p(off_fields(f_id%DIC)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%DSi) = EXP(state_p(off_fields(f_id%DSi)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%NanoN) = EXP(state_p(off_fields(f_id%NanoN)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%NanoC) = EXP(state_p(off_fields(f_id%NanoC)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%NanoChl) = EXP(state_p(off_fields(f_id%NanoChl)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%DiaN) = EXP(state_p(off_fields(f_id%DiaN)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%DiaC) = EXP(state_p(off_fields(f_id%DiaC)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%DiaChl) = EXP(state_p(off_fields(f_id%DiaChl)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%DiaSi) = EXP(state_p(off_fields(f_id%DiaSi)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%NanoCaCO3) = EXP(state_p(off_fields(f_id%NanoCaCO3)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%DON) = EXP(state_p(off_fields(f_id%DON)+1:dim_field_1d))

  tr_arr(1:dim_field_1d, tr_id%DOC) = EXP(state_p(off_fields(f_id%DOC)+1:dim_field_1d))

END SUBROUTINE distribute_state_pdaf
