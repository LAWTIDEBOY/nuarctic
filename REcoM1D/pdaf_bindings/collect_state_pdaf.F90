!$Id: collect_state_pdaf.F90 871 2021-11-22 16:44:34Z lnerger $
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
!! * 2004-11 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE collect_state_pdaf(dim_p, state_p)



  USE ocean_module, &
      ONLY: tr_arr


  USE mod_assimilation, &
      ONLY: off_fields, dim_fields, dim_field_1d, f_id, tr_id

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector

! *** local variables ***


! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  state_p(off_fields(f_id%DIN)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%DIN))

  state_p(off_fields(f_id%DIC)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%DIC))

  state_p(off_fields(f_id%DSi)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%DSi))

  state_p(off_fields(f_id%NanoN)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%NanoN))

  state_p(off_fields(f_id%NanoC)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%NanoC))

  state_p(off_fields(f_id%NanoChl)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%NanoChl))

  state_p(off_fields(f_id%DiaN)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%DiaN))

  state_p(off_fields(f_id%DIN)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%DIN))

  state_p(off_fields(f_id%DiaC)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%DiaC))

  state_p(off_fields(f_id%DiaChl)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%DiaChl))

  state_p(off_fields(f_id%DiaSi)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%DiaSi))

  state_p(off_fields(f_id%NanoCaCO3)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%NanoCaCO3))

  state_p(off_fields(f_id%DON)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%DON))

  state_p(off_fields(f_id%DOC)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%DOC))

  state_p(off_fields(f_id%TotChl)+1:dim_field_1d) = LOG(tr_arr(1:dim_field_1d, tr_id%NanoChl) + tr_arr(1:dim_field_1d, tr_id%DiaChl))

END SUBROUTINE collect_state_pdaf
