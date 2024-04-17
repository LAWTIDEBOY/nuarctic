!$Id: prepoststep_ens_pdaf.F90 870 2021-11-22 14:02:55Z lnerger $
!>  Used-defined Pre/Poststep routine for PDAF
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all ensemble filters.
!!
!! The routine is called for global filters (e.g. ESTKF)
!! before the analysis and after the ensemble transformation.
!! For local filters (e.g. LESTKF) the routine is called
!! before and after the loop over all local analysis
!! domains.
!!
!! The routine provides full access to the state
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep
!! operations can be performed here. For example
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by
!! computing the estimated variances.
!! For the offline mode, this routine is the place
!! in which the writing of the analysis ensemble
!! can be performed.
!!
!! If a user considers to perform adjustments to the
!! estimates (e.g. for balances), this routine is
!! the right place for it.
!!
!! Implementation for the 2D example
!! without model parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code based on offline_1D
!! * Later revisions - see repository log
!!
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mod_assimilation, &         ! Variables for assimilation
     ONLY: parameter_estimation, step_null, step_assim, out_dir, f_id, tr_id,   &
     off_fields, dim_field_1d, write_state_variable, write_ens

  
  USE mod_parallel_pdaf, &    ! Parallelization variables
       ONLY: mype_world

  USE mod_utils, &
      ONLY: file_exist, error_handler, get_unit,  e_warn, e_err

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step (negative for call after forecast)
  INTEGER, INTENT(in) :: dim_p       !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     !< Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   !< PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   !< PE-local dimension of observation vector
  REAL(kind=8), INTENT(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
  !< (The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.)
  REAL(kind=8), INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  REAL(kind=8), INTENT(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
  INTEGER, INTENT(in) :: flag        !< PDAF status flag


! *** local variables ***
  INTEGER                     :: i, j, member   ! Counters
  ! LOGICAL, SAVE               :: firsttime = .TRUE. ! Routine is called for first time?
  REAL(kind=8)                :: invdim_ens     ! Inverse ensemble size
  REAL(kind=8)                :: invdim_ensm1   ! Inverse of ensemble size minus 1
  REAL(kind=8)                :: rmserror_est   ! estimated RMS error
  REAL(kind=8), ALLOCATABLE   :: variance(:)    ! model state variances
  REAL(kind=8), ALLOCATABLE   :: field(:)       ! global model field
  CHARACTER(len=110)          :: file_name      ! String for ensemble member
  CHARACTER(len=4)            :: ensstr         ! String for ensemble member
  CHARACTER(len=3)            :: anastr         ! String for call type (initial, forecast, analysis)
  CHARACTER(len=256)          :: msgstring      ! String for error handling message 
  INTEGER :: io, iunit


! **********************
! *** INITIALIZATION ***
! **********************

  

  IF (step - step_null==0) THEN
    WRITE (*, '(8x, a)') 'Analyze initial state ensemble'
    anastr = 'ini'
  ELSE
    IF (step<0) THEN
      WRITE (*, '(8x, a)') 'Analyze and write forecasted state ensemble'
      anastr = 'for'
    ELSE
      WRITE (*, '(8x, a)') 'Analyze and write assimilated state ensemble'
      anastr = 'ana'
    END IF
  END IF

  invdim_ens    = 1.0 / REAL(dim_ens, 8)
  invdim_ensm1  = 1.0 / REAL(dim_ens - 1, 8)

! **************************************************************
! *** Perform prepoststep for SEIK with re-inititialization. ***
! *** The state and error information is completely in the   ***
! *** ensemble.                                              ***
! *** Also performed for SEIK without re-init at the initial ***
! *** time.                                                  ***
! **************************************************************

  ! *** Compute mean state
  WRITE (*, '(8x, a)') '--- compute ensemble mean'

  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)


  ! *** Compute sampled variances ***
    ! Allocate fields
  ALLOCATE(variance(dim_p))

  variance(:) = 0.0
  DO member = 1, dim_ens
     DO j = 1, dim_p
        variance(j) = variance(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
     END DO
  END DO
  variance(:) = invdim_ensm1 * variance(:)



! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************
  ! ! Initialize numbers
  ! rmserror_est  = 0.0
  ! ! total estimated RMS error
  ! DO i = 1, dim_p
  !    rmserror_est = rmserror_est + variance(i)
  ! ENDDO
  ! rmserror_est = SQRT(rmserror_est / dim_p)
  ! ! *****************
  ! ! *** Screen IO ***
  ! ! *****************
  ! ! Output RMS errors given by sampled covar matrix
  ! WRITE (*, '(12x, a, es12.4)') &
  !      'RMS error according to sampled variance: ', rmserror_est



  !*******************
  !*** File output ***
  !*******************


  out_dir = '../da_outputs/'

  IF (write_state_variable) THEN 

    IF (f_id%DIN /= 0) THEN

      file_name = TRIM(out_dir)//'DIN.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND',              &
            ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DIN  ', anastr, step, '=', &
            state_p(off_fields(f_id%DIN)+1:off_fields(f_id%DIN)+dim_field_1d)

      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%DIC /= 0) THEN

      file_name = TRIM(out_dir)//'DIC.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND',              &
              ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DIC  ', anastr, step, '=', &
            state_p(off_fields(f_id%DIC)+1:off_fields(f_id%DIC)+dim_field_1d)

      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%DSi /= 0) THEN

      file_name = TRIM(out_dir)//'DSi.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND',              &
              ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DSi  ', anastr, step, '=', &
            state_p(off_fields(f_id%DSi)+1:off_fields(f_id%DSi)+dim_field_1d)

      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%NanoN /= 0) THEN

      file_name = TRIM(out_dir)//'NanoN.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND',              &
                ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'NanoN  ', anastr, step, '=', &
            state_p(off_fields(f_id%NanoN)+1:off_fields(f_id%NanoN)+dim_field_1d)

      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%NanoC /= 0) THEN

      file_name = TRIM(out_dir)//'NanoC.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'NanoC  ', anastr, step, '=', &
            state_p(off_fields(f_id%NanoC)+1:off_fields(f_id%NanoC)+dim_field_1d)

      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%NanoChl /= 0) THEN

      file_name = TRIM(out_dir)//'NanoChl.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'NanoChl  ', anastr, step, '=', &
            state_p(off_fields(f_id%NanoChl)+1:off_fields(f_id%NanoChl)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%DiaN /= 0) THEN

      file_name = TRIM(out_dir)//'DiaN.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DiaN  ', anastr, step, '=', &
            state_p(off_fields(f_id%DiaN)+1:off_fields(f_id%DiaN)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%DiaC /= 0) THEN

      file_name = TRIM(out_dir)//'DiaC.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DiaC  ', anastr, step, '=', &
            state_p(off_fields(f_id%DiaC)+1:off_fields(f_id%DiaC)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%DiaChl /= 0) THEN

      file_name = TRIM(out_dir)//'DiaChl.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DiaChl  ', anastr, step, '=', &
            state_p(off_fields(f_id%DiaChl)+1:off_fields(f_id%DiaChl)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%DiaSi /= 0) THEN

      file_name = TRIM(out_dir)//'DiaSi.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DiaSi  ', anastr, step, '=', &
            state_p(off_fields(f_id%DiaSi)+1:off_fields(f_id%DiaSi)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%NanoCaCO3 /= 0) THEN

      file_name = TRIM(out_dir)//'NanoCaCO3.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'NanoCaCO3  ', anastr, step, '=', &
            state_p(off_fields(f_id%NanoCaCO3)+1:off_fields(f_id%NanoCaCO3)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%DON /= 0) THEN

      file_name = TRIM(out_dir)//'DON.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DON  ', anastr, step, '=', &
            state_p(off_fields(f_id%DON)+1:off_fields(f_id%DON)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%DOC /= 0) THEN

      file_name = TRIM(out_dir)//'DOC.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DOC  ', anastr, step, '=', &
            state_p(off_fields(f_id%DOC)+1:off_fields(f_id%DOC)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%DetN /= 0) THEN

      file_name = TRIM(out_dir)//'DetN.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DetN  ', anastr, step, '=', &
            state_p(off_fields(f_id%DetN)+1:off_fields(f_id%DetN)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 


    IF (f_id%DetC /= 0) THEN

      file_name = TRIM(out_dir)//'DetC.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'DetC  ', anastr, step, '=', &
            state_p(off_fields(f_id%DetC)+1:off_fields(f_id%DetC)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 




    IF (f_id%TotChl /= 0) THEN

      file_name = TRIM(out_dir)//'TotChl.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'TotChl  ', anastr, step, '=', &
            state_p(off_fields(f_id%TotChl)+1:off_fields(f_id%TotChl)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 





    IF (f_id%NPP /= 0) THEN

      file_name = TRIM(out_dir)//'NPP.dat'
      iunit = get_unit()

      IF (file_exist(file_name) ) THEN
        OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
      ELSE
        OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
      END IF

      IF (io /= 0) THEN
        WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
        CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
      ELSE
        WRITE(iunit, *) 'NPP  ', anastr, step, '=', &
            state_p(off_fields(f_id%NPP)+1:off_fields(f_id%NPP)+dim_field_1d)
      END IF

      CLOSE(UNIT=iunit)

    END IF 

  END IF 




  ! IF (parameter_estimation) THEN 


      IF (f_id%NCuptakeRatio /= 0) THEN

        file_name = TRIM(out_dir)//'NCuptakeRatio.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'NCuptakeRatio  ', anastr, step, '=', &
              state_p(off_fields(f_id%NCuptakeRatio)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%NCUptakeRatio_d /= 0) THEN

        file_name = TRIM(out_dir)//'NCUptakeRatio_d.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'NCUptakeRatio_d  ', anastr, step, '=', &
              state_p(off_fields(f_id%NCUptakeRatio_d)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%SiCUptakeRatio /= 0) THEN

        file_name = TRIM(out_dir)//'SiCUptakeRatio.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'SiCUptakeRatio  ', anastr, step, '=', &
              state_p(off_fields(f_id%SiCUptakeRatio)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%k_din /= 0) THEN

        file_name = TRIM(out_dir)//'k_din.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'k_din  ', anastr, step, '=', &
              state_p(off_fields(f_id%k_din)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%k_din_d /= 0) THEN

        file_name = TRIM(out_dir)//'k_din_d.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'k_din_d  ', anastr, step, '=', &
              state_p(off_fields(f_id%k_din_d)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%alfa /= 0) THEN

        file_name = TRIM(out_dir)//'alfa.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'alfa  ', anastr, step, '=', &
              state_p(off_fields(f_id%alfa)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%alfa_d /= 0) THEN

        file_name = TRIM(out_dir)//'alfa_d.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'alfa_d  ', anastr, step, '=', &
              state_p(off_fields(f_id%alfa_d)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%P_cm /= 0) THEN

        file_name = TRIM(out_dir)//'P_cm.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'P_cm  ', anastr, step, '=', &
              state_p(off_fields(f_id%P_cm)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%P_cm_d /= 0) THEN

        file_name = TRIM(out_dir)//'P_cm_d.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'P_cm_d  ', anastr, step, '=', &
              state_p(off_fields(f_id%P_cm_d)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%Chl2N_max /= 0) THEN

        file_name = TRIM(out_dir)//'Chl2N_max.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'Chl2N_max  ', anastr, step, '=', &
              state_p(off_fields(f_id%Chl2N_max)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%Chl2N_max_d /= 0) THEN

        file_name = TRIM(out_dir)//'Chl2N_max_d.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'Chl2N_max_d  ', anastr, step, '=', &
              state_p(off_fields(f_id%Chl2N_max_d)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%deg_Chl /= 0) THEN

        file_name = TRIM(out_dir)//'deg_Chl.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'deg_Chl  ', anastr, step, '=', &
              state_p(off_fields(f_id%deg_Chl)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%deg_Chl_d /= 0) THEN

        file_name = TRIM(out_dir)//'deg_Chl_d.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'deg_Chl_d  ', anastr, step, '=', &
              state_p(off_fields(f_id%deg_Chl_d)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%graz_max /= 0) THEN

        file_name = TRIM(out_dir)//'graz_max.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'graz_max  ', anastr, step, '=', &
              state_p(off_fields(f_id%graz_max)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%graz_max2 /= 0) THEN

        file_name = TRIM(out_dir)//'graz_max2.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'graz_max2  ', anastr, step, '=', &
              state_p(off_fields(f_id%graz_max2)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%grazEff /= 0) THEN

        file_name = TRIM(out_dir)//'grazEff.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'grazEff  ', anastr, step, '=', &
              state_p(off_fields(f_id%grazEff)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%grazEff2 /= 0) THEN

        file_name = TRIM(out_dir)//'grazEff2.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'grazEff2  ', anastr, step, '=', &
              state_p(off_fields(f_id%grazEff2)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%lossN /= 0) THEN

        file_name = TRIM(out_dir)//'lossN.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'lossN  ', anastr, step, '=', &
              state_p(off_fields(f_id%lossN)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%lossN_d /= 0) THEN

        file_name = TRIM(out_dir)//'lossN_d.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'lossN_d  ', anastr, step, '=', &
              state_p(off_fields(f_id%lossN_d)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%lossN_z /= 0) THEN

        file_name = TRIM(out_dir)//'lossN_z.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'lossN_z  ', anastr, step, '=', &
              state_p(off_fields(f_id%lossN_z)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%lossN_z2 /= 0) THEN

        file_name = TRIM(out_dir)//'lossN_z2.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'lossN_z2  ', anastr, step, '=', &
              state_p(off_fields(f_id%lossN_z2)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%lossC_z /= 0) THEN

        file_name = TRIM(out_dir)//'lossC_z.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'lossC_z  ', anastr, step, '=', &
              state_p(off_fields(f_id%lossC_z)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%lossC_z2 /= 0) THEN

        file_name = TRIM(out_dir)//'lossC_z2.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'lossC_z2  ', anastr, step, '=', &
              state_p(off_fields(f_id%lossC_z2)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%reminN /= 0) THEN

        file_name = TRIM(out_dir)//'reminN.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'reminN  ', anastr, step, '=', &
              state_p(off_fields(f_id%reminN)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 


      IF (f_id%reminC /= 0) THEN

        file_name = TRIM(out_dir)//'reminC.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'reminC  ', anastr, step, '=', &
              state_p(off_fields(f_id%reminC)+1)

        END IF

        CLOSE(UNIT=iunit)
      
      END IF 



      IF (f_id%VDet /= 0) THEN

        file_name = TRIM(out_dir)//'VDet.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'VDet  ', anastr, step, '=', &
              state_p(off_fields(f_id%VDet)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 




      IF (f_id%VDet_zoo2 /= 0) THEN

        file_name = TRIM(out_dir)//'VDet_zoo2.dat'
        iunit = get_unit()

        IF (file_exist(file_name) ) THEN
          OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
        ELSE
          OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
        END IF

        IF (io /= 0) THEN
          WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
          CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
        ELSE
          WRITE(iunit, *) 'VDet_zoo2  ', anastr, step, '=', &
              state_p(off_fields(f_id%VDet_zoo2)+1)

        END IF

        CLOSE(UNIT=iunit)

      END IF 









!!     write out ensemble of parameters 
  IF (write_ens) THEN

    DO member = 1, dim_ens
      ! formate string for ensemble members
      IF (member<10) THEN
        WRITE(ensstr,'(I1)') member
      ELSEIF (member>=10 .AND. member<100) THEN
        WRITE(ensstr,'(I2)') member
      ELSEIF (member>=100 .AND. member<1000) THEN
        WRITE(ensstr,'(I3)') member
      ELSE
        WRITE(ensstr,'(I4)') member
      ENDIF




      IF (write_state_variable) THEN 

        IF (f_id%DIN /= 0) THEN

          file_name = TRIM(out_dir)//'DIN_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND',              &
                ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DIN  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DIN)+1:off_fields(f_id%DIN)+dim_field_1d, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%DIC /= 0) THEN

          file_name = TRIM(out_dir)//'DIC_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND',              &
                  ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DIC  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DIC)+1:off_fields(f_id%DIC)+dim_field_1d, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%DSi /= 0) THEN

          file_name = TRIM(out_dir)//'DSi_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND',              &
                  ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DSi  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DSi)+1:off_fields(f_id%DSi)+dim_field_1d, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%NanoN /= 0) THEN

          file_name = TRIM(out_dir)//'NanoN_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND',              &
                    ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'NanoN  ', anastr, step, '=', &
                ens_p(off_fields(f_id%NanoN)+1:off_fields(f_id%NanoN)+dim_field_1d, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%NanoC /= 0) THEN

          file_name = TRIM(out_dir)//'NanoC_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'NanoC  ', anastr, step, '=', &
                ens_p(off_fields(f_id%NanoC)+1:off_fields(f_id%NanoC)+dim_field_1d, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%NanoChl /= 0) THEN

          file_name = TRIM(out_dir)//'NanoChl_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'NanoChl  ', anastr, step, '=', &
                ens_p(off_fields(f_id%NanoChl)+1:off_fields(f_id%NanoChl)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%DiaN /= 0) THEN

          file_name = TRIM(out_dir)//'DiaN_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DiaN  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DiaN)+1:off_fields(f_id%DiaN)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%DiaC /= 0) THEN

          file_name = TRIM(out_dir)//'DiaC_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DiaC  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DiaC)+1:off_fields(f_id%DiaC)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%DiaChl /= 0) THEN

          file_name = TRIM(out_dir)//'DiaChl_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DiaChl  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DiaChl)+1:off_fields(f_id%DiaChl)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%DiaSi /= 0) THEN

          file_name = TRIM(out_dir)//'DiaSi_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DiaSi  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DiaSi)+1:off_fields(f_id%DiaSi)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%NanoCaCO3 /= 0) THEN

          file_name = TRIM(out_dir)//'NanoCaCO3_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'NanoCaCO3  ', anastr, step, '=', &
                ens_p(off_fields(f_id%NanoCaCO3)+1:off_fields(f_id%NanoCaCO3)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%DON /= 0) THEN

          file_name = TRIM(out_dir)//'DON_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DON  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DON)+1:off_fields(f_id%DON)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%DOC /= 0) THEN

          file_name = TRIM(out_dir)//'DOC_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DOC  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DOC)+1:off_fields(f_id%DOC)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%DetN /= 0) THEN

          file_name = TRIM(out_dir)//'DetN_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DetN  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DetN)+1:off_fields(f_id%DetN)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%DetC /= 0) THEN

          file_name = TRIM(out_dir)//'DetC_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'DetC  ', anastr, step, '=', &
                ens_p(off_fields(f_id%DetC)+1:off_fields(f_id%DetC)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 




        IF (f_id%TotChl /= 0) THEN

          file_name = TRIM(out_dir)//'TotChl_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'TotChl  ', anastr, step, '=', &
                ens_p(off_fields(f_id%TotChl)+1:off_fields(f_id%TotChl)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 





        IF (f_id%NPP /= 0) THEN

          file_name = TRIM(out_dir)//'NPP_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'NPP  ', anastr, step, '=', &
                ens_p(off_fields(f_id%NPP)+1:off_fields(f_id%NPP)+dim_field_1d, member)
          END IF

          CLOSE(UNIT=iunit)

        END IF 

      END IF 

      !! Write out the parameters ensemble to the files
      !! only in analysis step because the the nextt forcast step has the same values. 

      IF (step>=0) THEN
        
        IF (f_id%NCuptakeRatio /= 0) THEN

          file_name = TRIM(out_dir)//'NCuptakeRatio_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'NCuptakeRatio  ', anastr, step, '=', &
                ens_p(off_fields(f_id%NCuptakeRatio)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%NCUptakeRatio_d /= 0) THEN

          file_name = TRIM(out_dir)//'NCUptakeRatio_d_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'NCUptakeRatio_d  ', anastr, step, '=', &
                ens_p(off_fields(f_id%NCUptakeRatio_d)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%SiCUptakeRatio /= 0) THEN

          file_name = TRIM(out_dir)//'SiCUptakeRatio_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'SiCUptakeRatio  ', anastr, step, '=', &
                ens_p(off_fields(f_id%SiCUptakeRatio)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%k_din /= 0) THEN

          file_name = TRIM(out_dir)//'k_din_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'k_din  ', anastr, step, '=', &
                ens_p(off_fields(f_id%k_din)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%k_din_d /= 0) THEN

          file_name = TRIM(out_dir)//'k_din_d_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'k_din_d  ', anastr, step, '=', &
                ens_p(off_fields(f_id%k_din_d)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%alfa /= 0) THEN

          file_name = TRIM(out_dir)//'alfa_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'alfa  ', anastr, step, '=', &
                ens_p(off_fields(f_id%alfa)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%alfa_d /= 0) THEN

          file_name = TRIM(out_dir)//'alfa_d_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'alfa_d  ', anastr, step, '=', &
                ens_p(off_fields(f_id%alfa_d)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%P_cm /= 0) THEN

          file_name = TRIM(out_dir)//'P_cm_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'P_cm  ', anastr, step, '=', &
                ens_p(off_fields(f_id%P_cm)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%P_cm_d /= 0) THEN

          file_name = TRIM(out_dir)//'P_cm_d_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'P_cm_d  ', anastr, step, '=', &
                ens_p(off_fields(f_id%P_cm_d)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%Chl2N_max /= 0) THEN

          file_name = TRIM(out_dir)//'Chl2N_max_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'Chl2N_max  ', anastr, step, '=', &
                ens_p(off_fields(f_id%Chl2N_max)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%Chl2N_max_d /= 0) THEN

          file_name = TRIM(out_dir)//'Chl2N_max_d_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'Chl2N_max_d  ', anastr, step, '=', &
                ens_p(off_fields(f_id%Chl2N_max_d)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%deg_Chl /= 0) THEN

          file_name = TRIM(out_dir)//'deg_Chl_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'deg_Chl  ', anastr, step, '=', &
                ens_p(off_fields(f_id%deg_Chl)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%deg_Chl_d /= 0) THEN

          file_name = TRIM(out_dir)//'deg_Chl_d_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'deg_Chl_d  ', anastr, step, '=', &
                ens_p(off_fields(f_id%deg_Chl_d)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%graz_max /= 0) THEN

          file_name = TRIM(out_dir)//'graz_max_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'graz_max  ', anastr, step, '=', &
                ens_p(off_fields(f_id%graz_max)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%graz_max2 /= 0) THEN

          file_name = TRIM(out_dir)//'graz_max2_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'graz_max2  ', anastr, step, '=', &
                ens_p(off_fields(f_id%graz_max2)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%grazEff /= 0) THEN

          file_name = TRIM(out_dir)//'grazEff_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'grazEff  ', anastr, step, '=', &
                ens_p(off_fields(f_id%grazEff)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%grazEff2 /= 0) THEN

          file_name = TRIM(out_dir)//'grazEff2_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'grazEff2  ', anastr, step, '=', &
                ens_p(off_fields(f_id%grazEff2)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%lossN /= 0) THEN

          file_name = TRIM(out_dir)//'lossN_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'lossN  ', anastr, step, '=', &
                ens_p(off_fields(f_id%lossN)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%lossN_d /= 0) THEN

          file_name = TRIM(out_dir)//'lossN_d_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'lossN_d  ', anastr, step, '=', &
                ens_p(off_fields(f_id%lossN_d)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%lossN_z /= 0) THEN

          file_name = TRIM(out_dir)//'lossN_z_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'lossN_z  ', anastr, step, '=', &
                ens_p(off_fields(f_id%lossN_z)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%lossN_z2 /= 0) THEN

          file_name = TRIM(out_dir)//'lossN_z2_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'lossN_z2  ', anastr, step, '=', &
                ens_p(off_fields(f_id%lossN_z2)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%lossC_z /= 0) THEN

          file_name = TRIM(out_dir)//'lossC_z_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'lossC_z  ', anastr, step, '=', &
                ens_p(off_fields(f_id%lossC_z)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%lossC_z2 /= 0) THEN

          file_name = TRIM(out_dir)//'lossC_z2_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'lossC_z2  ', anastr, step, '=', &
                ens_p(off_fields(f_id%lossC_z2)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%reminN /= 0) THEN

          file_name = TRIM(out_dir)//'reminN_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'reminN  ', anastr, step, '=', &
                ens_p(off_fields(f_id%reminN)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 


        IF (f_id%reminC /= 0) THEN

          file_name = TRIM(out_dir)//'reminC_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'reminC  ', anastr, step, '=', &
                ens_p(off_fields(f_id%reminC)+1, member)

          END IF

          CLOSE(UNIT=iunit)
        
        END IF 



        IF (f_id%VDet /= 0) THEN

          file_name = TRIM(out_dir)//'VDet_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'VDet  ', anastr, step, '=', &
                ens_p(off_fields(f_id%VDet)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 




        IF (f_id%VDet_zoo2 /= 0) THEN

          file_name = TRIM(out_dir)//'VDet_zoo2_'//TRIM(ensstr)//'.dat'
          iunit = get_unit()

          IF (file_exist(file_name) ) THEN
            OPEN(UNIT=iunit, FILE=file_name, STATUS='OLD', POSITION='APPEND', ACTION='WRITE', IOSTAT=io)
          ELSE
            OPEN(UNIT=iunit, FILE=file_name, STATUS='NEW', ACTION='WRITE', IOSTAT=io)
          END IF

          IF (io /= 0) THEN
            WRITE(msgstring, * ) "cannot open (",io,") file ", TRIM(file_name)
            CALL error_handler(e_warn, "prepoststep_ens_pdaf", msgstring )
          ELSE
            WRITE(iunit, *) 'VDet_zoo2  ', anastr, step, '=', &
                ens_p(off_fields(f_id%VDet_zoo2)+1, member)

          END IF

          CLOSE(UNIT=iunit)

        END IF 

      END IF 

    END DO 

  END IF 

! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(variance)

  ! firsttime = .FALSE.

END SUBROUTINE prepoststep_ens_pdaf
