!$Id: init_ens_pdaf.F90 883 2021-11-27 14:16:40Z lnerger $
!> Initialize ensemble
!!
!! User-supplied routine for PDAF.
!!
!! Used for all filters except SEEK, EnKF, LEnKF
!!
!! The routine is called when the filter is
!! initialized in PDAF_filter_init.
!!
!! More information on this initialization variant can be
!! found on the PDAF web site on ensemble initialization.
!!
!! The routine is called by all filter processes and
!! initializes the ensemble for the local domain.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! *Later revisions - see svn log
!!
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, ens_p, flag)

  USE mod_parallel_pdaf, &
      ONLY: mype_world, task_id, mype_filter

  USE mod_assimilation, &
      ONLY: dim_state_p, screen, n_fields_1d, n_fields_0d, n_fields,          &
      n_params, off_fields, dim_fields, dim_field_1d, f_id, tr_id,            &
      parameter_estimation,  perturb_scale

  USE recom_config, &             ! REcoM parameters
      ONLY: NCuptakeRatio, NCUptakeRatio_d, k_din, k_din_d, alfa, alfa_d, P_cm,     & 
      P_cm_d, Chl2N_max, Chl2N_max_d, deg_Chl, deg_Chl_d, graz_max, graz_max2,      &
      grazEff, grazEff2, lossN, lossN_d, lossN_z, lossC_z, lossN_z2, lossC_z2,      &
      reminN, reminC, VDet, VDet_zoo2

  USE mod_perturbation_pdaf, &
      ONLY:  perturb_lognorm, perturb_lognorm_ens, perturb_beta_ens

  IMPLICIT NONE

  ! *** Arguments ***
  INTEGER, INTENT(in)         :: filtertype              !< Type of filter to initialize
  INTEGER, INTENT(in)         :: dim_p                   !< PE-local state dimension
  INTEGER, INTENT(in)         :: dim_ens                 !< Size of ensemble
  REAL(kind=8), INTENT(out)   :: state_p(dim_p)          !< PE-local model state
  REAL(kind=8), INTENT(out)   :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for SEIK
  REAL(kind=8), INTENT(out)   :: ens_p(dim_p, dim_ens)   !< PE-local state ensemble
  INTEGER, INTENT(inout)      :: flag                 !< PDAF status flag


  ! *** local variables ***
  INTEGER       :: i, j, cnt, member ! counter 
  REAL(kind=8)  :: mean_state(dim_p)       ! Ensemble mean state
  REAL(kind=8)  :: varience_state(dim_p)   ! Varience of ensemble state 
  REAL(kind=8)  :: inv_value              ! Inverse value 
  INTEGER       :: iseed(4)               ! Seed for random number
  ! parameters ensemble       
  REAL(kind=8)  :: NCuptakeRatio_ens(dim_ens)
  REAL(kind=8)  :: NCUptakeRatio_d_ens(dim_ens)
  REAL(kind=8)  :: k_din_ens(dim_ens)
  REAL(kind=8)  :: k_din_d_ens(dim_ens)
  REAL(kind=8)  :: alfa_ens(dim_ens)
  REAL(kind=8)  :: alfa_d_ens(dim_ens)
  REAL(kind=8)  :: P_cm_ens(dim_ens)
  REAL(kind=8)  :: P_cm_d_ens(dim_ens)
  REAL(kind=8)  :: Chl2N_max_ens(dim_ens)
  REAL(kind=8)  :: Chl2N_max_d_ens(dim_ens)
  REAL(kind=8)  :: deg_Chl_ens(dim_ens)
  REAL(kind=8)  :: deg_Chl_d_ens(dim_ens)
  REAL(kind=8)  :: graz_max_ens(dim_ens)
  REAL(kind=8)  :: graz_max2_ens(dim_ens)
  REAL(kind=8)  :: grazEff_ens(dim_ens)
  REAL(kind=8)  :: grazEff2_ens(dim_ens)
  REAL(kind=8)  :: lossN_ens(dim_ens)
  REAL(kind=8)  :: lossN_d_ens(dim_ens)
  REAL(kind=8)  :: lossN_z_ens(dim_ens)
  REAL(kind=8)  :: lossN_z2_ens(dim_ens)
  REAL(kind=8)  :: lossC_z_ens(dim_ens)
  REAL(kind=8)  :: lossC_z2_ens(dim_ens)
  REAL(kind=8)  :: reminN_ens(dim_ens)
  REAL(kind=8)  :: reminC_ens(dim_ens)

  REAL(kind=8)  :: VDet_ens(dim_ens)
  REAL(kind=8)  :: VDet_zoo2_ens(dim_ens)


  ! **********************
  ! *** INITIALIZATION ***
  ! **********************

  screen = 1

  IF (mype_filter == 0) THEN
    WRITE(*, '(9x, a)') '-- initialize ensemble from model state.'
    WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
  ENDIF

    state_p = 0.0D+00

  ! collect the initial state vectors from the model
  CALL collect_state_pdaf(dim_p, state_p)


  IF (mype_filter == 0)  WRITE (*,'(9x, a)') '--- generate state ensemble with zero spread'

  ens_p = 0.0D+00

  DO member=1, dim_ens
    ens_p(1 : off_fields(n_fields+1), member) = state_p( 1 : off_fields(n_fields+1))
  ENDDO



 
  ! ***               Parameter Perturbation               ***
      !  NCuptakeRatio
    IF (f_id%NCuptakeRatio /= 0) THEN 
      iseed = (/478, 1927, 221, 3605/)
      NCuptakeRatio_ens = perturb_lognorm_ens(  NCuptakeRatio, &
                                                perturb_scale, &
                                                dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%NCuptakeRatio) + 1, member) = NCuptakeRatio_ens(member)
      END DO
      WRITE(*,*)"NCuptakeRatio_ens:", NCuptakeRatio_ens
    END IF 

    !  NCUptakeRatio_d
    IF (f_id%NCUptakeRatio_d /= 0) THEN 
      iseed = (/126, 293, 1870, 177/)
      NCUptakeRatio_d_ens = perturb_lognorm_ens(  NCUptakeRatio_d, &
                                                  perturb_scale, &
                                                  dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%NCUptakeRatio_d) + 1, member) = NCUptakeRatio_d_ens(member)
      END DO 
      WRITE(*,*)"NCUptakeRatio_d_ens:", NCUptakeRatio_d_ens
    END IF 

    !  k_din
    IF (f_id%k_din /= 0) THEN 
      iseed = (/3786, 1403, 1552, 1069/)
      k_din_ens = perturb_lognorm_ens(  k_din, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%k_din) + 1, member) = k_din_ens(member)
      END DO 
      WRITE(*,*)"k_din_ens:", k_din_ens
    END IF

    !  k_din_d
    IF (f_id%k_din_d /= 0) THEN 
      iseed = (/3035, 2794, 2464, 2505/)
      k_din_d_ens = perturb_lognorm_ens(  k_din_d, &
                                          perturb_scale, &
                                          dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%k_din_d) + 1, member) = k_din_d_ens(member)
      END DO 
      WRITE(*,*)"k_din_d_ens:", k_din_d_ens
    END IF

    !  alfa
    IF (f_id%alfa /= 0) THEN 
      iseed = (/2911, 2490, 1042, 1441/)
      alfa_ens = perturb_lognorm_ens(  alfa, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%alfa) + 1, member) = alfa_ens(member)
      END DO 
      WRITE(*,*)"alfa_ens:", alfa_ens
    END IF

    !  alfa_d
    IF (f_id%alfa_d /= 0) THEN 
      iseed = (/379, 2434, 13, 1343/)
      alfa_d_ens = perturb_lognorm_ens(  alfa_d, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%alfa_d) + 1, member) = alfa_d_ens(member)
      END DO
      WRITE(*,*)"alfa_d_ens:", alfa_d_ens
    END IF

    !   P_cm
    IF (f_id%P_cm /= 0) THEN 
      iseed = (/2133, 3284, 2572, 1797/)
      P_cm_ens = perturb_lognorm_ens(  P_cm, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%P_cm) + 1, member) = P_cm_ens(member)
      END DO
      WRITE(*,*)"P_cm_ens:", P_cm_ens
    END IF
    !   P_cm_d
    IF (f_id%P_cm_d /= 0) THEN 
      iseed = (/3769, 1543, 29, 671/)
      P_cm_d_ens = perturb_lognorm_ens(  P_cm_d, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%P_cm_d) + 1, member) = P_cm_d_ens(member)
      END DO
      WRITE(*,*)"P_cm_d_ens:", P_cm_d_ens
    END IF

    !  Chl2N_max
    IF (f_id%Chl2N_max /= 0) THEN 
      iseed = (/503, 3237, 1066, 2627/)
      Chl2N_max_ens = perturb_lognorm_ens(  Chl2N_max, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%Chl2N_max) + 1, member) = Chl2N_max_ens(member)
      END DO
      WRITE(*,*)"Chl2N_max_ens:", Chl2N_max_ens
    END IF

    !  Chl2N_max_d
    IF (f_id%Chl2N_max_d /= 0) THEN
      iseed = (/3056, 3948, 2634, 4015/)
      Chl2N_max_d_ens = perturb_lognorm_ens(  Chl2N_max_d, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%Chl2N_max_d) + 1, member) = Chl2N_max_d_ens(member)
      END DO
      WRITE(*,*)"Chl2N_max_d_ens:", Chl2N_max_d_ens
    END IF

    !  deg_Chl
    IF (f_id%deg_Chl /= 0) THEN
      iseed = (/1029, 686, 491, 3863/)
      deg_Chl_ens = perturb_lognorm_ens(  deg_Chl, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%deg_Chl) + 1, member) = deg_Chl_ens(member)
      END DO
      WRITE(*,*)"deg_Chl_ens:", deg_Chl_ens
    END IF

    !  deg_Chl_d
    IF (f_id%deg_Chl_d /= 0) THEN
      iseed = (/1084, 3015, 2325, 1877/)
      deg_Chl_d_ens = perturb_lognorm_ens(  deg_Chl_d, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%deg_Chl_d) + 1, member) = deg_Chl_d_ens(member)
      END DO
      WRITE(*,*)"deg_Chl_d_ens:", deg_Chl_d_ens
    END IF

    !  graz_max
    IF (f_id%graz_max /= 0) THEN
      iseed = (/2744, 2001, 2997, 2729/)
      graz_max_ens = perturb_lognorm_ens(  graz_max, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%graz_max) + 1, member) = graz_max_ens(member)
      END DO
      WRITE(*,*)"graz_max_ens:", graz_max_ens
    END IF

    !  graz_max2
    IF (f_id%graz_max2 /= 0) THEN
      iseed = (/451, 2858, 3472, 3871/)
      graz_max2_ens = perturb_lognorm_ens(  graz_max2, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%graz_max2) + 1, member) = graz_max2_ens(member)
      END DO
      WRITE(*,*)"graz_max2_ens:", graz_max2_ens
    END IF

    !  grazEff
    IF (f_id%grazEff /= 0) THEN
      iseed = (/2331, 2119, 3579, 909/)
      grazEff_ens = perturb_lognorm_ens(  grazEff, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%grazEff) + 1, member) = grazEff_ens(member)
      END DO
      WRITE(*,*)"grazEff_ens:", grazEff_ens
    END IF

    !  grazEff2
    IF (f_id%grazEff2 /= 0) THEN
      iseed = (/1758, 1975, 1378, 2519/)
      grazEff2_ens = perturb_lognorm_ens(grazEff2, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%grazEff2) + 1, member) = grazEff2_ens(member)
      END DO
      WRITE(*,*)"grazEff2_ens:", grazEff2_ens
    END IF

    !  lossN
    IF (f_id%lossN /= 0) THEN
      iseed = (/720, 3498, 2865, 3553/)
      lossN_ens = perturb_lognorm_ens(  lossN, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%lossN) + 1, member) = lossN_ens(member)
      END DO
      WRITE(*,*)"lossN_ens:", lossN_ens
    END IF

    !  lossN_d
    IF (f_id%lossN_d /= 0) THEN
      iseed = (/527, 3342, 46, 521/)
      lossN_d_ens = perturb_lognorm_ens(  lossN_d, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%lossN_d) + 1, member) = lossN_d_ens(member)
      END DO
      WRITE(*,*)"lossN_d_ens:", lossN_d_ens
    END IF

    !  lossN_z
    IF (f_id%lossN_z /= 0) THEN
      iseed = (/2209, 1656, 2423, 1055/)
      lossN_z_ens = perturb_lognorm_ens(  lossN_z, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%lossN_z) + 1, member) = lossN_z_ens(member)
      END DO
      WRITE(*,*)"lossN_z_ens:", lossN_z_ens
    END IF

    !  lossC_z
    IF (f_id%lossC_z /= 0) THEN
      iseed = (/66, 10, 89, 51/)
      lossC_z_ens = perturb_lognorm_ens(  lossC_z, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%lossC_z) + 1, member) = lossC_z_ens(member)
      END DO
      WRITE(*,*)"lossC_z_ens:", lossC_z_ens
    END IF

    !  lossN_z2
    IF (f_id%lossC_z /= 0) THEN
      iseed = (/4, 38, 85, 11/)
      lossC_z_ens = perturb_lognorm_ens(  lossC_z, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%lossC_z) + 1, member) = lossC_z_ens(member)
      END DO
      WRITE(*,*)"lossC_z_ens:", lossC_z_ens
    END IF

    !  lossC_z2
    IF (f_id%lossC_z2 /= 0) THEN
      iseed = (/21, 34, 36, 83/)
      lossC_z2_ens = perturb_lognorm_ens(  lossC_z2, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%lossC_z2) + 1, member) = lossC_z2_ens(member)
      END DO
      WRITE(*,*)"lossC_z2_ens:", lossC_z2_ens
    END IF

    !  reminN
    IF (f_id%reminN /= 0) THEN
      iseed = (/31, 71, 90, 73/)
      reminN_ens = perturb_lognorm_ens(  reminN, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%reminN) + 1, member) = reminN_ens(member)
      END DO
      WRITE(*,*)"reminN_ens:", reminN_ens
    END IF

    !  reminC
    IF (f_id%reminC /= 0) THEN
      iseed = (/68, 7, 74, 35/)
      reminC_ens = perturb_lognorm_ens(  reminC, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%reminC) + 1, member) = reminC_ens(member)
      END DO
      WRITE(*,*)"reminC_ens:", reminC_ens
    END IF






    !  VDet
    IF (f_id%VDet /= 0) THEN
      iseed = (/68, 7, 74, 33/)
      VDet_ens = perturb_lognorm_ens(  VDet, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%VDet) + 1, member) = VDet_ens(member)
      END DO
      WRITE(*,*)"VDet_ens:", VDet_ens
    END IF



    !  VDet_zoo2
    IF (f_id%VDet_zoo2 /= 0) THEN
      iseed = (/68, 7, 74, 31/)
      VDet_zoo2_ens = perturb_lognorm_ens(  VDet_zoo2, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%VDet_zoo2) + 1, member) = VDet_zoo2_ens(member)
      END DO
      WRITE(*,*)"VDet_zoo2_ens:", VDet_zoo2_ens
    END IF















  mean_state = 0.0D+00
  DO member = 1, dim_ens
     DO i = 1, dim_p
        mean_state(i) = mean_state(i) + ens_p(i, member)
     END DO
  END DO
  inv_value = 1.0D+00/REAL(dim_ens, 8)
  mean_state(:) = inv_value * mean_state(:)

  IF (f_id%DIN /= 0) WRITE (*,*) 'DIN:', &
    mean_state(off_fields(f_id%DIN)+1:off_fields(f_id%DIN)+dim_field_1d)

  IF (f_id%DIC /= 0) WRITE (*,*) 'DIC:', &
    mean_state(off_fields(f_id%DIC)+1:off_fields(f_id%DIC)+dim_field_1d)

  IF (f_id%DSi /= 0) WRITE (*,*) 'DSi:', &
    mean_state(off_fields(f_id%DSi)+1:off_fields(f_id%DSi)+dim_field_1d)

  IF (f_id%NanoN /= 0) WRITE (*,*) 'NanoN:', &
    mean_state(off_fields(f_id%NanoN)+1:off_fields(f_id%NanoN)+dim_field_1d)

  IF (f_id%NanoC /= 0) WRITE (*,*) 'NanoC:', &
    mean_state(off_fields(f_id%NanoC)+1:off_fields(f_id%NanoC)+dim_field_1d)

  IF (f_id%NanoChl /= 0) WRITE (*,*) 'NanoChl:', &
    mean_state(off_fields(f_id%NanoChl)+1:off_fields(f_id%NanoChl)+dim_field_1d)

  IF (f_id%DiaN /= 0) WRITE (*,*) 'DiaN:', &
    mean_state(off_fields(f_id%DiaN)+1:off_fields(f_id%DiaN)+dim_field_1d)

  IF (f_id%DiaC /= 0) WRITE (*,*) 'DiaC:', &
    mean_state(off_fields(f_id%DiaC)+1:off_fields(f_id%DiaC)+dim_field_1d)

  IF (f_id%DiaChl /= 0) WRITE (*,*) 'DiaChl:', &
    mean_state(off_fields(f_id%DiaChl)+1:off_fields(f_id%DiaChl)+dim_field_1d)

  IF (f_id%DiaSi /= 0) WRITE (*,*) 'DiaSi:', &
    mean_state(off_fields(f_id%DiaSi)+1:off_fields(f_id%DiaSi)+dim_field_1d)

  IF (f_id%NanoCaCO3 /= 0) WRITE (*,*) 'NanoCaCO3:', &
    mean_state(off_fields(f_id%NanoCaCO3)+1:off_fields(f_id%NanoCaCO3)+dim_field_1d)

  IF (f_id%DON /= 0) WRITE (*,*) 'DON:', &
    mean_state(off_fields(f_id%DON)+1:off_fields(f_id%DON)+dim_field_1d)

  IF (f_id%DOC /= 0) WRITE (*,*) 'DOC:', &
    mean_state(off_fields(f_id%DOC)+1:off_fields(f_id%DOC)+dim_field_1d)

  IF (f_id%DetN /= 0) WRITE (*,*) 'DetN:', &
    mean_state(off_fields(f_id%DetN)+1:off_fields(f_id%DetN)+dim_field_1d)

  IF (f_id%DetC /= 0) WRITE (*,*) 'DetC:', &
    mean_state(off_fields(f_id%DetC)+1:off_fields(f_id%DetC)+dim_field_1d)

  IF (f_id%TotChl /= 0) WRITE (*,*) 'TotChl:', &
    mean_state(off_fields(f_id%TotChl)+1:off_fields(f_id%TotChl)+dim_field_1d)

  IF (f_id%NPP /= 0) WRITE (*,*) 'NPP:', mean_state(off_fields(f_id%NPP)+1)




  varience_state = 0.0D+00
  DO member = 1, dim_ens
     DO j = 1, dim_p
        varience_state(j) = varience_state(j)       &
             + (ens_p(j, member) - mean_state(j))   &
             * (ens_p(j, member) - mean_state(j))
     END DO
  END DO
  inv_value = 1.0D+00/REAL(dim_ens - 1, 8)
  varience_state(:) = inv_value * varience_state(:)



  IF (f_id%NCuptakeRatio /= 0) &
    WRITE (*,*) 'NCuptakeRatio=', mean_state(off_fields(f_id%NCuptakeRatio)+1), &
                'Spread=', varience_state(off_fields(f_id%NCuptakeRatio)+1)

  IF (f_id%NCUptakeRatio_d /= 0) &
    WRITE (*,*) 'NCUptakeRatio_d=', mean_state(off_fields(f_id%NCUptakeRatio_d)+1), &
                'Spread=', varience_state(off_fields(f_id%NCUptakeRatio_d)+1)

  IF (f_id%k_din /= 0) &
    WRITE (*,*) 'k_din=', mean_state(off_fields(f_id%k_din)+1), &
                'Spread=', varience_state(off_fields(f_id%k_din)+1)

  IF (f_id%k_din_d /= 0) &
    WRITE (*,*) 'k_din_d=', mean_state(off_fields(f_id%k_din_d)+1), &
                'Spread=', varience_state(off_fields(f_id%k_din_d)+1)

  IF (f_id%alfa /= 0) &
    WRITE (*,*) 'alfa=', mean_state(off_fields(f_id%alfa)+1), &
                'Spread=', varience_state(off_fields(f_id%alfa)+1)

  IF (f_id%alfa_d /= 0) &
    WRITE (*,*) 'alfa_d=', mean_state(off_fields(f_id%alfa_d)+1), &
                'Spread=', varience_state(off_fields(f_id%alfa_d)+1)

  IF (f_id%P_cm /= 0) &
    WRITE (*,*) 'P_cm=', mean_state(off_fields(f_id%P_cm)+1), &
                'Spread=', varience_state(off_fields(f_id%P_cm)+1)

  IF (f_id%P_cm_d /= 0) &
    WRITE (*,*) 'P_cm_d=', mean_state(off_fields(f_id%P_cm_d)+1), &
                'Spread=', varience_state(off_fields(f_id%P_cm_d)+1)

  IF (f_id%Chl2N_max /= 0) &
    WRITE (*,*) 'Chl2N_max=', mean_state(off_fields(f_id%Chl2N_max)+1), &
                'Spread=', varience_state(off_fields(f_id%Chl2N_max)+1)

  IF (f_id%Chl2N_max_d /= 0) &
    WRITE (*,*) 'Chl2N_max_d=', mean_state(off_fields(f_id%Chl2N_max_d)+1), &
                'Spread=', varience_state(off_fields(f_id%Chl2N_max_d)+1)

  IF (f_id%deg_Chl /= 0) &
    WRITE (*,*) 'deg_Chl=', mean_state(off_fields(f_id%deg_Chl)+1), &
                'Spread=', varience_state(off_fields(f_id%deg_Chl)+1)

  IF (f_id%deg_Chl_d /= 0) &
    WRITE (*,*) 'deg_Chl_d=', mean_state(off_fields(f_id%deg_Chl_d)+1), &
                'Spread=', varience_state(off_fields(f_id%deg_Chl_d)+1)

  IF (f_id%graz_max /= 0) &
    WRITE (*,*) 'graz_max=', mean_state(off_fields(f_id%graz_max)+1), &
                'Spread=', varience_state(off_fields(f_id%graz_max)+1)

  IF (f_id%NCuptakeRatio /= 0) &
    WRITE (*,*) 'graz_max2=', mean_state(off_fields(f_id%graz_max2)+1), &
                'Spread=', varience_state(off_fields(f_id%graz_max2)+1)

  IF (f_id%grazEff /= 0) &
    WRITE (*,*) 'grazEff=', mean_state(off_fields(f_id%grazEff)+1), &
                'Spread=', varience_state(off_fields(f_id%grazEff)+1)

  IF (f_id%grazEff2 /= 0) &
    WRITE (*,*) 'grazEff2=', mean_state(off_fields(f_id%grazEff2)+1), &
                'Spread=', varience_state(off_fields(f_id%grazEff2)+1)

  IF (f_id%lossN /= 0) &
    WRITE (*,*) 'lossN=', mean_state(off_fields(f_id%lossN)+1), &
                'Spread=', varience_state(off_fields(f_id%lossN)+1)

  IF (f_id%lossN_d /= 0) &
    WRITE (*,*) 'lossN_d=', mean_state(off_fields(f_id%lossN_d)+1), &
                'Spread=', varience_state(off_fields(f_id%lossN_d)+1)

  IF (f_id%lossN_z /= 0) &
    WRITE (*,*) 'lossN_z=', mean_state(off_fields(f_id%lossN_z)+1), &
                'Spread=', varience_state(off_fields(f_id%lossN_z)+1)

  IF (f_id%lossN_z2 /= 0) &
    WRITE (*,*) 'lossN_z2=', mean_state(off_fields(f_id%lossN_z2)+1), &
                'Spread=', varience_state(off_fields(f_id%lossN_z2)+1)

  IF (f_id%lossC_z /= 0) &
    WRITE (*,*) 'lossC_z=', mean_state(off_fields(f_id%lossC_z)+1), &
                'Spread=', varience_state(off_fields(f_id%lossC_z)+1)

  IF (f_id%lossC_z2 /= 0) &
    WRITE (*,*) 'lossC_z2=', mean_state(off_fields(f_id%lossC_z2)+1), &
                'Spread=', varience_state(off_fields(f_id%lossC_z2)+1)


  IF (f_id%reminN /= 0) &
    WRITE (*,*) 'reminN=', mean_state(off_fields(f_id%reminN)+1), &
                'Spread=', varience_state(off_fields(f_id%reminN)+1)


  IF (f_id%reminC /= 0) &
    WRITE (*,*) 'reminC=', mean_state(off_fields(f_id%reminC)+1), &
                'Spread=', varience_state(off_fields(f_id%reminC)+1)

  
  IF (f_id%VDet /= 0) &
    WRITE (*,*) 'VDet=', mean_state(off_fields(f_id%VDet)+1), &
                'Spread=', varience_state(off_fields(f_id%VDet)+1)

  
  IF (f_id%VDet_zoo2 /= 0) &
    WRITE (*,*) 'VDet_zoo2=', mean_state(off_fields(f_id%VDet_zoo2)+1), &
                'Spread=', varience_state(off_fields(f_id%VDet_zoo2)+1)




END SUBROUTINE init_ens_pdaf
