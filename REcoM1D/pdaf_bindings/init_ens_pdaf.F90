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
      reminN, reminC

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
      iseed(1)=1
      iseed(2)=3
      iseed(3)=20
      iseed(4)=41
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
      iseed(1)=2
      iseed(2)=5
      iseed(3)=19
      iseed(4)=39
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
      iseed(1)=3
      iseed(2)=7
      iseed(3)=18
      iseed(4)=37
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
      iseed(1)=4
      iseed(2)=11
      iseed(3)=17
      iseed(4)=35
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
      iseed(1)=13
      iseed(2)=5
      iseed(3)=16
      iseed(4)=99
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
      iseed(1)=17
      iseed(2)=6
      iseed(3)=15
      iseed(4)=97
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
      iseed(1)=17
      iseed(2)=6
      iseed(3)=15
      iseed(4)=97
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
      iseed(1)=17
      iseed(2)=6
      iseed(3)=15
      iseed(4)=97
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
      iseed(1)=5
      iseed(2)=13
      iseed(3)=16*member + 5
      iseed(4)=33
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
      iseed(1)=6
      iseed(2)=17
      iseed(3)=15*member + 6
      iseed(4)=31
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
      iseed(1)=7
      iseed(2)=19
      iseed(3)=14*member + 7
      iseed(4)=29  
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
      iseed(1)=8
      iseed(2)=23
      iseed(3)=13*member + 8
      iseed(4)=27
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
      iseed(1)=9
      iseed(2)=29
      iseed(3)=12*member + 9
      iseed(4)=25
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
      iseed(1)=10
      iseed(2)=1
      iseed(3)=11*member + 10
      iseed(4)=23
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
      iseed(1)=11
      iseed(2)=31
      iseed(3)=10*member + 11
      iseed(4)=21
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
      iseed(1)=12
      iseed(2)=37
      iseed(3)=9*member + 12
      iseed(4)=19
      grazEff2_ens = perturb_lognorm_ens(  grazEff2, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%grazEff2) + 1, member) = grazEff2_ens(member)
      END DO
      WRITE(*,*)"grazEff2_ens:", grazEff2_ens
    END IF

    !  lossN
    IF (f_id%lossN /= 0) THEN
      iseed(1)=13
      iseed(2)=41
      iseed(3)=8*member + 13
      iseed(4)=17
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
      iseed(1)=14
      iseed(2)=43
      iseed(3)=7*member + 14
      iseed(4)=15
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
      iseed(1)=15
      iseed(2)=47
      iseed(3)=6*member + 15
      iseed(4)=13
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
      iseed(1)=16
      iseed(2)=53
      iseed(3)=5*member + 16
      iseed(4)=11
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
      iseed(1)=17
      iseed(2)=59
      iseed(3)=4*member + 17
      iseed(4)=9
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
      iseed(1)=18
      iseed(2)=61
      iseed(3)=3*member + 18
      iseed(4)=7
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
      iseed(1)=19
      iseed(2)=67
      iseed(3)=2*member + 19
      iseed(4)=5
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
      iseed(1)=20
      iseed(2)=71
      iseed(3)=1*member + 20
      iseed(4)=3
      reminC_ens = perturb_lognorm_ens(  reminC, &
                                        perturb_scale, &
                                        dim_ens, iseed)     
      DO member = 1, dim_ens
        ens_p(off_fields(f_id%reminC) + 1, member) = reminC_ens(member)
      END DO
      WRITE(*,*)"reminC_ens:", reminC_ens
    END IF




! C 
!   ! ***               Parameter Perturbation               ***
!   DO member = 1, dim_ens


!     !  NCuptakeRatio
!     IF (f_id%NCuptakeRatio /= 0) THEN 
!       iseed(1)=1
!       iseed(2)=3
!       iseed(3)=20*member + 1
!       iseed(4)=41
!       ens_p(off_fields(f_id%NCuptakeRatio) + 1, member)         &
!         = perturb_lognorm(NCuptakeRatio, perturb_scale, iseed)
!     END IF 

!     !  NCUptakeRatio_d
!     IF (f_id%NCUptakeRatio_d /= 0) THEN 
!       iseed(1)=2
!       iseed(2)=5
!       iseed(3)=19*member + 2
!       iseed(4)=39      
!       ens_p(off_fields(f_id%NCUptakeRatio_d) + 1, member)       &
!         = perturb_lognorm(NCUptakeRatio_d, perturb_scale, iseed)
!     END IF 

!     !  k_din
!     IF (f_id%k_din /= 0) THEN 
!       iseed(1)=3
!       iseed(2)=7
!       iseed(3)=18*member + 3
!       iseed(4)=37
!       ens_p(off_fields(f_id%k_din) + 1, member)                 &
!         = perturb_lognorm(k_din, perturb_scale, iseed)
!     END IF

!     !  k_din_d
!     IF (f_id%k_din_d /= 0) THEN 
!       iseed(1)=4
!       iseed(2)=11
!       iseed(3)=17*member + 4
!       iseed(4)=35
!       ens_p(off_fields(f_id%k_din_d) + 1, member)               &
!         = perturb_lognorm(k_din_d, perturb_scale, iseed)
!     END IF

!     !  alfa
!     IF (f_id%alfa /= 0) THEN 
!       iseed(1)=13
!       iseed(2)=5
!       iseed(3)=16*member + 5
!       iseed(4)=99
!       ens_p(off_fields(f_id%alfa) + 1, member)                  &
!         = perturb_lognorm(alfa, perturb_scale, iseed)
!     END IF

!     !  alfa_d
!     IF (f_id%alfa_d /= 0) THEN 
!       iseed(1)=17
!       iseed(2)=6
!       iseed(3)=15*member + 6
!       iseed(4)=97      
!       ens_p(off_fields(f_id%alfa_d) + 1, member)                &
!         = perturb_lognorm(alfa_d, perturb_scale, iseed)
!     END IF

!     !   P_cm
!     IF (f_id%P_cm /= 0) THEN 
!       iseed(1)=17
!       iseed(2)=6
!       iseed(3)=15*member + 6
!       iseed(4)=97      
!       ens_p(off_fields(f_id%P_cm) + 1, member)                  &
!         = perturb_lognorm(P_cm, perturb_scale, iseed)
!     END IF
!     !   P_cm_d
!     IF (f_id%P_cm_d /= 0) THEN 
!       iseed(1)=17
!       iseed(2)=6
!       iseed(3)=15*member + 6
!       iseed(4)=97      
!       ens_p(off_fields(f_id%P_cm_d) + 1, member)                &
!         = perturb_lognorm(P_cm_d, perturb_scale, iseed)
!     END IF

!     !  Chl2N_max
!     IF (f_id%Chl2N_max /= 0) THEN 
!       iseed(1)=5
!       iseed(2)=13
!       iseed(3)=16*member + 5
!       iseed(4)=33
!       ens_p(off_fields(f_id%Chl2N_max) + 1, member)             &
!         = perturb_lognorm(Chl2N_max, perturb_scale, iseed)
!     END IF

!     !  Chl2N_max_d
!     IF (f_id%Chl2N_max_d /= 0) THEN
!       iseed(1)=6
!       iseed(2)=17
!       iseed(3)=15*member + 6
!       iseed(4)=31
!       ens_p(off_fields(f_id%Chl2N_max_d) + 1, member)           &
!         = perturb_lognorm(Chl2N_max_d, perturb_scale, iseed)
!     END IF

!     !  deg_Chl
!     IF (f_id%deg_Chl /= 0) THEN
!       iseed(1)=7
!       iseed(2)=19
!       iseed(3)=14*member + 7
!       iseed(4)=29      
!       ens_p(off_fields(f_id%deg_Chl) + 1, member)               &
!         = perturb_lognorm(deg_Chl, perturb_scale, iseed)
!     END IF

!     !  deg_Chl_d
!     IF (f_id%deg_Chl_d /= 0) THEN
!       iseed(1)=8
!       iseed(2)=23
!       iseed(3)=13*member + 8
!       iseed(4)=27
!       ens_p(off_fields(f_id%deg_Chl_d) + 1, member)             &
!         = perturb_lognorm(deg_Chl_d, perturb_scale, iseed)
!     END IF

!     !  graz_max
!     IF (f_id%graz_max /= 0) THEN
!       iseed(1)=9
!       iseed(2)=29
!       iseed(3)=12*member + 9
!       iseed(4)=25
!       ens_p(off_fields(f_id%graz_max) + 1, member)              &
!         = perturb_lognorm(graz_max, perturb_scale, iseed)
!     END IF

!     !  graz_max2
!     IF (f_id%graz_max2 /= 0) THEN
!       iseed(1)=10
!       iseed(2)=1
!       iseed(3)=11*member + 10
!       iseed(4)=23
!       ens_p(off_fields(f_id%graz_max2) + 1, member)              &
!         = perturb_lognorm(graz_max2, perturb_scale, iseed)
!     END IF

!     !  grazEff
!     IF (f_id%grazEff /= 0) THEN
!       iseed(1)=11
!       iseed(2)=31
!       iseed(3)=10*member + 11
!       iseed(4)=21
!       ens_p(off_fields(f_id%grazEff) + 1, member)              &
!         = perturb_lognorm(grazEff, perturb_scale, iseed)
!     END IF

!     !  grazEff2
!     IF (f_id%grazEff2 /= 0) THEN
!       iseed(1)=12
!       iseed(2)=37
!       iseed(3)=9*member + 12
!       iseed(4)=19
!       ens_p(off_fields(f_id%grazEff2) + 1, member)              &
!         = perturb_lognorm(grazEff2, perturb_scale, iseed)
!     END IF

!     !  lossN
!     IF (f_id%lossN /= 0) THEN
!       iseed(1)=13
!       iseed(2)=41
!       iseed(3)=8*member + 13
!       iseed(4)=17
!       ens_p(off_fields(f_id%lossN) + 1, member)              &
!         = perturb_lognorm(lossN, perturb_scale, iseed)
!     END IF

!     !  lossN_d
!     IF (f_id%lossN_d /= 0) THEN
!       iseed(1)=14
!       iseed(2)=43
!       iseed(3)=7*member + 14
!       iseed(4)=15
!       ens_p(off_fields(f_id%lossN_d) + 1, member)              &
!         = perturb_lognorm(lossN_d, perturb_scale, iseed)
!     END IF

!     !  lossN_z
!     IF (f_id%lossN_z /= 0) THEN
!       iseed(1)=15
!       iseed(2)=47
!       iseed(3)=6*member + 15
!       iseed(4)=13
!       ens_p(off_fields(f_id%lossN_z) + 1, member)              &
!         = perturb_lognorm(lossN_z, perturb_scale, iseed)
!     END IF

!     !  lossC_z
!     IF (f_id%lossC_z /= 0) THEN
!       iseed(1)=16
!       iseed(2)=53
!       iseed(3)=5*member + 16
!       iseed(4)=11
!       lossC_z = perturb_lognorm(lossC_z, perturb_scale, iseed)
!       ens_p(off_fields(f_id%lossC_z) + 1, member)              &
!         = perturb_lognorm(lossC_z, perturb_scale, iseed)
!     END IF

!     !  lossN_z2
!     IF (f_id%lossC_z /= 0) THEN
!       iseed(1)=17
!       iseed(2)=59
!       iseed(3)=4*member + 17
!       iseed(4)=9
!       ens_p(off_fields(f_id%lossN_z2) + 1, member)              &
!         = perturb_lognorm(lossN_z2, perturb_scale, iseed)
!     END IF

!     !  lossC_z2
!     IF (f_id%lossC_z2 /= 0) THEN
!       iseed(1)=18
!       iseed(2)=61
!       iseed(3)=3*member + 18
!       iseed(4)=7
!       ens_p(off_fields(f_id%lossC_z2) + 1, member)              &
!         = perturb_lognorm(lossC_z2, perturb_scale, iseed)
!     END IF

!     !  reminN
!     IF (f_id%reminN /= 0) THEN
!       iseed(1)=19
!       iseed(2)=67
!       iseed(3)=2*member + 19
!       iseed(4)=5
!       ens_p(off_fields(f_id%reminN) + 1, member)              &
!         = perturb_lognorm(reminN, perturb_scale, iseed)
!     END IF

!     !  reminC
!     IF (f_id%reminC /= 0) THEN
!       iseed(1)=20
!       iseed(2)=71
!       iseed(3)=1*member + 20
!       iseed(4)=3
!       ens_p(off_fields(f_id%reminC) + 1, member)              &
!         = perturb_lognorm(reminC, perturb_scale, iseed)
!     END IF

!   END DO 









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





END SUBROUTINE init_ens_pdaf
