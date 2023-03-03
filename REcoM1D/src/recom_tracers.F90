module REcoM_tracers

 use general_config
 use recom_config
 use ocean_module
 
 ! remineralisation and sinking specific arrays
 !real(kind=8), allocatable    :: dtr_bf(:), str_bf(:) !  fluxes from benthos
 !real(kind=8), allocatable    :: vert_sink(:)           ! vertical sinking 

contains


!---------------------------------------------------------------------

subroutine REcoM_mixing

!--------------------------------------------------------------------
! routine to solve vertical mixing of tracers due to turbulence 
! The mixing  solves the equation: d/dz(kv*d(tracers)/dz (cf. MOPS)
!--------------------------------------------------------------------

  !use recom_config
  !use ocean_module
  Implicit none
  
  Integer                                       :: tr_num, k, id
  Integer         				:: nzmin, nzmax
  Real(kind=8)                                  :: deepscale  
  Real(kind=8), dimension(:), allocatable       :: dz_trr, vflux
  Real(kind=8), dimension(:), allocatable       :: deepflux, bc_bottom_tracers
  
  
  nzmin = ulevel
  nzmax = nlevel
  
  allocate (dz_trr(nzmax), vflux(nzmax), deepflux(num_tracers), bc_bottom_tracers(num_tracers))
  ! compute vertical discretization step (between center of each cell
  dz_trr                = 0.0d0
  dz_trr(nzmin+1:nzmax) = abs(Z(nzmin:nzmax-1) - Z(nzmin+1:nzmax))
  dz_trr(nzmin)         = hnode(nzmin)/2.0d0

  ! perform the mixing on each tracer  
  do tr_num = 3,num_tracers
    vflux=0.d0
    ! compute the mixing flux in the water column
    do k = max(2,nzmin), nzmax
    	if (isnan(Kz(k))) Kz(k)=1e-7
    	vflux(k) = Kz(k) * (tr_arr(k,tr_num) - tr_arr(k-1,tr_num))/dz_trr(k)
    	if (isnan(vflux(k))) vflux(k) = 0 
    enddo 

    ! rk: no flux at the surface
      
    ! update the tracers accordingly
    do k = max(2,nzmin), nzmax
    	tr_arr(k-1,tr_num) = tr_arr(k-1,tr_num) + vflux(k) / dz_trr(k-1)
        tr_arr(k,tr_num)   = tr_arr(k,tr_num)   - vflux(k) / dz_trr(k)
    enddo
  enddo
  
  ! compute the residual flux across the lower boundary (might be negligible for significant depth)
   bc_bottom_tracers = 0.d0
   ! set boundary conditions at the bottom
   do tr_num = 6, num_tracers
      id = tracer_id(tr_num)
      SELECT CASE (id)
        CASE(1004, 1013)
       ! if (tracer_id(tr_num) .eq. 1004 .or. tracer_id(tr_num) .eq. 1013) then
  	  bc_bottom_tracers(tr_num) = tiny_chl/chl2N_max
  	CASE (1005, 1014)
      !elseif (tracer_id(tr_num) .eq. 1005 .or. tracer_id(tr_num) .eq. 1014) then
  	  bc_bottom_tracers(tr_num) = tiny_chl/chl2N_max/NCmax
  	 CASE (1006, 1015)
      !elseif (tracer_id(tr_num) .eq. 1006 .or. tracer_id(tr_num) .eq. 1015) then
  	  bc_bottom_tracers(tr_num) = tiny_chl
      !elseif (tracer_id(tr_num) .eq. 1010 .or. tracer_id(tr_num) .eq. 1024) then
      	 CASE (1010, 1024)
      	  bc_bottom_tracers(tr_num) = tiny * Redfield  
      !elseif (tracer_id(tr_num) .eq. 1016) then
         CASE (1016)
  	  bc_bottom_tracers(tr_num) = tiny_chl/chl2N_max_d/NCmax_d/SiCmax
         CASE DEFAULT
          bc_bottom_tracers(tr_num) = tiny
      !endif
      END SELECT
   enddo 
   ! scale deep diffusion coefficient for cross-boundary flux, (0-> closed lower boundary) 
   deepscale = 0.d0
   ! compute flux at lower boundary   
   do tr_num = 3, num_tracers
        ! flux
        deepflux(tr_num) = deepscale * Kz(nzmax) *    &
                     (bc_bottom_tracers(tr_num) -    & 
                     tr_arr(nzmax,tr_num)) / dz_trr(nzmax) 
        if (isnan(deepflux(tr_num))) deepflux(tr_num)=0.
         ! tracer update at bottom boundary
         tr_arr(nzmax,tr_num) = tr_arr(nzmax,tr_num) + & 
         deepflux(tr_num)/dz_trr(nzmax)                    
    enddo
  
  deallocate (dz_trr, vflux, deepflux, bc_bottom_tracers)
  
end subroutine
!---------------------------------------------------------------------

subroutine update_tracers

 !use recom_config
 !use ocean_module
 
 implicit none
 
 real(kind=8), dimension(:), allocatable :: dtr_bf, str_bf, vert_sink 
 Integer				 :: tr_num, id
 Integer   				 :: nzmax, nzmin


  ! allocate sinking and remineralisation auxiliary arrays
   allocate(dtr_bf(nlevel-1), str_bf(nlevel-1), vert_sink(nlevel-1))

   do tr_num = 3, num_tracers
     
     ! initialization
     dtr_bf         = 0.d0
     str_bf         = 0.d0
     vert_sink      = 0.d0
     id = tracer_id(tr_num)

     SELECT CASE (id)
        CASE (1001:1003, 1018:1019, 1022)
           !ids: DIN, DIC, Alk, Si, Fe, O2
           ! 1) Remineralization from the benthos
           !    Nutrient fluxes come from the bottom boundary
           !    Unit [mmol/m2/s]
     !if (tracer_id(tr_num) == 1001 .or.    &   ! DIN
     ! 	 tracer_id(tr_num) == 1002 .or.    &   ! DIC
     ! 	 tracer_id(tr_num) == 1003 .or.    &   ! Alk
     !	 tracer_id(tr_num) == 1018 .or.    &   ! Si
     !	 tracer_id(tr_num) == 1019 .or.    &   ! Fe
     !    tracer_id(tr_num) == 1022     ) then  ! Oxy
	   
	   call diff_ver_recom_expl(tr_num, dtr_bf)

           ! update tracer fields
           nzmax=nlevel-1
           nzmin=ulevel
           tr_arr(nzmin:nzmax,tr_num)=tr_arr(nzmin:nzmax,tr_num)+ &
                                        dtr_bf(nzmin:nzmax)
    !end if
        CASE (1004:1008, 1013:1017, 1020:1021, 1025:1028)
        ! ids: DetN, DetC, DetSi, DetCalc
        !      PhyN, PhyC, PhyCalc, PhyChl
        !      DiaN, DiaC, DiaSi, DiaChl
        !      Detz2n, Detz2c, Detz2si, Detz2calc   
        ! 2) Sinking in water column 

    !if (tracer_id(tr_num) == 1007 .or.    &   ! idetn
    !    tracer_id(tr_num) == 1008 .or.    &   ! idetc
    !    tracer_id(tr_num) == 1017 .or.    &   ! idetsi
    !    tracer_id(tr_num) == 1021 .or.    &   ! idetcal
    !    tracer_id(tr_num) == 1004 .or.    &   !iphyn
    !    tracer_id(tr_num) == 1005 .or.    &   !iphyc
    !    tracer_id(tr_num) == 1020 .or.    &   !iphycal
    !    tracer_id(tr_num) == 1006 .or.    &   !ipchl
    !    tracer_id(tr_num) == 1013 .or.    &   !idian
    !    tracer_id(tr_num) == 1014 .or.    &   !idiac
    !    tracer_id(tr_num) == 1016 .or.    &   !idiasi
    !    tracer_id(tr_num) == 1015 .or.    &   !idchl
    !    tracer_id(tr_num) == 1025 .or.    &   !idetz2n
    !    tracer_id(tr_num) == 1026 .or.    &   !idetz2c
    !    tracer_id(tr_num) == 1027 .or.    &   !idetz2si
    !    tracer_id(tr_num) == 1028 ) then      !idetz2calc

          ! sinking
          call recom_sinking(tr_num, vert_sink) !vert_sink

          ! sinking into the benthos
          call ver_sinking_recom_benthos(tr_num, str_bf) !str_bf
       
          ! update tracer fields
           nzmax=nlevel-1
           nzmin=ulevel
	   tr_arr(nzmin:nzmax,tr_num)=tr_arr(nzmin:nzmax,tr_num)+ &
                                        vert_sink(nzmin:nzmax)
           tr_arr(nzmin:nzmax,tr_num)=tr_arr(nzmin:nzmax,tr_num)+ &
                                        str_bf(nzmin:nzmax)                           
     END SELECT
  enddo
deallocate (dtr_bf, str_bf, vert_sink)

end subroutine

!---------------------------------------------------------------------

subroutine diff_ver_recom_expl(tr_num,dtr_bf)
   !use ocean_module
   use REcoM_GloVar
   !use recom_config

    implicit none
    
    integer, intent(in)                  		:: tr_num
    real(kind=8),dimension(:), intent(inout)		:: dtr_bf
    integer						:: nzmin,nzmax
    integer                  				:: nz,id
    integer                  				:: nlevels_minimum
    real(kind=8)            				:: bottom_flux
    real(kind=8)                                        :: area, areasvol
    real(kind=8), dimension(:), allocatable            	:: vd_flux
   
    ! initialization and allocation
    allocate(vd_flux(nlevel))
    bottom_flux = 0.d0
    vd_flux = 0.d0
    id = tracer_id(tr_num)
    area = 1.d0 ! we consider unit area for our column model
    areasvol = 1.d0
    nzmax=nlevel-1
    nzmin=ulevel	

    ! compute bottom fluxes according to tracers
    SELECT CASE (id)
      CASE (1001)
         bottom_flux = GlodecayBenthos(1) ! DIN [mmolN/m^2/s]
      CASE (1002)
      	 bottom_flux = GlodecayBenthos(2) + GlodecayBenthos(4) ! DIC + calcification
      CASE (1003)
         bottom_flux = GlodecayBenthos(4) * 2.d0 ! Alk
      CASE (1018)
      	 bottom_flux = GlodecayBenthos(3) ! Si
      CASE (1019)
      	 if(use_Fe2N) then 
	    bottom_flux = GlodecayBenthos(1) * Fe2N_benthos ! Fe
      	 else
	    bottom_flux = GlodecayBenthos(2) * Fe2C_benthos
	 end if
      CASE (1022)
      	 bottom_flux = -GlodecayBenthos(2) * redO2C ! Oxy

      END SELECT

      !_______________________________________________________________________
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  this part needs to be rechecked (we consider unit areas)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
      ! Bottom flux
      vd_flux(nzmax + 1) = bottom_flux
      !_______________________________________________________________________
      ! writing flux into rhs
      do nz=nzmin,nzmax
         ! flux contribute only the cell through its bottom !!!
         dtr_bf(nz) = dtr_bf(nz) + vd_flux(nz+1)*dt/areasvol/hnode(nz)
     end do
     
     deallocate (vd_flux)
     
end subroutine diff_ver_recom_expl

!---------------------------------------------------------------------
subroutine recom_sinking(tr_num, vert_sink)

  implicit none

  Integer, intent(in)		    	:: tr_num
  real(kind=8),dimension(:), intent(inout)	:: vert_sink
  Integer 			    	:: id, nz, nzmin,nzmax
  Real(kind=8)                      	:: Vsink, dt_sink, area, areasvol
  Real(kind=8), dimension(:), allocatable :: dz_trr, Wvel_flux, wflux, vd_flux 
  
!< Constant sinking velocities (we prescribe under namelist recom)
!< This hardcoded part is temporary 
!< .OG. 07.07.2021
    area = 1.d0 ! computation are done per unit area
    areasvol = 1.d0
    
    id = tracer_id(tr_num)
    SELECT CASE (id)
    	  CASE (1001, 1008, 1017, 1021) ! ids: DetN, DetC, DetSi, DetCalc    
               Vsink = VDet 
          CASE (1004, 1005, 1006, 1020) ! ids: PhyN, PhyC, PhyChl, PhyCalc
       	       Vsink = VPhy
	  CASE (1013:1016)              ! ids: DiaN, DiaC, DiaChl, DiaSi
	       Vsink = VDia
	  CASE (1025:1026)
	       Vsink = VDet_zoo2
    END SELECT

    if (Vsink .gt. 0.1) then
	nzmin = ulevel
  	nzmax = nlevel-1
  
        allocate (dz_trr(nzmax+1),Wvel_flux(nzmax+1), wflux(nzmax+1), vd_flux(nzmax+1))
  	! compute vertical discretization step (between center of each cell
  	dz_trr                = 0.0d0
  	dz_trr(nzmin+1:nzmax) = abs(Z(nzmin:nzmax-1) - Z(nzmin+1:nzmax))
  	dz_trr(nzmin)         = hnode(nzmin)/2.0d0
	dz_trr(nzmax+1)       = hnode(nzmax)/2.0d0
        ! initialization
        Wvel_flux(nzmin:nzmax+1)= 0.d0  ! Vertical velocity for BCG tracers ( it can be variable)
                                       
        do nz=nzmin,nzmax+1

            if (allow_var_sinking .and. id<1025) then 
                Wvel_flux(nz) = -((Vdet_a * abs(zbar(nz))/SecondsPerDay) + Vsink/SecondsPerDay)
            elseif(.not.allow_var_sinking .and. id<1025) then
                Wvel_flux(nz) = -Vsink/SecondsPerDay
            elseif (id>=1025) then ! 
	     ! We assume constant sinking for second detritus (Detz2n, Detz2c, Detz2si. Detz2Calc)    
               Wvel_flux(nz) = -VSink/SecondsPerDay 
            endif
!            endif
        end do

        wflux = 0.d0	
        dt_sink = dt
        vd_flux = 0.d0

	! compute flux with the appropriate numerical scheme
	!!!!!! check again
        if (Sinking_scheme.eq.1) then 
		call Sinking_upwind_scheme(nzmin, nzmax, tr_num, tr_arr, Wvel_flux, area, vd_flux)
        elseif (Sinking_scheme.eq.3) then
        	call Sinking_DST_3rd_order_scheme(nzmin, nzmax, tr_num, tr_arr, dt_sink, dz_trr, Wvel_flux, area, vd_flux) 
        endif

        ! estimate sinking
        do nz=nzmin,nzmax
   		vert_sink(nz) = vert_sink(nz) + (vd_flux(nz)-vd_flux(nz+1))	&
   			                      * dt_sink/areasvol	        &  
			                      / (zbar(nz)-zbar(nz+1))
	end do

	deallocate(dz_trr, Wvel_flux, wflux, vd_flux)
     endif 
     
end subroutine
!-----------------------------------     
subroutine Sinking_DST_3rd_order_scheme(nzm, nzp, trnum, trarr, &
					dtsink, dz, Wvelflux, area, vdflux)       
   
   implicit none

   ! global variables   
   integer, intent(in) 			    :: nzm, nzp, trnum	! nzmin, nzmax, tr_num
   real(kind=8), dimension(:,:), intent(in) :: trarr
   real(kind=8), dimension(:), intent(in)   :: dz, Wvelflux
   real(kind=8), intent(in)		    :: dtsink, area
   real(kind=8), dimension(:), intent(inout) :: vdflux 
   ! local variables
   integer			:: nz
   Real(kind=8)			:: Rjp,Rj,Rjm
   Real(kind=8)			:: cfl, wPs, wM, d0, d1
   Real(kind=8)			:: onesixth
   Real(kind=8)			:: thetaP, thetaM, psiP, psiM
   Real(kind=8)			:: tv

   ! initialisation
   onesixth = 1.d0/6.d0
   vdflux(nzm:nzp+1)= 0.d0
   do nz=nzp, nzm+1,-1
	Rjp = trarr(nz,trnum) - trarr(min(nz+1,nzp),trnum)
        Rj  = trarr(max(nzm,nz-1),trnum) - trarr(nz,trnum) 
        Rjm = trarr(max(nzm,nz-2),trnum) - tr_arr(max(nzm,nz-1),trnum)
        cfl = abs(Wvelflux(nz) * dtsink / dz(nz)) ! [m/day] * [day] * [1/m]
  
	wPs = Wvelflux(nz) + abs(Wvelflux(nz)) ! --> Positive vertical velocity
        wM  = Wvelflux(nz) - abs(Wvelflux(nz)) ! --> Negative vertical velocity
	d0 = (2.d0 - cfl)*(1.d0 - cfl)*onesixth
	d1 = (1.d0 - cfl*cfl)*onesixth
        thetaP = Rjm/(1.d-20+Rj)
        psiP = d0 + d1*thetaP
        psiP = max(0.d0, min(min(1.d0,psiP), &
               (1.d0-cfl)/(1.d-20+cfl)*thetaP))
        thetaM = Rjp/(1.d-20 + Rj)	
        psiM = d0 + d1*thetaM
        psiM = max(0.d0, min(min(1.d0,psiM), &
               (1.d0-cfl)/(1.d-20-cfl)*thetaM))
        tv = (0.5 * wPs * (trarr(nz,trnum) + psiM * Rj)		&
	    + 0.5 * wM  * (trarr(max(nzm,nz-1),trnum) 	&
	    + psiP * Rj))
        vdflux(nz) = - tv * area
     end do

end subroutine
!-----------------------------------     

subroutine Sinking_upwind_scheme(nzm, nzp, trnum, trarr, Wvelflux, area, vdflux)       
   
   implicit none

   ! global variables   
   integer, intent(in) 			    :: nzm, nzp, trnum	! nzmin, nzmax, tr_num
   real(kind=8), dimension(:,:), intent(in) :: trarr
   real(kind=8), dimension(:), intent(in)   :: Wvelflux
   real(kind=8), intent(in)		    :: area
   real(kind=8), dimension(:), intent(inout):: vdflux 
   
   ! local variables
   integer			:: nz
   real(kind=8)			:: tv, wPs, wM
   
   ! initialisation
   vdflux(nzm:nzp+1)= 0.d0
   
   !!!! to be recheck here
   do nz=nzm+1,nzp-1
        tv = - 0.5 * (trarr(nz-1,trnum) * (Wvelflux(nz) - abs(Wvelflux(nz))) &
             + trarr(nz,trnum) * (Wvelflux(nz) + abs(Wvelflux(nz))))
        vdflux(nz)= tv * area
    end do
   ! flux down are equal to zero (cf. original code in FESOM)
end subroutine 

!---------------------------------------------------------------------

subroutine ver_sinking_recom_benthos(tr_num, str_bf)
   
    use REcoM_GloVar
    
    implicit none

    integer, intent (in)			:: tr_num
    real(kind=8), dimension(:),intent(inout)    :: str_bf
    integer					:: nzmin, nzmax, id, nz
    real(kind=8)				:: add_benthos, Vsink, area, areasvol
    real(kind=8), dimension(:), allocatable 	:: aux, Vben

    ! initialization and allocation 
    nzmin = ulevel
    nzmax = nlevel-1
    
    allocate (aux(nlevel-1), Vben(nlevel))
    aux  = 0.d0
    Vben = 0.d0
    add_benthos = 0.d0
    
    id = tracer_id(tr_num)
    
    area = 1.d0 ! unit area
    areasvol = 1.d0
    
! 1) Calculate sinking velociy for vertical sinking case
! ****************************************************** 
    SELECT CASE (id)
   	CASE(1007:1008,1017,1021)	! ids: DetN, DetC, DetSi, DetCalc
		Vsink = VDet
	CASE(1004:1006,1020)		! ids: PhyN, PhyC, PhyChl, PhyCalc 		
		Vsink = VPhy
	CASE(1013:1016)			! ids: DiaN. DiaC, DiaChl, DiaCalc
		Vsink = VDia
	CASE(1025:1028)			! ids: Detz2n, Detz2c, Detz2si, Detz2Calc 
		Vsink = VDet_zoo2
     END SELECT
     
     if (allow_var_sinking .and. id<25) then 
        ! allow for depth-dependent sinking velocity
        Vben = Vdet_a * abs(zbar(:)) + VSink
     else
        ! constant vertical sinking (second detritus class for sure)
        Vben = Vsink
     endif
        
     ! conversion [m/d] --> [m/s] (vertical velocity, note that it is positive here)
     Vben= Vben/SecondsPerDay 
     
     !!!! this part need to be checked again
     ! compute addition to benthos (last layer, the rest is 0)
     aux(nzmax)= - tr_arr(nzmax,tr_num) * Vben(nzmax) * area
        
     do nz=nzmin,nzmax
     	str_bf(nz) = str_bf(nz) + (aux(nz))* dt / areasvol / (zbar(nz) - zbar(nz+1))
     	!!!! pb here aux(nz+1) out of bounds aux(nz) instead
        add_benthos = add_benthos - (aux(nz)) * dt
     end do

     ! update benthos add tracer to benthos
     SELECT CASE (id)
        CASE(1004, 1007, 1013, 1025) ! Nitrogen, ids: PhyN, DetN, DiaN, Detz2n 
     		Benthos(1)= Benthos(1) +  add_benthos ![mmol]
        CASE(1005, 1008, 1014, 1026) ! Carbon, ids: PhyC, DetC, DiaC, Detz2c         
                Benthos(2)= Benthos(2) + add_benthos
	CASE(1016:1017,1027)	     ! Silicon, ids: DiaSi, DetSi, Detz2si
		Benthos(3)= Benthos(3) + add_benthos
	CASE(1020:1021, 1028)	     ! Calcite, ids: PhyCalc, DetCalc, Detz2Calc
		Benthos(4)= Benthos(4) + add_benthos 
     END SELECT
     
     deallocate(aux, Vben)
     
end subroutine ver_sinking_recom_benthos

end module
