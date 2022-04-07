!> \file eos.f90
!! \BRIEF 
!> Module with EOS-80/TEOS-10 conversion subroutines
MODULE meos
CONTAINS

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

SUBROUTINE sa2sp_chem(SA, TA, DIC, NO3, SIOH4, SP)
! sa2sp_chem:             From Absolute to Practical Salinity
! 
! Description:
!      Convert from Absolute (SA) to Practical Salinity (SP), as 1D arrays,
!       based on Total Alkalinity, Dissolved Inorganic Carbon and Nitrate and
!      Silicate concentrations.
! 
! Input:
!       SA: Absolute Salinity in g/kg
!       TA: Total Alkalinity, in mol/kg
!      DIC: Dissolved Inorganic Carbone concentration in mol/kg
!      NO3: Total Nitrate concentration in mol/kg
!    SIOH4: Total Silicate concentration in mol/kg
!        N: Size of input and output 1D arrays
! 
! Details:
!      Convert from Absolute (SA) to Practical Salinity (SP) from carbon
!      system parameters and ion concentration which most affect water
!      density anomalies (nitrate and silicate).
! 
! Output:
!       SP: Practical Salinity (in psu)
! 
! Author(s):
!      Jean-Marie Epitalon
! 
! References:
!      TEOS-10 web site: http://www.teos-10.org/
! 
!      What every oceanographer needs to know about TEOS-10 (The TEOS-10
!      Primer) by Rich Pawlowicz (on TEOS-10 web site)
! 
!      R. Pawlowicz, D. G. Wright, and F. J. Millero, 2011: The
!      effects of biogeochemical processes on oceanic conductivity/
!      salinity/density relationships and the characterization of real
!      seawater
! 
!      T. J. McDougall, D. R. Jackett, F. J. Millero, R. Pawlowicz, 
!      and P. M. Barker, 2012: Algorithm for estimating
!      Absolute Salinity
! 
! See Also:
!      sp2sa_chem does the reverse
!      sa2sp_geo
! 
! Examples:
!         ! Calculate the practical salinity of a sample with Absolute Salinity of 35 g/kg,
!         ! Total Alkalinity of 0.00234 mol/kg, DIC of 0.00202 mol/kg, zero Nitrate and Silicate
!         CALL sa2sp_chem(35, 2340, 0.2020, 0, 0, 1, SP)
!      

USE msingledouble
USE gsw_mod_toolbox, only: gsw_sp_from_sa, gsw_sa_from_sp, gsw_ct_from_t, gsw_t_from_ct

IMPLICIT NONE

! Input variables:
  !>  SA: Absolute Salinity (g/kg)
  REAL(kind=rx), INTENT(in) :: SA
  !>  TA: Total Alkalinity, in mol/kg 
  REAL(kind=rx), INTENT(in) :: TA
  !> DIC: Dissolved Inorganic Carbone concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: DIC
  !> NO3: Total Nitrate concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: NO3
  !> SIOH4: Total Silicate concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: SIOH4

! Output variables:
  !>  SP: Practical Salinity on the practical salinity scale
  REAL(kind=rx), INTENT(out) :: SP

  
  ! Reverse conversion (SA from SP) follows this equation :
  !  SP = 1/r * (SA - sa_anomaly)
  !    with
  !        r = 35.16504 / 35
  !  SP = 1/r * (SA - 55.6 * (TA - x*SP) - 4.7 * (DIC - y*SP) - nutrients)
  !    with
  !       x  = 0.00230 / 35
  !       y  = 0.00208 / 35
  !       nutrients = 38.9 * NO3 + 50.7 * SIOH4
  !
  !  SP * (1 - 55.6 * x/r - 4.7 * y/r) = 1/r * (SA - 55.6 * TA - 4.7 * DIC - nutrients)
  !
  !  SP * (1 - 55.6 * 0.00230 / 35.16504 - 4.7 * 0.00208 / 35.16504) = ....
  !
  !  SP * 0.9960854303 = 1/r * (SA - 55.6 * TA - 4.7 * DIC - nutrients)
  !
  !  SP = 0.999218211671 * (SA - 55.6 * TA - 4.7 * DIC - nutrients)


   SP = 0.999218211671 * (SA - 55.6 * TA - 4.7 * DIC - 38.9 * NO3 - 50.7 * SIOH4)

  
END SUBROUTINE

SUBROUTINE sa2sp_geo(SA, SP, P, lon, lat)
! sa2sp_geo:         From Absolute to Practical Salinity
! 
! Description:
!      Convert from Absolute (SA) to Practical Salinity (SP) based on
!      depth (pressure) and geographic location.
! 
! Input:
!       SA: Absolute Salinity in g/kg
!        N: Size of input and output 1D arrays
!        P: Sea water pressure in dbar (optional, default = 0)
!      lon: Longitude, optional, in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
!      lat: Latitude, optional, in decimal degrees [-90 ... 90]
! 
! Details:
!      This subroutine is almost an alias of subroutine gsw_sp_from_sa
!      from gsw package on which it relies. The difference is in
!      that pressure and location are optional and that input/output are
!      vectors instead of scalars.
!      If location is not given or incomplete (either longitude or latitude missing),
!      an arbitrary location is chosen: the mid equatorial atlantic ocean. Note that
!      this implies an error on computed SP up to 0.02 psu
! 
! Output:
!       SP: Practical Salinity (psu)
! 
! Author(s):
!      Jean-Marie Epitalon
! 
! References:
!      TEOS-10 web site: http://www.teos-10.org/
! 
!      What every oceanographer needs to know about TEOS-10 (The TEOS-10
!      Primer) by Rich Pawlowicz (on TEOS-10 web site)
! 
!      T. J. McDougall, D. R. Jackett, F. J. Millero, R. Pawlowicz, 
!      and P. M. Barker, 2012: Algorithm for estimating
!      Absolute Salinity
! 
! See Also:
!      sp2sa_geo does the reverse
!      sa2sp_chem
! 
! Examples:
!         ! Calculate the practical salinity of a sample whose absolute Salinity is 35,
!         ! depth is 10 dbar and location is 188 degrees East and 4 degrees North.
!         CALL sa2sp_geo(35, 1, SP, 10, 188, 4)
!      

USE msingledouble
USE gsw_mod_toolbox, only: gsw_sp_from_sa, gsw_sa_from_sp, gsw_ct_from_t, gsw_t_from_ct

IMPLICIT NONE

! Input variables:

  !> SA: Absolute Salinity (g/kg)
  REAL(kind=rx), INTENT(in) :: SA
  !> P: Sea water pressure in dbar (optional, default = 0)
!f2py real(8) intent(in), optional, dimension(n) :: p = 0.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: P
  !> lon: Longitude, optional, in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
!f2py real(8) intent(in), optional, dimension(n) :: lon = -25.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: lon
  !> lat: Latitude, optional, in decimal degrees [-90 ... 90]
!f2py real(8) intent(in), optional, dimension(n) :: lat = 0.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: lat

! Output variables:
  !> SP: Practical Salinity on the practical salinity scale
  REAL(kind=rx), INTENT(out) :: SP

  ! Default latitude and longitude
  REAL(kind=r8) :: def_lat, def_lon

  IF (PRESENT(lat) .AND. PRESENT(lon)) THEN
     IF (PRESENT(P)) THEN
        SP = SGLE(gsw_sp_from_sa(DBLE(SA),DBLE(P),DBLE(lon),DBLE(lat)))
     ELSE
        SP = SGLE(gsw_sp_from_sa(DBLE(SA), 0.0d0, DBLE(lon),DBLE(lat)))
     ENDIF
  ELSE
     def_lat = 0.0
     def_lon = -25
     IF (PRESENT(P)) THEN
        SP = SGLE(gsw_sp_from_sa(DBLE(SA), DBLE(P), def_lon, def_lat))
     ELSE
        SP = SGLE(gsw_sp_from_sa(DBLE(SA), 0.0d0, def_lon, def_lat))
     ENDIF
  ENDIF

END SUBROUTINE

SUBROUTINE sp2sa_chem(SP, TA, DIC, NO3, SIOH4, SA)
! sp2sa_chem :           From Practical to Absolute Salinity
! 
! Description:
!      Convert from Practical (SP) to Absolute Salinity (SA) based on
!      Total Alkalinity, Dissolved Inorganic Carbon and Nitrate and
!      Silicate concentrations.
!      
! Input:
!       SP: Practical Salinity on the practical salinity scale
!       TA: Total Alkalinity, in mol/kg
!      DIC: Dissolved Inorganic Carbone concentration in mol/kg
!      NO3: Total Nitrate concentration in mol/kg
!    SIOH4: Total Silicate concentration in mol/kg
! 
! Details:
!      Convert from Practical (SP) to Absolute Salinity (SA) from carbon
!      system parameters and ion concentration which most affect water
!      density anomalies (nitrate and silicate).
! 
! Output:
!       SA: Absolute Salinity (g/kg)
! 
! Author(s):
!      Jean-Marie Epitalon
! 
! References:
!      TEOS-10 web site: http://www.teos-10.org/
! 
!      What every oceanographer needs to know about TEOS-10 (The TEOS-10
!      Primer) by Rich Pawlowicz (on TEOS-10 web site)
! 
!      R. Pawlowicz, D. G. Wright, and F. J. Millero, 2011: The
!      effects of biogeochemical processes on oceanic conductivity/
!      salinity/density relationships and the characterization of real
!      seawater
! 
!      T. J. McDougall, D. R. Jackett, F. J. Millero, R. Pawlowicz, 
!      and P. M. Barker, 2012: Algorithm for estimating
!      Absolute Salinity
! 
! See Also:
!      sa2sp_chem does the reverse
!      sp2sa_geo
! 
! Examples:
!         ! Calculate the absolute salinity of a sample with practical Salinity of 35,
!         ! Total Alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
!         sp2sa_chem(35, 0.00234, 0.00202, 1, SA)
!

USE gsw_mod_toolbox, only: gsw_sp_from_sa, gsw_sa_from_sp, gsw_ct_from_t, gsw_t_from_ct
USE msingledouble

IMPLICIT NONE

! Input variables:
  !>  SP: Practical Salinity on the practical salinity scale
  REAL(kind=rx), INTENT(in) :: SP
  !>  TA: Total Alkalinity, in mol/kg 
  REAL(kind=rx), INTENT(in):: TA
  !> DIC: Dissolved Inorganic Carbone concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: DIC
  !> NO3: Total Nitrate concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: NO3
  !> SIOH4: Total Silicate concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: SIOH4

! Output variables:
  !>  SA: Absolute Salinity (g/kg)
  REAL(kind=rx), INTENT(out) :: SA

  REAL(kind=r8) :: DIC_stdsw, TA_stdsw
  REAL(kind=r8) :: SR, sal_anomaly
  INTEGER :: i

  ! This function is made of two parts :
  !   1) conversion from Practical to Reference Salinity (major part)
  !   2) conversion from Reference to Absolute Salinity (Absolute Anomaly, minor part)
  !
  ! DIC and TA values of standard sea water
  TA_stdsw  = 2300.e-6 * SP /35
  DIC_stdsw = 2080.e-6 * SP /35
    
  ! Absolute salinity anomaly (g/kg)
  sal_anomaly = 55.6 * (TA - TA_stdsw) + 4.7 * (DIC - DIC_stdsw) + 38.9 * NO3 + 50.7 * SIOH4
  SR = SP * 35.16504 / 35
  SA = SR + sal_anomaly

END SUBROUTINE

SUBROUTINE sp2sa_geo(SP, SA, P, lon, lat)
! sp2sa_geo:         From Practical to Absolute Salinity
! 
! Description:
!      Convert from Practical (SP) to Absolute Salinity (SA) based on
!      depth (pressure) and geographic location.
!      
! Input:
!       SP: Practical Salinity on the practical salinity scale
!        P: Sea water pressure in dbar, optional, default = 0
!      lon: Longitude, optional, in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
!      lat: Latitude, optional, in decimal degrees [-90 ... 90]
! 
! Details:
!      This subroutine is almost an alias of subroutine gsw_sa_from_sp of
!      from gsw package on which it relies. The difference is in
!      that pressure and location are optional and that input/output are
!      vectors instead of scalars.
!      If location is not given or incomplete (either longitude or latitude missing),
!      an arbitrary location is chosen: the mid equatorial atlantic ocean. Note that
!      this implies an error on computed SA up to 0.02 g/kg
! 
! Output:
!       SA: Absolute Salinity (g/kg)
! 
! Author(s):
!      Jean-Marie Epitalon
! 
! References:
!      TEOS-10 web site: http://www.teos-10.org/
! 
!      What every oceanographer needs to know about TEOS-10 (The TEOS-10
!      Primer) by Rich Pawlowicz (on TEOS-10 web site)
! 
!      T. J. McDougall, D. R. Jackett, F. J. Millero, R. Pawlowicz, 
!      and P. M. Barker, 2012: Algorithm for estimating
!      Absolute Salinity
! 
! See Also:
!      sa2sp_geo does the reverse
!      sp2sa_chem
! 
! Examples:
!         ! Calculate the absolute salinity of a sample whose practical Salinity is 35,
!         ! depth is 10 dbar and location is 188 degrees East and 4 degrees North.
!         CALL sp2sa_geo(35, 1, SA, 10, 188, 4)

USE gsw_mod_toolbox, only: gsw_sp_from_sa, gsw_sa_from_sp, gsw_ct_from_t, gsw_t_from_ct
USE msingledouble

IMPLICIT NONE

! Input variables:
  !> SP: Practical Salinity on the practical salinity scale
  REAL(kind=rx), INTENT(in) :: SP
  !> P: Sea water pressure in dbar (optional, default = 0)
!f2py real(8) intent(in), optional, dimension(n) :: p = 0.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: P
  !> lon: Longitude, optional, in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
!f2py real(8) intent(in), optional, dimension(n) :: lon = -25.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: lon
  !> lat: Latitude, optional, in decimal degrees [-90 ... 90]
!f2py real(8) intent(in), optional, dimension(n) :: lat = 0.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: lat

! Output variables:
  !>  SA: Absolute Salinity (g/kg)
  REAL(kind=rx), INTENT(out) :: SA


  ! Default latitude and longitude
  REAL(kind=r8) :: def_lat, def_lon


  IF (PRESENT(lat) .AND. PRESENT(lon)) THEN
     IF (PRESENT(P)) THEN
           SA = SGLE(gsw_sa_from_sp(DBLE(SP), DBLE(P), DBLE(lon), DBLE(lat)))
     ELSE
           SA = SGLE(gsw_sa_from_sp(DBLE(SP), 0.0d0, DBLE(lon), DBLE(lat)))
     ENDIF
  ELSE
     def_lat = 0
     def_lon = -25
     IF (PRESENT(P)) THEN
           SA = SGLE(gsw_sa_from_sp(DBLE(SP), DBLE(P), def_lon, def_lat))
     ELSE
           SA = SGLE(gsw_sa_from_sp(DBLE(SP), 0.0d0, def_lon, def_lat))
     ENDIF
  ENDIF 
 
END SUBROUTINE

SUBROUTINE teos2eos_chem(SA, CVT, P, TA, DIC, NO3, SIOH4, T, SP)
! teos2eos_chem:      Convert temperature and salinity from TEOS-10 to EOS-80
! 
! Description:
!      Convert conservative temperature to in-situ temperature and
!      Absolute Salinity (SA) to Practical Salinity (SP). Salinity
!      conversion depends on Total Alkalinity, Dissolved Inorganic Carbon
!      and Nitrate and Silicate concentrations.
! 
! Input:
!       SA: Absolute Salinity in g/kg
!      CvT: Conservative Temperature in degrees C
!        P: Sea water pressure in dbar
!       TA: Total Alkalinity, in µmol/kg
!      DIC: Dissolved Inorganic Carbone concentration in µmol/kg
!      NO3: Total Nitrate concentration in µmol/kg
!    SIOH4: Total Silicate concentration in µmol/kg
! 
! Details:
!      Conversion from Absolute (SA) to Practical Salinity (SP) depends
!      on carbon system parameters and ion concentration which most
!      affect water density anomalies.
! 
! Output:
!        T: in situ Temperature (deg C)
!       SP: Practical Salinity (psu)
! 
! Author(s):
!      Jean-Marie Epitalon
! 
! References:
!      TEOS-10 web site: http://www.teos-10.org/
! 
!      What every oceanographer needs to know about TEOS-10 (The TEOS-10
!      Primer) by Rich Pawlowicz (on TEOS-10 web site)
! 
!      R. Pawlowicz, D. G. Wright, and F. J. Millero, 2011: The
!      effects of biogeochemical processes on oceanic conductivity/
!      salinity/density relationships and the characterization of real
!      seawater
! 
!      T. J. McDougall, D. R. Jackett, F. J. Millero, R. Pawlowicz, 
!      and P. M. Barker, 2012: Algorithm for estimating
!      Absolute Salinity
! 
! See Also:
!      eos2teos_chem does the reverse, teos2eos_geo, sa2sp_cham
!      package GSW for Fortran
! 
! Examples:
!         ! Calculate in-situ Temperature and practical salinity of a sample with 
!         ! Absolute Salinity of 35 g/kg, Conservative Temperature of 18 deg C,
!         ! at 0 dbar and Total Alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg, 
!         ! zero Nitrate and Silicate
!         CALL teos2eos_chem(35, 18, 0, 0.00234, 0.00202, 0, 0, 1, T, SP)
!      

USE gsw_mod_toolbox, only: gsw_sp_from_sa, gsw_sa_from_sp, gsw_ct_from_t, gsw_t_from_ct
USE msingledouble

IMPLICIT NONE

! Input variables:
  !> CvT: Conservative Temperature (deg C)
  REAL(kind=rx), INTENT(in) :: CvT
  !> SA: Absolute Salinity (g/kg)
  REAL(kind=rx), INTENT(in) :: SA
  !>  P: Sea water pressure in dbar
  REAL(kind=rx), INTENT(in) :: P
  !> TA: Total Alkalinity, in mol/kg 
  REAL(kind=rx), INTENT(in) :: TA
  !> DIC: Dissolved Inorganic Carbone concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: DIC
  !> NO3: Total Nitrate concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: NO3
  !> SIOH4: Total Silicate concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: SIOH4

! Output variables:
  !> SP: Practical Salinity on the practical salinity scale
  REAL(kind=rx), INTENT(out) :: SP
  !>  T: in-situ temperature in deg. C
  REAL(kind=rx), INTENT(out) :: T


  ! convert temperature
  T = SGLE(gsw_t_from_ct (DBLE(SA), DBLE(CvT), DBLE(P)))

    
  ! convert salinity
  CALL sa2sp_chem (SA, TA, DIC, NO3, SIOH4, SP)

END SUBROUTINE

SUBROUTINE teos2eos_geo(SA, CvT, P, T, SP, lon, lat)
! teos2eos_geo:       Convert temperature and salinity from TEOS-10 to EOS-80
! 
! Description:
!      Convert conservative temperature to in-situ temperature and
!      Absolute Salinity (SA) to Practical Salinity (SP). Salinity
!      conversion depends on depth and geographic location.
! 
! Input:
!       SA: Absolute Salinity in g/kg
!      CvT: Conservative Temperature in degrees C
!        P: Sea water pressure in dbar
!      lon: Longitude, optional, in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
!      lat: Latitude, optional, in decimal degrees [-90 ... 90]
! 
! Details:
!      Conversion from Absolute (SA) to Practical Salinity (SP) depends
!      on water density anomaly which is correlated with Silicate
!      concentration. This function relies on silicate concentration
!      taken from WOA (World Ocean Atlas) to evaluate density anomaly.
! 
! Output:
!        T: in situ Temperature (deg C)
!       SP: Practical Salinity (psu)
! 
! Author(s):
!      Jean-Marie Epitalon
! 
! References:
!      TEOS-10 web site: http://www.teos-10.org/
! 
!      What every oceanographer needs to know about TEOS-10 (The TEOS-10
!      Primer) by Rich Pawlowicz (on TEOS-10 web site)
! 
!      T. J. McDougall, D. R. Jackett, F. J. Millero, R. Pawlowicz, 
!      and P. M. Barker, 2012: Algorithm for estimating
!      Absolute Salinity
! 
! See Also:
!      eos2teos_geo does the reverse, teos2eos_chem, sa2sp_geo
!      package GSW for Fortran
! 
! Examples:
!         ! Calculate insitu Temperature and practical salinity of a sample with 
!         ! Absolute Salinity of 35 g/kg, Conservative Temperature of 18 deg C,
!         ! depth is 10 dbar and location is 188 degrees East and 4 degrees North.
!         CALL teos2eos_geo(35, 18, 10, 1, T, SP, 188, 4)
!      

USE gsw_mod_toolbox, only: gsw_sp_from_sa, gsw_sa_from_sp, gsw_ct_from_t, gsw_t_from_ct
USE msingledouble

IMPLICIT NONE

! Input variables:

  !> CvT: Conservative Temperature (deg C)
  REAL(kind=rx), INTENT(in) :: CvT
  !> SA: Absolute Salinity (g/kg)
  REAL(kind=rx), INTENT(in) :: SA
  !> P: Sea water pressure in dbar (optional, default = 0)
!f2py real(8) intent(in), optional, dimension(n) :: p = 0.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: P
  !> lon: Longitude, optional, in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
!f2py real(8) intent(in), optional, dimension(n) :: lon = -25.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: lon
  !> lat: Latitude, optional, in decimal degrees [-90 ... 90]
!f2py real(8) intent(in), optional, dimension(n) :: lat = 0.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: lat

! Output variables:
  !> SP: Practical Salinity on the practical salinity scale
  REAL(kind=rx), INTENT(out) :: SP
  !>  T: in-situ temperature in deg. C
  REAL(kind=rx), INTENT(out) :: T

  
  ! convert temperature

  T = SGLE(gsw_t_from_ct (DBLE(SA), DBLE(CvT), DBLE(P)))
  ! convert salinity
  IF (PRESENT(lat) .AND. PRESENT(lon)) THEN
     CALL sa2sp_geo(SA, SP, P, lon, lat)
  ELSE    
     CALL sa2sp_geo(SA, SP, P)
  ENDIF
  
END SUBROUTINE

SUBROUTINE eos2teos_geo(SP, T, P, CvT, SA, lon, lat)
! eos2teos_geo:        Convert temperature and salinity from EOS-80 to TEOS-10
! 
! Description:
!      Convert in-situ to Conservative temperature and Practical (SP) to
!      Absolute Salinity (SA), AS 1d ARRAYS.
!      Salinity conversion depends on depth and geographic location.
! 
! Input:
!       SP: Practical Salinity on the practical salinity scale
!        T: in-situ temperature in deg. C
!        P: Sea water pressure in dbar
!      lon: Longitude, optional, in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
!      lat: Latitude, optional, in decimal degrees [-90 ... 90]
!        N: Size of input and output 1D arrays
! 
! Details:
!      Conversion from Practical (SP) to Absolute Salinity (SA) depends
!      on water density anomaly which is correlated with Silicate
!      concentration. This function relies on silicate concentration
!      taken from WOA (World Ocean Atlas) to evaluate density anomaly.
! 
! Output:
!       CVT: Conservative Temperature (deg C)
!       SA: Absolute Salinity (g/kg)
! 
! Author(s):
!      Jean-Marie Epitalon
! 
! References:
!      TEOS-10 web site: http://www.teos-10.org/
! 
!      What every oceanographer needs to know about TEOS-10 (The TEOS-10
!      Primer) by Rich Pawlowicz (on TEOS-10 web site)
! 
!      T. J. McDougall, D. R. Jackett, F. J. Millero, R. Pawlowicz, 
!      and P. M. Barker, 2012: Algorithm for estimating
!      Absolute Salinity
! 
! See Also:
!      teos2eos_geo does the reverse, eos2teos_chem, sp2sa_geo
!      package GSW for Fortran
! 
! Examples:
!         ! Calculate Conservative Temperature and Absolute Salinity of a sample with 
!         ! Practical Salinity of 35 psu, in-situ Temperature of 18 deg C,
!         ! depth is 10 dbar and location is 188 degrees East and 4 degrees North.
!         call eos2teos_geo(35, 18, 10, 188, 4, 1, CvT, SA)
! 
USE gsw_mod_toolbox, only: gsw_sp_from_sa, gsw_sa_from_sp, gsw_ct_from_t, gsw_t_from_ct
USE msingledouble

IMPLICIT NONE

! Input variables:

  !> SP: Practical Salinity on the practical salinity scale
  REAL(kind=rx), INTENT(in) :: SP
  !> T: in-situ temperature in deg. C
  REAL(kind=rx), INTENT(in) :: T
  !> P: Sea water pressure in dbar (optional, default = 0)
!f2py real(8) intent(in), optional, dimension(n) :: p = 0.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: P
  !> lon: Longitude, optional, in decimal degrees [ 0 ... +360 ] or [ -180 ... +180 ]
!f2py real(8) intent(in), optional, dimension(n) :: lon = -25.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: lon
  !> lat: Latitude, optional, in decimal degrees [-90 ... 90]
!f2py real(8) intent(in), optional, dimension(n) :: lat = 0.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: lat

! Output variables:
  !> CvT: Conservative Temperature (deg C)
  REAL(kind=rx), INTENT(out) :: CvT
  !> SA: Absolute Salinity (g/kg)
  REAL(kind=rx), INTENT(out) :: SA

  INTEGER :: i

  ! convert salinity
  IF (PRESENT(lat) .AND. PRESENT(lon)) THEN
     CALL sp2sa_geo(SP, SA, P, lon, lat)
  ELSE    
     CALL sp2sa_geo(SP, SA, P)
  ENDIF
  ! convert temperature
  CvT = SGLE(gsw_ct_from_t (DBLE(SA), DBLE(T), DBLE(P)))
END SUBROUTINE

SUBROUTINE eos2teos_chem(SP, T, P, TA, DIC, NO3, SIOH4, CvT, SA)
! eos2teos_chem.Rd
! Convert temperature and salinity from EOS-80 to TEOS-10
! 
! Description:
!      Convert in-situ to Conservative temperature and Practical (SP) 
!      to Absolute Salinity (SA), as 1D arrays
!      Salinity conversion depends on Total Alkalinity, Dissolved Inorganic Carbon 
!      and Nitrate and Silicate concentrations.
! 
! Usage:
!      call eos2teos_chem(SP, T, P, TA, DIC, NO3, SIOH4, N, CvT, SA)
!      
! Input: all 1D arrays of length N
!       SP: Practical Salinity on the practical salinity scale
!        T: in-situ temperature in deg. C
!        P: Sea water pressure in dbar
!       TA: Total Alkalinity, in mol/kg 
!      DIC: Dissolved Inorganic Carbone concentration in mol/kg
!      NO3: Total Nitrate concentration in mol/kg
!    SIOH4: Total Silicate concentration in mol/kg
!        N: Size of input and output 1D arrays
! 
! Details:
!      Conversion from Practical (SP) to Absolute Salinity (SA) depends
!      on carbon system parameters and ion concentration which most
!      affect water density anomalies.
! 
! Output: all 1D arrays of length N
!       CvT: Conservative Temperature (deg C)
!       SA: Absolute Salinity (g/kg)
! 
! Author(s):
!      Jean-Marie Epitalon
! 
! References:
!      TEOS-10 web site: http://www.teos-10.org/
! 
!      What every oceanographer needs to know about TEOS-10 (The TEOS-10
!      Primer) by Rich Pawlowicz (on TEOS-10 web site)
! 
!      R. Pawlowicz, D. G. Wright, and F. J. Millero, 2011: The
!      effects of biogeochemical processes on oceanic conductivity/
!      salinity/density relationships and the characterization of real
!      seawater
! 
!      T. J. McDougall, D. R. Jackett, F. J. Millero, R. Pawlowicz, 
!      and P. M. Barker, 2012: Algorithm for estimating
!      Absolute Salinity
! 
! See Also:
!      teos2eos_chem does the reverse, eos2teos_geo, sp2sa_chem
!      package GSW for Fortran
! 
! Examples:
!         ! Calculate Conservative Temperature and Absolute Salinity of a sample with 
!         ! Practical Salinity of 35 psu, in-situ Temperature of 18 deg C,
!         ! at 0 dbar and Total Alkalinity of 0.00234 mol/kg and DIC of 0.00202 mol/kg
!         ! zero Nitrate and Silicate
!         call eos2teos_chem(35, 18, 0, 0.00234, 0.00202, 0, 0, 1, CvT, SA)
!      

USE gsw_mod_toolbox, only: gsw_sp_from_sa, gsw_sa_from_sp, gsw_ct_from_t, gsw_t_from_ct
USE msingledouble

IMPLICIT NONE

! Input variables:
  !>  SP: Practical Salinity on the practical salinity scale
  REAL(kind=rx), INTENT(in):: SP
  !>   T: in-situ temperature in deg. C
  REAL(kind=rx), INTENT(in) :: T
  !>   P: Sea water pressure in dbar
  REAL(kind=rx), INTENT(in) :: P
  !>  TA: Total Alkalinity, in mol/kg 
  REAL(kind=rx), INTENT(in) :: TA
  !> DIC: Dissolved Inorganic Carbone concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: DIC
  !> NO3: Total Nitrate concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: NO3
  !> SIOH4: Total Silicate concentration in mol/kg
  REAL(kind=rx), INTENT(in) :: SIOH4

! Output variables:
  !> CvT: Conservative Temperature (deg C)
  REAL(kind=rx), INTENT(out) :: CvT
  !> SA: Absolute Salinity (g/kg)
  REAL(kind=rx), INTENT(out) :: SA


  ! convert salinity
  CALL sp2sa_chem (SP, TA, DIC, NO3, SIOH4, SA)
  ! convert temperature
  CvT = SGLE(gsw_ct_from_t (DBLE(SA), DBLE(T), DBLE(P)))
  
END SUBROUTINE

END MODULE meos
