!$Id: mod_assimilation.F90 1444 2013-10-04 10:54:08Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_constants

! !DESCRIPTION:
        
!
! !USES:

  IMPLICIT NONE
  SAVE
!EOP
  ! define the constants
  INTEGER, PARAMETER        :: rdk = kind ( 1.0D+00 )                 !kind of double precision real number
  REAL(kind=rdk), PARAMETER :: pi = 4.0D+00*DATAN(1.0D+00)
  REAL(kind=rdk), PARAMETER :: deg2rad = pi/1.80D+02                  ! Degrees to radians conversion factor
  REAL(kind=rdk), PARAMETER :: rad2deg = 1.0D+00/deg2rad              ! Radians to degrees conversion factor
  REAL(kind=rdk), PARAMETER :: rEarth = 6.371D+06                     ! Radious of Earth in meter 
  REAL(kind=rdk), PARAMETER :: temp_C2K = 2.7315D+02
  REAL(kind=rdk), PARAMETER :: sec_per_day = 2.4D+01*6.0D+01*6.0D+01  
  REAL(kind=rdk), PARAMETER :: sec_per_year = 3.6525D+02*sec_per_day
  REAL(kind=rdk), PARAMETER :: mmol2mg_C = 1.201D+01                  ! molar mass of carbon (C) is approximately 12.01 grams per
  REAL(kind=rdk), PARAMETER :: mg2mmol_C = 1.0D+00/mmol2mg_C
  REAL(kind=rdk), PARAMETER :: mmol2mg_Chl = 8.9349D+02               ! The molar mass of Chlorophyll-a (C55H72O5N4Mg) is approximately 893.49 grams per mole (g/mol)
  REAL(kind=rdk), PARAMETER :: mg2mmol_Chl = 1.0D+00/mmol2mg_Chl  
  REAL(kind=rdk), PARAMETER :: mmol2mg_N = 2.802D+01                  ! molar mass of nitrogen (N₂) is approximately 28.02 grams per mole (g/mol).
  REAL(kind=rdk), PARAMETER :: mg2mmol_N = 1.0D+00/mmol2mg_N
  REAL(kind=rdk), PARAMETER :: universal_gas_constant = 8.314D+00     ! Universal gas constant R 
  REAL(kind=rdk), PARAMETER :: atomic_mass_unit = 1.66054D-27         ! Atomic mass unit u
  REAL(kind=rdk), PARAMETER :: molar_mass_o2 = 3.2D+01
  REAL(kind=rdk), PARAMETER :: redfield_ratio = 1.06D+02/1.6D+01      ! The stoichiometric ratio of carbon and nitrogen
  REAL(kind=rdk), PARAMETER :: cubic_m2liter = 1.0D+03
  REAL(kind=rdk), PARAMETER :: liter2cubic_m = 1.0D+01/cubic_m2liter

END MODULE mod_constants
