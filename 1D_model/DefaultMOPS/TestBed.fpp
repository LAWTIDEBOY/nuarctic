! Iris Kriest (ikriest@geomar.de), 03 November 2016
! Iris Kriest (ikriest@geomar.de), 03 February 2017
!
! This program serves as a testbed for MOPS coupled to the TMM.
! It runs BGC_MODEL with initialization carried out as in BGC_INI over a specified time.
! To test any new development of MOPS, you just have to do the following things:
!
! (1) Copy or link the following source code files into the current directory:
!     BGC_MODEL.F, BGC_INI.F, BGC_PARAMS.h, BGC_CONTROL.h, BGC_BUDGET.h
!     BGC_MISFIT.h, BGC_DIAGNOSTICS.h, CAR_PARAMS.h
!     (The latter 3 are just needed for consistency)
!
! (2) Change - if required - the following parameters:
!
! (2.1) Choose vertical domain over which to run the model:
! 
! If you want to run the model over the full vertical domain, choose
! #define DEEP 1
! If you only want to simulate the upper 220 m, choose
! #define SHALLOW 1
! Default is to only simulate the upper 1000 m
!
! (2.2) Choose length of simulation:
! 
! nyears:    Number of years you want to simulate. 
!            If you only want to simulate a few days, set this to 1 and change ndays accordingly.
! 
! (2.3) Choose mixing coefficient
! 
! kv:        Mixing coefficient (m2/d)
! 
! (2.4) Set temporal and spatial setup (these are the parameters that will be passed during runtime from TMM)
! 
! bgc_timesteps: Number of timesteps for biogeochemistry (loop inside BGC_MODEL.F) 
! bgc_keuph:     Number of layers in the euphotic zone
!
! (2.5) Set physical parameters (these are the parameters that will be passed during runtime from TMM)
!
! bgc_swr:       Surface radiation (W/m2)
! bgc_tau:       Daylength (days)
! bgc_seaice:    Sea ice cover (fraction of surface)
! bgc_wind:      Wind speed (m/s)?
! bgc_atmosp:    Atmospheric pressure at sea surface (dbar)
! bgc_theta(k):  Temperature (deg C)
! bgc_salt(k):   Salinity (PSU)
!
! (2.6) If you have changed the number of biogeochemical tracers:
!
! adapt initialization statements
! adapt budgeting statements
! adapt output statements
!
! (3) Compile the model via compile script "compilemodel" or from the command line via
!
!     ifort -o testbed TestBed.fpp BGC_INI.F BGC_MODEL.F
!     or
!     gfortran -cpp  -o testbed TestBed.fpp BGC_INI.F BGC_MODEL.F
!
! (4) Start the simulation via
!    
!     ./testbed
! 
! (5) If model runs properly, it will write the following output:
!
! testbed.zax      File containing vertical grid information
! testbed.txt      File with traver concentrations at the end of each year
! testbed.log      File with budget information at the end of a simulation
!
! You can analyse the output for budget by looking at testbed.log 
! You can analyse the output for tracer fields by starting ferret and calling
! script "readtestbed":
!
! go readtestbed <time>
!
! where <time> is the number of output times (e.g., "10" for a 10 year simulation)
!

      PROGRAM TestBed

      implicit none

! communication between modules

#include "BGC_PARAMS.h"
#include "BGC_CONTROL.h"

! Needed for calculation of mass conservation of P, N, O

#include "BGC_BUDGET.h"
 
      real*8 psum0,psum,pgain,pmasscons
      real*8 nsum0,nsum,ngain,nmasscons
      real*8 osum0,osum,ogain,omasscons
      real*8 co2gain,h2ogain

! The vertical grid definition

! (2.1) select the vertical domain over which to run the model

#define DEEP 1
!#define SHALLOW 1

#ifdef DEEP
      integer*8  bgc_kloc
      parameter(bgc_kloc=15) ! the full vertical domain of the MIT grid
      real*8 mit_dz(bgc_kloc)
      data mit_dz/50.0d0,70.0d0,100.0d0,140.0d0,190.0d0,240.0d0,
     &            290.0d0,340.0d0,390.0d0,440.0d0,490.0d0,540.0d0,
     &            590.0d0,640d0,690.0d0/
#elif SHALLOW 
      integer*8  bgc_kloc
      parameter(bgc_kloc=7) ! only the upper 7 layers, down to 1080m
      real*8 mit_dz(bgc_kloc)
      data mit_dz/50.0d0,70.0d0,100.0d0,140.0d0,190.0d0,240.0d0,290.0d0/
#else
      integer*8  bgc_kloc
      parameter(bgc_kloc=7) ! only the upper 7 layers, down to 1080m
      real*8 mit_dz(bgc_kloc)
      data mit_dz/50.0d0,70.0d0,100.0d0,140.0d0,190.0d0,240.0d0,
     &            290.0d0/
#endif 

! Note: bgc_kloc is defined by the model geometry, and may vary (in the TMM) from profile to profile
!       In TMM mops_biogeochem_model.F defines this from Nrloc, and passes it to BGC_MODEL
!       bgc_ktotal is defined in BGC_PARAMS.h as 100, and served as a dummy vector length

      real*8 bgc_dz(bgc_ktotal),bgc_zu(bgc_ktotal+1)
      real*8 bgc_center(bgc_ktotal),bgc_dzc(bgc_ktotal)

! Some local variables. 
! When biogeochemistry is run in TMM, these will be provided by the main driver routine.

      integer*8 i,j,k,n,id,iy,ndays,nyears
      real*8 bgc_runoffvol,bgc_globalrunoff,bgc_swr,bgc_tau,
     &       bgc_seaice,bgc_wind,bgc_atmosp
      real*8 bgc_theta(bgc_ktotal),bgc_salt(bgc_ktotal)
      real*8 kv(bgc_kloc),flux(bgc_kloc)

#ifdef SHALLOW
! Introduce a potential lower boundary condition when simuzlated for a "shallow" domain       
      real*8 deepscale
      real*8 deep(bgc_ntracer),deepflux(bgc_ntracer)
      real*8 deepflux_p,deepflux_n,deepflux_o,sumdeepflux
#endif            

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C INITIALIZATION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! OPEN OUTPUT FILES

      open(99,file='testbed.txt',status='unknown',form='formatted')
      open(98,file='testbed.log',status='unknown',form='formatted')
      open(97,file='testbed.zax',status='unknown',form='formatted')

! INITIALIZE THE VERTICAL GRID

      bgc_kmax = bgc_kloc !bgc_kmax is used by BGC_INI.F; it gets this parameter from BGC_CONTROL.h

! initially, set all values which are define over the dummy vector length to 0
      do k=1,bgc_ktotal
        bgc_zu(k)=0.0d0
        bgc_dz(k)=0.0d0
        bgc_center(k)=0.0d0
      enddo

! thickness of each layer (= distance between upper and lower boundary)
      do k=1,bgc_kloc
         bgc_dz(k) = mit_dz(k)
      enddo

! distance between centers of each layer (distance from center of layer k to center of layer k-1)
      do k=2,bgc_kloc
         bgc_dzc(k) = (bgc_dz(k)+bgc_dz(k-1))/2.0d0
      enddo

! depth of upper boundary of each layer
      bgc_zu(1)=0.0d0    
      do k=2,bgc_kloc
         bgc_zu(k) = bgc_zu(k-1)+bgc_dz(k-1)
      enddo
      bgc_zu(bgc_kloc+1) = bgc_zu(bgc_kloc)+bgc_dz(bgc_kloc)

! center of each layer
      do k=1,bgc_kloc
         bgc_center(k) = bgc_zu(k)+bgc_dz(k)/2.0d0
      enddo

! write vertical grid definition for ferret input
      do k=1,bgc_kloc
         write(97,'(3f8.1)') bgc_zu(k),bgc_dz(k),bgc_center(k)
      enddo

! (2.2) CHOOSE LENGTH OF SIMULATION

      nyears = 3000
      ndays = 360

! (2.3) MIXING PARAMETERIZATION

      do k=1,bgc_kloc
        if(k.le.bgc_keuph) then
           kv(k) = 1.0d0 !m2/d, Brandt et al., 2015, OMZ off Mauritania had about 10 max. value; but generally 1 
        else
           kv(k) = 1.0d-1 !m2/d, Brandt et al., 2015, OMZ off Mauritania had about 10 max. value; but generally 1 
        endif
      enddo

! (2.4) INITIALIZE THE TIME STEPPING
! In TMM these values will be defined during runtime by mops_biogeochem_ini.F

      bgc_timesteps = 16 
      bgc_dt = 1.0d0/DBLE(bgc_timesteps)

! (2.5) SET PHYSICAL PARAMETERS
! In TMM these values will be passed from the main driver,
! and handed over during runtime by mops_biogeochem_model.F

      bgc_keuph = 2
      bgc_runoffvol = 0.0d0
      bgc_globalrunoff = 0.0d0   
      bgc_swr =200.0d0 
      bgc_tau = 0.5d0  
      bgc_seaice = 0.0d0 
      bgc_wind = 10.0d0 
      bgc_atmosp = 1.0d0

      do k=1,bgc_ktotal
        bgc_theta(k) = 20.0d0
        bgc_salt(k) = 37.0d0
      enddo

! INITIALIZE BIOGEOCHEMICAL PARAMETERS

      call BGC_INI(bgc_zu,1)

! (2.6) INITIALIZE THE TRACERS FIELDS

      do k=1,bgc_kloc
        do j=1,bgc_ntracer
           bgc_tracer(k,j) = 0.001d0
        enddo
        bgc_tracer(k,ipo4) = 1.00d0
        bgc_tracer(k,idin) = 30.00d0
        bgc_tracer(k,ioxy) = 100.0d0
      enddo

! (2.6) PHOSPHORUS, NITROGEN, AND OXYGEN - BUDGET (IN NUMBER OF ATOMS)
! Total inventory of P, N, O of the water column at start of simulation
! rnp is the N:P ratio in organic matter (set in BGC_INI.F)
! rop is the O:P ratio in organic matter (set in BGC_INI.F)

      psum0 = 0.0d0
      nsum0 = 0.0d0
      osum0 = 0.0d0
      do k=1,bgc_kloc
        psum0 = psum0 + bgc_dz(k)*(
     &    bgc_tracer(k,ipo4)
     &   +bgc_tracer(k,idop)
     &   +bgc_tracer(k,iphy)
     &   +bgc_tracer(k,izoo)
     &   +bgc_tracer(k,idet))
        nsum0 = nsum0 + bgc_dz(k)*(
     &    bgc_tracer(k,idin)
     &   +bgc_tracer(k,idop)*rnp
     &   +bgc_tracer(k,iphy)*rnp
     &   +bgc_tracer(k,izoo)*rnp
     &   +bgc_tracer(k,idet)*rnp)
        osum0 = osum0 + bgc_dz(k)*(
     &    bgc_tracer(k,ioxy)*2.0d0
     &   +bgc_tracer(k,ipo4)*4.0d0
     &   +bgc_tracer(k,idin)*3.0d0
     &   +bgc_tracer(k,idop)*rop
     &   +bgc_tracer(k,iphy)*rop
     &   +bgc_tracer(k,izoo)*rop
     &   +bgc_tracer(k,idet)*rop)
      enddo
      
! Variables for sources and sinks of P, N and O. 
! Variable for (implicit) water and CO2 (for accounting of oxygen) 
      pgain = 0.0d0
      ngain = 0.0d0
      ogain = 0.0d0
      co2gain = 0.0d0
      h2ogain = 0.0d0
      
#ifdef SHALLOW

! scale deep diffusion coefficient for cross-boundary flux, or
! set to 0 if you want to have a closed lower boundary 
      deepscale = 0.0d0
      
! the lower boundary conditions
      deep(ipo4) = 1.0d0
      deep(idop) = 2.0d-4
      deep(iphy) = 2.0d-4
      deep(izoo) = 2.0d-4
      deep(idet) = 2.0d-4
      deep(idin) = 30.0d0
      deep(ioxy) = 100.0d0

! account for cross-boundary fluxes when calculating the budget
      deepflux_p = 0.0d0
      deepflux_n = 0.0d0
      deepflux_o = 0.0d0
      
#endif
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C RUN THE MODEL FORWARD IN TIME
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C OUTER TIME LOOP FOR BGC

      flux_bury = 0.0d0

      do iy=1,nyears ! cycle over years 
      do id=1,ndays  ! cycle over days

        bgc_globalrunoff = flux_bury

        call BGC_MODEL(bgc_kloc,bgc_dz,     
     &            bgc_runoffvol,bgc_globalrunoff,
     &            bgc_swr,bgc_tau,bgc_seaice,bgc_wind,bgc_atmosp,
     &            bgc_theta,bgc_salt)

! (2.6) NITROGEN BUDGET (IN NUMBER OF ATOMS)
! Nitrogen fixation releases HNO3. Denitrification reduces HNO3 to N2.
! Loss (gain) of nitrogen in denitrification (nitrogen fixation) is accounted for in budget_n2.

         ngain = ngain + budget_n2          

! (2.6) OXYGEN BUDGET (IN NUMBER OF ATOMS)
! Air-sea gas exchange of oxygen creates a flux at the sea surface. This is accounted for in budget_o2. 
! For oxygen there is a loss at the sea floor via burial of POM.
! The only compensating supply is via PO4 and DIN supply (runoff) at the surface.
! Remineralization and denitrification produce CO2 and water (see r.h.s. of Eqns. 10 and 18 of Paulmier et al., 2009).
! Photosynthesis consumes CO2 and water. Nitrogen fixation consumes oxygen and water.
! The oxygen contained in water and CO2 is accounted for in budget_h2o and budget_co2.

         ogain = ogain + budget_o2*2.0d0
         ogain = ogain - flux_bury*rop
         ogain = ogain + bgc_globalrunoff*(rnp*3.0d0+4.0d0)
         co2gain = co2gain + budget_co2*2.0d0
         h2ogain = h2ogain + budget_h2o

! MIX EVERYTHING

        do j=1,bgc_ntracer
          do k=2,bgc_kloc
          flux(k) = kv(k)*(bgc_tracer(k,j)-bgc_tracer(k-1,j))/bgc_dzc(k)
          enddo              
          do k=2,bgc_kloc
            bgc_tracer(k-1,j) = bgc_tracer(k-1,j)+flux(k)/bgc_dz(k-1)
            bgc_tracer(k,j)   = bgc_tracer(k,j)  -flux(k)/bgc_dz(k)
          enddo
        enddo

#ifdef SHALLOW

! compute the flux across the lower boundary 
        do j=1,bgc_ntracer
          deepflux(j) = deepscale*kv(bgc_kloc)*
     &               (deep(j)-bgc_tracer(bgc_kloc,j))/bgc_dzc(bgc_kloc) 
          bgc_tracer(bgc_kloc,j) = bgc_tracer(bgc_kloc,j) + 
     &               deepflux(j)/bgc_dz(bgc_kloc)                    
        enddo

! sum cross boundary fluxes for budget 
        deepflux_p = deepflux_p  
     &               +deepflux(ipo4)
     &               +deepflux(idop)
     &               +deepflux(iphy)
     &               +deepflux(izoo)
     &               +deepflux(idet)
        deepflux_n = deepflux_n  
     &               +deepflux(idop)*rnp
     &               +deepflux(iphy)*rnp
     &               +deepflux(izoo)*rnp
     &               +deepflux(idet)*rnp
     &               +deepflux(idin)
        deepflux_o = deepflux_o  
     &               +deepflux(ioxy)*2.0d0
     &               +deepflux(ipo4)*4.0d0
     &               +deepflux(idin)*3.0d0
     &               +deepflux(idop)*rop
     &               +deepflux(iphy)*rop
     &               +deepflux(izoo)*rop
     &               +deepflux(idet)*rop

#endif

      enddo ! cycle over days

      write(6,*) iy

! (2.6) WRITE TRACER CONCENTRATIONS EVERY YEAR

      do k=1,bgc_kloc
         write(99,'(i8,8e13.5)') iy, 
     &                 bgc_tracer(k,ipo4),
     &                 bgc_tracer(k,idop),
     &                 bgc_tracer(k,iphy),
     &                 bgc_tracer(k,izoo),
     &                 bgc_tracer(k,idet),
     &                 bgc_tracer(k,ioxy),
     &                 bgc_tracer(k,idin)
      enddo

      enddo ! cycle over years

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C END OF SIMULATION: WRITE OUTPUT FOR MASS BUDGET
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! (2.6) PHOSPHORUS, NITROGEN, AND OXYGEN - BUDGET (IN NUMBER OF ATOMS)
! Total inventory of P, N, O of the water column at end of simulation
! rnp is the N:P ratio in organic matter (set in BGC_INI.F)
! rop is the O:P ratio in organic matter (set in BGC_INI.F)
! There is one last burial flux of P and N that is not taken care of in the above calulations:

        pgain = pgain - flux_bury
        ngain = ngain - flux_bury*rnp

#ifdef SHALLOW

! add cross boundary fluxes to overall P, N and O gain
        pgain = pgain + deepflux_p
        ngain = ngain + deepflux_n
        ogain = ogain + deepflux_o
        
#endif

        psum = 0.0d0
        nsum = 0.0d0
        osum = 0.0d0
        do k=1,bgc_kloc
          psum = psum + bgc_dz(k)*(
     &      bgc_tracer(k,ipo4)
     &     +bgc_tracer(k,idop)
     &     +bgc_tracer(k,iphy)
     &     +bgc_tracer(k,izoo)
     &     +bgc_tracer(k,idet))
          nsum = nsum + bgc_dz(k)*(
     &      bgc_tracer(k,idin)
     &     +bgc_tracer(k,idop)*rnp
     &     +bgc_tracer(k,iphy)*rnp
     &     +bgc_tracer(k,izoo)*rnp
     &     +bgc_tracer(k,idet)*rnp)
          osum = osum + bgc_dz(k)*(
     &      bgc_tracer(k,ioxy)*2.0d0
     &     +bgc_tracer(k,ipo4)*4.0d0
     &     +bgc_tracer(k,idin)*3.0d0
     &     +bgc_tracer(k,idop)*rop
     &     +bgc_tracer(k,iphy)*rop
     &     +bgc_tracer(k,izoo)*rop
     &     +bgc_tracer(k,idet)*rop)
        enddo
        
! (2.6) Mass budget including sources and sinks, expressed as % of initial mass

      pmasscons = (psum-psum0-pgain)/psum0*100.0d0
      nmasscons = (nsum-nsum0-ngain)/nsum0*100.0d0
      omasscons = (osum-osum0-ogain+co2gain+h2ogain)/osum0*100.0d0

      write(98,'(A, i8)') 'last year: ',iy-1
      write(98,'(A, i8)') 'biogeochemical time steps: ',bgc_timesteps
      write(98,'(A, i8)') 'number of layers: ',bgc_kloc
      write(98,'(A, f8.3, A)') 'surface radiation: ',bgc_swr, ' W/m2'
      write(98,'(A, f8.3, A)') 'day length: ',bgc_tau, ' days'
      write(98,'(A)') '  '

      write(98,'(A)') 'Phosphorus budget'
      write(98,'(A, g12.4)') 'total P inventory at start: ',psum0
      write(98,'(A, g12.4)') 'total P inventory at end: ',psum
      write(98,'(A, g12.4)') 'total P gain: ',pgain
      write(98,'(A, g12.4)') 'percent numerical change: ',pmasscons
      write(98,'(A)') '  '

      write(98,'(A)') 'Nitrogen budget'
      write(98,'(A, g12.4)') 'total N inventory at start: ',nsum0
      write(98,'(A, g12.4)') 'total N inventory at end: ',nsum
      write(98,'(A, g12.4)') 'total N gain: ',ngain
      write(98,'(A, g12.4)') 'percent numerical change: ',nmasscons
      write(98,'(A)') '  '

      write(98,'(A)') 'Oxygen budget'
      write(98,'(A, g12.4)') 'total O inventory at start: ',osum0
      write(98,'(A, g12.4)') 'total O inventory at end: ',osum
      write(98,'(A, g12.4)') 'total O gain: ',ogain
      write(98,'(A, g12.4)') 'O in H2O: ',h2ogain
      write(98,'(A, g12.4)') 'O in CO2: ',co2gain
      write(98,'(A, g12.4)') 'percent numerical change: ',omasscons

#ifdef SHALLOW
      write(98,'(A, g12.4)') 'P flux at lower boundary: ',deepflux_p
      write(98,'(A, g12.4)') 'N flux at lower boundary: ',deepflux_n
      write(98,'(A, g12.4)') 'O2 flux at lower boundary: ',deepflux_o
# endif


      END

