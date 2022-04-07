############################################
## Makefile pour la compilation de REcoM1D #
## 2021/12/01 | F.Birrien | Creation       #
############################################
## make ou make all : compilation 
## make clean       : nettoyage des objets, modules et executables


## Compilateur 
COMPILATEUR=gfortran 

## Options de compilation
#OPTIONS= -O2 -IPF-flt-eval-method0
#OPTIONS= -O2 -g -traceback -fp-stack-check
#OPTIONS= -O2 -g -traceback -check bounds
#OPTIONS= -O2 -ffree-line-length-none
#OPTIONS= -O2 -openmp -parallel -fpp
OPTIONS= -O0 -g -fbacktrace -ffree-line-length-none
#OPTIONS= -O0 -fbounds-check -ffree-line-length-none
#OPTIONS= -fbacktrace 
#OPTIONS= -ffpe-trap=zero
## Nom de l'executable final
EXECUTABLE=recom1D.exe

## Liste des sources
SOURCES=recom_modules.F90        \
	recom_clock.F90          \
	forcing_modules.F90      \
	recom_setup.F90          \
	recom_init.F90	         \
	recom_extra.F90          \
	recom_sms.F90		 \
	recom_forcing.F90        \
	recom_main.F90           \
					REcoM1D_main.F90 
	
## Liste des objets a partir des sources (changement .f90 en .o)               
OBJETS=$(SOURCES:.f90=.o)

## Librairies
LIB_NETCDF_PATH=/usr/include
LIB_RECOM_PATH=/home/fbirrien/NuArctic/nuarctic/REcoM1D/src
#RLIBS=-L$(LIB_RECOM_PATH) -lrecomlib.so
NLIBS=-L$(LIB_NETCDF_PATH) -lnetcdf -lnetcdff

## Regle 0 de make, tout compiler 
all: $(EXECUTABLE)
## Regle 1 de make, liens et creation de l'executable
$(EXECUTABLE): $(OBJETS)
	$(COMPILATEUR) $(OPTIONS) $^ -o $@ $(NLIBS) 
#$(RLIBS)
 
      
## Regle 2 de make, compilation des objets separement
%.o: %.f90
	$(COMPILATEUR) -c $(OPTIONS) $^ -I$(LIB_NETCDF_PATH) 
#$(RLIBS)
               

clean:
	@rm -f *.o *.mod *~ $(EXECUTABLE)
	@echo " "
	@echo "Nettoyage OK."
	@echo "-------------"
