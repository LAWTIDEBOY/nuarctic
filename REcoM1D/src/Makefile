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
	ocean_ice_modules.F90        \
	recom_setup.F90          \
	recom_init.F90	         \
					REcoM1D_main.F90 

## Liste des objets a partir des sources (changement .f90 en .o)               
OBJETS=$(SOURCES:.f90=.o)

## Librairies
LIB_NETCDF_PATH=/usr/include
LIBS=-L$(LIB_NETCDF_PATH) -lnetcdf -lnetcdff
#LIBS = -L/usr/lib -llapack -lblas
#LIBS = /usr/lib64/atlas/liblapack.so.3 /usr/lib64/libblas.so.3
#LIBS = -L/usr/lib -llapack.so.3 -lblas.so.3
## Regle 0 de make, tout compiler 
all: $(EXECUTABLE)
## Regle 1 de make, liens et creation de l'executable
$(EXECUTABLE): $(OBJETS)
	$(COMPILATEUR) $(OPTIONS) $^ -o $@ $(LIBS)
 
      
## Regle 2 de make, compilation des objets separement
%.o: %.f90
	$(COMPILATEUR) -c $(OPTIONS) $^ -I$(LIB_NETCDF_PATH)
               

clean:
	@rm -f *.o *.mod *~ $(EXECUTABLE)
	@echo " "
	@echo "Nettoyage OK."
	@echo "-------------"
