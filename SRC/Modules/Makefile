# Makefile for SRC/Modules/
#

F90  =  gfortran

#OPT  = -O0 -fdefault-double-8 -fdefault-real-8 -fdefault-integer-8
OPT  = -O3  -fdefault-double-8 -fdefault-real-8 -fdefault-integer-8


MODULEDIR = .
F90MOD = -I

F90FLAGS = $(F90MOD) $(MODULEDIR)

OBJ := $(patsubst %.f90,%.o,$(wildcard *.f90) )

all : $(OBJ)

clean :
	rm -f *.o *.mod *~

%.o: %.f90
	$(F90) $(F90FLAGS) $(OPT) -c $<


include Makefile.depend
#
