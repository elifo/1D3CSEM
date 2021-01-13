# Makefile for SRC/
#
EXEC = ../SEM1D3C

# For gfortran compiler
FC   = gfortran

#OPT  = -O0  -fdefault-double-8 -fdefault-real-8 -fdefault-integer-8
OPT  = -O3  -fdefault-double-8 -fdefault-real-8 -fdefault-integer-8

MMOD= -I
MODULEDIR =./Modules/

F90FLAGS =$(MMOD) $(MODULEDIR)

#
OBJ = $(patsubst %.f90,%.o, $(wildcard *.f90) )
#
OBJ_MOD = $(patsubst %.f90,%.o,$(wildcard $(MODULEDIR)*.f90) )


all :  $(EXEC)

$(EXEC) : make_module $(OBJ) 
	$(FC) -o $@ $(OBJ) $(OBJ_MOD) 

make_module : 
	$(MAKE) all -C $(MODULEDIR) "OPT=$(OPT)" "F90MOD=$(MMOD)"

%.o : %.f90 
	$(FC) $(F90FLAGS) $(OPT) -c $< 
clean:
	rm -f *.o *~
#
