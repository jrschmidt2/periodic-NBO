# Makefile created by mkmf $Id: mkmf,v 13.0 2006/04/10 21:20:01 fms Exp $ 

MKLROOT=/home/schmidt/intel/Compiler/11.1/069/mkl/lib/em64t


LD=ifort
FC=ifort
FFLAGS= -I$(MKLROOT) -qopenmp
LDFLAGS=-L$(MKLROOT) -qopenmp -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lpthread -lm -I$(MKLROOT) -mkl=sequential



.DEFAULT:
	-touch $@
all: a.out
rd_basis.o: ./rd_basis.f90 rd_wavefunction.o projection_shared.o
	$(FC) $(FFLAGS) -c	./rd_basis.f90
projection_main.o: ./projection_main.f90 rd_wavefunction.o projection_shared.o rd_basis.o bloch_overlap.o PAW.o
	$(FC) $(FFLAGS) -c	./projection_main.f90
PAW.o: ./PAW.f90 projection_shared.o rd_wavefunction.o
	$(FC) $(FFLAGS) -c	./PAW.f90
rd_wavefunction.o: ./rd_wavefunction.f90
	$(FC) $(FFLAGS) -c	./rd_wavefunction.f90
projection_shared.o: ./projection_shared.f90 rd_wavefunction.o
	$(FC) $(FFLAGS) -c	./projection_shared.f90
bloch_overlap.o: ./bloch_overlap.f90 rd_wavefunction.o projection_shared.o
	$(FC) $(FFLAGS) -c	./bloch_overlap.f90
SRC = ./rd_basis.f90 ./projection_main.f90 ./PAW.f90 ./rd_wavefunction.f90 ./projection_shared.f90 ./bloch_overlap.f90
OBJ = rd_basis.o projection_main.o PAW.o rd_wavefunction.o projection_shared.o bloch_overlap.o
clean: neat
	-rm -f .cppdefs $(OBJ) projection.exe
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
a.out: $(OBJ)
	$(LD) $(OBJ) -o projection.exe $(LDFLAGS)
