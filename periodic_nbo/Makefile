# Makefile created by mkmf $Id: mkmf,v 13.0 2006/04/10 21:20:01 fms Exp $ 

include Makefile.template


.DEFAULT:
	-touch $@
all: nbo.exe
matutil.o: ./matutil.f90
	$(FC) $(FFLAGS) -c	./matutil.f90
periodic_matutil.o: ./periodic_matutil.f90 matutil.o
	$(FC) $(FFLAGS) -c	./periodic_matutil.f90
nbo.o: ./nbo.f90 matutil.o nbo_shared.o 
	$(FC) $(FFLAGS) -c	./nbo.f90
nao.o: ./nao.f90 matutil.o nbo_shared.o sortutil.o
	$(FC) $(FFLAGS) -c	./nao.f90
nbo_shared.o: ./nbo_shared.f90 matutil.o
	$(FC) $(FFLAGS) -c	./nbo_shared.f90
pre_nao.o: ./pre_nao.f90 matutil.o periodic_matutil.o nbo_shared.o sortutil.o
	$(FC) $(FFLAGS) -c	./pre_nao.f90
sortutil.o: ./sortutil.f90
	$(FC) $(FFLAGS) -c	./sortutil.f90
diag.o: ./diag.f
	$(FC) $(FFLAGS) -c	./diag.f
visual.o: ./visual.f90 nbo_shared.o
	$(FC) $(FFLAGS) -c	./visual.f90
nbo_main.o: ./nbo_main.f90 periodic_matutil.o nbo_shared.o pre_nao.o nao.o nbo.o visual.o
	$(FC) $(FFLAGS) -c	./nbo_main.f90
SRC = ./matutil.f90 ./periodic_matutil.f90 ./nbo.f90 ./nao.f90 ./nbo_shared.f90 ./pre_nao.f90 ./sortutil.f90 ./diag.f ./visual.f90 ./nbo_main.f90
OBJ = matutil.o periodic_matutil.o nbo.o nao.o nbo_shared.o pre_nao.o sortutil.o diag.o visual.o nbo_main.o
clean: neat
	-rm -f .cppdefs $(OBJ) nbo.exe
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
nbo.exe: $(OBJ)
	$(LD) $(OBJ) -o nbo.exe $(LDFLAGS)
