base=process_crystal.exe
compiler=ifort
MKROOT=/opt/intel/mkl/10.0.3.020

crystal_nbo:  crystal_nbo.f90
	$(compiler)  shared.f90  bloch_overlap.f90 crystal_nbo.f90  -o $(base) -I$(MKLROOT)/lib/em64t -L$(MKLROOT)/lib/em64t -lmkl_lapack95 -lmkl_solver_lp64_sequential -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread -lmkl_blas95 -openmp

