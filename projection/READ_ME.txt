***********************************************************
***********************************************************

Projection Code Instructions and guide lines
By BDD 2/25/15

A description of the projection algorithm can be found in:
Dunnington and Schmidt, JCTC; 8, 1902-1911 (2012)
***********************************************************
***********************************************************

This projection code is capable of representing any plane-wave bands output by VASP 
into any user specified Gaussian type orbital basis set.
The code is set up to then output all information necessary for periodic 
NBO analysis using the code of JRS and BDD.


*************************************
*************************************
************BASIC RECIPE*************
*************************************
*************************************


1. Insert projection_output.F into directory containing VASP source code and apply 
proj_VASP.patch

2. Run VASP calculation with LNBO=.TRUE. in INCAR (creates wavefunction.dat)

3. Command to run the projection with optional arguments in brackets
    ./projection.exe [basis.inp] [wavefunction.dat] [NBO.out]
    NOTE: basis.inp must obey Gaussian94 format

*************************************
*************************************
***********COMMON PROBLEMS***********
*************************************
*************************************


If the code crashes and prints "Segmentation fault." you need to increase 
OMP_STACKSIZE.

Make sure the order of atom types in the basis.inp files matches that of the 
POTCAR file used for the VAsP calculation. 

Using diffuse AO-basis functions will cause problems due to the creation of linear 
dependencies. These functions are not useful for bulk systems and are also 
incompatible with NBO analysis. A good rule of thumb is no exponents much lower 
than 0.1 bohr^-2. 

Balance is the most important characteristic when determining the quality of a 
basis set for projecting into. For instance, using a quadruple-zeta basis for atoms
in an adosrbate with a double-zeta basis for surface atoms will give artificially
polarized results. 
This means you should be careful when using basis sets for different atom types 
from different sources. 


*************************************
*************************************
**********VASP Modification**********
*************************************
*************************************


VASP must be modified to output the necessary information that can then be
read by the projection code. There are two aspects of this modification process:
1. Inclusion of module, 'projection_output.F'
   This module is included in the tar and needs to be copied into the directory
   containing the VASP source code. 
   For VASP 4 and VASP 5.3 you will also need to modify the makefile. In the
   list labelled 'SOURCE' you will need to include 'projection_output.o'at the end
   of the list.
   No modification of the makefile is required for VASP 5.4.
2. Apply proj_VASP.patch
   You will need to apply the appropriate patch based on the version of VASP you are
   using. The patch is called 'proj_VASP_{version #}.patch'. It can be applied
   to the VASP source code using the command 
   "patch < proj_VASP_{version #}.patch"
   This will modify main.F to call our output subroutine as well as create a 
   logical variable that controls printing of the output.  


*************************************
*************************************
************Running VASP*************
*************************************
*************************************


For VASP to print the necessary output, the line LNBO=.TRUE. must be included in 
the INCAR. This will produce an unformatted file called wavefunction.dat that can
be read by the projection algorithm. 
There are a few limiations on other parameters of the VASP run as well:

1. ISYM=0 Turning off k-point symmetry is not required for the projection algorithm,
   as the process is performed uniquely at each k-point. However, if you are
   interested in performing a real-space analysis, i.e. NBO, k-point symmetry will 
   need to be turned off in the VASP calculatio as our implementation of inverse
   Fourier Transforms is not compatible with symmetry.  
   NOTE: ISYM=-1 is not recommended.

2. Gamma-point containing k-point mesh. To take advantage of the inversion symmetry
   of the Brillouin zone, the first k-point is assumed to be the gamma point.

3. NPAR must be left at default values. 

4. NO ULTRASOFT PSEUDOPOTENTIALS.  There is no functional form associated with this 
   pseudopotential and it is therefore impossible to project a correction for the 
   valence electrons' interaction with the core electrons.
   This means the valence bands will each reprensent a (different) incorrect number
   of electrons.
   Norm conserving and PAW type pseudopotentials are compatible with the code.
 
The algorithm should be compatible with all other VASP options, including spin 
polarized calculations.


*************************************
*************************************
******Compiling the Projection*******
*************************************
*************************************


The makefile contained with this file will produce an executable called 
'projection.exe'.
The code uses algorithms from BLAS and LAPACK. We have used the MKL versions. 
Modifications to both linking and the subroutine calls will need to be made if
MKL libraries are not available.


*************************************
*************************************
*******Running the Projection********
*************************************
*************************************


This executable then runs with three (optional) inputs.
They are listed below in order they will be read in by the program, along with the 
default file name the program will look for in the working directory.

1. Basis set file. 'basis.inp'
   This should contain all information on the atomic orbital basis set to be used, 
   in the format of Gaussian94. This means any header can be used, but each atom 
   type must begin and end with a line of '****'.
   Basis set orbital types can include everything up to f-type, including 
   'sp'-labeling.
   Order of the atom types must match that of the POTCAR used in the VASP run. 

2. VASP output file. 'wavefunction.dat'
   This file is output from our customized version of VASP as wavefunction.dat and
   contains all necessary information about the plane-wave output for use in the 
   projection.

3. NBO input file. 'NBO.out'
   This file contains information, in a readable format, necessary for the periodic
   NBO code of JRS and BDD.
   The header of the file contains information on the basis set as well as atoms of 
   the system.
   It also contains the k-points utilized in the VASP calculation (where the 
   projection has been performed) as well as the indices of the real space unit 
   cells that will be used in the NBO computation. This set of unit cells is only 
   those which are nearset neighbors to the central unit cell.
   The actual matrices are stored in 'NBO_mat.out'. For each k-point there is an 
   overlap, density and fock matrix (two density and fock for spin polarized 
   calculations). This file is unformatted for memory considerations as well as
   efficient input/output. The file name of the matrix storage is placed at the end
   of 'NBO.out' and read from there by the NBO code, so it will have to be changed
   if the file name is changed.

A few steps in the algorithm parallelized using Open MP.
Some of the default environmental variables related to this are not optimal and 
should be changed:

1. OMP_NUM_THREADS controls how many threads are parallelized over. If this is not
   set, the code will try to use as many as possible. This should be set to how
   many processors are actually available. Not a functionality concern, just 
   courteous.

2. OMP_STACKSIZE controls the amount of memory given to each thread of
   parallelization. The default is way too low. If this is too low, the code will
   crash and simply display 'Segmentaion fault.' I have found that setting this to
   800mb is sufficient for even surface supercells. 


*************************************
*************************************
************Output Files*************
*************************************
*************************************


The projection code outputs a couple of additional files.

1.  'spillover.out' Contains quantitative information on the quality of the
     projection. SPILLOVER is defined as in Eq 10 of the periodic NBO paper. For 
     norm-conserving pseudopotentials, this value is rigorously bound by 0 and 1.
     However, for PAW-type pseudopotentials, approximations are made in calculating
     atomic orbital-augmenter overlap and this constraint is lifted. As a result the
     spillover sum includes both positive and negative terms. To account for this we
     have defined the SPREAD parameter, which is defined the same as spillover, 
     except using the absolute value for all terms within the summation.
     Both spillover and spread are calculated across all bands, as well as only 
     those that are occupied, and therefore contribute to the density matrix.
     Additionally, NO norm should be above 1.01 and a tally of occupied bands who 
     breach this limit is included in this file.

     Additionally, atomic weighted spread and spillover are provided based on 
     occupied bands. These are calculated similarly to Eq 10, except there is now a
     weighting (and normalization factor) calculated as the sum of the squares of
     the coefficients in the projected band of all basis functions centered on a 
     particular atom.

     Finally is a listing of the occupied band that the AO basis does the worst job
     in representing, gauged by the spread of that band, for each k-point. These 
     give an upward bound on the error incurred using a projection


2.  'band_spillover.out' This contains the norm of each projected band, for each 
     spin at each k-point. This just gives a more detailed picture of the 
     information summarized in the spillover.out file.

The program also outputs information to the screen.
This mainly consists of system information from the VASP wavefunction file.
Additionally after the projection the density and fock matrices (in the projected 
matrices) are checked by calculating the associated observable, number of electrons
and sum of orbital eigenvalues respectively. These are also calculated from the VASP
input at the beginning of the calculation.


*************************************
*************************************
*********Sample Calculations*********
*************************************
*************************************


Included are three sample systems to run; bulk silicon, magnesium oxide, and nickel.

First I have included a script to run the projection program. 
Note that two environmental parameters related to Open_MP paralllelization have been
set. Both should be set when running the code, especially OMP_STACKSIZE which 
controls the memory for each thread included. The default is low and even if only 
one processor is used, the code can crash here.
If the code crashes and just prints 'Segmentation fault.' this is the problem. 
I have found that export OMP_STACKSIZE=800mb is sufficient for most systems. 

For each system, I have included 3 of the 4 VASP input files. 
There are also AO basis sets included for each system.

Si: Standard 3-21G(d) basis set

MgO: 6-311G(d) basis set with outermost valence functions removed.  
     Linear dependencies arise very easily in periodic systems and any basis 
     function with an exponent below 0.1 bohr^-2 is probably going to cause problem.
     While the projection code is enabled to deal with these linear dependencies,
     NBO is not. In general, diffuse AO-functions should always be avoided for 
     periodic systems. 

Ni: Modified 6-31G basis set. Details of the modification process can be found 
    in the supporting information of: Journal of Catalysis, 324, 50-58, (2015).
    The overall strategy is to modulate the sp-functions' exponent to limit linear
    dependencies, then modify the d-function's exponent to account for the bulk
    environment. For main group or molecular solids, gas-phase basis sets can often
    be applied (perhaps requiring trimming). However, bulk transition metals are 
    not represented as well with gas-phase basis sets and requie more involved basis
    set development.

One final note about basis sets:
It should be noted that even though the VASP results only explicitly include valence
electrons, the core basis functions remain in all sets.  These functions help in 
capturing the PAW effect's on the orbitals.  The core region of the valence orbitals
is oscillatory and can in general not be well represented by only the smoothly 
varying valence like AO-basis functions.
