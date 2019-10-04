*************************************************
*************************************************

Periodic NBO instructions and guidelines
by BDD 2/26/15

A description of the algorithm can be found in:
Dunnington and Schmidt, JCTC; 8, 1902-1911 (2012)
*************************************************
*************************************************

This code will calculate the complete set of Natural Bond Orbitals for a periodic
system.


*************************************
*************************************
************BASIC RECIPE*************
*************************************
*************************************


1. Run an interface code to get NBO.out and NBO_mat.out

2. Command to run NBO analysis with optional arguements in brackets
    ./nbo.exe NBO.out [nbo.chk]

3. Control calculation parameters via nbo.config file. 


*************************************
*************************************
***********COMMON PROBLEMS***********
*************************************
*************************************


If the code crashes and prints "Segmentation fault." you need to increase 
OMP_STACKSIZE.

Diffuse AO-basis functions are not compatible with NPA and NBO analysis. The 
definition of 'too diffuse' is more restrictive for bulk systems due to the 
increased density and number of nearest neighbors. 

If you are planning to write a checkpoint file, the name you supply must not match
that of an existing file in the working directory. Otherwise, the code will try to 
read the already existing file as a checkpoint.


*************************************
*************************************
********Compiling Periodic NBO*******
*************************************
*************************************


The makefile contained with this file will produce an executable called 
'nbo.exe'.
The code uses algorithms from BLAS and LAPACK. We have used the MKL versions. 
Modifications to both linking and the subroutine calls will need to be modified 
if the MKL libraries are not available.


*************************************
*************************************
********Running Periodic NBO*********
*************************************
*************************************


The executable, 'nbo.exe,' runs with two inputs.
They are listed below in order they will be read in by the program, along with 
sample values.

1. System parameters file 'NBO.out'
   This is a mandatory input and should be the name of the file containing all the
   necessary information for NBO analysis. 
   This file is produced from the projection code or CRYSTAL interface. 
   The code assumes a specific format and ordering of information in the file.
   This file does not contain any matrices, which should be stored in a separate
   unformatted file. The name of the unformatted file should be included here. 

2. Checkpoint file 'nbo.chk'
   The second input is optional and is the name of a checkpoint file.
   Since performing the transformation into the NAO basis is the rate limiting step
   of the calculation, the density matrix in this basis can be stored in a
   checkpoint file so that this step need only be performed once. 
   The code will check to see if the name given matches an existing file in the 
   working directory. 
   If a file with that name does not exist, the code will run the calculation to
   obtain the NAO basis, then write a checkpoint file with the given name. 
   If a file with that name does exist, the code will try to read it as a checkpoint
   file instead of calculating the NAO basis.


*************************************
*************************************
***************Output****************
*************************************
*************************************


Most code output is sent to the screen. This includes:

  Natural Population Analysis
  Natural Atomic Orbital list for all atoms, with occupancy
  Natural Bond Orbitals
    -Lone Pairs
    -Bonds & Anti-Bonds
    -Rydberg Orbitals w/ Occupancy above 10^-3

For each NBO, the following information is given:

  Occupancy
  Hybridization
  Coefficients of Natural Hybrid Orbitals in NAO basis. The order is the same as 
    in NAO occupancy listing.

The only output file is nbo_vis.out that contains information that can be used to
generate .cube files of each NBO obtained in the analysis. 
The program capable of 



*************************************
*************************************
*********Configuration File**********
*************************************
*************************************


The nbo.config file is the primary method for controlling parameters in of the 
calculation. Upon execution of the program, it will look to see if if such a file
exists.
If the nbo.config file does not exist, default parameters will be used and a 
new nbo.config file will be produced. 
If the nbo.config file does exist, parameters will be read in. A given 
format/ordering of parameters is expected by the code, so modification of an 
automatically generated nbo.config file is the safest method of creation. 

As of now this file controls two things:
  1. Occupancy cutoffs used in the NBO search algorithm.
     The code is setup to perform an NBO search with a hard occupancy cutoff. 
     Different cutoffs can be set for lone pairs and bonds. 
     Use of a checkpoint file is encouraged if you will be varying these cutoffs.
  2. Visualization output of NBOs
     The code can output .cube files containing gridded representations of any 
     desired NBOs. 
     Default is to output nothing, and it is recommended to identify relevant NBOs 
     before creating any .cube files.
     Visualization will be explained more in its own section.

*************************************
*************************************
********Visualization Output*********
*************************************
*************************************


If visualization of the NBOs is desired, the code can generate .cube files of the 
gridded electron density, or wavefunction of each orbital. 
Due to periodic boundary conditions the volume of a given orbital need not reside
solely in the central unit cell and thus any grid must extend beyond the central 
unit cell. However an explicit representation of the central unit cell is useful 
for visualization purpose. TO address this, two types of files are created 

1. lattice_vec.cube - This is a cube file containing only atoms, where the a,b,c 
   vectors are the same as the central unit cell.  The idea is to have a set of 
   vectors for a single unit cell, but include some atoms from neighboring unit 
   cells.

2. nbo_xxx.cube - This contains an isosurface of the xxxth nbo.  The a,b,c vectors
   used here are defined by the user and may be larger than a standard unit cell
   and the 'origin' is shifted off of the unit cell's.
   The ordering of the NBO's is:
     i. Lone pairs in the order they were printed to the screen.
    ii. Paired bonds and antibonds in the order they were printed to the screen.
        Bonds are first.
   iii. Rydberg orbitals in order of increasing occupancy by atom. 
          The ones printed to the screen had the highest occupancy and are therefore
          last in a given atom's list.

In addition to whether or not to visualize, a variety of parameters are controlled
in the nbo.config file. In order they are:
  1. density - [T or F] Whether to grid out the density or the wavefunction.
  2. vis_start vis_end - [Two integers] The beginning and end of the range of NBO's
     to plot. 
  3. mesh - [Three integers] Resolution of grid point along each lattice direction.
     Note this is of the larger box.
  4. box_int - [Three integers] How many total lattice cell vector lengths to
     include along each side of the large box for orbital visualization.
  5. origin_fact - [One real] How to shift the origin of the large box with respect to the 
     central unit cell origin.

Some useful values for box size:
  -To include the central unit cell along with half of each surrounding unit cell
   for orbital visualization:
   bulk_int = 2 2 2
   origin_fact = -0.5
  -To only incloud the central unit cell for orbital visualization (useful for 
   VASP supercells):
   bulk_int = 1 1 1
   origin_fact = 0.0




