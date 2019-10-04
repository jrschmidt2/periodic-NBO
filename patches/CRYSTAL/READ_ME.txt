******************
Interface for CRYSTAL09 with Periodic NBO code of JRS and BDD

5/28/14

******************

The purpose of this code is to take the output of a CRYSTAL09 calculation and produce an input file for the periodic NBO code of BDD and JRS.
Since, the density is already represented by an atom-centered basis in the CRYSTAL program, the main purpose of this code is simply to collect the necessary data from various CRYSTAL output files, and convert it to the appropriate format.

******************
Running Crystal Software
******************

No modification are necessary to the CRYSTAL source code to obtain all the necessary data, but you will have to run the associated program 'properties' to get all the data in a readable format.  There are a few requirements for running the codes to produce all the required output.

1. Running CRYSTAL
        -Gamma centered k-point meshed must be used (odd number of k-points in each direction), but Brillouin Zone symmetry can be employed.
        -No printing statements must be made
2. Running PROPERTIES
        -Must be run after the CRYSTAL program has completed.
        -This will be used to generate the ${system name}_${system name}_dat.KRED file, which contains the coefficients of each basis function in all the bands.  This is used to generate a density matrix and Fock matrix, for use in NBO calculations.  
        -The input file is a *.d3 file.  An example of which can be found in the sample/ folder.  The parameters in that sample file are general to any system, except for the second line.  The second line should contain the same k-point mesh description used in the *.d12 file used for the CRYSTAL calculation.  

******************
Running Interface
******************

There are then two codes that must be run for the creation of the NBO.out file: read_crystal.py and process_crystal.exe.  The second is from FORTRAN90 source code that must be compiled, for which a basic makefile is included.  An example shell script (crys_2_nbo.sh) to run these is included in the sample/ directory. 

1. read_crystal.py
        -This should be run first.  It simply finds relevant information in the CRYSTAL *.out file and prints out the relevant information to the screen.  Thus, the output should be forwarded to a file to be read by the process_crystal.exe program.
        -The name of the *.out file must be specified when the program is called.

2. process_crystal.exe
        -This should be run second.
        -In addition to taking in the formatted information provided by read_crystal.py, this will also read in the formatted information for the *_*_dat.KRED file produced by the properties program.  
        -This also performs some additional data manipulation to get all the information needed for an NBO calculation.  
        -Upon calling this code, the output of read_crytal.py and the *.KRED file must be given by name in order.  
        -This will output a NBO.out file that can then be read by the periodic NBO executable.

******************
Sample Systems
******************

A sample bulk fcc Cu system is included, with both the crystal and properties input files.  The basis set included is a customized bulk fcc Cu basis set.  

