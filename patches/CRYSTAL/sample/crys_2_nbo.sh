#!/bin/bash
#PBS -l nodes=1:ppn=1,mem=1800mb


export OMP_NUM_THREADS=1
export OMP_STACKSIZE=800m




NAME=Cu

      #Parse the relevant information from teh crystal output into a readable format
      ./read_crystal.py $NAME.out  > cry2nbo.out

      #Then execute fortran processing to generate NBO.out file
      ./process_crystal.exe cry2nbo.out ${NAME}_${NAME}_dat.KRED





