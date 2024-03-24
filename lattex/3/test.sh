#!/bin/bash

#module load intel > /dev/null 2>&1
module load intel/oneAPI_2021.04 2>&1
module load nvhpc > /dev/null 2>&1
module load gcc/9.2.0 > /dev/null 2>&1

echo "GCC"
\rm a.out ; gfortran -fopenmp test.F90 ; OMP_NUM_THREADS=4 ./a.out

echo "Intel"
\rm a.out ; ifort   -qopenmp  test.F90 ; OMP_NUM_THREADS=4 ./a.out

echo "PGI"
\rm a.out ; pgf90   -mp       test.F90 ; OMP_NUM_THREADS=4 ./a.out


