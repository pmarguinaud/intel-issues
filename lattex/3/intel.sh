#!/bin/bash

module load intel > /dev/null 2>&1

\rm a.out ; ifort   -qopenmp  test.F90 ; OMP_NUM_THREADS=4 ./a.out

