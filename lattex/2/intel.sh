#!/bin/bash

ulimit -s unlimited
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=4

for v in intel/2018.5.274 intel/2019.5.281 intel/2020.2.254 intel/oneAPI_2021.04                 
do
  module load $v 2> /dev/null
  echo "$v"
  set -x
  type ifort
  ifort -qopenmp -g -traceback master.F90
  ./a.out 1
  ./a.out 0
  set +x
  \rm -f a.out
  module unload $v 2> /dev/null
done

