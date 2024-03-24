#!/bin/bash

#module load intel/2018.5.274 > /dev/null 2>&1
module load intel/2020.2.254 > /dev/null 2>&1

FRTFLAGS="-convert big_endian -assume byterecl -align array64byte,all -traceback -fpic -qopenmp -qopenmp-threadprivate compat -fp-model source -ftz -diag-disable=remark,cpu-dispatch -qopt-report=5 -qopt-report-phase=vec"
OPT_FRTFLAGS="-g -O2 -march=core-avx2 -finline-functions -finline-limit=500 -Winline -qopt-prefetch=4 -fast-transcendentals -fimf-use-svml -no-fma"

set -x
ifort $FRTFLAGS $OPT_FRTFLAGS lattex.F90 lattex_dnt.F90 sc2prg.F90

ulimit -s unlimited

./a.out 0
./a.out 1
