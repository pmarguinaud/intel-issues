#!/bin/bash

module () 
{ 
    eval `/usr/bin/modulecmd bash $*`
}

INTELONEAPI="intel/oneapi/2023.2"
MPIPACKAGE="mpi/2021.10.0"
MKL="mkl/2023.2.0"
COMPILER="compiler/2023.2.0"
module purge                   2>/dev/null
module load $INTELONEAPI       2>/dev/null
module load $COMPILER          2>/dev/null
module load $MKL               2>/dev/null
module load $MPIPACKAGE        2>/dev/null
module load gcc/9.2.0          2>/dev/null


FRTFLAGS="-convert big_endian -assume byterecl -align array64byte,all -traceback -fpic -qopenmp -qopenmp-threadprivate compat -fp-model source -ftz -diag-disable=remark,cpu-dispatch -qopt-report=5 -qopt-report-phase=vec"
OPT_FRTFLAGS="-g -O2 -march=core-avx2 -finline-functions -finline-limit=500 -Winline -qopt-prefetch=4 -fast-transcendentals -fimf-use-svml -no-fma"

for f in yomrip1.F90 ss1.F90 ss2.F90 
do
/home/gmap/mrpm/khatib/public/bin/mpiifort_wrapper $FRTFLAGS $OPT_FRTFLAGS  -c $f
done

/home/gmap/mrpm/khatib/public/bin/mpiifort_wrapper $FRTFLAGS $OPT_FRTFLAGS  -o ss.x ss.F90 *.o

./ss.x




