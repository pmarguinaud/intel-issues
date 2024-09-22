#!/bin/bash


INTELONEAPI="intel/oneapi/2023.2"
COMPILER="compiler/2023.2.0"
module purge                   
module load $INTELONEAPI       
module load $COMPILER          
module load gcc/5.3.0          


set -x

cmake ../src -DCMAKE_Fortran_COMPILER=ifort

perl -i -pe 's/XX = 123/XX = 456/; s/XTT/XTT,YTT/' ../src/b.F90
cat ../src/b.F90

make 

./main.x

perl -i -pe 's/XX = 456/XX = 123/; s/XTT,YTT/XTT/' ../src/b.F90
cat ../src/b.F90

make 

./main.x

