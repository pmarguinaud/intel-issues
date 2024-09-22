#!/bin/bash

module purge                       2>/dev/null
module load intel/oneapi/2023.2    2>/dev/null
module load compiler/2023.2.0      2>/dev/null
module load gcc/9.2.0              2>/dev/null
module load nvhpc                  2>/dev/null



for cc in gcc icx icc pgcc
do

for opt in 0 2
do

echo "==> $cc -O$opt <=="

for n in 32 64
do

$cc -O$opt -g maxval.c

./a.out $n

done

done

done

for cxx in g++ icpx icpc pgc++ 
do

for opt in 0 2
do

echo "==> $cxx -O$opt <=="

for n in 32 64
do

$cxx -O$opt -g maxval.cc

./a.out $n

done

done

done

