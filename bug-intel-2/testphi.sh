#!/bin/bash

#module load intel/2018.5.274
#module load intel/2019.5.281
module load intel/2020.0.166

ulimit -s unlimited
export OMP_STACKSIZE=1G

set -x

\rm *.x *.out

echo -e "\n### COMPILER ENVIRONMENT ###"
icpc --version
g++ --version

echo -e "\n### CODE COMPILE ###"
icpc -std=c++11 -qopenmp testphi.cc -o testphi.icpc.openmp.x 
g++  -std=c++11 -fopenmp testphi.cc -o testphi.g++.openmp.x 
icpc -std=c++11          testphi.cc -o testphi.icpc.x
g++  -std=c++11          testphi.cc -o testphi.g++.x

echo -e "\n### RUN ###"
for prog in testphi.*.x
do
  OMP_NUM_THREADS=1 ./$prog > $prog.out
done

echo -e "\n### COMPARE RESULTS ###" 
for prog in testphi.*.x
do
  echo -e "\n==> $prog <=="
  diff testphi.g++.x.out $prog.out
done



# adding a simple printf line in the code change the execution behavior and permit to retrieve the expected result (comparison order is completely modified)

#echo -e "\n### TEST WITH THE ADDITION OF A SINGLE PRINTF ###"
#
#icpc -std=c++11 -qopenmp testphi_printf_OK.cc -o testphi_printf_OK.x 
#icpc -std=c++11 -qopenmp testphi_printf_KO.cc -o testphi_printf_KO.x 
#
#OMP_NUM_THREADS=1 ./testphi_printf_OK.x > testphi_printf_OK.out
#OMP_NUM_THREADS=1 ./testphi_printf_KO.x > testphi_printf_KO.out
#
#diff -y testphi_printf_OK.out testphi_printf_KO.out

