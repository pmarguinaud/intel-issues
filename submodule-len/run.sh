#!/bin/bash

module purge                   
module load intel/oneapi/2024.1
module load compiler-rt/2024.1.0
module load ifort/2024.1.0
module load gcc/9.2.0          

ifx -c a.F90
ifx -c sub.F90

