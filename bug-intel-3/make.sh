#!/bin/bash

source /opt/softs/intel/2020.02/compilers_and_libraries_2020.2.254/linux/bin/compilervars.sh intel64

ifort -c tt_owned_mod.F90 
ifort -qopenmp -O2 -c toto.F90 

