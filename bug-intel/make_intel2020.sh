#!/bin/bash

\rm *.o *.mod

module purge
module load intel/2020.0.166

FC=ifort
$FC --version

$FC -c fckit_log.F90
$FC -c atlas_Grid_module.F90
$FC -c atlas_module.F90
$FC -c atlas-arpege_module.F90
$FC -c atlas-arpege_f.F90

