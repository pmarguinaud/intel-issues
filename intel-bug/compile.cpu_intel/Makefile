FC=gfortran -g -O2
FC=pgf90 -g -O2
FC=ifort -g -O2

all: wrap_cpg_dia_flux.x

cpg_dia_flu_parallel.o: cpg_dia_flu_parallel.F90 mf_phys_type_mod.o
	$(FC) -c cpg_dia_flu_parallel.F90

field_module.o: field_module.F90 
	$(FC) -c field_module.F90

mf_phys_type_mod.o: mf_phys_type_mod.F90 field_module.o
	$(FC) -c mf_phys_type_mod.F90

wrap_cpg_dia_flux.o: wrap_cpg_dia_flux.F90 mf_phys_type_mod.o cpg_dia_flu_parallel.intfb.h
	$(FC) -c wrap_cpg_dia_flux.F90

wrap_cpg_dia_flux.x: wrap_cpg_dia_flux.o cpg_dia_flu_parallel.o field_module.o mf_phys_type_mod.o
	$(FC) -o wrap_cpg_dia_flux.x wrap_cpg_dia_flux.o cpg_dia_flu_parallel.o field_module.o mf_phys_type_mod.o 


clean: 
	\rm -f *.o *.x *.mod *.optrpt 

