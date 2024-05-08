include Makefile.inc

all: apl_arpege.o

./apl_arpege.intfb.ok: $(SRC)/./apl_arpege.intfb.h geometry_mod.o mf_phys_type_mod.o cpg_type_mod.o cpg_sl2_type_mod.o cpg_opts_type_mod.o mf_phys_surface_type_mod.o field_variables_mod.o mf_phys_base_state_type_mod.o mf_phys_next_state_type_mod.o type_model.o parkind1.o yomhook.o ddh_mix.o yoe_phys_mwave.o yomcfu.o
	touch ./apl_arpege.intfb.ok

ecrad_config.ok: $(SRC)/ecrad_config.h
	touch ecrad_config.ok
netcdf3_file.ok: $(SRC)/netcdf3_file.h
	touch netcdf3_file.ok
netcdf_attributes.ok: $(SRC)/netcdf_attributes.h
	touch netcdf_attributes.ok
netcdf_constants.ok: $(SRC)/netcdf_constants.h
	touch netcdf_constants.ok
netcdf_dims.ok: $(SRC)/netcdf_dims.h
	touch netcdf_dims.ok
netcdf_eightbyte.ok: $(SRC)/netcdf_eightbyte.h
	touch netcdf_eightbyte.ok
netcdf_expanded.ok: $(SRC)/netcdf_expanded.h
	touch netcdf_expanded.ok
netcdf_externals.ok: $(SRC)/netcdf_externals.h
	touch netcdf_externals.ok
netcdf_file.ok: $(SRC)/netcdf_file.h
	touch netcdf_file.ok
netcdf_overloads.ok: $(SRC)/netcdf_overloads.h
	touch netcdf_overloads.ok
netcdf_text_variables.ok: $(SRC)/netcdf_text_variables.h
	touch netcdf_text_variables.ok
netcdf_variables.ok: $(SRC)/netcdf_variables.h
	touch netcdf_variables.ok
netcdf_visibility.ok: $(SRC)/netcdf_visibility.h
	touch netcdf_visibility.ok
radiation_delta_eddington.ok: $(SRC)/radiation_delta_eddington.h
	touch radiation_delta_eddington.ok
radiation_optical_depth_scaling.ok: $(SRC)/radiation_optical_depth_scaling.h
	touch radiation_optical_depth_scaling.ok
algorithm_state_mod.o: $(SRC)/algorithm_state_mod.F90 parkind1.o
	$(FC) -c $(SRC)/algorithm_state_mod.F90

apl_arpege.o: $(SRC)/apl_arpege.F90 geometry_mod.o mf_phys_type_mod.o cpg_type_mod.o cpg_sl2_type_mod.o cpg_opts_type_mod.o mf_phys_surface_type_mod.o field_variables_mod.o mf_phys_base_state_type_mod.o mf_phys_next_state_type_mod.o type_model.o parkind1.o yomhook.o ddh_mix.o yoe_phys_mwave.o yomcfu.o apl_arpege_part0.intfb.ok accldia.intfb.ok acclph.intfb.ok acdnshf.intfb.ok acdrag.intfb.ok acdrme.intfb.ok acevadcape.intfb.ok achmt.intfb.ok achmtls.intfb.ok acpluis.intfb.ok acsol.intfb.ok actqsat.intfb.ok acvisih.intfb.ok cpqsol.intfb.ok mf_phys_bayrad.intfb.ok mf_phys_corwat.intfb.ok mf_phys_cvv.intfb.ok mf_phys_fpl_part2.intfb.ok mf_phys_mocon.intfb.ok mf_phys_precips.intfb.ok mf_phys_save_phsurf_part2.intfb.ok mf_phys_transfer.intfb.ok ppwetpoint.intfb.ok qngcor.intfb.ok aplpar_flexdia.intfb.ok apl_arpege_oceanic_fluxes.intfb.ok apl_wind_gust.intfb.ok apl_arpege_shallow_convection_and_turbulence.intfb.ok apl_arpege_albedo_computation.intfb.ok apl_arpege_aerosols_for_radiation.intfb.ok apl_arpege_cloudiness.intfb.ok apl_arpege_radiation.intfb.ok apl_arpege_soil_hydro.intfb.ok apl_arpege_surface.intfb.ok apl_arpege_deep_convection.intfb.ok apl_arpege_turbulence_part2.intfb.ok apl_arpege_precipitation.intfb.ok apl_arpege_hydro_budget.intfb.ok apl_arpege_surface_update.intfb.ok apl_arpege_atmosphere_update.intfb.ok apl_arpege_init.intfb.ok apl_arpege_init_surfex.intfb.ok apl_arpege_dprecips.intfb.ok apl_arpege_mwave.intfb.ok
	$(FC) -c $(SRC)/apl_arpege.F90

compensated_summation_mod.o: $(SRC)/compensated_summation_mod.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/compensated_summation_mod.F90

control_vectors_base_mix.o: $(SRC)/control_vectors_base_mix.F90 parkind1.o yomhook.o spectral_fields_mod.o control_vectors_data_mix.o jb_control_vectors_base_mod.o yomct0.o
	$(FC) -c $(SRC)/control_vectors_base_mix.F90

control_vectors_data_mix.o: $(SRC)/control_vectors_data_mix.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/control_vectors_data_mix.F90

cpg_opts_type_mod.o: $(SRC)/cpg_opts_type_mod.F90 parkind1.o geometry_mod.o yomcli.o surface_fields_mix.o fields_mod.o type_model.o yomct0.o yomct3.o yomsta.o yomvert.o yomdyncore.o yomlsforc.o yomlun.o yomnud.o yomrip0.o yomsnu.o yomini.o yommp0.o field_variables_mod.o abor1.intfb.ok
	$(FC) -c $(SRC)/cpg_opts_type_mod.F90

cpg_sl1_type_mod.o: $(SRC)/cpg_sl1_type_mod.F90 field_module.o field_factory_module.o field_registry_mod.o parkind1.o yom_ygfl.o type_model.o ptrslb1.o eint_mod.o yomhook.o field_access_module.o cpg_opts_type_mod.o abor1.intfb.ok
	$(FC) -c $(SRC)/cpg_sl1_type_mod.F90

cpg_sl2_type_mod.o: $(SRC)/cpg_sl2_type_mod.F90 field_module.o field_factory_module.o field_registry_mod.o parkind1.o ptrslb2.o
	$(FC) -c $(SRC)/cpg_sl2_type_mod.F90

cpg_type_mod.o: $(SRC)/cpg_type_mod.F90 field_module.o field_array_module.o field_factory_module.o field_registry_mod.o parkind1.o yom_ygfl.o intdyn_mod.o yomparar.o cpg_opts_type_mod.o yomdyna.o type_model.o ddh_mix.o
	$(FC) -c $(SRC)/cpg_type_mod.F90

cplng_types_mod.o: $(SRC)/cplng_types_mod.F90 parkind1.o parkind_ocean.o
	$(FC) -c $(SRC)/cplng_types_mod.F90

ddh_mix.o: $(SRC)/ddh_mix.F90 parkind1.o yomhook.o yomlun_ifsaux.o yomlddh.o yommddh.o yomtddh.o abor1.intfb.ok
	$(FC) -c $(SRC)/ddh_mix.F90

deallocate_if_allocated_mod.o: $(SRC)/deallocate_if_allocated_mod.F90 ec_parkind.o yomhook.o
	$(FC) -c $(SRC)/deallocate_if_allocated_mod.F90

deallocate_if_associated_mod.o: $(SRC)/deallocate_if_associated_mod.F90 ec_parkind.o yomhook.o
	$(FC) -c $(SRC)/deallocate_if_associated_mod.F90

dev_alloc_module.o: $(SRC)/dev_alloc_module.F90 parkind1.o
	$(FC) -c $(SRC)/dev_alloc_module.F90

easy_netcdf.o: $(SRC)/easy_netcdf.F90 parkind1.o radiation_io.o netcdf.o
	$(FC) -c $(SRC)/easy_netcdf.F90

easy_netcdf_read_mpi.o: $(SRC)/easy_netcdf_read_mpi.F90 easy_netcdf.o parkind1.o radiation_io.o
	$(FC) -c $(SRC)/easy_netcdf_read_mpi.F90

ec_lun.o: $(SRC)/ec_lun.F90 ec_parkind.o
	$(FC) -c $(SRC)/ec_lun.F90

ec_parkind.o: $(SRC)/ec_parkind.F90 
	$(FC) -c $(SRC)/ec_parkind.F90

ec_phys_fields_mod.o: $(SRC)/ec_phys_fields_mod.F90 yoe_tile_prop.o yoe_phys_mwave.o geometry_mod.o yomdphy.o
	$(FC) -c $(SRC)/ec_phys_fields_mod.F90

eint_mod.o: $(SRC)/eint_mod.F90 parkind1.o yomhook.o yomlun.o yommp0.o abor1.intfb.ok
	$(FC) -c $(SRC)/eint_mod.F90

fdbsubs_mod.o: $(SRC)/fdbsubs_mod.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/fdbsubs_mod.F90

field_1im_access_module.o: $(SRC)/field_1im_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_1im_access_module.F90

field_1im_array_util_module.o: $(SRC)/field_1im_array_util_module.F90 field_array_module.o parkind1.o
	$(FC) -c $(SRC)/field_1im_array_util_module.F90

field_1im_data_module.o: $(SRC)/field_1im_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_1im_data_module.F90

field_1im_factory_module.o: $(SRC)/field_1im_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_1im_factory_module.F90

field_1im_gang_module.o: $(SRC)/field_1im_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o
	$(FC) -c $(SRC)/field_1im_gang_module.F90

field_1im_gather_module.o: $(SRC)/field_1im_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o
	$(FC) -c $(SRC)/field_1im_gather_module.F90

field_1im_module.o: $(SRC)/field_1im_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_1im_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_1im_module.F90

field_1im_shuffle_module.o: $(SRC)/field_1im_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_1im_shuffle_module.F90

field_1im_util_module.o: $(SRC)/field_1im_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_1im_util_module.F90

field_1lm_access_module.o: $(SRC)/field_1lm_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_1lm_access_module.F90

field_1lm_array_util_module.o: $(SRC)/field_1lm_array_util_module.F90 field_array_module.o parkind1.o
	$(FC) -c $(SRC)/field_1lm_array_util_module.F90

field_1lm_data_module.o: $(SRC)/field_1lm_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_1lm_data_module.F90

field_1lm_factory_module.o: $(SRC)/field_1lm_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_1lm_factory_module.F90

field_1lm_gang_module.o: $(SRC)/field_1lm_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o
	$(FC) -c $(SRC)/field_1lm_gang_module.F90

field_1lm_gather_module.o: $(SRC)/field_1lm_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o
	$(FC) -c $(SRC)/field_1lm_gather_module.F90

field_1lm_module.o: $(SRC)/field_1lm_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_1lm_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_1lm_module.F90

field_1lm_shuffle_module.o: $(SRC)/field_1lm_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_1lm_shuffle_module.F90

field_1lm_util_module.o: $(SRC)/field_1lm_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_1lm_util_module.F90

field_1rb_access_module.o: $(SRC)/field_1rb_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rb_access_module.F90

field_1rb_array_util_module.o: $(SRC)/field_1rb_array_util_module.F90 field_array_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rb_array_util_module.F90

field_1rb_data_module.o: $(SRC)/field_1rb_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_1rb_data_module.F90

field_1rb_factory_module.o: $(SRC)/field_1rb_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rb_factory_module.F90

field_1rb_gang_module.o: $(SRC)/field_1rb_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o
	$(FC) -c $(SRC)/field_1rb_gang_module.F90

field_1rb_gather_module.o: $(SRC)/field_1rb_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rb_gather_module.F90

field_1rb_module.o: $(SRC)/field_1rb_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_1rb_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_1rb_module.F90

field_1rb_shuffle_module.o: $(SRC)/field_1rb_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rb_shuffle_module.F90

field_1rb_util_module.o: $(SRC)/field_1rb_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_1rb_util_module.F90

field_1rd_access_module.o: $(SRC)/field_1rd_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rd_access_module.F90

field_1rd_array_util_module.o: $(SRC)/field_1rd_array_util_module.F90 field_array_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rd_array_util_module.F90

field_1rd_data_module.o: $(SRC)/field_1rd_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_1rd_data_module.F90

field_1rd_factory_module.o: $(SRC)/field_1rd_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rd_factory_module.F90

field_1rd_gang_module.o: $(SRC)/field_1rd_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o
	$(FC) -c $(SRC)/field_1rd_gang_module.F90

field_1rd_gather_module.o: $(SRC)/field_1rd_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rd_gather_module.F90

field_1rd_module.o: $(SRC)/field_1rd_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_1rd_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_1rd_module.F90

field_1rd_shuffle_module.o: $(SRC)/field_1rd_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rd_shuffle_module.F90

field_1rd_util_module.o: $(SRC)/field_1rd_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_1rd_util_module.F90

field_1rm_access_module.o: $(SRC)/field_1rm_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rm_access_module.F90

field_1rm_array_util_module.o: $(SRC)/field_1rm_array_util_module.F90 field_array_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rm_array_util_module.F90

field_1rm_data_module.o: $(SRC)/field_1rm_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_1rm_data_module.F90

field_1rm_factory_module.o: $(SRC)/field_1rm_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rm_factory_module.F90

field_1rm_gang_module.o: $(SRC)/field_1rm_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o
	$(FC) -c $(SRC)/field_1rm_gang_module.F90

field_1rm_gather_module.o: $(SRC)/field_1rm_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rm_gather_module.F90

field_1rm_module.o: $(SRC)/field_1rm_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_1rm_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_1rm_module.F90

field_1rm_shuffle_module.o: $(SRC)/field_1rm_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_1rm_shuffle_module.F90

field_1rm_util_module.o: $(SRC)/field_1rm_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_1rm_util_module.F90

field_2im_access_module.o: $(SRC)/field_2im_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_2im_access_module.F90

field_2im_array_util_module.o: $(SRC)/field_2im_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_2im_array_util_module.F90

field_2im_data_module.o: $(SRC)/field_2im_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_2im_data_module.F90

field_2im_factory_module.o: $(SRC)/field_2im_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_2im_factory_module.F90

field_2im_gang_module.o: $(SRC)/field_2im_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o
	$(FC) -c $(SRC)/field_2im_gang_module.F90

field_2im_gather_module.o: $(SRC)/field_2im_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_2im_shuffle_module.o
	$(FC) -c $(SRC)/field_2im_gather_module.F90

field_2im_module.o: $(SRC)/field_2im_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_2im_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_2im_module.F90

field_2im_shuffle_module.o: $(SRC)/field_2im_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_2im_shuffle_module.F90

field_2im_util_module.o: $(SRC)/field_2im_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_2im_util_module.F90

field_2lm_access_module.o: $(SRC)/field_2lm_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_2lm_access_module.F90

field_2lm_array_util_module.o: $(SRC)/field_2lm_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_2lm_array_util_module.F90

field_2lm_data_module.o: $(SRC)/field_2lm_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_2lm_data_module.F90

field_2lm_factory_module.o: $(SRC)/field_2lm_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_2lm_factory_module.F90

field_2lm_gang_module.o: $(SRC)/field_2lm_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o
	$(FC) -c $(SRC)/field_2lm_gang_module.F90

field_2lm_gather_module.o: $(SRC)/field_2lm_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_2lm_shuffle_module.o
	$(FC) -c $(SRC)/field_2lm_gather_module.F90

field_2lm_module.o: $(SRC)/field_2lm_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_2lm_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_2lm_module.F90

field_2lm_shuffle_module.o: $(SRC)/field_2lm_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_2lm_shuffle_module.F90

field_2lm_util_module.o: $(SRC)/field_2lm_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_2lm_util_module.F90

field_2rb_access_module.o: $(SRC)/field_2rb_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_2rb_access_module.F90

field_2rb_array_util_module.o: $(SRC)/field_2rb_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_2rb_array_util_module.F90

field_2rb_data_module.o: $(SRC)/field_2rb_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_2rb_data_module.F90

field_2rb_factory_module.o: $(SRC)/field_2rb_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_2rb_factory_module.F90

field_2rb_gang_module.o: $(SRC)/field_2rb_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o
	$(FC) -c $(SRC)/field_2rb_gang_module.F90

field_2rb_gather_module.o: $(SRC)/field_2rb_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_2rb_shuffle_module.o
	$(FC) -c $(SRC)/field_2rb_gather_module.F90

field_2rb_module.o: $(SRC)/field_2rb_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_2rb_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_2rb_module.F90

field_2rb_shuffle_module.o: $(SRC)/field_2rb_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_2rb_shuffle_module.F90

field_2rb_util_module.o: $(SRC)/field_2rb_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_2rb_util_module.F90

field_2rd_access_module.o: $(SRC)/field_2rd_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_2rd_access_module.F90

field_2rd_array_util_module.o: $(SRC)/field_2rd_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_2rd_array_util_module.F90

field_2rd_data_module.o: $(SRC)/field_2rd_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_2rd_data_module.F90

field_2rd_factory_module.o: $(SRC)/field_2rd_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_2rd_factory_module.F90

field_2rd_gang_module.o: $(SRC)/field_2rd_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o
	$(FC) -c $(SRC)/field_2rd_gang_module.F90

field_2rd_gather_module.o: $(SRC)/field_2rd_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_2rd_shuffle_module.o
	$(FC) -c $(SRC)/field_2rd_gather_module.F90

field_2rd_module.o: $(SRC)/field_2rd_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_2rd_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_2rd_module.F90

field_2rd_shuffle_module.o: $(SRC)/field_2rd_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_2rd_shuffle_module.F90

field_2rd_util_module.o: $(SRC)/field_2rd_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_2rd_util_module.F90

field_2rm_access_module.o: $(SRC)/field_2rm_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_2rm_access_module.F90

field_2rm_array_util_module.o: $(SRC)/field_2rm_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_2rm_array_util_module.F90

field_2rm_data_module.o: $(SRC)/field_2rm_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_2rm_data_module.F90

field_2rm_factory_module.o: $(SRC)/field_2rm_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_2rm_factory_module.F90

field_2rm_gang_module.o: $(SRC)/field_2rm_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o
	$(FC) -c $(SRC)/field_2rm_gang_module.F90

field_2rm_gather_module.o: $(SRC)/field_2rm_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_2rm_shuffle_module.o
	$(FC) -c $(SRC)/field_2rm_gather_module.F90

field_2rm_module.o: $(SRC)/field_2rm_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_2rm_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_2rm_module.F90

field_2rm_shuffle_module.o: $(SRC)/field_2rm_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_2rm_shuffle_module.F90

field_2rm_util_module.o: $(SRC)/field_2rm_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_2rm_util_module.F90

field_3im_access_module.o: $(SRC)/field_3im_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_3im_access_module.F90

field_3im_array_util_module.o: $(SRC)/field_3im_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_3im_array_util_module.F90

field_3im_data_module.o: $(SRC)/field_3im_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_3im_data_module.F90

field_3im_factory_module.o: $(SRC)/field_3im_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_3im_factory_module.F90

field_3im_gang_module.o: $(SRC)/field_3im_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_3im_module.o field_2im_module.o parkind1.o
	$(FC) -c $(SRC)/field_3im_gang_module.F90

field_3im_gather_module.o: $(SRC)/field_3im_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_3im_shuffle_module.o
	$(FC) -c $(SRC)/field_3im_gather_module.F90

field_3im_module.o: $(SRC)/field_3im_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_3im_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_3im_module.F90

field_3im_shuffle_module.o: $(SRC)/field_3im_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_3im_shuffle_module.F90

field_3im_util_module.o: $(SRC)/field_3im_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_3im_util_module.F90

field_3lm_access_module.o: $(SRC)/field_3lm_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_3lm_access_module.F90

field_3lm_array_util_module.o: $(SRC)/field_3lm_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_3lm_array_util_module.F90

field_3lm_data_module.o: $(SRC)/field_3lm_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_3lm_data_module.F90

field_3lm_factory_module.o: $(SRC)/field_3lm_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_3lm_factory_module.F90

field_3lm_gang_module.o: $(SRC)/field_3lm_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_3lm_module.o field_2lm_module.o parkind1.o
	$(FC) -c $(SRC)/field_3lm_gang_module.F90

field_3lm_gather_module.o: $(SRC)/field_3lm_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_3lm_shuffle_module.o
	$(FC) -c $(SRC)/field_3lm_gather_module.F90

field_3lm_module.o: $(SRC)/field_3lm_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_3lm_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_3lm_module.F90

field_3lm_shuffle_module.o: $(SRC)/field_3lm_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_3lm_shuffle_module.F90

field_3lm_util_module.o: $(SRC)/field_3lm_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_3lm_util_module.F90

field_3rb_access_module.o: $(SRC)/field_3rb_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rb_access_module.F90

field_3rb_array_util_module.o: $(SRC)/field_3rb_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_3rb_array_util_module.F90

field_3rb_data_module.o: $(SRC)/field_3rb_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_3rb_data_module.F90

field_3rb_factory_module.o: $(SRC)/field_3rb_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rb_factory_module.F90

field_3rb_gang_module.o: $(SRC)/field_3rb_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_3rb_module.o field_2rb_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rb_gang_module.F90

field_3rb_gather_module.o: $(SRC)/field_3rb_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_3rb_shuffle_module.o
	$(FC) -c $(SRC)/field_3rb_gather_module.F90

field_3rb_module.o: $(SRC)/field_3rb_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_3rb_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_3rb_module.F90

field_3rb_shuffle_module.o: $(SRC)/field_3rb_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rb_shuffle_module.F90

field_3rb_util_module.o: $(SRC)/field_3rb_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_3rb_util_module.F90

field_3rd_access_module.o: $(SRC)/field_3rd_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rd_access_module.F90

field_3rd_array_util_module.o: $(SRC)/field_3rd_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_3rd_array_util_module.F90

field_3rd_data_module.o: $(SRC)/field_3rd_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_3rd_data_module.F90

field_3rd_factory_module.o: $(SRC)/field_3rd_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rd_factory_module.F90

field_3rd_gang_module.o: $(SRC)/field_3rd_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_3rd_module.o field_2rd_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rd_gang_module.F90

field_3rd_gather_module.o: $(SRC)/field_3rd_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_3rd_shuffle_module.o
	$(FC) -c $(SRC)/field_3rd_gather_module.F90

field_3rd_module.o: $(SRC)/field_3rd_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_3rd_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_3rd_module.F90

field_3rd_shuffle_module.o: $(SRC)/field_3rd_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rd_shuffle_module.F90

field_3rd_util_module.o: $(SRC)/field_3rd_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_3rd_util_module.F90

field_3rm_access_module.o: $(SRC)/field_3rm_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rm_access_module.F90

field_3rm_array_util_module.o: $(SRC)/field_3rm_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_3rm_array_util_module.F90

field_3rm_data_module.o: $(SRC)/field_3rm_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_3rm_data_module.F90

field_3rm_factory_module.o: $(SRC)/field_3rm_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rm_factory_module.F90

field_3rm_gang_module.o: $(SRC)/field_3rm_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_3rm_module.o field_2rm_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rm_gang_module.F90

field_3rm_gather_module.o: $(SRC)/field_3rm_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_3rm_shuffle_module.o
	$(FC) -c $(SRC)/field_3rm_gather_module.F90

field_3rm_module.o: $(SRC)/field_3rm_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_3rm_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_3rm_module.F90

field_3rm_shuffle_module.o: $(SRC)/field_3rm_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_3rm_shuffle_module.F90

field_3rm_util_module.o: $(SRC)/field_3rm_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_3rm_util_module.F90

field_4im_access_module.o: $(SRC)/field_4im_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_4im_access_module.F90

field_4im_array_util_module.o: $(SRC)/field_4im_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_4im_array_util_module.F90

field_4im_data_module.o: $(SRC)/field_4im_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_4im_data_module.F90

field_4im_factory_module.o: $(SRC)/field_4im_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_4im_factory_module.F90

field_4im_gang_module.o: $(SRC)/field_4im_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_4im_module.o field_3im_module.o parkind1.o
	$(FC) -c $(SRC)/field_4im_gang_module.F90

field_4im_gather_module.o: $(SRC)/field_4im_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_4im_shuffle_module.o
	$(FC) -c $(SRC)/field_4im_gather_module.F90

field_4im_module.o: $(SRC)/field_4im_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_4im_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_4im_module.F90

field_4im_shuffle_module.o: $(SRC)/field_4im_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_4im_shuffle_module.F90

field_4im_util_module.o: $(SRC)/field_4im_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_4im_util_module.F90

field_4lm_access_module.o: $(SRC)/field_4lm_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_4lm_access_module.F90

field_4lm_array_util_module.o: $(SRC)/field_4lm_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_4lm_array_util_module.F90

field_4lm_data_module.o: $(SRC)/field_4lm_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_4lm_data_module.F90

field_4lm_factory_module.o: $(SRC)/field_4lm_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_4lm_factory_module.F90

field_4lm_gang_module.o: $(SRC)/field_4lm_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_4lm_module.o field_3lm_module.o parkind1.o
	$(FC) -c $(SRC)/field_4lm_gang_module.F90

field_4lm_gather_module.o: $(SRC)/field_4lm_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_4lm_shuffle_module.o
	$(FC) -c $(SRC)/field_4lm_gather_module.F90

field_4lm_module.o: $(SRC)/field_4lm_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_4lm_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_4lm_module.F90

field_4lm_shuffle_module.o: $(SRC)/field_4lm_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_4lm_shuffle_module.F90

field_4lm_util_module.o: $(SRC)/field_4lm_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_4lm_util_module.F90

field_4rb_access_module.o: $(SRC)/field_4rb_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rb_access_module.F90

field_4rb_array_util_module.o: $(SRC)/field_4rb_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_4rb_array_util_module.F90

field_4rb_data_module.o: $(SRC)/field_4rb_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_4rb_data_module.F90

field_4rb_factory_module.o: $(SRC)/field_4rb_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rb_factory_module.F90

field_4rb_gang_module.o: $(SRC)/field_4rb_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_4rb_module.o field_3rb_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rb_gang_module.F90

field_4rb_gather_module.o: $(SRC)/field_4rb_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_4rb_shuffle_module.o
	$(FC) -c $(SRC)/field_4rb_gather_module.F90

field_4rb_module.o: $(SRC)/field_4rb_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_4rb_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_4rb_module.F90

field_4rb_shuffle_module.o: $(SRC)/field_4rb_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rb_shuffle_module.F90

field_4rb_util_module.o: $(SRC)/field_4rb_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_4rb_util_module.F90

field_4rd_access_module.o: $(SRC)/field_4rd_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rd_access_module.F90

field_4rd_array_util_module.o: $(SRC)/field_4rd_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_4rd_array_util_module.F90

field_4rd_data_module.o: $(SRC)/field_4rd_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_4rd_data_module.F90

field_4rd_factory_module.o: $(SRC)/field_4rd_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rd_factory_module.F90

field_4rd_gang_module.o: $(SRC)/field_4rd_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_4rd_module.o field_3rd_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rd_gang_module.F90

field_4rd_gather_module.o: $(SRC)/field_4rd_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_4rd_shuffle_module.o
	$(FC) -c $(SRC)/field_4rd_gather_module.F90

field_4rd_module.o: $(SRC)/field_4rd_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_4rd_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_4rd_module.F90

field_4rd_shuffle_module.o: $(SRC)/field_4rd_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rd_shuffle_module.F90

field_4rd_util_module.o: $(SRC)/field_4rd_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_4rd_util_module.F90

field_4rm_access_module.o: $(SRC)/field_4rm_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rm_access_module.F90

field_4rm_array_util_module.o: $(SRC)/field_4rm_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_4rm_array_util_module.F90

field_4rm_data_module.o: $(SRC)/field_4rm_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_4rm_data_module.F90

field_4rm_factory_module.o: $(SRC)/field_4rm_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rm_factory_module.F90

field_4rm_gang_module.o: $(SRC)/field_4rm_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_4rm_module.o field_3rm_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rm_gang_module.F90

field_4rm_gather_module.o: $(SRC)/field_4rm_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_4rm_shuffle_module.o
	$(FC) -c $(SRC)/field_4rm_gather_module.F90

field_4rm_module.o: $(SRC)/field_4rm_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_4rm_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_4rm_module.F90

field_4rm_shuffle_module.o: $(SRC)/field_4rm_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_4rm_shuffle_module.F90

field_4rm_util_module.o: $(SRC)/field_4rm_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_4rm_util_module.F90

field_5im_access_module.o: $(SRC)/field_5im_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_5im_access_module.F90

field_5im_array_util_module.o: $(SRC)/field_5im_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_5im_array_util_module.F90

field_5im_data_module.o: $(SRC)/field_5im_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_5im_data_module.F90

field_5im_factory_module.o: $(SRC)/field_5im_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_5im_factory_module.F90

field_5im_gang_module.o: $(SRC)/field_5im_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_5im_module.o field_4im_module.o parkind1.o
	$(FC) -c $(SRC)/field_5im_gang_module.F90

field_5im_gather_module.o: $(SRC)/field_5im_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_5im_shuffle_module.o
	$(FC) -c $(SRC)/field_5im_gather_module.F90

field_5im_module.o: $(SRC)/field_5im_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_5im_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_5im_module.F90

field_5im_shuffle_module.o: $(SRC)/field_5im_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_5im_shuffle_module.F90

field_5im_util_module.o: $(SRC)/field_5im_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_5im_util_module.F90

field_5lm_access_module.o: $(SRC)/field_5lm_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_5lm_access_module.F90

field_5lm_array_util_module.o: $(SRC)/field_5lm_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_5lm_array_util_module.F90

field_5lm_data_module.o: $(SRC)/field_5lm_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_5lm_data_module.F90

field_5lm_factory_module.o: $(SRC)/field_5lm_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_5lm_factory_module.F90

field_5lm_gang_module.o: $(SRC)/field_5lm_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_5lm_module.o field_4lm_module.o parkind1.o
	$(FC) -c $(SRC)/field_5lm_gang_module.F90

field_5lm_gather_module.o: $(SRC)/field_5lm_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_5lm_shuffle_module.o
	$(FC) -c $(SRC)/field_5lm_gather_module.F90

field_5lm_module.o: $(SRC)/field_5lm_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_5lm_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_5lm_module.F90

field_5lm_shuffle_module.o: $(SRC)/field_5lm_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_5lm_shuffle_module.F90

field_5lm_util_module.o: $(SRC)/field_5lm_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_5lm_util_module.F90

field_5rb_access_module.o: $(SRC)/field_5rb_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rb_access_module.F90

field_5rb_array_util_module.o: $(SRC)/field_5rb_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_5rb_array_util_module.F90

field_5rb_data_module.o: $(SRC)/field_5rb_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_5rb_data_module.F90

field_5rb_factory_module.o: $(SRC)/field_5rb_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rb_factory_module.F90

field_5rb_gang_module.o: $(SRC)/field_5rb_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_5rb_module.o field_4rb_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rb_gang_module.F90

field_5rb_gather_module.o: $(SRC)/field_5rb_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_5rb_shuffle_module.o
	$(FC) -c $(SRC)/field_5rb_gather_module.F90

field_5rb_module.o: $(SRC)/field_5rb_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_5rb_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_5rb_module.F90

field_5rb_shuffle_module.o: $(SRC)/field_5rb_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rb_shuffle_module.F90

field_5rb_util_module.o: $(SRC)/field_5rb_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_5rb_util_module.F90

field_5rd_access_module.o: $(SRC)/field_5rd_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rd_access_module.F90

field_5rd_array_util_module.o: $(SRC)/field_5rd_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_5rd_array_util_module.F90

field_5rd_data_module.o: $(SRC)/field_5rd_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_5rd_data_module.F90

field_5rd_factory_module.o: $(SRC)/field_5rd_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rd_factory_module.F90

field_5rd_gang_module.o: $(SRC)/field_5rd_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_5rd_module.o field_4rd_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rd_gang_module.F90

field_5rd_gather_module.o: $(SRC)/field_5rd_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_5rd_shuffle_module.o
	$(FC) -c $(SRC)/field_5rd_gather_module.F90

field_5rd_module.o: $(SRC)/field_5rd_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_5rd_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_5rd_module.F90

field_5rd_shuffle_module.o: $(SRC)/field_5rd_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rd_shuffle_module.F90

field_5rd_util_module.o: $(SRC)/field_5rd_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_5rd_util_module.F90

field_5rm_access_module.o: $(SRC)/field_5rm_access_module.F90 field_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rm_access_module.F90

field_5rm_array_util_module.o: $(SRC)/field_5rm_array_util_module.F90 field_array_module.o parkind1.o field_util_module.o
	$(FC) -c $(SRC)/field_5rm_array_util_module.F90

field_5rm_data_module.o: $(SRC)/field_5rm_data_module.F90 field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_5rm_data_module.F90

field_5rm_factory_module.o: $(SRC)/field_5rm_factory_module.F90 field_module.o field_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rm_factory_module.F90

field_5rm_gang_module.o: $(SRC)/field_5rm_gang_module.F90 field_basic_module.o field_constants_module.o dev_alloc_module.o field_5rm_module.o field_4rm_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rm_gang_module.F90

field_5rm_gather_module.o: $(SRC)/field_5rm_gather_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_5rm_shuffle_module.o
	$(FC) -c $(SRC)/field_5rm_gather_module.F90

field_5rm_module.o: $(SRC)/field_5rm_module.F90 dev_alloc_module.o host_alloc_module.o field_basic_module.o field_constants_module.o field_defaults_module.o parkind1.o field_5rm_data_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_5rm_module.F90

field_5rm_shuffle_module.o: $(SRC)/field_5rm_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_5rm_shuffle_module.F90

field_5rm_util_module.o: $(SRC)/field_5rm_util_module.F90 field_module.o field_access_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_5rm_util_module.F90

field_abort_module.o: $(SRC)/field_abort_module.F90 
	$(FC) -c $(SRC)/field_abort_module.F90

field_access_module.o: $(SRC)/field_access_module.F90 field_module.o parkind1.o field_1rm_access_module.o field_2rm_access_module.o field_3rm_access_module.o field_4rm_access_module.o field_5rm_access_module.o field_1rb_access_module.o field_2rb_access_module.o field_3rb_access_module.o field_4rb_access_module.o field_5rb_access_module.o field_1rd_access_module.o field_2rd_access_module.o field_3rd_access_module.o field_4rd_access_module.o field_5rd_access_module.o field_1im_access_module.o field_2im_access_module.o field_3im_access_module.o field_4im_access_module.o field_5im_access_module.o field_1lm_access_module.o field_2lm_access_module.o field_3lm_access_module.o field_4lm_access_module.o field_5lm_access_module.o
	$(FC) -c $(SRC)/field_access_module.F90

field_array_module.o: $(SRC)/field_array_module.F90 field_module.o field_factory_module.o parkind1.o
	$(FC) -c $(SRC)/field_array_module.F90

field_array_util_module.o: $(SRC)/field_array_util_module.F90 field_1rm_array_util_module.o field_2rm_array_util_module.o field_3rm_array_util_module.o field_4rm_array_util_module.o field_5rm_array_util_module.o field_1rb_array_util_module.o field_2rb_array_util_module.o field_3rb_array_util_module.o field_4rb_array_util_module.o field_5rb_array_util_module.o field_1rd_array_util_module.o field_2rd_array_util_module.o field_3rd_array_util_module.o field_4rd_array_util_module.o field_5rd_array_util_module.o field_1im_array_util_module.o field_2im_array_util_module.o field_3im_array_util_module.o field_4im_array_util_module.o field_5im_array_util_module.o field_1lm_array_util_module.o field_2lm_array_util_module.o field_3lm_array_util_module.o field_4lm_array_util_module.o field_5lm_array_util_module.o
	$(FC) -c $(SRC)/field_array_util_module.F90

field_basic_module.o: $(SRC)/field_basic_module.F90 dev_alloc_module.o field_constants_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_basic_module.F90

field_constants_module.o: $(SRC)/field_constants_module.F90 parkind1.o
	$(FC) -c $(SRC)/field_constants_module.F90

field_container_base_mod.o: $(SRC)/field_container_base_mod.F90 parkind1.o yomhook.o yomlun.o field_definitions_base.o field_definitions.o field_gfl_wrapper.o
	$(FC) -c $(SRC)/field_container_base_mod.F90

field_container_gp_mod.o: $(SRC)/field_container_gp_mod.F90 parkind1.o yomhook.o yomlun.o yomgfl.o yomgmv.o surface_fields_mix.o field_definitions_base.o field_definitions.o field_gfl_wrapper.o field_container_base_mod.o abor1.intfb.ok
	$(FC) -c $(SRC)/field_container_gp_mod.F90

field_container_oper_mod.o: $(SRC)/field_container_oper_mod.F90 parkind1.o yomlun.o field_container_gp_mod.o field_container_sp_mod.o geometry_mod.o field_definitions.o yomhook.o yommp0.o yom_grid_biconserv.o yomdim.o yomlap.o iostream_mix.o abor1.intfb.ok
	$(FC) -c $(SRC)/field_container_oper_mod.F90

field_container_sp_mod.o: $(SRC)/field_container_sp_mod.F90 parkind1.o yomhook.o field_container_base_mod.o field_definitions_base.o
	$(FC) -c $(SRC)/field_container_sp_mod.F90

field_defaults_module.o: $(SRC)/field_defaults_module.F90 parkind1.o
	$(FC) -c $(SRC)/field_defaults_module.F90

field_definitions.o: $(SRC)/field_definitions.F90 field_definitions_base.o parkind1.o yomhook.o yomcst.o yomlun.o yom_grib_codes.o abor1.intfb.ok
	$(FC) -c $(SRC)/field_definitions.F90

field_definitions_base.o: $(SRC)/field_definitions_base.F90 parkind1.o yomlun.o yomhook.o abor1.intfb.ok
	$(FC) -c $(SRC)/field_definitions_base.F90

field_factory_module.o: $(SRC)/field_factory_module.F90 field_1rm_factory_module.o field_2rm_factory_module.o field_3rm_factory_module.o field_4rm_factory_module.o field_5rm_factory_module.o field_1rb_factory_module.o field_2rb_factory_module.o field_3rb_factory_module.o field_4rb_factory_module.o field_5rb_factory_module.o field_1rd_factory_module.o field_2rd_factory_module.o field_3rd_factory_module.o field_4rd_factory_module.o field_5rd_factory_module.o field_1im_factory_module.o field_2im_factory_module.o field_3im_factory_module.o field_4im_factory_module.o field_5im_factory_module.o field_1lm_factory_module.o field_2lm_factory_module.o field_3lm_factory_module.o field_4lm_factory_module.o field_5lm_factory_module.o
	$(FC) -c $(SRC)/field_factory_module.F90

field_gang_module.o: $(SRC)/field_gang_module.F90 field_3rm_gang_module.o field_4rm_gang_module.o field_5rm_gang_module.o field_3rb_gang_module.o field_4rb_gang_module.o field_5rb_gang_module.o field_3rd_gang_module.o field_4rd_gang_module.o field_5rd_gang_module.o field_3im_gang_module.o field_4im_gang_module.o field_5im_gang_module.o field_3lm_gang_module.o field_4lm_gang_module.o field_5lm_gang_module.o parkind1.o
	$(FC) -c $(SRC)/field_gang_module.F90

field_gathscat_module.o: $(SRC)/field_gathscat_module.F90 field_module.o field_access_module.o field_factory_module.o field_gathscat_type_module.o parkind1.o field_1rm_gather_module.o field_2rm_gather_module.o field_3rm_gather_module.o field_4rm_gather_module.o field_5rm_gather_module.o field_1rb_gather_module.o field_2rb_gather_module.o field_3rb_gather_module.o field_4rb_gather_module.o field_5rb_gather_module.o field_1rd_gather_module.o field_2rd_gather_module.o field_3rd_gather_module.o field_4rd_gather_module.o field_5rd_gather_module.o field_1im_gather_module.o field_2im_gather_module.o field_3im_gather_module.o field_4im_gather_module.o field_5im_gather_module.o field_1lm_gather_module.o field_2lm_gather_module.o field_3lm_gather_module.o field_4lm_gather_module.o field_5lm_gather_module.o
	$(FC) -c $(SRC)/field_gathscat_module.F90

field_gathscat_type_module.o: $(SRC)/field_gathscat_type_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/field_gathscat_type_module.F90

field_gfl_wrapper.o: $(SRC)/field_gfl_wrapper.F90 parkind1.o yomhook.o yomct0.o surface_fields_mix.o yomgmv.o yomgfl.o yom_ygfl.o yomlun.o field_definitions_base.o field_definitions.o abor1.intfb.ok
	$(FC) -c $(SRC)/field_gfl_wrapper.F90

field_module.o: $(SRC)/field_module.F90 dev_alloc_module.o field_1rm_module.o field_2rm_module.o field_3rm_module.o field_4rm_module.o field_5rm_module.o field_1rb_module.o field_2rb_module.o field_3rb_module.o field_4rb_module.o field_5rb_module.o field_1rd_module.o field_2rd_module.o field_3rd_module.o field_4rd_module.o field_5rd_module.o field_1im_module.o field_2im_module.o field_3im_module.o field_4im_module.o field_5im_module.o field_1lm_module.o field_2lm_module.o field_3lm_module.o field_4lm_module.o field_5lm_module.o parkind1.o
	$(FC) -c $(SRC)/field_module.F90

field_registry_mod.o: $(SRC)/field_registry_mod.F90 parkind1.o field_module.o geometry_mod.o variable_module.o field_factory_module.o field_variables_mod.o surface_variables_mod.o yomgmv.o yom_ygfl.o yomgfl.o surface_fields_mix.o ec_phys_fields_mod.o model_physics_radiation_mod.o type_model.o
	$(FC) -c $(SRC)/field_registry_mod.F90

field_shuffle_module.o: $(SRC)/field_shuffle_module.F90 field_module.o field_access_module.o field_factory_module.o field_shuffle_type_module.o parkind1.o field_1rm_shuffle_module.o field_1rm_gather_module.o field_2rm_shuffle_module.o field_2rm_gather_module.o field_3rm_shuffle_module.o field_3rm_gather_module.o field_4rm_shuffle_module.o field_4rm_gather_module.o field_5rm_shuffle_module.o field_5rm_gather_module.o field_1rb_shuffle_module.o field_1rb_gather_module.o field_2rb_shuffle_module.o field_2rb_gather_module.o field_3rb_shuffle_module.o field_3rb_gather_module.o field_4rb_shuffle_module.o field_4rb_gather_module.o field_5rb_shuffle_module.o field_5rb_gather_module.o field_1rd_shuffle_module.o field_1rd_gather_module.o field_2rd_shuffle_module.o field_2rd_gather_module.o field_3rd_shuffle_module.o field_3rd_gather_module.o field_4rd_shuffle_module.o field_4rd_gather_module.o field_5rd_shuffle_module.o field_5rd_gather_module.o field_1im_shuffle_module.o field_1im_gather_module.o field_2im_shuffle_module.o field_2im_gather_module.o field_3im_shuffle_module.o field_3im_gather_module.o field_4im_shuffle_module.o field_4im_gather_module.o field_5im_shuffle_module.o field_5im_gather_module.o field_1lm_shuffle_module.o field_1lm_gather_module.o field_2lm_shuffle_module.o field_2lm_gather_module.o field_3lm_shuffle_module.o field_3lm_gather_module.o field_4lm_shuffle_module.o field_4lm_gather_module.o field_5lm_shuffle_module.o field_5lm_gather_module.o
	$(FC) -c $(SRC)/field_shuffle_module.F90

field_shuffle_type_module.o: $(SRC)/field_shuffle_type_module.F90 field_module.o field_access_module.o field_factory_module.o parkind1.o field_2rm_shuffle_module.o field_3rm_shuffle_module.o field_4rm_shuffle_module.o field_5rm_shuffle_module.o field_2rb_shuffle_module.o field_3rb_shuffle_module.o field_4rb_shuffle_module.o field_5rb_shuffle_module.o field_2rd_shuffle_module.o field_3rd_shuffle_module.o field_4rd_shuffle_module.o field_5rd_shuffle_module.o field_2im_shuffle_module.o field_3im_shuffle_module.o field_4im_shuffle_module.o field_5im_shuffle_module.o field_2lm_shuffle_module.o field_3lm_shuffle_module.o field_4lm_shuffle_module.o field_5lm_shuffle_module.o field_abort_module.o
	$(FC) -c $(SRC)/field_shuffle_type_module.F90

field_util_module.o: $(SRC)/field_util_module.F90 field_2rm_util_module.o field_3rm_util_module.o field_4rm_util_module.o field_5rm_util_module.o field_2rb_util_module.o field_3rb_util_module.o field_4rb_util_module.o field_5rb_util_module.o field_2rd_util_module.o field_3rd_util_module.o field_4rd_util_module.o field_5rd_util_module.o field_2im_util_module.o field_3im_util_module.o field_4im_util_module.o field_5im_util_module.o field_2lm_util_module.o field_3lm_util_module.o field_4lm_util_module.o field_5lm_util_module.o
	$(FC) -c $(SRC)/field_util_module.F90

field_variables_mod.o: $(SRC)/field_variables_mod.F90 parkind1.o variable_module.o yomgmv.o yomgfl.o
	$(FC) -c $(SRC)/field_variables_mod.F90

fields_base_mod.o: $(SRC)/fields_base_mod.F90 parkind1.o yomlun.o yomct0.o geometry_mod.o type_model.o field_container_gp_mod.o variables_mod.o order_independent_summation_mod.o yommp0.o field_definitions.o iostream_mix.o yomgrib.o grib_handles_mod.o yomhook.o field_container_oper_mod.o field_container_sp_mod.o abor1.intfb.ok
	$(FC) -c $(SRC)/fields_base_mod.F90

fields_mod.o: $(SRC)/fields_mod.F90 parkind1.o yomlun.o fields_base_mod.o geometry_mod.o variables_mod.o type_model.o yomgfl.o yomgmv.o surface_fields_mix.o ec_phys_fields_mod.o spectral_fields_mod.o field_container_gp_mod.o field_variables_mod.o surface_variables_mod.o field_registry_mod.o yomcfu.o yomxfu.o yommcuf.o fullpos.o yemlbc_fields.o yomct0.o yommp0.o yomhook.o type_fluxes.o yomct3.o field_container_oper_mod.o field_definitions.o field_definitions_base.o abor1.intfb.ok
	$(FC) -c $(SRC)/fields_mod.F90

fitspectrum_mod.o: $(SRC)/fitspectrum_mod.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/fitspectrum_mod.F90

fullpos.o: $(SRC)/fullpos.F90 yomfpcnt.o yomfpgeometry.o yomvert.o yomfpfilters.o eint_mod.o yomwfpb.o yomafn.o yomfpop.o yomfpc.o type_fposbuf.o
	$(FC) -c $(SRC)/fullpos.F90

fullpos_mix.o: $(SRC)/fullpos_mix.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/fullpos_mix.F90

geometry_mod.o: $(SRC)/geometry_mod.F90 parkind1.o type_geometry.o yomlun.o yomhook.o yomcsgeom.o yomgsgeom.o yomvert.o abor1.intfb.ok
	$(FC) -c $(SRC)/geometry_mod.F90

grib_api_interface.o: $(SRC)/grib_api_interface.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/grib_api_interface.F90

grib_handles_mod.o: $(SRC)/grib_handles_mod.F90 parkind1.o grib_api_interface.o yomhook.o yommp0.o
	$(FC) -c $(SRC)/grib_handles_mod.F90

grib_utils_mod.o: $(SRC)/grib_utils_mod.F90 parkind1.o grib_api_interface.o yom_grib_codes.o fitspectrum_mod.o algorithm_state_mod.o yomvar.o yomct0.o yomgrib.o yomhook.o abor1.intfb.ok
	$(FC) -c $(SRC)/grib_utils_mod.F90

gridpoint_fields_mix.o: $(SRC)/gridpoint_fields_mix.F90 parkind1.o yomhook.o geometry_mod.o yom_grib_codes.o indexfind_mod.o deallocate_if_associated_mod.o order_independent_summation_mod.o abor1.intfb.ok
	$(FC) -c $(SRC)/gridpoint_fields_mix.F90

host_alloc_module.o: $(SRC)/host_alloc_module.F90 parkind1.o field_abort_module.o
	$(FC) -c $(SRC)/host_alloc_module.F90

indexfind_mod.o: $(SRC)/indexfind_mod.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/indexfind_mod.F90

intdyn_mod.o: $(SRC)/intdyn_mod.F90 parkind1.o yomhook.o geometry_mod.o field_module.o field_factory_module.o
	$(FC) -c $(SRC)/intdyn_mod.F90

intdynsl_mod.o: $(SRC)/intdynsl_mod.F90 parkind1.o yomhook.o yomdyn.o yomdyna.o yom_ygfl.o
	$(FC) -c $(SRC)/intdynsl_mod.F90

iostream_mix.o: $(SRC)/iostream_mix.F90 parkind1.o pardim.o fdbsubs_mod.o yomct0.o yomhook.o grib_handles_mod.o yom_grib_codes.o grib_utils_mod.o grib_api_interface.o deallocate_if_allocated_mod.o
	$(FC) -c $(SRC)/iostream_mix.F90

jb_control_vectors_base_mod.o: $(SRC)/jb_control_vectors_base_mod.F90 parkind1.o yomlun.o yommp0.o yomhook.o geometry_mod.o spectral_fields_mod.o control_vectors_data_mix.o
	$(FC) -c $(SRC)/jb_control_vectors_base_mod.F90

mf_phys_base_state_type_mod.o: $(SRC)/mf_phys_base_state_type_mod.F90 parkind1.o cpg_type_mod.o field_module.o mf_phys_surface_type_mod.o field_variables_mod.o sc2prg_mod.o type_model.o
	$(FC) -c $(SRC)/mf_phys_base_state_type_mod.F90

mf_phys_next_state_type_mod.o: $(SRC)/mf_phys_next_state_type_mod.F90 parkind1.o field_module.o cpg_opts_type_mod.o geometry_mod.o field_variables_mod.o cpg_sl1_type_mod.o type_model.o variable_module.o
	$(FC) -c $(SRC)/mf_phys_next_state_type_mod.F90

mf_phys_surface_type_mod.o: $(SRC)/mf_phys_surface_type_mod.F90 parkind1.o surface_variables_mod.o surface_views_module.o yomdyn.o cpg_opts_type_mod.o
	$(FC) -c $(SRC)/mf_phys_surface_type_mod.F90

mf_phys_type_mod.o: $(SRC)/mf_phys_type_mod.F90 parkind1.o field_module.o field_registry_mod.o cpg_opts_type_mod.o
	$(FC) -c $(SRC)/mf_phys_type_mod.F90

mo_gas_concentrations.o: $(SRC)/mo_gas_concentrations.F90 mo_rte_kind.o mo_rte_rrtmgp_config.o mo_rrtmgp_util_string.o mo_rte_util_array.o
	$(FC) -c $(SRC)/mo_gas_concentrations.F90

mo_gas_optics_kernels.o: $(SRC)/mo_gas_optics_kernels.F90 mo_rte_kind.o mod_mlp_rrtmgp.o
	$(FC) -c $(SRC)/mo_gas_optics_kernels.F90

mo_gas_optics_rrtmgp.o: $(SRC)/mo_gas_optics_rrtmgp.F90 mo_rte_kind.o mo_rte_rrtmgp_config.o mo_rte_util_array.o mo_gas_optics_kernels.o mo_rrtmgp_constants.o mo_rrtmgp_util_string.o mo_gas_concentrations.o mo_gas_ref_concentrations.o mod_mlp_rrtmgp.o
	$(FC) -c $(SRC)/mo_gas_optics_rrtmgp.F90

mo_gas_ref_concentrations.o: $(SRC)/mo_gas_ref_concentrations.F90 mo_rte_kind.o mo_rrtmgp_util_string.o
	$(FC) -c $(SRC)/mo_gas_ref_concentrations.F90

mo_load_coefficients.o: $(SRC)/mo_load_coefficients.F90 mo_rte_kind.o mo_gas_concentrations.o mo_gas_optics_rrtmgp.o mo_simple_netcdf.o netcdf.o
	$(FC) -c $(SRC)/mo_load_coefficients.F90

mo_rrtmgp_constants.o: $(SRC)/mo_rrtmgp_constants.F90 mo_rte_kind.o
	$(FC) -c $(SRC)/mo_rrtmgp_constants.F90

mo_rrtmgp_util_string.o: $(SRC)/mo_rrtmgp_util_string.F90 
	$(FC) -c $(SRC)/mo_rrtmgp_util_string.F90

mo_rte_kind.o: $(SRC)/mo_rte_kind.F90 parkind1.o
	$(FC) -c $(SRC)/mo_rte_kind.F90

mo_rte_rrtmgp_config.o: $(SRC)/mo_rte_rrtmgp_config.F90 mo_rte_kind.o
	$(FC) -c $(SRC)/mo_rte_rrtmgp_config.F90

mo_rte_util_array.o: $(SRC)/mo_rte_util_array.F90 mo_rte_kind.o
	$(FC) -c $(SRC)/mo_rte_util_array.F90

mo_simple_netcdf.o: $(SRC)/mo_simple_netcdf.F90 mo_rte_kind.o netcdf.o
	$(FC) -c $(SRC)/mo_simple_netcdf.F90

mod_activation.o: $(SRC)/mod_activation.F90 mo_rte_kind.o
	$(FC) -c $(SRC)/mod_activation.F90

mod_layer.o: $(SRC)/mod_layer.F90 mod_activation.o mo_rte_kind.o
	$(FC) -c $(SRC)/mod_layer.F90

mod_mlp.o: $(SRC)/mod_mlp.F90 mo_rte_kind.o mod_mlp_layer.o
	$(FC) -c $(SRC)/mod_mlp.F90

mod_mlp_activation.o: $(SRC)/mod_mlp_activation.F90 mo_rte_kind.o
	$(FC) -c $(SRC)/mod_mlp_activation.F90

mod_mlp_layer.o: $(SRC)/mod_mlp_layer.F90 mod_mlp_activation.o mo_rte_kind.o
	$(FC) -c $(SRC)/mod_mlp_layer.F90

mod_mlp_rrtmgp.o: $(SRC)/mod_mlp_rrtmgp.F90 mo_rte_kind.o mod_mlp.o mo_simple_netcdf.o netcdf.o radiation_io.o
	$(FC) -c $(SRC)/mod_mlp_rrtmgp.F90

mod_network_rrtmgp.o: $(SRC)/mod_network_rrtmgp.F90 mo_rte_kind.o mod_neural_network.o mo_simple_netcdf.o netcdf.o
	$(FC) -c $(SRC)/mod_network_rrtmgp.F90

mod_neural_network.o: $(SRC)/mod_neural_network.F90 mo_rte_kind.o mod_layer.o
	$(FC) -c $(SRC)/mod_neural_network.F90

modd_budget.o: $(SRC)/modd_budget.F90 ddh_mix.o yomlddh.o yommddh.o
	$(FC) -c $(SRC)/modd_budget.F90

modd_cloudparn.o: $(SRC)/modd_cloudparn.F90 modd_parameters.o
	$(FC) -c $(SRC)/modd_cloudparn.F90

modd_convpar.o: $(SRC)/modd_convpar.F90 yomhook.o
	$(FC) -c $(SRC)/modd_convpar.F90

modd_convpar_shal.o: $(SRC)/modd_convpar_shal.F90 
	$(FC) -c $(SRC)/modd_convpar_shal.F90

modd_convparext.o: $(SRC)/modd_convparext.F90 
	$(FC) -c $(SRC)/modd_convparext.F90

modd_cst.o: $(SRC)/modd_cst.F90 
	$(FC) -c $(SRC)/modd_cst.F90

modd_cturb.o: $(SRC)/modd_cturb.F90 
	$(FC) -c $(SRC)/modd_cturb.F90

modd_field.o: $(SRC)/modd_field.F90 modd_parameters.o
	$(FC) -c $(SRC)/modd_field.F90

modd_io.o: $(SRC)/modd_io.F90 modd_parameters.o
	$(FC) -c $(SRC)/modd_io.F90

modd_misc.o: $(SRC)/modd_misc.F90 modd_budget.o
	$(FC) -c $(SRC)/modd_misc.F90

modd_nebn.o: $(SRC)/modd_nebn.F90 modd_parameters.o mode_posnam_phy.o mode_check_nam_val.o
	$(FC) -c $(SRC)/modd_nebn.F90

modd_nsv.o: $(SRC)/modd_nsv.F90 modd_field.o modd_parameters.o
	$(FC) -c $(SRC)/modd_nsv.F90

modd_param_icen.o: $(SRC)/modd_param_icen.F90 modd_parameters.o mode_posnam_phy.o mode_msg.o mode_check_nam_val.o
	$(FC) -c $(SRC)/modd_param_icen.F90

modd_param_lima.o: $(SRC)/modd_param_lima.F90 modd_parameters.o mode_posnam_phy.o mode_msg.o mode_check_nam_val.o
	$(FC) -c $(SRC)/modd_param_lima.F90

modd_param_lima_cold.o: $(SRC)/modd_param_lima_cold.F90 modd_parameters.o
	$(FC) -c $(SRC)/modd_param_lima_cold.F90

modd_param_lima_mixed.o: $(SRC)/modd_param_lima_mixed.F90 
	$(FC) -c $(SRC)/modd_param_lima_mixed.F90

modd_param_lima_warm.o: $(SRC)/modd_param_lima_warm.F90 modd_parameters.o
	$(FC) -c $(SRC)/modd_param_lima_warm.F90

modd_param_mfshalln.o: $(SRC)/modd_param_mfshalln.F90 modd_parameters.o mode_posnam_phy.o mode_check_nam_val.o
	$(FC) -c $(SRC)/modd_param_mfshalln.F90

modd_parameters.o: $(SRC)/modd_parameters.F90 
	$(FC) -c $(SRC)/modd_parameters.F90

modd_phyex.o: $(SRC)/modd_phyex.F90 modd_cst.o modd_param_icen.o modd_rain_ice_descrn.o modd_rain_ice_paramn.o modd_cloudparn.o modd_param_mfshalln.o modd_turbn.o modd_cturb.o modd_nebn.o modd_param_lima.o modd_param_lima_warm.o modd_param_lima_cold.o modd_param_lima_mixed.o modd_nsv.o modd_misc.o
	$(FC) -c $(SRC)/modd_phyex.F90

modd_rain_ice_descrn.o: $(SRC)/modd_rain_ice_descrn.F90 modd_parameters.o
	$(FC) -c $(SRC)/modd_rain_ice_descrn.F90

modd_rain_ice_paramn.o: $(SRC)/modd_rain_ice_paramn.F90 modd_parameters.o
	$(FC) -c $(SRC)/modd_rain_ice_paramn.F90

modd_turbn.o: $(SRC)/modd_turbn.F90 modd_parameters.o mode_posnam_phy.o mode_check_nam_val.o
	$(FC) -c $(SRC)/modd_turbn.F90

mode_check_nam_val.o: $(SRC)/mode_check_nam_val.F90 mode_msg.o
	$(FC) -c $(SRC)/mode_check_nam_val.F90

mode_msg.o: $(SRC)/mode_msg.F90 modd_io.o yomlun.o abor1.intfb.ok
	$(FC) -c $(SRC)/mode_msg.F90

mode_posnam_phy.o: $(SRC)/mode_posnam_phy.F90 
	$(FC) -c $(SRC)/mode_posnam_phy.F90

model_atmos_ocean_coupling_mod.o: $(SRC)/model_atmos_ocean_coupling_mod.F90 yommcc.o yomcom.o yomcou.o
	$(FC) -c $(SRC)/model_atmos_ocean_coupling_mod.F90

model_chem_mod.o: $(SRC)/model_chem_mod.F90 yomozo.o yomchem.o yomcompo.o
	$(FC) -c $(SRC)/model_chem_mod.F90

model_diagnostics_mod.o: $(SRC)/model_diagnostics_mod.F90 yomcddh.o yomlddh.o yommddh.o yomsddh.o yomtddh.o yomgpddh.o yompaddh.o yomspddh.o
	$(FC) -c $(SRC)/model_diagnostics_mod.F90

model_dynamics_mod.o: $(SRC)/model_dynamics_mod.F90 yomdyna.o yomdyn.o yemdyn.o spng_mod.o ptrgppc.o intdynsl_mod.o yomslint.o yomslrep.o ptrslb1.o ptrslb2.o ptrslb15.o yomtnh.o eint_mod.o
	$(FC) -c $(SRC)/model_dynamics_mod.F90

model_general_conf_mod.o: $(SRC)/model_general_conf_mod.F90 type_geometry.o yomdimf.o yom_ygfl.o yomrip.o yommoderrmod.o yomspsdt.o spp_mod.o type_ecv.o yomhook.o
	$(FC) -c $(SRC)/model_general_conf_mod.F90

model_physics_aerosol_mod.o: $(SRC)/model_physics_aerosol_mod.F90 yoeaerlid.o yoeaermap.o yoeaersnk.o yoeaersrc.o yoeaervol.o yoedbug.o
	$(FC) -c $(SRC)/model_physics_aerosol_mod.F90

model_physics_ecmwf_mod.o: $(SRC)/model_physics_ecmwf_mod.F90 yoephy.o yoecld.o yoecldp.o yoecnd.o yoecumf.o yoe_cuconvca.o yoegwd.o yoegwwms.o yoethf.o
	$(FC) -c $(SRC)/model_physics_ecmwf_mod.F90

model_physics_general_mod.o: $(SRC)/model_physics_general_mod.F90 yomdphy.o yomslphy.o yoevdf.o
	$(FC) -c $(SRC)/model_physics_general_mod.F90

model_physics_mf_mod.o: $(SRC)/model_physics_mf_mod.F90 yomphy.o yomphy0.o yomphy1.o yomphy2.o yomphy3.o yomphyds.o yomcvmnh.o yomtoph.o yomvdoz.o yomsimphl.o yomarphy.o yomparar.o yommse.o yomlouis.o yomnorgwd.o yomnrtaer.o yoaiop.o eint_mod.o yomdgradient.o parkind1.o
	$(FC) -c $(SRC)/model_physics_mf_mod.F90

model_physics_radiation_mod.o: $(SRC)/model_physics_radiation_mod.F90 yomradf.o yoerad.o yoesw.o yoeovlp.o yoeneur.o yoelwrad.o yoeaerd.o yoeaeratm.o yoe_uvrad.o yoerdi.o yoe_mcica.o yomrcoef.o yomtrc.o yomprad.o yomgsgeom.o yoerip.o radiation_setup.o eint_mod.o parkind1.o yomhook.o
	$(FC) -c $(SRC)/model_physics_radiation_mod.F90

model_physics_simplinear_mod.o: $(SRC)/model_physics_simplinear_mod.F90 yoephli.o yomcumfs.o yoegwdwms.o yoecumf2.o yophlc.o yophnc.o yomncl.o yomsrftlad.o yomsphyhist.o
	$(FC) -c $(SRC)/model_physics_simplinear_mod.F90

model_physics_stochast_mod.o: $(SRC)/model_physics_stochast_mod.F90 stoph_mix.o yomrandom_streams.o yompertpar.o
	$(FC) -c $(SRC)/model_physics_stochast_mod.F90

netcdf.o: $(SRC)/netcdf.F90 typeSizes.o
	$(FC) -c $(SRC)/netcdf.F90

order_independent_summation_mod.o: $(SRC)/order_independent_summation_mod.F90 parkind1.o compensated_summation_mod.o yomhook.o
	$(FC) -c $(SRC)/order_independent_summation_mod.F90

par_gfl.o: $(SRC)/par_gfl.F90 parkind1.o
	$(FC) -c $(SRC)/par_gfl.F90

pardim.o: $(SRC)/pardim.F90 parkind1.o
	$(FC) -c $(SRC)/pardim.F90

parfpos.o: $(SRC)/parfpos.F90 parkind1.o pardim.o
	$(FC) -c $(SRC)/parfpos.F90

parkind1.o: $(SRC)/parkind1.F90 
	$(FC) -c $(SRC)/parkind1.F90

parkind_ocean.o: $(SRC)/parkind_ocean.F90 
	$(FC) -c $(SRC)/parkind_ocean.F90

parrrtm.o: $(SRC)/parrrtm.F90 parkind1.o
	$(FC) -c $(SRC)/parrrtm.F90

parsrtm.o: $(SRC)/parsrtm.F90 parkind1.o
	$(FC) -c $(SRC)/parsrtm.F90

ptrgfu.o: $(SRC)/ptrgfu.F90 ptrgfu_type.o
	$(FC) -c $(SRC)/ptrgfu.F90

ptrgfu_type.o: $(SRC)/ptrgfu_type.F90 parkind1.o
	$(FC) -c $(SRC)/ptrgfu_type.F90

ptrgppc.o: $(SRC)/ptrgppc.F90 parkind1.o
	$(FC) -c $(SRC)/ptrgppc.F90

ptrslb1.o: $(SRC)/ptrslb1.F90 parkind1.o
	$(FC) -c $(SRC)/ptrslb1.F90

ptrslb15.o: $(SRC)/ptrslb15.F90 parkind1.o
	$(FC) -c $(SRC)/ptrslb15.F90

ptrslb2.o: $(SRC)/ptrslb2.F90 parkind1.o
	$(FC) -c $(SRC)/ptrslb2.F90

ptrxfu.o: $(SRC)/ptrxfu.F90 ptrxfu_type.o
	$(FC) -c $(SRC)/ptrxfu.F90

ptrxfu_type.o: $(SRC)/ptrxfu_type.F90 parkind1.o
	$(FC) -c $(SRC)/ptrxfu_type.F90

radiation_adding_ica_lw.o: $(SRC)/radiation_adding_ica_lw.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/radiation_adding_ica_lw.F90

radiation_adding_ica_sw.o: $(SRC)/radiation_adding_ica_sw.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/radiation_adding_ica_sw.F90

radiation_aerosol.o: $(SRC)/radiation_aerosol.F90 parkind1.o radiation_io.o yomhook.o radiation_config.o radiation_check.o
	$(FC) -c $(SRC)/radiation_aerosol.F90

radiation_aerosol_optics.o: $(SRC)/radiation_aerosol_optics.F90 parkind1.o yomhook.o radiation_config.o radiation_aerosol_optics_data.o radiation_io.o easy_netcdf.o radiation_spectral_definition.o radiation_thermodynamics.o radiation_gas.o radiation_aerosol.o radiation_constants.o ecrad_config.ok radiation_delta_eddington.ok
	$(FC) -c $(SRC)/radiation_aerosol_optics.F90

radiation_aerosol_optics_data.o: $(SRC)/radiation_aerosol_optics_data.F90 parkind1.o radiation_io.o yomhook.o easy_netcdf.o ecrad_config.ok
	$(FC) -c $(SRC)/radiation_aerosol_optics_data.F90

radiation_check.o: $(SRC)/radiation_check.F90 parkind1.o radiation_io.o
	$(FC) -c $(SRC)/radiation_check.F90

radiation_cloud.o: $(SRC)/radiation_cloud.F90 parkind1.o yomhook.o radiation_thermodynamics.o radiation_constants.o radiation_io.o radiation_check.o
	$(FC) -c $(SRC)/radiation_cloud.F90

radiation_cloud_cover.o: $(SRC)/radiation_cloud_cover.F90 parkind1.o yomhook.o ecrad_config.ok
	$(FC) -c $(SRC)/radiation_cloud_cover.F90

radiation_cloud_generator.o: $(SRC)/radiation_cloud_generator.F90 parkind1.o yomhook.o radiation_io.o random_numbers_mix.o radiation_pdf_sampler.o radiation_cloud_cover.o radiation_random_numbers.o
	$(FC) -c $(SRC)/radiation_cloud_generator.F90

radiation_cloud_generator_acc.o: $(SRC)/radiation_cloud_generator_acc.F90 parkind1.o radiation_random_numbers.o
	$(FC) -c $(SRC)/radiation_cloud_generator_acc.F90

radiation_cloud_optics.o: $(SRC)/radiation_cloud_optics.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_ice_optics_fu.o radiation_ice_optics_baran.o radiation_ice_optics_baran2017.o radiation_ice_optics_yi.o radiation_liquid_optics_socrates.o radiation_liquid_optics_slingo.o radiation_thermodynamics.o radiation_cloud.o radiation_constants.o radiation_delta_eddington.ok
	$(FC) -c $(SRC)/radiation_cloud_optics.F90

radiation_cloud_optics_data.o: $(SRC)/radiation_cloud_optics_data.F90 parkind1.o yomhook.o easy_netcdf.o ecrad_config.ok
	$(FC) -c $(SRC)/radiation_cloud_optics_data.F90

radiation_cloudless_lw.o: $(SRC)/radiation_cloudless_lw.F90 parkind1.o yomhook.o radiation_config.o radiation_flux.o radiation_two_stream.o radiation_adding_ica_lw.o radiation_lw_derivatives.o
	$(FC) -c $(SRC)/radiation_cloudless_lw.F90

radiation_cloudless_sw.o: $(SRC)/radiation_cloudless_sw.F90 parkind1.o yomhook.o radiation_config.o radiation_single_level.o radiation_flux.o radiation_two_stream.o radiation_constants.o radiation_adding_ica_sw.o radiation_delta_eddington.ok
	$(FC) -c $(SRC)/radiation_cloudless_sw.F90

radiation_config.o: $(SRC)/radiation_config.F90 parkind1.o radiation_cloud_optics_data.o radiation_general_cloud_optics_data.o radiation_aerosol_optics_data.o radiation_pdf_sampler.o radiation_cloud_cover.o radiation_ecckd.o mo_gas_optics_rrtmgp.o mod_network_rrtmgp.o yomhook.o radiation_io.o radiation_spectral_definition.o ecrad_config.ok
	$(FC) -c $(SRC)/radiation_config.F90

radiation_constants.o: $(SRC)/radiation_constants.F90 parkind1.o
	$(FC) -c $(SRC)/radiation_constants.F90

radiation_ecckd.o: $(SRC)/radiation_ecckd.F90 parkind1.o radiation_gas_constants.o radiation_ecckd_gas.o radiation_spectral_definition.o easy_netcdf.o yomhook.o radiation_io.o radiation_constants.o ecrad_config.ok
	$(FC) -c $(SRC)/radiation_ecckd.F90

radiation_ecckd_gas.o: $(SRC)/radiation_ecckd_gas.F90 parkind1.o radiation_gas_constants.o easy_netcdf.o ecrad_config.ok
	$(FC) -c $(SRC)/radiation_ecckd_gas.F90

radiation_ecckd_interface.o: $(SRC)/radiation_ecckd_interface.F90 parkind1.o radiation_config.o yomhook.o radiation_gas.o radiation_thermodynamics.o radiation_single_level.o radiation_gas_constants.o
	$(FC) -c $(SRC)/radiation_ecckd_interface.F90

radiation_flux.o: $(SRC)/radiation_flux.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_check.o radiation_constants.o ecrad_config.ok
	$(FC) -c $(SRC)/radiation_flux.F90

radiation_gas.o: $(SRC)/radiation_gas.F90 parkind1.o radiation_gas_constants.o yomhook.o radiation_io.o radiation_check.o
	$(FC) -c $(SRC)/radiation_gas.F90

radiation_gas_constants.o: $(SRC)/radiation_gas_constants.F90 parkind1.o
	$(FC) -c $(SRC)/radiation_gas_constants.F90

radiation_general_cloud_optics.o: $(SRC)/radiation_general_cloud_optics.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_spectral_definition.o radiation_thermodynamics.o radiation_cloud.o radiation_constants.o easy_netcdf.o radiation_delta_eddington.ok
	$(FC) -c $(SRC)/radiation_general_cloud_optics.F90

radiation_general_cloud_optics_data.o: $(SRC)/radiation_general_cloud_optics_data.F90 parkind1.o yomhook.o easy_netcdf.o radiation_spectral_definition.o radiation_io.o radiation_constants.o ecrad_config.ok radiation_delta_eddington.ok
	$(FC) -c $(SRC)/radiation_general_cloud_optics_data.F90

radiation_homogeneous_lw.o: $(SRC)/radiation_homogeneous_lw.F90 parkind1.o yomhook.o radiation_config.o radiation_cloud.o radiation_flux.o radiation_two_stream.o radiation_adding_ica_lw.o radiation_lw_derivatives.o
	$(FC) -c $(SRC)/radiation_homogeneous_lw.F90

radiation_homogeneous_sw.o: $(SRC)/radiation_homogeneous_sw.F90 parkind1.o yomhook.o radiation_config.o radiation_single_level.o radiation_cloud.o radiation_flux.o radiation_two_stream.o radiation_constants.o radiation_adding_ica_sw.o radiation_delta_eddington.ok
	$(FC) -c $(SRC)/radiation_homogeneous_sw.F90

radiation_ice_optics_baran.o: $(SRC)/radiation_ice_optics_baran.F90 parkind1.o
	$(FC) -c $(SRC)/radiation_ice_optics_baran.F90

radiation_ice_optics_baran2017.o: $(SRC)/radiation_ice_optics_baran2017.F90 parkind1.o
	$(FC) -c $(SRC)/radiation_ice_optics_baran2017.F90

radiation_ice_optics_fu.o: $(SRC)/radiation_ice_optics_fu.F90 parkind1.o
	$(FC) -c $(SRC)/radiation_ice_optics_fu.F90

radiation_ice_optics_yi.o: $(SRC)/radiation_ice_optics_yi.F90 parkind1.o
	$(FC) -c $(SRC)/radiation_ice_optics_yi.F90

radiation_ifs_rrtm.o: $(SRC)/radiation_ifs_rrtm.F90 yoerrtm.o yoesrtm.o yoerrtftr.o yomhook.o radiation_config.o radiation_spectral_definition.o radiation_gas.o parkind1.o radiation_io.o parrrtm.o radiation_thermodynamics.o radiation_single_level.o yoerrtwn.o
	$(FC) -c $(SRC)/radiation_ifs_rrtm.F90

radiation_ifs_rrtmgp.o: $(SRC)/radiation_ifs_rrtmgp.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_spectral_definition.o mo_load_coefficients.o mo_gas_concentrations.o mo_gas_optics_rrtmgp.o radiation_gas.o radiation_thermodynamics.o radiation_single_level.o radiation_gas_constants.o mo_rte_util_array.o
	$(FC) -c $(SRC)/radiation_ifs_rrtmgp.F90

radiation_interface.o: $(SRC)/radiation_interface.F90 radiation_interface_gpu.o radiation_interface_cpu.o
	$(FC) -c $(SRC)/radiation_interface.F90

radiation_interface_cpu.o: $(SRC)/radiation_interface_cpu.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_spectral_definition.o radiation_monochromatic.o radiation_ifs_rrtm.o radiation_ecckd_interface.o radiation_ifs_rrtmgp.o radiation_cloud_optics.o radiation_general_cloud_optics.o radiation_aerosol_optics.o yoerdi.o radiation_gas.o radiation_single_level.o radiation_thermodynamics.o radiation_cloud.o radiation_aerosol.o radiation_flux.o radiation_spartacus_sw.o radiation_spartacus_lw.o radiation_tripleclouds_sw.o radiation_tripleclouds_lw.o radiation_mcica_sw.o radiation_mcica_lw.o radiation_mcica_acc_sw.o radiation_mcica_acc_lw.o radiation_cloudless_sw.o radiation_cloudless_lw.o radiation_homogeneous_sw.o radiation_homogeneous_lw.o radiation_save.o
	$(FC) -c $(SRC)/radiation_interface_cpu.F90

radiation_interface_gpu.o: $(SRC)/radiation_interface_gpu.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_spectral_definition.o radiation_monochromatic.o radiation_ifs_rrtm.o radiation_ecckd_interface.o radiation_ifs_rrtmgp.o radiation_cloud_optics.o radiation_general_cloud_optics.o radiation_aerosol_optics.o yoerdi.o radiation_gas.o radiation_single_level.o radiation_thermodynamics.o radiation_cloud.o radiation_aerosol.o radiation_flux.o radiation_spartacus_sw.o radiation_spartacus_lw.o radiation_tripleclouds_sw.o radiation_tripleclouds_lw.o radiation_mcica_sw.o radiation_mcica_lw.o radiation_mcica_acc_sw.o radiation_mcica_acc_lw.o radiation_cloudless_sw.o radiation_cloudless_lw.o radiation_homogeneous_sw.o radiation_homogeneous_lw.o radiation_save.o
	$(FC) -c $(SRC)/radiation_interface_gpu.F90

radiation_io.o: $(SRC)/radiation_io.F90 yomlun_ifsaux.o
	$(FC) -c $(SRC)/radiation_io.F90

radiation_liquid_optics_slingo.o: $(SRC)/radiation_liquid_optics_slingo.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/radiation_liquid_optics_slingo.F90

radiation_liquid_optics_socrates.o: $(SRC)/radiation_liquid_optics_socrates.F90 parkind1.o
	$(FC) -c $(SRC)/radiation_liquid_optics_socrates.F90

radiation_lw_derivatives.o: $(SRC)/radiation_lw_derivatives.F90 parkind1.o yomhook.o radiation_matrix.o
	$(FC) -c $(SRC)/radiation_lw_derivatives.F90

radiation_matrix.o: $(SRC)/radiation_matrix.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/radiation_matrix.F90

radiation_mcica_acc_lw.o: $(SRC)/radiation_mcica_acc_lw.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_single_level.o radiation_cloud.o radiation_flux.o radiation_two_stream.o radiation_adding_ica_lw.o radiation_lw_derivatives.o radiation_cloud_generator_acc.o radiation_cloud_cover.o
	$(FC) -c $(SRC)/radiation_mcica_acc_lw.F90

radiation_mcica_acc_sw.o: $(SRC)/radiation_mcica_acc_sw.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_single_level.o radiation_cloud.o radiation_flux.o radiation_two_stream.o radiation_adding_ica_sw.o radiation_cloud_generator_acc.o radiation_cloud_cover.o radiation_delta_eddington.ok
	$(FC) -c $(SRC)/radiation_mcica_acc_sw.F90

radiation_mcica_lw.o: $(SRC)/radiation_mcica_lw.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_single_level.o radiation_cloud.o radiation_flux.o radiation_two_stream.o radiation_adding_ica_lw.o radiation_lw_derivatives.o radiation_cloud_generator.o ecrad_config.ok
	$(FC) -c $(SRC)/radiation_mcica_lw.F90

radiation_mcica_sw.o: $(SRC)/radiation_mcica_sw.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_single_level.o radiation_cloud.o radiation_flux.o radiation_two_stream.o radiation_adding_ica_sw.o radiation_cloud_generator.o ecrad_config.ok radiation_delta_eddington.ok
	$(FC) -c $(SRC)/radiation_mcica_sw.F90

radiation_monochromatic.o: $(SRC)/radiation_monochromatic.F90 radiation_config.o radiation_gas.o parkind1.o radiation_thermodynamics.o radiation_single_level.o radiation_constants.o radiation_cloud.o radiation_aerosol.o radiation_delta_eddington.ok
	$(FC) -c $(SRC)/radiation_monochromatic.F90

radiation_overlap.o: $(SRC)/radiation_overlap.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/radiation_overlap.F90

radiation_pdf_sampler.o: $(SRC)/radiation_pdf_sampler.F90 parkind1.o yomhook.o easy_netcdf.o
	$(FC) -c $(SRC)/radiation_pdf_sampler.F90

radiation_random_numbers.o: $(SRC)/radiation_random_numbers.F90 parkind1.o
	$(FC) -c $(SRC)/radiation_random_numbers.F90

radiation_regions.o: $(SRC)/radiation_regions.F90 parkind1.o yomhook.o radiation_io.o
	$(FC) -c $(SRC)/radiation_regions.F90

radiation_save.o: $(SRC)/radiation_save.F90 parkind1.o yomhook.o easy_netcdf.o radiation_io.o radiation_config.o radiation_thermodynamics.o radiation_flux.o radiation_single_level.o radiation_cloud.o radiation_gas.o radiation_aerosol.o radiation_matrix.o
	$(FC) -c $(SRC)/radiation_save.F90

radiation_setup.o: $(SRC)/radiation_setup.F90 parkind1.o radiation_config.o yoeaerop.o yoe_spectral_planck.o yomhook.o yomlun.o yoerad.o yoephy.o yoeaeratm.o yomcompo.o yoerdi.o radiation_interface.o radiation_aerosol_optics.o radiation_aerosol_optics_data.o
	$(FC) -c $(SRC)/radiation_setup.F90

radiation_single_level.o: $(SRC)/radiation_single_level.F90 parkind1.o yomhook.o radiation_config.o radiation_io.o radiation_check.o
	$(FC) -c $(SRC)/radiation_single_level.F90

radiation_spartacus_lw.o: $(SRC)/radiation_spartacus_lw.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_thermodynamics.o radiation_cloud.o radiation_regions.o radiation_overlap.o radiation_flux.o radiation_matrix.o radiation_two_stream.o radiation_lw_derivatives.o radiation_constants.o radiation_optical_depth_scaling.ok
	$(FC) -c $(SRC)/radiation_spartacus_lw.F90

radiation_spartacus_sw.o: $(SRC)/radiation_spartacus_sw.F90 parkind1.o yomhook.o radiation_io.o radiation_config.o radiation_single_level.o radiation_thermodynamics.o radiation_cloud.o radiation_regions.o radiation_overlap.o radiation_flux.o radiation_matrix.o radiation_two_stream.o radiation_constants.o radiation_optical_depth_scaling.ok
	$(FC) -c $(SRC)/radiation_spartacus_sw.F90

radiation_spectral_definition.o: $(SRC)/radiation_spectral_definition.F90 parkind1.o easy_netcdf.o yomhook.o radiation_io.o radiation_constants.o ecrad_config.ok
	$(FC) -c $(SRC)/radiation_spectral_definition.F90

radiation_thermodynamics.o: $(SRC)/radiation_thermodynamics.F90 parkind1.o yomhook.o radiation_constants.o radiation_check.o
	$(FC) -c $(SRC)/radiation_thermodynamics.F90

radiation_tripleclouds_lw.o: $(SRC)/radiation_tripleclouds_lw.F90 parkind1.o yomhook.o radiation_config.o radiation_cloud.o radiation_regions.o radiation_overlap.o radiation_flux.o radiation_matrix.o radiation_two_stream.o radiation_adding_ica_lw.o radiation_lw_derivatives.o radiation_optical_depth_scaling.ok
	$(FC) -c $(SRC)/radiation_tripleclouds_lw.F90

radiation_tripleclouds_sw.o: $(SRC)/radiation_tripleclouds_sw.F90 parkind1.o yomhook.o radiation_config.o radiation_single_level.o radiation_cloud.o radiation_regions.o radiation_overlap.o radiation_flux.o radiation_matrix.o radiation_two_stream.o radiation_delta_eddington.ok radiation_optical_depth_scaling.ok
	$(FC) -c $(SRC)/radiation_tripleclouds_sw.F90

radiation_two_stream.o: $(SRC)/radiation_two_stream.F90 parkind1.o ecrad_config.ok
	$(FC) -c $(SRC)/radiation_two_stream.F90

random_numbers_mix.o: $(SRC)/random_numbers_mix.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/random_numbers_mix.F90

reglatlon_field_mix.o: $(SRC)/reglatlon_field_mix.F90 parkind1.o yomlun.o yomhook.o deallocate_if_associated_mod.o
	$(FC) -c $(SRC)/reglatlon_field_mix.F90

sc2prg_mod.o: $(SRC)/sc2prg_mod.F90 parkind1.o yomhook.o yomtraj.o
	$(FC) -c $(SRC)/sc2prg_mod.F90

spectral_arp_mod.o: $(SRC)/spectral_arp_mod.F90 geometry_mod.o parkind1.o yomhook.o spectral_fields_mod.o spectral_fields_data.o random_numbers_mix.o
	$(FC) -c $(SRC)/spectral_arp_mod.F90

spectral_fields_data.o: $(SRC)/spectral_fields_data.F90 parkind1.o
	$(FC) -c $(SRC)/spectral_fields_data.F90

spectral_fields_mod.o: $(SRC)/spectral_fields_mod.F90 parkind1.o yomhook.o geometry_mod.o spectral_variables_mod.o spectral_fields_data.o spectral_fields_oper_mod.o spectral_fields_para_mod.o random_numbers_mix.o
	$(FC) -c $(SRC)/spectral_fields_mod.F90

spectral_fields_oper_mod.o: $(SRC)/spectral_fields_oper_mod.F90 parkind1.o yomhook.o spectral_fields_data.o yomlun.o
	$(FC) -c $(SRC)/spectral_fields_oper_mod.F90

spectral_fields_para_mod.o: $(SRC)/spectral_fields_para_mod.F90 parkind1.o yomhook.o spectral_fields_data.o
	$(FC) -c $(SRC)/spectral_fields_para_mod.F90

spectral_variables_mod.o: $(SRC)/spectral_variables_mod.F90 parkind1.o
	$(FC) -c $(SRC)/spectral_variables_mod.F90

spng_mod.o: $(SRC)/spng_mod.F90 parkind1.o yomhook.o yomcst.o yomlun.o yomrip.o yomdyna.o
	$(FC) -c $(SRC)/spng_mod.F90

spp_def_mod.o: $(SRC)/spp_def_mod.F90 parkind1.o spp_gen_mod.o abor1.intfb.ok
	$(FC) -c $(SRC)/spp_def_mod.F90

spp_gen_mod.o: $(SRC)/spp_gen_mod.F90 parkind1.o
	$(FC) -c $(SRC)/spp_gen_mod.F90

spp_mod.o: $(SRC)/spp_mod.F90 parkind1.o yomhook.o spp_gen_mod.o spp_def_mod.o spectral_arp_mod.o gridpoint_fields_mix.o yomlun.o abor1.intfb.ok
	$(FC) -c $(SRC)/spp_mod.F90

stoph_mix.o: $(SRC)/stoph_mix.F90 parkind1.o yomrandom_streams.o random_numbers_mix.o yomlun.o yomhook.o yomcst.o geometry_mod.o yomct3.o yommp0.o yomrip.o yomlap.o yomdim.o yomtag.o abor1.intfb.ok
	$(FC) -c $(SRC)/stoph_mix.F90

surface_fields_mix.o: $(SRC)/surface_fields_mix.F90 parkind1.o yomhook.o yomdim.o yomlun.o yomct0.o yomdyn.o geometry_mod.o abor1.intfb.ok
	$(FC) -c $(SRC)/surface_fields_mix.F90

surface_variables_mod.o: $(SRC)/surface_variables_mod.F90 parkind1.o variable_module.o field_module.o field_factory_module.o
	$(FC) -c $(SRC)/surface_variables_mod.F90

surface_views_diagnostic_module.o: $(SRC)/surface_views_diagnostic_module.F90 parkind1.o field_module.o surface_variables_mod.o
	$(FC) -c $(SRC)/surface_views_diagnostic_module.F90

surface_views_module.o: $(SRC)/surface_views_module.F90 surface_views_diagnostic_module.o surface_views_prognostic_module.o
	$(FC) -c $(SRC)/surface_views_module.F90

surface_views_prognostic_module.o: $(SRC)/surface_views_prognostic_module.F90 parkind1.o field_module.o surface_variables_mod.o yomdyn.o cpg_opts_type_mod.o field_access_module.o
	$(FC) -c $(SRC)/surface_views_prognostic_module.F90

type_ecv.o: $(SRC)/type_ecv.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/type_ecv.F90

type_faoph.o: $(SRC)/type_faoph.F90 parkind1.o
	$(FC) -c $(SRC)/type_faoph.F90

type_fluxes.o: $(SRC)/type_fluxes.F90 parkind1.o pardim.o yomhook.o abor1.intfb.ok
	$(FC) -c $(SRC)/type_fluxes.F90

type_fpdsphys.o: $(SRC)/type_fpdsphys.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/type_fpdsphys.F90

type_fpofn.o: $(SRC)/type_fpofn.F90 
	$(FC) -c $(SRC)/type_fpofn.F90

type_fposbuf.o: $(SRC)/type_fposbuf.F90 parkind1.o yomfp4l.o
	$(FC) -c $(SRC)/type_fposbuf.F90

type_fpusergeo.o: $(SRC)/type_fpusergeo.F90 parkind1.o yomhook.o abor1.intfb.ok
	$(FC) -c $(SRC)/type_fpusergeo.F90

type_geometry.o: $(SRC)/type_geometry.F90 yomvert.o yomsta.o yomlap.o yomleg.o yomdim.o yomdimv.o yommp.o yomgem.o yomcsgeom.o yomgsgeom.o yomorog.o type_spgeom.o yemdim.o yemgeo.o yemmp.o yemlap.o yemgsl.o yemlbc_geo.o yomcver.o
	$(FC) -c $(SRC)/type_geometry.F90

type_gmvs.o: $(SRC)/type_gmvs.F90 parkind1.o
	$(FC) -c $(SRC)/type_gmvs.F90

type_model.o: $(SRC)/type_model.F90 parkind1.o model_general_conf_mod.o model_atmos_ocean_coupling_mod.o model_dynamics_mod.o model_physics_general_mod.o model_physics_ecmwf_mod.o model_physics_simplinear_mod.o model_physics_aerosol_mod.o model_physics_radiation_mod.o model_physics_stochast_mod.o model_physics_mf_mod.o model_chem_mod.o model_diagnostics_mod.o yoewcou.o yomcst.o yomspsdt.o spp_mod.o yemlbc_model.o yomlun.o yomhook.o
	$(FC) -c $(SRC)/type_model.F90

type_spgeom.o: $(SRC)/type_spgeom.F90 parkind1.o
	$(FC) -c $(SRC)/type_spgeom.F90

types_fpcat.o: $(SRC)/types_fpcat.F90 parkind1.o yomhook.o yomlun.o
	$(FC) -c $(SRC)/types_fpcat.F90

typeSizes.o: $(SRC)/typeSizes.F90 
	$(FC) -c $(SRC)/typeSizes.F90

variable_module.o: $(SRC)/variable_module.F90 parkind1.o field_module.o yom_ygfl.o
	$(FC) -c $(SRC)/variable_module.F90

variables_mod.o: $(SRC)/variables_mod.F90 parkind1.o abor1.intfb.ok
	$(FC) -c $(SRC)/variables_mod.F90

yemdim.o: $(SRC)/yemdim.F90 parkind1.o
	$(FC) -c $(SRC)/yemdim.F90

yemdyn.o: $(SRC)/yemdyn.F90 parkind1.o
	$(FC) -c $(SRC)/yemdyn.F90

yemgeo.o: $(SRC)/yemgeo.F90 parkind1.o
	$(FC) -c $(SRC)/yemgeo.F90

yemgsl.o: $(SRC)/yemgsl.F90 parkind1.o
	$(FC) -c $(SRC)/yemgsl.F90

yemlap.o: $(SRC)/yemlap.F90 parkind1.o
	$(FC) -c $(SRC)/yemlap.F90

yemlbc_fields.o: $(SRC)/yemlbc_fields.F90 parkind1.o yomhook.o yomct0.o yomlun.o yomgmv.o yemlbc_model.o geometry_mod.o yom_ygfl.o yomdyna.o
	$(FC) -c $(SRC)/yemlbc_fields.F90

yemlbc_geo.o: $(SRC)/yemlbc_geo.F90 parkind1.o
	$(FC) -c $(SRC)/yemlbc_geo.F90

yemlbc_model.o: $(SRC)/yemlbc_model.F90 parkind1.o model_general_conf_mod.o yomhook.o yomdfi.o yomct0.o yomgmv.o yomini.o yomlun.o yommp0.o geometry_mod.o yom_ygfl.o yomdyna.o abor1.intfb.ok
	$(FC) -c $(SRC)/yemlbc_model.F90

yemmp.o: $(SRC)/yemmp.F90 parkind1.o
	$(FC) -c $(SRC)/yemmp.F90

yoaiop.o: $(SRC)/yoaiop.F90 parkind1.o
	$(FC) -c $(SRC)/yoaiop.F90

yoe_aerodiag.o: $(SRC)/yoe_aerodiag.F90 parkind1.o
	$(FC) -c $(SRC)/yoe_aerodiag.F90

yoe_cuconvca.o: $(SRC)/yoe_cuconvca.F90 parkind1.o random_numbers_mix.o geometry_mod.o yomdyna.o yomhook.o yomct0.o yomcst.o yomlun.o yomgrib.o eint_mod.o yommp0.o yoecumf.o yomct3.o order_independent_summation_mod.o
	$(FC) -c $(SRC)/yoe_cuconvca.F90

yoe_mcica.o: $(SRC)/yoe_mcica.F90 parkind1.o
	$(FC) -c $(SRC)/yoe_mcica.F90

yoe_phys_mwave.o: $(SRC)/yoe_phys_mwave.F90 parkind1.o field_module.o field_factory_module.o geometry_mod.o
	$(FC) -c $(SRC)/yoe_phys_mwave.F90

yoe_spectral_planck.o: $(SRC)/yoe_spectral_planck.F90 parkind1.o yomhook.o yomcst.o yomlun.o abor1.intfb.ok
	$(FC) -c $(SRC)/yoe_spectral_planck.F90

yoe_tile_prop.o: $(SRC)/yoe_tile_prop.F90 parkind1.o geometry_mod.o yomdphy.o
	$(FC) -c $(SRC)/yoe_tile_prop.F90

yoe_uvrad.o: $(SRC)/yoe_uvrad.F90 parkind1.o
	$(FC) -c $(SRC)/yoe_uvrad.F90

yoeaeratm.o: $(SRC)/yoeaeratm.F90 parkind1.o yoe_aerodiag.o yomhook.o
	$(FC) -c $(SRC)/yoeaeratm.F90

yoeaerc.o: $(SRC)/yoeaerc.F90 parkind1.o yomhook.o yomlun.o easy_netcdf_read_mpi.o yommp0.o abor1.intfb.ok
	$(FC) -c $(SRC)/yoeaerc.F90

yoeaerd.o: $(SRC)/yoeaerd.F90 parkind1.o
	$(FC) -c $(SRC)/yoeaerd.F90

yoeaerlid.o: $(SRC)/yoeaerlid.F90 parkind1.o
	$(FC) -c $(SRC)/yoeaerlid.F90

yoeaermap.o: $(SRC)/yoeaermap.F90 parkind1.o
	$(FC) -c $(SRC)/yoeaermap.F90

yoeaerop.o: $(SRC)/yoeaerop.F90 parkind1.o
	$(FC) -c $(SRC)/yoeaerop.F90

yoeaersnk.o: $(SRC)/yoeaersnk.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/yoeaersnk.F90

yoeaersrc.o: $(SRC)/yoeaersrc.F90 parkind1.o
	$(FC) -c $(SRC)/yoeaersrc.F90

yoeaervol.o: $(SRC)/yoeaervol.F90 parkind1.o
	$(FC) -c $(SRC)/yoeaervol.F90

yoecld.o: $(SRC)/yoecld.F90 parkind1.o
	$(FC) -c $(SRC)/yoecld.F90

yoecldp.o: $(SRC)/yoecldp.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/yoecldp.F90

yoecmip.o: $(SRC)/yoecmip.F90 parkind1.o
	$(FC) -c $(SRC)/yoecmip.F90

yoecnd.o: $(SRC)/yoecnd.F90 parkind1.o
	$(FC) -c $(SRC)/yoecnd.F90

yoecumf.o: $(SRC)/yoecumf.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/yoecumf.F90

yoecumf2.o: $(SRC)/yoecumf2.F90 parkind1.o
	$(FC) -c $(SRC)/yoecumf2.F90

yoedbug.o: $(SRC)/yoedbug.F90 parkind1.o
	$(FC) -c $(SRC)/yoedbug.F90

yoegwd.o: $(SRC)/yoegwd.F90 parkind1.o
	$(FC) -c $(SRC)/yoegwd.F90

yoegwdwms.o: $(SRC)/yoegwdwms.F90 parkind1.o
	$(FC) -c $(SRC)/yoegwdwms.F90

yoegwwms.o: $(SRC)/yoegwwms.F90 parkind1.o
	$(FC) -c $(SRC)/yoegwwms.F90

yoelwrad.o: $(SRC)/yoelwrad.F90 parkind1.o
	$(FC) -c $(SRC)/yoelwrad.F90

yoeneur.o: $(SRC)/yoeneur.F90 parkind1.o
	$(FC) -c $(SRC)/yoeneur.F90

yoeovlp.o: $(SRC)/yoeovlp.F90 parkind1.o
	$(FC) -c $(SRC)/yoeovlp.F90

yoeozoc.o: $(SRC)/yoeozoc.F90 parkind1.o
	$(FC) -c $(SRC)/yoeozoc.F90

yoephli.o: $(SRC)/yoephli.F90 parkind1.o
	$(FC) -c $(SRC)/yoephli.F90

yoephy.o: $(SRC)/yoephy.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/yoephy.F90

yoerad.o: $(SRC)/yoerad.F90 parkind1.o yoe_spectral_planck.o
	$(FC) -c $(SRC)/yoerad.F90

yoeradghg.o: $(SRC)/yoeradghg.F90 parkind1.o yomhook.o yomcst.o
	$(FC) -c $(SRC)/yoeradghg.F90

yoerdi.o: $(SRC)/yoerdi.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/yoerdi.F90

yoerip.o: $(SRC)/yoerip.F90 parkind1.o
	$(FC) -c $(SRC)/yoerip.F90

yoerrtftr.o: $(SRC)/yoerrtftr.F90 parkind1.o
	$(FC) -c $(SRC)/yoerrtftr.F90

yoerrtm.o: $(SRC)/yoerrtm.F90 parkind1.o parrrtm.o
	$(FC) -c $(SRC)/yoerrtm.F90

yoerrtwn.o: $(SRC)/yoerrtwn.F90 parkind1.o
	$(FC) -c $(SRC)/yoerrtwn.F90

yoesrtm.o: $(SRC)/yoesrtm.F90 parkind1.o parsrtm.o
	$(FC) -c $(SRC)/yoesrtm.F90

yoesw.o: $(SRC)/yoesw.F90 parkind1.o
	$(FC) -c $(SRC)/yoesw.F90

yoethf.o: $(SRC)/yoethf.F90 parkind1.o
	$(FC) -c $(SRC)/yoethf.F90

yoevdf.o: $(SRC)/yoevdf.F90 parkind1.o
	$(FC) -c $(SRC)/yoevdf.F90

yoewcou.o: $(SRC)/yoewcou.F90 parkind1.o
	$(FC) -c $(SRC)/yoewcou.F90

yom_grib_codes.o: $(SRC)/yom_grib_codes.F90 parkind1.o par_gfl.o
	$(FC) -c $(SRC)/yom_grib_codes.F90

yom_grid_biconserv.o: $(SRC)/yom_grid_biconserv.F90 parkind1.o
	$(FC) -c $(SRC)/yom_grid_biconserv.F90

yom_ygfl.o: $(SRC)/yom_ygfl.F90 parkind1.o yoe_aerodiag.o par_gfl.o
	$(FC) -c $(SRC)/yom_ygfl.F90

yomafn.o: $(SRC)/yomafn.F90 parkind1.o parfpos.o fullpos_mix.o type_fpdsphys.o
	$(FC) -c $(SRC)/yomafn.F90

yomarphy.o: $(SRC)/yomarphy.F90 parkind1.o
	$(FC) -c $(SRC)/yomarphy.F90

yomcape.o: $(SRC)/yomcape.F90 parkind1.o
	$(FC) -c $(SRC)/yomcape.F90

yomcddh.o: $(SRC)/yomcddh.F90 
	$(FC) -c $(SRC)/yomcddh.F90

yomcfu.o: $(SRC)/yomcfu.F90 yomcfu_type.o
	$(FC) -c $(SRC)/yomcfu.F90

yomcfu_type.o: $(SRC)/yomcfu_type.F90 parkind1.o yomhook.o type_fluxes.o ptrgfu.o field_module.o field_factory_module.o
	$(FC) -c $(SRC)/yomcfu_type.F90

yomchem.o: $(SRC)/yomchem.F90 parkind1.o
	$(FC) -c $(SRC)/yomchem.F90

yomcli.o: $(SRC)/yomcli.F90 parkind1.o
	$(FC) -c $(SRC)/yomcli.F90

yomcom.o: $(SRC)/yomcom.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/yomcom.F90

yomcompo.o: $(SRC)/yomcompo.F90 parkind1.o yomlun.o
	$(FC) -c $(SRC)/yomcompo.F90

yomcosjb.o: $(SRC)/yomcosjb.F90 parkind1.o
	$(FC) -c $(SRC)/yomcosjb.F90

yomcou.o: $(SRC)/yomcou.F90 parkind1.o
	$(FC) -c $(SRC)/yomcou.F90

yomcsgeom.o: $(SRC)/yomcsgeom.F90 parkind1.o yomhook.o deallocate_if_associated_mod.o
	$(FC) -c $(SRC)/yomcsgeom.F90

yomcst.o: $(SRC)/yomcst.F90 parkind1.o yomhook.o yomlun.o yomdyncore.o
	$(FC) -c $(SRC)/yomcst.F90

yomct0.o: $(SRC)/yomct0.F90 parkind1.o
	$(FC) -c $(SRC)/yomct0.F90

yomct1.o: $(SRC)/yomct1.F90 parkind1.o
	$(FC) -c $(SRC)/yomct1.F90

yomct3.o: $(SRC)/yomct3.F90 parkind1.o
	$(FC) -c $(SRC)/yomct3.F90

yomcumfs.o: $(SRC)/yomcumfs.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/yomcumfs.F90

yomcver.o: $(SRC)/yomcver.F90 parkind1.o yomhook.o yomlun.o yomct0.o abor1.intfb.ok
	$(FC) -c $(SRC)/yomcver.F90

yomcvmnh.o: $(SRC)/yomcvmnh.F90 parkind1.o modd_nsv.o modd_convpar.o modd_cst.o modd_convparext.o modd_convpar_shal.o
	$(FC) -c $(SRC)/yomcvmnh.F90

yomdfi.o: $(SRC)/yomdfi.F90 parkind1.o
	$(FC) -c $(SRC)/yomdfi.F90

yomdgradient.o: $(SRC)/yomdgradient.F90 yomhook.o
	$(FC) -c $(SRC)/yomdgradient.F90

yomdiagvar.o: $(SRC)/yomdiagvar.F90 parkind1.o
	$(FC) -c $(SRC)/yomdiagvar.F90

yomdim.o: $(SRC)/yomdim.F90 parkind1.o
	$(FC) -c $(SRC)/yomdim.F90

yomdimf.o: $(SRC)/yomdimf.F90 parkind1.o
	$(FC) -c $(SRC)/yomdimf.F90

yomdimv.o: $(SRC)/yomdimv.F90 parkind1.o
	$(FC) -c $(SRC)/yomdimv.F90

yomdphy.o: $(SRC)/yomdphy.F90 parkind1.o
	$(FC) -c $(SRC)/yomdphy.F90

yomdprecips.o: $(SRC)/yomdprecips.F90 parkind1.o
	$(FC) -c $(SRC)/yomdprecips.F90

yomdvisi.o: $(SRC)/yomdvisi.F90 parkind1.o
	$(FC) -c $(SRC)/yomdvisi.F90

yomdwet.o: $(SRC)/yomdwet.F90 parkind1.o
	$(FC) -c $(SRC)/yomdwet.F90

yomdyn.o: $(SRC)/yomdyn.F90 parkind1.o yomswe.o yomhook.o
	$(FC) -c $(SRC)/yomdyn.F90

yomdyna.o: $(SRC)/yomdyna.F90 parkind1.o intdyn_mod.o
	$(FC) -c $(SRC)/yomdyna.F90

yomdyncore.o: $(SRC)/yomdyncore.F90 parkind1.o
	$(FC) -c $(SRC)/yomdyncore.F90

yomfp4l.o: $(SRC)/yomfp4l.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/yomfp4l.F90

yomfpc.o: $(SRC)/yomfpc.F90 parkind1.o yomct0.o types_fpcat.o parfpos.o
	$(FC) -c $(SRC)/yomfpc.F90

yomfpcnt.o: $(SRC)/yomfpcnt.F90 parkind1.o
	$(FC) -c $(SRC)/yomfpcnt.F90

yomfpfilters.o: $(SRC)/yomfpfilters.F90 parkind1.o
	$(FC) -c $(SRC)/yomfpfilters.F90

yomfpgeo.o: $(SRC)/yomfpgeo.F90 parkind1.o
	$(FC) -c $(SRC)/yomfpgeo.F90

yomfpgeometry.o: $(SRC)/yomfpgeometry.F90 parkind1.o type_fpusergeo.o yomfpgeo.o yomfpgind.o abor1.intfb.ok
	$(FC) -c $(SRC)/yomfpgeometry.F90

yomfpgind.o: $(SRC)/yomfpgind.F90 parkind1.o
	$(FC) -c $(SRC)/yomfpgind.F90

yomfpios.o: $(SRC)/yomfpios.F90 parkind1.o yomct0.o
	$(FC) -c $(SRC)/yomfpios.F90

yomfpop.o: $(SRC)/yomfpop.F90 yomfpios.o type_faoph.o type_fpofn.o
	$(FC) -c $(SRC)/yomfpop.F90

yomgem.o: $(SRC)/yomgem.F90 parkind1.o
	$(FC) -c $(SRC)/yomgem.F90

yomgfl.o: $(SRC)/yomgfl.F90 parkind1.o yom_ygfl.o yomhook.o geometry_mod.o abor1.intfb.ok
	$(FC) -c $(SRC)/yomgfl.F90

yomgmv.o: $(SRC)/yomgmv.F90 parkind1.o type_gmvs.o yomlun.o yomhook.o geometry_mod.o abor1.intfb.ok
	$(FC) -c $(SRC)/yomgmv.F90

yomgpddh.o: $(SRC)/yomgpddh.F90 parkind1.o
	$(FC) -c $(SRC)/yomgpddh.F90

yomgrib.o: $(SRC)/yomgrib.F90 parkind1.o
	$(FC) -c $(SRC)/yomgrib.F90

yomgsgeom.o: $(SRC)/yomgsgeom.F90 parkind1.o yomhook.o deallocate_if_associated_mod.o
	$(FC) -c $(SRC)/yomgsgeom.F90

yomhook.o: $(SRC)/yomhook.F90 ec_parkind.o
	$(FC) -c $(SRC)/yomhook.F90

yomhslmer.o: $(SRC)/yomhslmer.F90 parkind1.o
	$(FC) -c $(SRC)/yomhslmer.F90

yomini.o: $(SRC)/yomini.F90 
	$(FC) -c $(SRC)/yomini.F90

yomjq.o: $(SRC)/yomjq.F90 parkind1.o spectral_fields_mod.o
	$(FC) -c $(SRC)/yomjq.F90

yomlap.o: $(SRC)/yomlap.F90 parkind1.o
	$(FC) -c $(SRC)/yomlap.F90

yomlddh.o: $(SRC)/yomlddh.F90 
	$(FC) -c $(SRC)/yomlddh.F90

yomleg.o: $(SRC)/yomleg.F90 parkind1.o
	$(FC) -c $(SRC)/yomleg.F90

yomlouis.o: $(SRC)/yomlouis.F90 parkind1.o
	$(FC) -c $(SRC)/yomlouis.F90

yomlsforc.o: $(SRC)/yomlsforc.F90 parkind1.o
	$(FC) -c $(SRC)/yomlsforc.F90

yomlun.o: $(SRC)/yomlun.F90 parkind1.o yomlun_ifsaux.o yomhook.o
	$(FC) -c $(SRC)/yomlun.F90

yomlun_ifsaux.o: $(SRC)/yomlun_ifsaux.F90 ec_lun.o
	$(FC) -c $(SRC)/yomlun_ifsaux.F90

yommcc.o: $(SRC)/yommcc.F90 parkind1.o cplng_types_mod.o par_gfl.o yomcompo.o yomhook.o
	$(FC) -c $(SRC)/yommcc.F90

yommcuf.o: $(SRC)/yommcuf.F90 parkind1.o yomhook.o abor1.intfb.ok
	$(FC) -c $(SRC)/yommcuf.F90

yommddh.o: $(SRC)/yommddh.F90 parkind1.o
	$(FC) -c $(SRC)/yommddh.F90

yommoderr.o: $(SRC)/yommoderr.F90 parkind1.o gridpoint_fields_mix.o spectral_fields_mod.o yommoderrconf.o
	$(FC) -c $(SRC)/yommoderr.F90

yommoderrconf.o: $(SRC)/yommoderrconf.F90 parkind1.o
	$(FC) -c $(SRC)/yommoderrconf.F90

yommoderrmod.o: $(SRC)/yommoderrmod.F90 parkind1.o
	$(FC) -c $(SRC)/yommoderrmod.F90

yommp.o: $(SRC)/yommp.F90 parkind1.o
	$(FC) -c $(SRC)/yommp.F90

yommp0.o: $(SRC)/yommp0.F90 parkind1.o
	$(FC) -c $(SRC)/yommp0.F90

yommse.o: $(SRC)/yommse.F90 parkind1.o yomct0.o
	$(FC) -c $(SRC)/yommse.F90

yomncl.o: $(SRC)/yomncl.F90 parkind1.o
	$(FC) -c $(SRC)/yomncl.F90

yomnorgwd.o: $(SRC)/yomnorgwd.F90 parkind1.o
	$(FC) -c $(SRC)/yomnorgwd.F90

yomnrtaer.o: $(SRC)/yomnrtaer.F90 parkind1.o
	$(FC) -c $(SRC)/yomnrtaer.F90

yomnud.o: $(SRC)/yomnud.F90 parkind1.o pardim.o
	$(FC) -c $(SRC)/yomnud.F90

yomorog.o: $(SRC)/yomorog.F90 parkind1.o
	$(FC) -c $(SRC)/yomorog.F90

yomozo.o: $(SRC)/yomozo.F90 parkind1.o
	$(FC) -c $(SRC)/yomozo.F90

yompaddh.o: $(SRC)/yompaddh.F90 parkind1.o
	$(FC) -c $(SRC)/yompaddh.F90

yomparar.o: $(SRC)/yomparar.F90 parkind1.o modd_phyex.o
	$(FC) -c $(SRC)/yomparar.F90

yompertpar.o: $(SRC)/yompertpar.F90 parkind1.o
	$(FC) -c $(SRC)/yompertpar.F90

yomphy.o: $(SRC)/yomphy.F90 parkind1.o yomdprecips.o yomdvisi.o yomcape.o yomdwet.o
	$(FC) -c $(SRC)/yomphy.F90

yomphy0.o: $(SRC)/yomphy0.F90 parkind1.o
	$(FC) -c $(SRC)/yomphy0.F90

yomphy1.o: $(SRC)/yomphy1.F90 parkind1.o
	$(FC) -c $(SRC)/yomphy1.F90

yomphy2.o: $(SRC)/yomphy2.F90 parkind1.o
	$(FC) -c $(SRC)/yomphy2.F90

yomphy3.o: $(SRC)/yomphy3.F90 parkind1.o
	$(FC) -c $(SRC)/yomphy3.F90

yomphyds.o: $(SRC)/yomphyds.F90 parkind1.o
	$(FC) -c $(SRC)/yomphyds.F90

yomprad.o: $(SRC)/yomprad.F90 parkind1.o
	$(FC) -c $(SRC)/yomprad.F90

yomradf.o: $(SRC)/yomradf.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/yomradf.F90

yomrandom_streams.o: $(SRC)/yomrandom_streams.F90 random_numbers_mix.o
	$(FC) -c $(SRC)/yomrandom_streams.F90

yomrcoef.o: $(SRC)/yomrcoef.F90 parkind1.o
	$(FC) -c $(SRC)/yomrcoef.F90

yomrip.o: $(SRC)/yomrip.F90 parkind1.o yoeozoc.o yoecmip.o yoeradghg.o yoeaerc.o reglatlon_field_mix.o
	$(FC) -c $(SRC)/yomrip.F90

yomrip0.o: $(SRC)/yomrip0.F90 parkind1.o
	$(FC) -c $(SRC)/yomrip0.F90

yomsddh.o: $(SRC)/yomsddh.F90 parkind1.o
	$(FC) -c $(SRC)/yomsddh.F90

yomsens.o: $(SRC)/yomsens.F90 parkind1.o
	$(FC) -c $(SRC)/yomsens.F90

yomsimphl.o: $(SRC)/yomsimphl.F90 parkind1.o
	$(FC) -c $(SRC)/yomsimphl.F90

yomslint.o: $(SRC)/yomslint.F90 parkind1.o yomvsplip.o yomvsleta.o yomhslmer.o
	$(FC) -c $(SRC)/yomslint.F90

yomslphy.o: $(SRC)/yomslphy.F90 parkind1.o yomhook.o
	$(FC) -c $(SRC)/yomslphy.F90

yomslrep.o: $(SRC)/yomslrep.F90 parkind1.o
	$(FC) -c $(SRC)/yomslrep.F90

yomsnu.o: $(SRC)/yomsnu.F90 parkind1.o
	$(FC) -c $(SRC)/yomsnu.F90

yomspddh.o: $(SRC)/yomspddh.F90 parkind1.o
	$(FC) -c $(SRC)/yomspddh.F90

yomsphyhist.o: $(SRC)/yomsphyhist.F90 parkind1.o
	$(FC) -c $(SRC)/yomsphyhist.F90

yomspsdt.o: $(SRC)/yomspsdt.F90 parkind1.o spectral_arp_mod.o gridpoint_fields_mix.o
	$(FC) -c $(SRC)/yomspsdt.F90

yomsrftlad.o: $(SRC)/yomsrftlad.F90 parkind1.o
	$(FC) -c $(SRC)/yomsrftlad.F90

yomsta.o: $(SRC)/yomsta.F90 parkind1.o
	$(FC) -c $(SRC)/yomsta.F90

yomswe.o: $(SRC)/yomswe.F90 parkind1.o
	$(FC) -c $(SRC)/yomswe.F90

yomtag.o: $(SRC)/yomtag.F90 parkind1.o
	$(FC) -c $(SRC)/yomtag.F90

yomtddh.o: $(SRC)/yomtddh.F90 parkind1.o
	$(FC) -c $(SRC)/yomtddh.F90

yomtnh.o: $(SRC)/yomtnh.F90 parkind1.o
	$(FC) -c $(SRC)/yomtnh.F90

yomtoph.o: $(SRC)/yomtoph.F90 parkind1.o
	$(FC) -c $(SRC)/yomtoph.F90

yomtraj.o: $(SRC)/yomtraj.F90 parkind1.o spectral_fields_mod.o
	$(FC) -c $(SRC)/yomtraj.F90

yomtrc.o: $(SRC)/yomtrc.F90 parkind1.o
	$(FC) -c $(SRC)/yomtrc.F90

yomvar.o: $(SRC)/yomvar.F90 parkind1.o yomhook.o yomlun.o yomct0.o yommp0.o yomdiagvar.o yomvcgl.o yommoderr.o yomjq.o yomsens.o yomvrtl.o yomcosjb.o algorithm_state_mod.o yomvrtlx.o abor1.intfb.ok
	$(FC) -c $(SRC)/yomvar.F90

yomvcgl.o: $(SRC)/yomvcgl.F90 parkind1.o control_vectors_base_mix.o
	$(FC) -c $(SRC)/yomvcgl.F90

yomvdoz.o: $(SRC)/yomvdoz.F90 parkind1.o
	$(FC) -c $(SRC)/yomvdoz.F90

yomvert.o: $(SRC)/yomvert.F90 parkind1.o yomhook.o yomcver.o
	$(FC) -c $(SRC)/yomvert.F90

yomvrtl.o: $(SRC)/yomvrtl.F90 
	$(FC) -c $(SRC)/yomvrtl.F90

yomvrtlx.o: $(SRC)/yomvrtlx.F90 parkind1.o
	$(FC) -c $(SRC)/yomvrtlx.F90

yomvsleta.o: $(SRC)/yomvsleta.F90 parkind1.o
	$(FC) -c $(SRC)/yomvsleta.F90

yomvsplip.o: $(SRC)/yomvsplip.F90 parkind1.o
	$(FC) -c $(SRC)/yomvsplip.F90

yomwfpb.o: $(SRC)/yomwfpb.F90 parkind1.o
	$(FC) -c $(SRC)/yomwfpb.F90

yomxfu.o: $(SRC)/yomxfu.F90 yomxfu_type.o
	$(FC) -c $(SRC)/yomxfu.F90

yomxfu_type.o: $(SRC)/yomxfu_type.F90 parkind1.o yomhook.o type_fluxes.o ptrxfu.o field_module.o field_factory_module.o yomct1.o yomct0.o yomct3.o yomrip.o yomphy.o yommcc.o
	$(FC) -c $(SRC)/yomxfu_type.F90

yophlc.o: $(SRC)/yophlc.F90 parkind1.o
	$(FC) -c $(SRC)/yophlc.F90

yophnc.o: $(SRC)/yophnc.F90 
	$(FC) -c $(SRC)/yophnc.F90


subclean:
	\rm -f apl_arpege.o

clean: 
	\rm -f *.o *.xml *.a *.x *.mod *.optrpt 

tidy:
	\rm -f *.xml *.optrpt
