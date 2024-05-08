! radiation_save.F90 - Save data to NetCDF files
!
! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Modifications
!   2017-04-22  R. Hogan  Adapt for new way of describing longwave properties
!   2019-01-02  R. Hogan  Only save cloud properties if do_clouds==.true.

module radiation_save

  use parkind1

  implicit none

  ! Save final fluxes and save intermediate radiative properties
  
  ! Save net fluxes IFS style, where upwelling fluxes are actually net down
  

contains

  !---------------------------------------------------------------------
  ! Save fluxes in "flux" to NetCDF file_name, plus pressure from the
  ! thermodynamics object
  


  !---------------------------------------------------------------------
  ! Save IFS-style net fluxes in "flux" to NetCDF file_name, plus
  ! pressure from the thermodynamics object
  


  !---------------------------------------------------------------------
  ! Save intermediate radiative properties, specifically the
  ! scattering and absorption properties at each g-point/band
  


  !---------------------------------------------------------------------
  ! Save inputs to the radiation scheme
  


  !---------------------------------------------------------------------
  ! Save shortwave diagnostics computed from "flux" to NetCDF
  ! file_name.  The "mapping" matrix maps from fluxes in bands or
  ! g-points to user specified spectral intervals, and should have
  ! been produced by config%get_sw_mapping. See the example in
  ! ecrad_driver.F90.
  

  subroutine save_fluxes_GPU(file_name, config, thermodynamics, flux, &
&                 iverbose, is_hdf5_file, experiment_name, &
&                 is_double_precision, lacc)
use yomhook
use easy_netcdf
use radiation_io
use radiation_config
use radiation_thermodynamics
use radiation_flux
character(len=*),           intent(in) :: file_name
type(config_type),          intent(in) :: config
type(thermodynamics_type),  intent(in) :: thermodynamics
type(flux_type),            intent(in) :: flux
integer,          optional, intent(in) :: iverbose
logical,          optional, intent(in) :: is_hdf5_file
logical,          optional, intent(in) :: is_double_precision
character(len=*), optional, intent(in) :: experiment_name






logical, intent (in) :: lacc









































end subroutine save_fluxes_GPU

  subroutine save_net_fluxes_GPU(file_name, config, thermodynamics, flux, &
&                     iverbose, is_hdf5_file, experiment_name, &
&                     is_double_precision, lacc)
use yomhook
use easy_netcdf
use radiation_io
use radiation_config
use radiation_thermodynamics
use radiation_flux
character(len=*),           intent(in) :: file_name
type(config_type),          intent(in) :: config
type(thermodynamics_type),  intent(in) :: thermodynamics
type(flux_type),            intent(in) :: flux
integer,          optional, intent(in) :: iverbose
logical,          optional, intent(in) :: is_hdf5_file
logical,          optional, intent(in) :: is_double_precision
character(len=*), optional, intent(in) :: experiment_name






logical, intent (in) :: lacc




































end subroutine save_net_fluxes_GPU

  subroutine save_radiative_properties_GPU(file_name, nlev, &
&  istartcol, iendcol, &
&  config, single_level, thermodynamics, cloud, &
&  planck_hl, lw_emission, lw_albedo, &
&  sw_albedo_direct, sw_albedo_diffuse, &
&  incoming_sw, &
&  od_lw, ssa_lw, g_lw, &
&  od_sw, ssa_sw, g_sw, &
&  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
&  od_sw_cloud, ssa_sw_cloud, g_sw_cloud, lacc)
use radiation_config
use radiation_single_level
use radiation_thermodynamics
use radiation_cloud
use easy_netcdf
character(len=*),         intent(in) :: file_name
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(cloud_type),         intent(in) :: cloud
integer, intent(in) :: nlev, istartcol, iendcol



real(jprb), intent(in), dimension(config%n_g_sw,nlev,istartcol:iendcol) :: od_sw, ssa_sw, g_sw


real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
&  od_sw_cloud, ssa_sw_cloud, g_sw_cloud



real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) &
&  :: sw_albedo_direct, sw_albedo_diffuse, incoming_sw




real(jprb), intent(in), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od_lw
real(jprb), intent(in), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: &
&  ssa_lw, g_lw




real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol) :: od_lw_cloud
real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol) :: &
&  ssa_lw_cloud, g_lw_cloud


real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: planck_hl


real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: lw_emission, lw_albedo

 


logical, intent (in) :: lacc













 












 













end subroutine save_radiative_properties_GPU

  subroutine save_inputs_GPU(file_name, config, single_level, thermodynamics, &
&                 gas, cloud, aerosol, lat, lon, iverbose, lacc)
use yomhook
use radiation_config
use radiation_single_level
use radiation_thermodynamics
use radiation_gas
use radiation_cloud
use radiation_aerosol
use easy_netcdf
character(len=*),             intent(in)   :: file_name
type(config_type),            intent(in)   :: config
type(single_level_type),      intent(in)   :: single_level
type(thermodynamics_type),    intent(in)   :: thermodynamics
type(gas_type),               intent(inout):: gas
type(cloud_type),             intent(in)   :: cloud
type(aerosol_type), optional, intent(in)   :: aerosol
real(jprb),         optional, intent(in)   :: lat(:), lon(:)
integer,            optional, intent(in)   :: iverbose











logical, intent (in) :: lacc
















 


























 



























end subroutine save_inputs_GPU

  subroutine save_sw_diagnostics_GPU(file_name, config, wavelength_bound, mapping, flux, &
&                         iverbose, is_hdf5_file, experiment_name, &
&                         is_double_precision, lacc)
use yomhook
use easy_netcdf
use radiation_io
use radiation_flux
use radiation_config
use radiation_matrix
character(len=*),           intent(in) :: file_name
type(config_type),          intent(in) :: config
real(jprb),                 intent(in) :: wavelength_bound(:) 
real(jprb),                 intent(in) :: mapping(:,:)
type(flux_type),            intent(in) :: flux
integer,          optional, intent(in) :: iverbose
logical,          optional, intent(in) :: is_hdf5_file
logical,          optional, intent(in) :: is_double_precision
character(len=*), optional, intent(in) :: experiment_name
 





logical, intent (in) :: lacc










































end subroutine save_sw_diagnostics_GPU

  subroutine save_fluxes_CPU(file_name, config, thermodynamics, flux, &
&                 iverbose, is_hdf5_file, experiment_name, &
&                 is_double_precision)
use yomhook
use easy_netcdf
use radiation_io
use radiation_config
use radiation_thermodynamics
use radiation_flux
character(len=*),           intent(in) :: file_name
type(config_type),          intent(in) :: config
type(thermodynamics_type),  intent(in) :: thermodynamics
type(flux_type),            intent(in) :: flux
integer,          optional, intent(in) :: iverbose
logical,          optional, intent(in) :: is_hdf5_file
logical,          optional, intent(in) :: is_double_precision
character(len=*), optional, intent(in) :: experiment_name















































end subroutine save_fluxes_CPU

  subroutine save_net_fluxes_CPU(file_name, config, thermodynamics, flux, &
&                     iverbose, is_hdf5_file, experiment_name, &
&                     is_double_precision)
use yomhook
use easy_netcdf
use radiation_io
use radiation_config
use radiation_thermodynamics
use radiation_flux
character(len=*),           intent(in) :: file_name
type(config_type),          intent(in) :: config
type(thermodynamics_type),  intent(in) :: thermodynamics
type(flux_type),            intent(in) :: flux
integer,          optional, intent(in) :: iverbose
logical,          optional, intent(in) :: is_hdf5_file
logical,          optional, intent(in) :: is_double_precision
character(len=*), optional, intent(in) :: experiment_name










































end subroutine save_net_fluxes_CPU

  subroutine save_radiative_properties_CPU(file_name, nlev, &
&  istartcol, iendcol, &
&  config, single_level, thermodynamics, cloud, &
&  planck_hl, lw_emission, lw_albedo, &
&  sw_albedo_direct, sw_albedo_diffuse, &
&  incoming_sw, &
&  od_lw, ssa_lw, g_lw, &
&  od_sw, ssa_sw, g_sw, &
&  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
&  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
use radiation_config
use radiation_single_level
use radiation_thermodynamics
use radiation_cloud
use easy_netcdf
character(len=*),         intent(in) :: file_name
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(cloud_type),         intent(in) :: cloud
integer, intent(in) :: nlev, istartcol, iendcol



real(jprb), intent(in), dimension(config%n_g_sw,nlev,istartcol:iendcol) :: od_sw, ssa_sw, g_sw


real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
&  od_sw_cloud, ssa_sw_cloud, g_sw_cloud



real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) &
&  :: sw_albedo_direct, sw_albedo_diffuse, incoming_sw




real(jprb), intent(in), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od_lw
real(jprb), intent(in), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: &
&  ssa_lw, g_lw




real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol) :: od_lw_cloud
real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol) :: &
&  ssa_lw_cloud, g_lw_cloud


real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: planck_hl


real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: lw_emission, lw_albedo

 















 












 













end subroutine save_radiative_properties_CPU

  subroutine save_inputs_CPU(file_name, config, single_level, thermodynamics, &
&                 gas, cloud, aerosol, lat, lon, iverbose)
use yomhook
use radiation_config
use radiation_single_level
use radiation_thermodynamics
use radiation_gas
use radiation_cloud
use radiation_aerosol
use easy_netcdf
character(len=*),             intent(in)   :: file_name
type(config_type),            intent(in)   :: config
type(single_level_type),      intent(in)   :: single_level
type(thermodynamics_type),    intent(in)   :: thermodynamics
type(gas_type),               intent(inout):: gas
type(cloud_type),             intent(in)   :: cloud
type(aerosol_type), optional, intent(in)   :: aerosol
real(jprb),         optional, intent(in)   :: lat(:), lon(:)
integer,            optional, intent(in)   :: iverbose



























 


























 



























end subroutine save_inputs_CPU

  subroutine save_sw_diagnostics_CPU(file_name, config, wavelength_bound, mapping, flux, &
&                         iverbose, is_hdf5_file, experiment_name, &
&                         is_double_precision)
use yomhook
use easy_netcdf
use radiation_io
use radiation_flux
use radiation_config
use radiation_matrix
character(len=*),           intent(in) :: file_name
type(config_type),          intent(in) :: config
real(jprb),                 intent(in) :: wavelength_bound(:) 
real(jprb),                 intent(in) :: mapping(:,:)
type(flux_type),            intent(in) :: flux
integer,          optional, intent(in) :: iverbose
logical,          optional, intent(in) :: is_hdf5_file
logical,          optional, intent(in) :: is_double_precision
character(len=*), optional, intent(in) :: experiment_name
 















































end subroutine save_sw_diagnostics_CPU

end module radiation_save

