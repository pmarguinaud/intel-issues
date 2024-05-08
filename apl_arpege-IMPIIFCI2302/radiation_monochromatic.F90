! radiation_interface.F90 - Monochromatic gas/cloud optics for testing
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
!   2017-04-11  R. Hogan  Receive "surface" dummy argument
!   2017-09-13  R. Hogan  Revert
!   2018-08-29  R. Hogan  Particulate single-scattering albedo / asymmetry from namelist

module radiation_monochromatic

  implicit none

  

contains

  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Setup the arrays in the config object corresponding to the
  ! monochromatic gas optics model.  The directory argument is not
  ! used, since no look-up tables need to be loaded.
  


  !---------------------------------------------------------------------
  ! Dummy routine for scaling gas mixing ratios
  


  !---------------------------------------------------------------------
  ! Dummy setup routine for cloud optics: in fact, no setup is
  ! required for monochromatic case
  


  !---------------------------------------------------------------------
  ! Dummy subroutine since no aerosols are represented in
  ! monochromatic case
  


  !---------------------------------------------------------------------
  ! Compute gas optical depths, shortwave scattering, Planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  


  !---------------------------------------------------------------------
  ! Compute cloud optical depth, single-scattering albedo and
  ! g factor in the longwave and shortwave
  


  !---------------------------------------------------------------------
  ! Dummy subroutine since no aerosols are represented in
  ! monochromatic case
  

  !---------------------------------------------------------------------
  ! Planck function in terms of wavelength
  

  subroutine setup_gas_optics_mono_GPU(config, directory, lacc)
use radiation_config
type(config_type), intent(inout) :: config
character(len=*),  intent(in)    :: directory
logical, intent (in) :: lacc




















end subroutine setup_gas_optics_mono_GPU

  subroutine set_gas_units_mono_GPU(gas, lacc)
use radiation_gas
type(gas_type),    intent(inout) :: gas
logical, intent (in) :: lacc
end subroutine set_gas_units_mono_GPU

  subroutine setup_cloud_optics_mono_GPU(config, lacc)
use radiation_config
type(config_type), intent(inout) :: config
logical, intent (in) :: lacc
end subroutine setup_cloud_optics_mono_GPU

  subroutine setup_aerosol_optics_mono_GPU(config, lacc)
use radiation_config
type(config_type), intent(inout) :: config
logical, intent (in) :: lacc
end subroutine setup_aerosol_optics_mono_GPU

  subroutine gas_optics_mono_GPU(ncol,nlev,istartcol,iendcol, &
config, single_level, thermodynamics, gas, lw_albedo, &
od_lw, od_sw, ssa_sw, planck_hl, lw_emission, &
incoming_sw, lacc)
use parkind1
use radiation_config
use radiation_thermodynamics
use radiation_single_level
use radiation_gas
use radiation_constants

integer, intent(in) :: ncol               
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(gas_type),           intent(in) :: gas

real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&  intent(in) :: lw_albedo




real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
&   od_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
&   od_sw, ssa_sw


real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out) :: &
&   planck_hl
real(jprb), dimension(config%n_g_lw,istartcol:iendcol), intent(out) :: &
&   lw_emission



real(jprb), dimension(config%n_g_sw,istartcol:iendcol), intent(out) :: &
&   incoming_sw










logical, intent (in) :: lacc







end subroutine gas_optics_mono_GPU

  subroutine cloud_optics_mono_GPU(nlev,istartcol,iendcol, &
&   config, thermodynamics, cloud, &
&   od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
&   od_sw_cloud, ssa_sw_cloud, g_sw_cloud, lacc)
use parkind1
use radiation_config
use radiation_thermodynamics
use radiation_cloud
use radiation_constants

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in) :: config
type(thermodynamics_type),intent(in) :: thermodynamics
type(cloud_type),   intent(in) :: cloud





real(jprb), dimension(config%n_bands_lw,nlev,istartcol:iendcol), intent(out) :: &
&   od_lw_cloud
real(jprb), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol), &
&   intent(out) :: ssa_lw_cloud, g_lw_cloud


real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
&   od_sw_cloud, ssa_sw_cloud, g_sw_cloud






logical, intent (in) :: lacc













end subroutine cloud_optics_mono_GPU

  subroutine add_aerosol_optics_mono_GPU(nlev,istartcol,iendcol, &
&  config, thermodynamics, gas, aerosol, &
&  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw, lacc)
use parkind1
use radiation_config
use radiation_thermodynamics
use radiation_gas
use radiation_aerosol
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in), target :: config
type(thermodynamics_type),intent(in)  :: thermodynamics
type(gas_type),           intent(in)  :: gas
type(aerosol_type),       intent(in)  :: aerosol





real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(inout) :: od_lw
real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
&  intent(out) :: ssa_lw, g_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(inout) &
&  :: od_sw, ssa_sw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: g_sw
logical, intent (in) :: lacc


end subroutine add_aerosol_optics_mono_GPU

  elemental function planck_function_GPU(wavelength, temperature, lacc)
use parkind1
use radiation_constants
real(jprb), intent(in) :: wavelength  
real(jprb), intent(in) :: temperature 

real(jprb)             :: planck_function_GPU
logical, intent (in) :: lacc

end function planck_function_GPU

  subroutine setup_gas_optics_mono_CPU(config, directory)
use radiation_config
type(config_type), intent(inout) :: config
character(len=*),  intent(in)    :: directory




















end subroutine setup_gas_optics_mono_CPU

  subroutine set_gas_units_mono_CPU(gas)
use radiation_gas
type(gas_type),    intent(inout) :: gas
end subroutine set_gas_units_mono_CPU

  subroutine setup_cloud_optics_mono_CPU(config)
use radiation_config
type(config_type), intent(inout) :: config
end subroutine setup_cloud_optics_mono_CPU

  subroutine setup_aerosol_optics_mono_CPU(config)
use radiation_config
type(config_type), intent(inout) :: config
end subroutine setup_aerosol_optics_mono_CPU

  subroutine gas_optics_mono_CPU(ncol,nlev,istartcol,iendcol, &
config, single_level, thermodynamics, gas, lw_albedo, &
od_lw, od_sw, ssa_sw, planck_hl, lw_emission, &
incoming_sw)
use parkind1
use radiation_config
use radiation_thermodynamics
use radiation_single_level
use radiation_gas
use radiation_constants

integer, intent(in) :: ncol               
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(gas_type),           intent(in) :: gas

real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&  intent(in) :: lw_albedo




real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
&   od_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
&   od_sw, ssa_sw


real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out) :: &
&   planck_hl
real(jprb), dimension(config%n_g_lw,istartcol:iendcol), intent(out) :: &
&   lw_emission



real(jprb), dimension(config%n_g_sw,istartcol:iendcol), intent(out) :: &
&   incoming_sw

















end subroutine gas_optics_mono_CPU

  subroutine cloud_optics_mono_CPU(nlev,istartcol,iendcol, &
&   config, thermodynamics, cloud, &
&   od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
&   od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
use parkind1
use radiation_config
use radiation_thermodynamics
use radiation_cloud
use radiation_constants

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in) :: config
type(thermodynamics_type),intent(in) :: thermodynamics
type(cloud_type),   intent(in) :: cloud





real(jprb), dimension(config%n_bands_lw,nlev,istartcol:iendcol), intent(out) :: &
&   od_lw_cloud
real(jprb), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol), &
&   intent(out) :: ssa_lw_cloud, g_lw_cloud


real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
&   od_sw_cloud, ssa_sw_cloud, g_sw_cloud



















end subroutine cloud_optics_mono_CPU

  subroutine add_aerosol_optics_mono_CPU(nlev,istartcol,iendcol, &
&  config, thermodynamics, gas, aerosol, &
&  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)
use parkind1
use radiation_config
use radiation_thermodynamics
use radiation_gas
use radiation_aerosol
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in), target :: config
type(thermodynamics_type),intent(in)  :: thermodynamics
type(gas_type),           intent(in)  :: gas
type(aerosol_type),       intent(in)  :: aerosol





real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(inout) :: od_lw
real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
&  intent(out) :: ssa_lw, g_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(inout) &
&  :: od_sw, ssa_sw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: g_sw


end subroutine add_aerosol_optics_mono_CPU

  elemental function planck_function_CPU(wavelength, temperature)
use parkind1
use radiation_constants
real(jprb), intent(in) :: wavelength  
real(jprb), intent(in) :: temperature 

real(jprb)             :: planck_function_CPU

end function planck_function_CPU

end module radiation_monochromatic

