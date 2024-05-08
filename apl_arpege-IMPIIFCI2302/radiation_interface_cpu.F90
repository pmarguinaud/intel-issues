! radiation_interface.F90 - Public interface to radiation scheme
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
!   2017-04-11  R. Hogan  Changes to enable generalized surface description
!   2017-09-08  R. Hogan  Reverted some changes
!   2022-01-18  P. Ukkonen Added support for RRTMGP gas optics
!
! To use the radiation scheme, create a configuration_type object,
! call "setup_radiation" on it once to load the look-up-tables and
! data describing how gas and hydrometeor absorption/scattering are to
! be represented, and call "radiation" multiple times on different
! input profiles.

module radiation_interface_CPU

  implicit none

  
  

contains

  !---------------------------------------------------------------------
  ! Load the look-up-tables and data describing how gas and
  ! hydrometeor absorption/scattering are to be represented
  


  !---------------------------------------------------------------------
  ! Scale the gas mixing ratios so that they have the units (and
  ! possibly scale factors) required by the specific gas absorption
  ! model.  This subroutine simply passes the gas object on to the
  ! module of the currently active gas model.
  


  !---------------------------------------------------------------------
  ! Run the radiation scheme according to the configuration in the
  ! config object. There are ncol profiles of which only istartcol to
  ! iendcol are to be processed, and there are nlev model levels.  The
  ! output fluxes are written to the flux object, and all other
  ! objects contain the input variables.  The variables may be defined
  ! either in order of increasing or decreasing pressure, but if in
  ! order of decreasing pressure then radiation_reverse will be called
  ! to reverse the order for the computation and then reverse the
  ! order of the output fluxes to match the inputs.
  


  !---------------------------------------------------------------------
  ! If the input arrays are arranged in order of decreasing pressure /
  ! increasing height then this subroutine reverses them, calls the
  ! radiation scheme and then reverses the returned fluxes. Since this
  ! subroutine calls, and is called by "radiation", it must be in this
  ! module to avoid circular dependencies.
  

  subroutine setup_radiation_CPU(yderdi,config)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_spectral_definition


use radiation_monochromatic
use radiation_ifs_rrtm
use radiation_ecckd_interface
use radiation_ifs_rrtmgp
use radiation_cloud_optics
use radiation_general_cloud_optics
use radiation_aerosol_optics
use yoerdi,                   only : terdi
type(terdi), intent(inout) :: yderdi
type(config_type), intent(inout) :: config
































end subroutine setup_radiation_CPU

  subroutine set_gas_units_CPU(config, gas, lacc)
use radiation_config
use radiation_gas
use radiation_monochromatic
use radiation_ifs_rrtm
use radiation_ifs_rrtmgp
use radiation_ecckd_interface
type(config_type), intent(in)    :: config
type(gas_type),    intent(inout) :: gas
logical, optional, intent(in) :: lacc

end subroutine set_gas_units_CPU

  subroutine radiation_CPU(ncol, nlev, istartcol, iendcol, config, &
&  single_level, thermodynamics, gas, cloud, aerosol, flux)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_single_level
use radiation_thermodynamics
use radiation_gas
use radiation_cloud
use radiation_aerosol
use radiation_flux
use radiation_spartacus_sw
use radiation_spartacus_lw
use radiation_tripleclouds_sw
use radiation_tripleclouds_lw
use radiation_mcica_sw
use radiation_mcica_lw
use radiation_mcica_acc_sw
use radiation_mcica_acc_lw
use radiation_cloudless_sw
use radiation_cloudless_lw
use radiation_homogeneous_sw
use radiation_homogeneous_lw
use radiation_save

use radiation_monochromatic
use radiation_ifs_rrtm
use radiation_ecckd_interface
use radiation_ifs_rrtmgp
use radiation_cloud_optics
use radiation_general_cloud_optics
use radiation_aerosol_optics

integer, intent(in) :: ncol               
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in)   :: config
type(single_level_type),  intent(in)   :: single_level
type(thermodynamics_type),intent(in)   :: thermodynamics
type(gas_type),           intent(in)   :: gas
type(cloud_type),         intent(inout):: cloud
type(aerosol_type),       intent(in)   :: aerosol

type(flux_type),          intent(inout):: flux












































end subroutine radiation_CPU

  subroutine radiation_reverse_CPU(ncol, nlev, istartcol, iendcol, config, &
&  single_level, thermodynamics, gas, cloud, aerosol, flux)
use parkind1
use radiation_io
use radiation_config
use radiation_single_level
use radiation_thermodynamics
use radiation_gas
use radiation_cloud
use radiation_aerosol
use radiation_flux

integer, intent(in) :: ncol               
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(gas_type),           intent(in) :: gas
type(cloud_type),         intent(in) :: cloud
type(aerosol_type),       intent(in) :: aerosol

type(flux_type),          intent(inout):: flux
































end subroutine radiation_reverse_CPU

end module radiation_interface_CPU

