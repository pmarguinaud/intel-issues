! radiation_ecckd_interface.F90 - Interface to ecCKD gas optics model
!
! (C) Copyright 2020- ECMWF.
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
! License: see the COPYING file for details
!

module radiation_ecckd_interface

  implicit none

   !, planck_function

contains

  !---------------------------------------------------------------------
  ! Setup the ecCKD generalized gas optics model
  


  !---------------------------------------------------------------------
  ! Scale gas mixing ratios according to required units
  


  !---------------------------------------------------------------------
  ! Compute gas optical depths, shortwave scattering, Planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  

  ! !---------------------------------------------------------------------
  ! ! Externally facing function for computing the Planck function
  ! ! without reference to any gas profile; typically this would be used
  ! ! for computing the emission by a surface.
  ! subroutine planck_function(config, temperature, planck_surf)

  !   use parkind1,                 only : jprb
  !   use radiation_config,         only : config_type

  !   type(config_type), intent(in) :: config
  !   real(jprb),        intent(in) :: temperature

  !   ! Planck function of the surface (W m-2)
  !   real(jprb), dimension(config%n_g_lw), intent(out) :: planck_surf

  ! end subroutine planck_function

  subroutine setup_gas_optics_ecckd_GPU(config, lacc)
use parkind1
use radiation_config
use yomhook
type(config_type), intent(inout), target :: config


logical, intent (in) :: lacc







end subroutine setup_gas_optics_ecckd_GPU

  subroutine set_gas_units_ecckd_GPU(gas, lacc)
use radiation_gas
use yomhook
type(gas_type),    intent(inout) :: gas
logical, intent(in), optional :: lacc




end subroutine set_gas_units_ecckd_GPU

  subroutine gas_optics_ecckd_GPU(ncol,nlev,istartcol,iendcol, &
&  config, single_level, thermodynamics, gas, &
&  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
&  incoming_sw, lacc)
use parkind1
use yomhook
use radiation_config
use radiation_thermodynamics
use radiation_single_level
use radiation_gas_constants
use radiation_gas
integer,                  intent(in) :: ncol               
integer,                  intent(in) :: nlev               
integer,                  intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(gas_type),           intent(in) :: gas

real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&  intent(in), optional :: lw_albedo



real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
&   od_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
&   od_sw, ssa_sw


real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), &
&   intent(out), optional :: planck_hl

real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&   intent(out), optional :: lw_emission



real(jprb), dimension(config%n_g_sw,istartcol:iendcol), &
&   intent(out), optional :: incoming_sw






logical, intent (in) :: lacc












end subroutine gas_optics_ecckd_GPU

  subroutine setup_gas_optics_ecckd_CPU(config)
use parkind1
use radiation_config
use yomhook
type(config_type), intent(inout), target :: config









end subroutine setup_gas_optics_ecckd_CPU

  subroutine set_gas_units_ecckd_CPU(gas, lacc)
use radiation_gas
use yomhook
type(gas_type),    intent(inout) :: gas
logical, intent(in), optional :: lacc




end subroutine set_gas_units_ecckd_CPU

  subroutine gas_optics_ecckd_CPU(ncol,nlev,istartcol,iendcol, &
&  config, single_level, thermodynamics, gas, &
&  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
&  incoming_sw)
use parkind1
use yomhook
use radiation_config
use radiation_thermodynamics
use radiation_single_level
use radiation_gas_constants
use radiation_gas
integer,                  intent(in) :: ncol               
integer,                  intent(in) :: nlev               
integer,                  intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(gas_type),           intent(in) :: gas

real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&  intent(in), optional :: lw_albedo



real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
&   od_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
&   od_sw, ssa_sw


real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), &
&   intent(out), optional :: planck_hl

real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&   intent(out), optional :: lw_emission



real(jprb), dimension(config%n_g_sw,istartcol:iendcol), &
&   intent(out), optional :: incoming_sw


















end subroutine gas_optics_ecckd_CPU

end module radiation_ecckd_interface

