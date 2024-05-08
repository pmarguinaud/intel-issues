! radiation_cloudless_sw.F90 - Shortwave homogeneous cloudless solver
!
! (C) Copyright 2019- ECMWF.
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

module radiation_cloudless_sw



contains

  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Shortwave homogeneous solver containing no clouds
  

  subroutine solver_cloudless_sw_GPU(nlev,istartcol,iendcol, &
&  config, single_level, &
&  od, ssa, g, albedo_direct, albedo_diffuse, incoming_sw, &
&  flux, lacc)
use parkind1
use yomhook
use radiation_config
use radiation_single_level
use radiation_flux
use radiation_two_stream
use radiation_constants
use radiation_adding_ica_sw
implicit none

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level


real(jprb), intent(in), dimension(config%n_g_sw, nlev, istartcol:iendcol) :: &
&  od, ssa, g



real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) :: &
&  albedo_direct, albedo_diffuse, incoming_sw

type(flux_type), intent(inout):: flux

























logical, intent (in) :: lacc





end subroutine solver_cloudless_sw_GPU

  subroutine solver_cloudless_sw_CPU(nlev,istartcol,iendcol, &
&  config, single_level, &
&  od, ssa, g, albedo_direct, albedo_diffuse, incoming_sw, &
&  flux)
use parkind1
use yomhook
use radiation_config
use radiation_single_level
use radiation_flux
use radiation_two_stream
use radiation_constants
use radiation_adding_ica_sw
implicit none

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level


real(jprb), intent(in), dimension(config%n_g_sw, nlev, istartcol:iendcol) :: &
&  od, ssa, g



real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) :: &
&  albedo_direct, albedo_diffuse, incoming_sw

type(flux_type), intent(inout):: flux






























end subroutine solver_cloudless_sw_CPU

end module radiation_cloudless_sw

