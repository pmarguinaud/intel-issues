! radiation_homogeneous_sw.F90 - Shortwave homogeneous-column (no cloud fraction) solver
!
! (C) Copyright 2016- ECMWF.
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
!   2017-04-11  R. Hogan  Receive albedos at g-points
!   2017-04-22  R. Hogan  Store surface fluxes at all g points
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2019-01-14  R. Hogan  Save spectral flux profile if required

module radiation_homogeneous_sw

  

contains

  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Shortwave homogeneous solver, in which clouds are assumed to fill
  ! the gridbox horizontally
  

  subroutine solver_homogeneous_sw_GPU(nlev,istartcol,iendcol, &
&  config, single_level, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
&  albedo_direct, albedo_diffuse, incoming_sw, &
&  flux, lacc)
use parkind1
use yomhook
use radiation_config
use radiation_single_level
use radiation_cloud
use radiation_flux
use radiation_two_stream
use radiation_constants
use radiation_adding_ica_sw
implicit none

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(cloud_type),         intent(in) :: cloud


real(jprb), intent(in), dimension(config%n_g_sw, nlev, istartcol:iendcol) :: &
&  od, ssa, g


real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
&  od_cloud, ssa_cloud, g_cloud



real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) :: &
&  albedo_direct, albedo_diffuse, incoming_sw

type(flux_type), intent(inout):: flux





























logical, intent (in) :: lacc





end subroutine solver_homogeneous_sw_GPU

  subroutine solver_homogeneous_sw_CPU(nlev,istartcol,iendcol, &
&  config, single_level, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
&  albedo_direct, albedo_diffuse, incoming_sw, &
&  flux)
use parkind1
use yomhook
use radiation_config
use radiation_single_level
use radiation_cloud
use radiation_flux
use radiation_two_stream
use radiation_constants
use radiation_adding_ica_sw
implicit none

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(cloud_type),         intent(in) :: cloud


real(jprb), intent(in), dimension(config%n_g_sw, nlev, istartcol:iendcol) :: &
&  od, ssa, g


real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
&  od_cloud, ssa_cloud, g_cloud



real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) :: &
&  albedo_direct, albedo_diffuse, incoming_sw

type(flux_type), intent(inout):: flux


































end subroutine solver_homogeneous_sw_CPU

end module radiation_homogeneous_sw

