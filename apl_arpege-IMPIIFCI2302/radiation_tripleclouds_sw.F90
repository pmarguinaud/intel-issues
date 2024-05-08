! radiation_tripleclouds_sw.F90 - Shortwave "Tripleclouds" solver
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
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-10-08  R. Hogan  Call calc_region_properties
!   2019-01-02  R. Hogan  Fixed problem of do_save_spectral_flux .and. .not. do_sw_direct
!   2020-09-18  R. Hogan  Replaced some array expressions with loops for speed
!   2021-10-01  P. Ukkonen Performance optimizations: batched computations

module radiation_tripleclouds_sw

  

contains
  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  ! Small routine for scaling cloud optical depth in the cloudy
  ! regions
#include "radiation_optical_depth_scaling.h"

  !---------------------------------------------------------------------
  ! This module contains just one subroutine, the shortwave
  ! "Tripleclouds" solver in which cloud inhomogeneity is treated by
  ! dividing each model level into three regions, one clear and two
  ! cloudy (with differing optical depth). This approach was described
  ! by Shonk and Hogan (2008).
  

  subroutine solver_tripleclouds_sw_GPU(ng_sw_in, nlev,istartcol,iendcol, &
&  config, single_level, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
&  albedo_direct, albedo_diffuse, incoming_sw, &
&  flux, lacc)
use parkind1
use yomhook

use radiation_config
use radiation_single_level
use radiation_cloud
use radiation_regions
use radiation_overlap
use radiation_flux
use radiation_matrix
use radiation_two_stream
implicit none




integer, intent(in) :: ng_sw_in           
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(cloud_type),         intent(in) :: cloud


real(jprb), intent(in), dimension(ng_sw_in,nlev,istartcol:iendcol) :: &
&  od, ssa, g


real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol) &
&  :: od_cloud, ssa_cloud, g_cloud










real(jprb), intent(in), dimension(ng_sw_in,istartcol:iendcol) :: &
&  albedo_direct, albedo_diffuse, incoming_sw

type(flux_type), intent(inout):: flux


















































logical, intent (in) :: lacc












 

end subroutine solver_tripleclouds_sw_GPU

  subroutine solver_tripleclouds_sw_CPU(ng_sw_in, nlev,istartcol,iendcol, &
&  config, single_level, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
&  albedo_direct, albedo_diffuse, incoming_sw, &
&  flux)
use parkind1
use yomhook

use radiation_config
use radiation_single_level
use radiation_cloud
use radiation_regions
use radiation_overlap
use radiation_flux
use radiation_matrix
use radiation_two_stream
implicit none




integer, intent(in) :: ng_sw_in           
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(cloud_type),         intent(in) :: cloud


real(jprb), intent(in), dimension(ng_sw_in,nlev,istartcol:iendcol) :: &
&  od, ssa, g


real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol) &
&  :: od_cloud, ssa_cloud, g_cloud










real(jprb), intent(in), dimension(ng_sw_in,istartcol:iendcol) :: &
&  albedo_direct, albedo_diffuse, incoming_sw

type(flux_type), intent(inout):: flux






























































 

end subroutine solver_tripleclouds_sw_CPU

end module radiation_tripleclouds_sw

