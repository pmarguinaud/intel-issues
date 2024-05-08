! radiation_tripleclouds_lw.F90 - Longwave "Tripleclouds" solver
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
!   2017-04-28  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-10-08  R. Hogan  Call calc_region_properties
!   2020-09-18  R. Hogan  Replaced some array expressions with loops
!   2020-09-19  R. Hogan  Implement the cloud-only-scattering optimization
!   2022-09-01  P. Ukkonen  Optimizations for much better performance with ECCKD, including:
!                           batching section 3 computations, faster kernels, ng can be defined at compile time

module radiation_tripleclouds_lw

  

contains
  ! Small routine for scaling cloud optical depth in the cloudy
  ! regions
#include "radiation_optical_depth_scaling.h"

  !---------------------------------------------------------------------
  ! This module contains just one subroutine, the longwave
  ! "Tripleclouds" solver in which cloud inhomogeneity is treated by
  ! dividing each model level into three regions, one clear and two
  ! cloudy (with differing optical depth). This approach was described
  ! by Shonk and Hogan (2008).
  

  subroutine solver_tripleclouds_lw_GPU(ng_lw_in, nlev,istartcol,iendcol, &
&  config, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
&  emission, albedo, &
&  flux, lacc)
use parkind1
use yomhook

use radiation_config
use radiation_cloud
use radiation_regions
use radiation_overlap
use radiation_flux
use radiation_matrix
use radiation_two_stream
use radiation_adding_ica_lw
use radiation_lw_derivatives
implicit none


integer, intent(in) :: ng_lw_in           
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(cloud_type),         intent(in) :: cloud


real(jprb), intent(in), dimension(ng_lw_in,nlev,istartcol:iendcol) :: od


real(jprb), intent(in), &
&  dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: ssa, g


real(jprb), intent(in) :: od_cloud(config%n_bands_lw,nlev,istartcol:iendcol)



real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
&                            nlev,istartcol:iendcol) :: ssa_cloud, g_cloud


real(jprb), intent(in), dimension(ng_lw_in,nlev+1,istartcol:iendcol) :: planck_hl


real(jprb), intent(in), dimension(ng_lw_in, istartcol:iendcol) :: emission, albedo










type(flux_type), intent(inout):: flux





















































logical, intent (in) :: lacc










 

end subroutine solver_tripleclouds_lw_GPU

  subroutine solver_tripleclouds_lw_CPU(ng_lw_in, nlev,istartcol,iendcol, &
&  config, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
&  emission, albedo, &
&  flux)
use parkind1
use yomhook

use radiation_config
use radiation_cloud
use radiation_regions
use radiation_overlap
use radiation_flux
use radiation_matrix
use radiation_two_stream
use radiation_adding_ica_lw
use radiation_lw_derivatives
implicit none


integer, intent(in) :: ng_lw_in           
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(cloud_type),         intent(in) :: cloud


real(jprb), intent(in), dimension(ng_lw_in,nlev,istartcol:iendcol) :: od


real(jprb), intent(in), &
&  dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: ssa, g


real(jprb), intent(in) :: od_cloud(config%n_bands_lw,nlev,istartcol:iendcol)



real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
&                            nlev,istartcol:iendcol) :: ssa_cloud, g_cloud


real(jprb), intent(in), dimension(ng_lw_in,nlev+1,istartcol:iendcol) :: planck_hl


real(jprb), intent(in), dimension(ng_lw_in, istartcol:iendcol) :: emission, albedo










type(flux_type), intent(inout):: flux































































 

end subroutine solver_tripleclouds_lw_CPU

end module radiation_tripleclouds_lw

