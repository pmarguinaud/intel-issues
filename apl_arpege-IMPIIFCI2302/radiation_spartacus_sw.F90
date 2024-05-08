! radiation_spartacus_sw.F90 - SPARTACUS shortwave solver
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
!   2017-04-11  R. Hogan  Receive albedos at g-points
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-06-30  R. Hogan  Reformulate to use total_albedo_direct not total_source
!   2017-07-03  R. Hogan  Explicit calculation of encroachment
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-02-20  R. Hogan  Corrected "computed" encroachment
!   2018-03-09  R. Hogan  Security on computed encroachment, transmittance and reflectance
!   2018-08-29  R. Hogan  Reformulate horizontal migration distances in step_migrations
!   2018-09-02  R. Hogan  Bug fix in x_direct in step_migrations
!   2018-09-03  R. Hogan  Security via min_cloud_effective_size
!   2018-09-04  R. Hogan  Use encroachment_scaling and read encroachment edge_length from upper layer
!   2018-09-13  R. Hogan  Added "Fractal" encroachment option
!   2018-09-14  R. Hogan  Added "Zero" encroachment option
!   2018-10-08  R. Hogan  Call calc_region_properties
!   2018-10-15  R. Hogan  Added call to fast_expm_exchange instead of expm
!   2019-01-12  R. Hogan  Use inv_inhom_effective_size if allocated
!   2019-02-10  R. Hogan  Renamed "encroachment" to "entrapment"
!   2022-09-01  P. Ukkonen  Optimizations for much better performance with ECCKD, including:
!                           batching section 3 computations, faster kernels, ng can be defined at compile time


module radiation_spartacus_sw

  

! Allow size of inner dimension (number of g-points) to be known at compile time if NG_SW is defined

contains

  ! Small routine for scaling cloud optical depth in the cloudy
  ! regions
#include "radiation_optical_depth_scaling.h"

  ! This module contains just one exported subroutine, the shortwave
  ! solver using the Speedy Algorithm for Radiative Transfer through
  ! Cloud Sides (SPARTACUS), which can represent 3D effects using a
  ! matrix form of the two-stream equations.
  !
  ! Sections:
  !   1: Prepare general variables and arrays
  !   2: Prepare column-specific variables and arrays
  !   3: First loop over layers
  !     3.1: Layer-specific initialization
  !     3.2: Compute gamma variables
  !       3.2a: Clear-sky case
  !       3.2b: Cloudy case
  !     3.3: Compute reflection, transmission and sources
  !       3.3a: g-points with 3D effects
  !       3.3b: g-points without 3D effects
  !   4: Compute total albedos
  !     4.1: Adding method
  !     4.2: Overlap and entrapment
  !   5: Compute fluxes
  


  ! Step the horizontal migration distances from the base of a layer
  ! to the top, accounting for the extra distance travelled within the
  ! layer
  

  

  subroutine solver_spartacus_sw_GPU(ng_sw_in, nlev,istartcol,iendcol, &
&  config, single_level, thermodynamics, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
&  albedo_direct, albedo_diffuse, incoming_sw, &
&  flux, lacc)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_single_level
use radiation_thermodynamics
use radiation_cloud
use radiation_regions
use radiation_overlap
use radiation_flux
use radiation_matrix
use radiation_two_stream
use radiation_constants
implicit none

integer, intent(in) :: ng_sw_in           
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(cloud_type),         intent(in) :: cloud


real(jprb), intent(in), dimension(ng_sw_in,nlev,istartcol:iendcol) :: &
&  od, ssa, g


real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
&  od_cloud, ssa_cloud, g_cloud



real(jprb), intent(in), dimension(ng_sw_in,istartcol:iendcol) :: &
&  albedo_direct, albedo_diffuse, incoming_sw

type(flux_type), intent(inout):: flux



 

 































































































































logical, intent (in) :: lacc














  










 





end subroutine solver_spartacus_sw_GPU

  pure subroutine step_migrations_GPU(ng_sw_in, nreg, cloud_frac, &
&  layer_depth, tan_diffuse_angle_3d, tan_sza, &
&  reflectance, transmittance, ref_dir, trans_dir_dir, &
&  trans_dir_diff, total_albedo_diff, total_albedo_dir, &
&  x_diffuse, x_direct, lacc)
use parkind1
implicit none


integer, intent(in) :: ng_sw_in, nreg

real(jprb), intent(in) :: cloud_frac


real(jprb), intent(in) :: layer_depth, tan_diffuse_angle_3d, tan_sza

real(jprb), intent(in), dimension(ng_sw_in, nreg, nreg) :: reflectance, transmittance

real(jprb), intent(in), dimension(ng_sw_in, nreg, nreg) :: ref_dir, trans_dir_dir


real(jprb), intent(in), dimension(ng_sw_in, nreg, nreg) :: trans_dir_diff


real(jprb), intent(in), dimension(ng_sw_in, nreg, nreg) &
&  :: total_albedo_diff, total_albedo_dir


real(jprb), intent(inout), dimension(ng_sw_in, nreg) :: x_diffuse, x_direct










logical, intent (in) :: lacc













end subroutine step_migrations_GPU

  pure subroutine write_gamma_diag_GPU(n, nreg, jreg, od_region, &
&   gamma1, gamma2, gamma3, ssa, one_over_mu0, Gamma_z1, lacc)
use parkind1


integer, intent(in) :: n

integer, intent(in) :: nreg, jreg
real(jprb), intent(in), dimension(n) :: od_region, gamma1, gamma2, gamma3, ssa
real(jprb), intent(in) :: one_over_mu0
real(jprb), intent(inout), dimension(n, 3*nreg, 3*nreg) :: Gamma_z1

logical, intent (in) :: lacc

end subroutine write_gamma_diag_GPU

  subroutine solver_spartacus_sw_CPU(ng_sw_in, nlev,istartcol,iendcol, &
&  config, single_level, thermodynamics, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
&  albedo_direct, albedo_diffuse, incoming_sw, &
&  flux)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_single_level
use radiation_thermodynamics
use radiation_cloud
use radiation_regions
use radiation_overlap
use radiation_flux
use radiation_matrix
use radiation_two_stream
use radiation_constants
implicit none

integer, intent(in) :: ng_sw_in           
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(cloud_type),         intent(in) :: cloud


real(jprb), intent(in), dimension(ng_sw_in,nlev,istartcol:iendcol) :: &
&  od, ssa, g


real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
&  od_cloud, ssa_cloud, g_cloud



real(jprb), intent(in), dimension(ng_sw_in,istartcol:iendcol) :: &
&  albedo_direct, albedo_diffuse, incoming_sw

type(flux_type), intent(inout):: flux



 

 













































































































































  










 





end subroutine solver_spartacus_sw_CPU

  pure subroutine step_migrations_CPU(ng_sw_in, nreg, cloud_frac, &
&  layer_depth, tan_diffuse_angle_3d, tan_sza, &
&  reflectance, transmittance, ref_dir, trans_dir_dir, &
&  trans_dir_diff, total_albedo_diff, total_albedo_dir, &
&  x_diffuse, x_direct)
use parkind1
implicit none


integer, intent(in) :: ng_sw_in, nreg

real(jprb), intent(in) :: cloud_frac


real(jprb), intent(in) :: layer_depth, tan_diffuse_angle_3d, tan_sza

real(jprb), intent(in), dimension(ng_sw_in, nreg, nreg) :: reflectance, transmittance

real(jprb), intent(in), dimension(ng_sw_in, nreg, nreg) :: ref_dir, trans_dir_dir


real(jprb), intent(in), dimension(ng_sw_in, nreg, nreg) :: trans_dir_diff


real(jprb), intent(in), dimension(ng_sw_in, nreg, nreg) &
&  :: total_albedo_diff, total_albedo_dir


real(jprb), intent(inout), dimension(ng_sw_in, nreg) :: x_diffuse, x_direct























end subroutine step_migrations_CPU

  pure subroutine write_gamma_diag_CPU(n, nreg, jreg, od_region, &
&   gamma1, gamma2, gamma3, ssa, one_over_mu0, Gamma_z1)
use parkind1


integer, intent(in) :: n

integer, intent(in) :: nreg, jreg
real(jprb), intent(in), dimension(n) :: od_region, gamma1, gamma2, gamma3, ssa
real(jprb), intent(in) :: one_over_mu0
real(jprb), intent(inout), dimension(n, 3*nreg, 3*nreg) :: Gamma_z1


end subroutine write_gamma_diag_CPU

end module radiation_spartacus_sw

