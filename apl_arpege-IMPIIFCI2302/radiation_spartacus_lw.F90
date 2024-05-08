! radiation_spartacus_lw.F90 - SPARTACUS longwave solver
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
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2018-09-03  R. Hogan  Security via min_cloud_effective_size
!   2018-10-08  R. Hogan  Call calc_region_properties
!   2019-01-12  R. Hogan  Use inv_inhom_effective_size if allocated
!   2022-09-01  P. Ukkonen  Optimizations for much better performance with ECCKD, including:
!                           batching section 3 computations, faster kernels, ng can be defined at compile time
module radiation_spartacus_lw

  

! Allow size of inner dimension (number of g-points) to be known at compile time if NG_LW is defined

contains

  ! Small routine for scaling cloud optical depth in the cloudy
  ! regions
#include "radiation_optical_depth_scaling.h"

  ! This module contains just one subroutine, the longwave solver
  ! using the Speedy Algorithm for Radiative Transfer through Cloud
  ! Sides (SPARTACUS), which can represent 3D effects using a matrix
  ! form of the two-stream equations.
  !
  ! Sections:
  !   1: Prepare general variables and arrays
  !   2: Prepare column-specific variables and arrays
  !   3: First loop over layers
  !     3.1: Layer-specific initialization
  !     3.2: Compute gamma variables
  !       3.2a: Clear-sky case
  !       3.2b: Cloudy case
  !     3.3: Compute reflection, transmission and emission
  !       3.3a: g-points with 3D effects
  !       3.3b: g-points without 3D effects
  !   4: Compute total sources and albedos
  !   5: Compute fluxes
  

  subroutine solver_spartacus_lw_GPU(ng_lw_in, nlev,istartcol,iendcol, &
&  config, thermodynamics, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
&  emission, albedo, &
&  flux, lacc)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_thermodynamics
use radiation_cloud
use radiation_regions
use radiation_overlap
use radiation_flux
use radiation_matrix
use radiation_two_stream
use radiation_lw_derivatives
use radiation_constants
implicit none

integer, intent(in) :: ng_lw_in           
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in)        :: config
type(thermodynamics_type),intent(in) :: thermodynamics
type(cloud_type),   intent(in)       :: cloud


real(jprb), intent(in), dimension(ng_lw_in,nlev,istartcol:iendcol) :: od


real(jprb), intent(in), target, &
&  dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: ssa, g


real(jprb), intent(in) :: od_cloud(config%n_bands_lw,nlev,istartcol:iendcol)



real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
&                            nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

real(jprb), intent(in), dimension(ng_lw_in,nlev+1,istartcol:iendcol) &
&  :: planck_hl


real(jprb), intent(in), dimension(ng_lw_in, istartcol:iendcol) &
&  :: emission, albedo

type(flux_type),          intent(inout):: flux


 

































































































































logical, intent (in) :: lacc







 

















 





end subroutine solver_spartacus_lw_GPU

  subroutine solver_spartacus_lw_CPU(ng_lw_in, nlev,istartcol,iendcol, &
&  config, thermodynamics, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
&  emission, albedo, &
&  flux)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_thermodynamics
use radiation_cloud
use radiation_regions
use radiation_overlap
use radiation_flux
use radiation_matrix
use radiation_two_stream
use radiation_lw_derivatives
use radiation_constants
implicit none

integer, intent(in) :: ng_lw_in           
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in)        :: config
type(thermodynamics_type),intent(in) :: thermodynamics
type(cloud_type),   intent(in)       :: cloud


real(jprb), intent(in), dimension(ng_lw_in,nlev,istartcol:iendcol) :: od


real(jprb), intent(in), target, &
&  dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: ssa, g


real(jprb), intent(in) :: od_cloud(config%n_bands_lw,nlev,istartcol:iendcol)



real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
&                            nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

real(jprb), intent(in), dimension(ng_lw_in,nlev+1,istartcol:iendcol) &
&  :: planck_hl


real(jprb), intent(in), dimension(ng_lw_in, istartcol:iendcol) &
&  :: emission, albedo

type(flux_type),          intent(inout):: flux


 








































































































































 

















 





end subroutine solver_spartacus_lw_CPU

end module radiation_spartacus_lw

