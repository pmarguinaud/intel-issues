! radiation_mcica_acc_lw.F90 - Monte-Carlo Independent Column Approximation longtwave solver
!
! (C) Copyright 2015- ECMWF.
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
!   2017-07-12  R. Hogan  Call fast adding method if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_mcica_acc_lw

  

contains

  !---------------------------------------------------------------------
  ! Longwave Monte Carlo Independent Column Approximation
  ! (McICA). This implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. This method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. The cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
  

  subroutine solver_mcica_acc_lw_GPU(nlev,istartcol,iendcol, &
&  config, single_level, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
&  emission, albedo, &
&  flux, lacc)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_single_level
use radiation_cloud
use radiation_flux
use radiation_two_stream
use radiation_adding_ica_lw
use radiation_lw_derivatives
use radiation_cloud_generator_acc
use radiation_cloud_cover
implicit none

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(cloud_type),         intent(in) :: cloud


real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
&  od
real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
&  ssa, g


real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
&  od_cloud
real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
&  nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
&  planck_hl


real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

type(flux_type), intent(inout):: flux

















































logical, intent (in) :: lacc


















































end subroutine solver_mcica_acc_lw_GPU

  subroutine solver_mcica_acc_lw_CPU(nlev,istartcol,iendcol, &
&  config, single_level, cloud, &
&  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
&  emission, albedo, &
&  flux)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_single_level
use radiation_cloud
use radiation_flux
use radiation_two_stream
use radiation_adding_ica_lw
use radiation_lw_derivatives
use radiation_cloud_generator_acc
use radiation_cloud_cover
implicit none

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(cloud_type),         intent(in) :: cloud


real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
&  od
real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
&  ssa, g


real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
&  od_cloud
real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
&  nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
&  planck_hl


real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

type(flux_type), intent(inout):: flux































































end subroutine solver_mcica_acc_lw_CPU

end module radiation_mcica_acc_lw

