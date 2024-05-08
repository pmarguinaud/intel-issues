! radiation_cloudless_lw.F90 - Longwave homogeneous cloudless solver
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

module radiation_cloudless_lw



contains

  !---------------------------------------------------------------------
  ! Longwave homogeneous solver containing no clouds
  

  subroutine solver_cloudless_lw_GPU(nlev,istartcol,iendcol, &
&  config, od, ssa, g, planck_hl, emission, albedo, flux, lacc)
use parkind1
use yomhook
use radiation_config
use radiation_flux
use radiation_two_stream
use radiation_adding_ica_lw
use radiation_lw_derivatives
implicit none

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config


real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
&  od
real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
&  ssa, g

real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
&  planck_hl


real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) &
&  :: emission, albedo

type(flux_type), intent(inout):: flux



















logical, intent (in) :: lacc





end subroutine solver_cloudless_lw_GPU

  subroutine solver_cloudless_lw_CPU(nlev,istartcol,iendcol, &
&  config, od, ssa, g, planck_hl, emission, albedo, flux)
use parkind1
use yomhook
use radiation_config
use radiation_flux
use radiation_two_stream
use radiation_adding_ica_lw
use radiation_lw_derivatives
implicit none

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in) :: config


real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
&  od
real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
&  ssa, g

real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
&  planck_hl


real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) &
&  :: emission, albedo

type(flux_type), intent(inout):: flux
























end subroutine solver_cloudless_lw_CPU

end module radiation_cloudless_lw

