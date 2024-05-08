! radiation_liquid_optics_slingo.F90 - Slingo SW & Lindner-Li LW parameterization of liquid droplet optics
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

module radiation_liquid_optics_slingo

  implicit none
  

  integer, parameter :: NLiqOpticsCoeffsSlingoSW = 6
  integer, parameter :: NLiqOpticsCoeffsLindnerLiLW = 13

contains

  !---------------------------------------------------------------------
  ! Compute liquid-droplet scattering properties in the shortwave from
  ! Slingo (1989). WARNING: this parameterization is known not to be
  ! very accurate: see Nielsen et al. (GMD 2014).
  


  !---------------------------------------------------------------------
  ! Compute liquid-droplet scattering properties in the longwave from
  ! Lindner & Li (2000)
  

  subroutine calc_liq_optics_slingo_GPU(nb, coeff, lwp, re, od, scat_od, g)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: lwp, re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)














end subroutine calc_liq_optics_slingo_GPU

  subroutine calc_liq_optics_lindner_li_GPU(nb, coeff, lwp, re, od, scat_od, g)
use parkind1
use yomhook

integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: lwp, re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)











end subroutine calc_liq_optics_lindner_li_GPU

  subroutine calc_liq_optics_slingo_CPU(nb, coeff, lwp, re, od, scat_od, g)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: lwp, re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)













end subroutine calc_liq_optics_slingo_CPU

  subroutine calc_liq_optics_lindner_li_CPU(nb, coeff, lwp, re, od, scat_od, g)
use parkind1
use yomhook

integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: lwp, re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)













end subroutine calc_liq_optics_lindner_li_CPU

end module radiation_liquid_optics_slingo

