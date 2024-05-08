! radiation_liquid_optics_socrates.F90 - SOCRATES method for parameterizing liquid droplet optics
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
!   2020-08-10  R. Hogan  Bounded re to be >=1.2um and <=50um

module radiation_liquid_optics_socrates

  use parkind1

  implicit none
  

  ! SOCRATES (Edwards-Slingo) parameterizes info on the dependence of
  ! the scattering properties in each band on effective radius in
  ! terms of 16 coefficients
  integer, parameter :: NLiqOpticsCoeffsSOCRATES = 16

  ! Range of valid input effective radius, in microns
  real(jprb), parameter :: MinEffectiveRadius = 1.2e-6
  real(jprb), parameter :: MaxEffectiveRadius = 50.0e-6

contains

  !---------------------------------------------------------------------
  ! Compute liquid-droplet scattering properties using a
  ! parameterization consisting of Pade approximants from the
  ! SOCRATES (Edwards-Slingo) code
  

  subroutine calc_liq_optics_socrates_GPU(nb, coeff, lwp, re_in, od, scat_od, g)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: lwp, re_in

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)












end subroutine calc_liq_optics_socrates_GPU

  subroutine calc_liq_optics_socrates_CPU(nb, coeff, lwp, re_in, od, scat_od, g)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: lwp, re_in

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)











end subroutine calc_liq_optics_socrates_CPU

end module radiation_liquid_optics_socrates

