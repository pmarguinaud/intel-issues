! radiation_ice_optics_fu.F90 - Fu's scheme for ice optical properties
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
!   2020-08-10  R. Hogan  Bounded re to be <= 100um and g to be < 1.0

module radiation_ice_optics_fu

  use parkind1

  implicit none
  

  ! The number of ice coefficients depends on the parameterization
  integer, parameter :: NIceOpticsCoeffsFuSW  = 10
  integer, parameter :: NIceOpticsCoeffsFuLW  = 11

  ! Limits based on the range of validity of the parameterizations
  real(jprb), parameter :: MaxAsymmetryFactor = 1.0_jprb - 10.0_jprb*epsilon(1.0_jprb)
  real(jprb), parameter :: MaxEffectiveRadius = 100.0e-6_jprb ! metres

contains

  !---------------------------------------------------------------------
  ! Compute shortwave ice-particle scattering properties using Fu
  ! (1996) parameterization.  The asymmetry factor in band 14 goes
  ! larger than one for re > 100.8 um, so we cap re at 100 um.
  ! Asymmetry factor is capped at just less than 1 because if it is
  ! exactly 1 then delta-Eddington scaling leads to a zero scattering
  ! optical depth and then division by zero.
  


  !---------------------------------------------------------------------
  ! Compute longwave ice-particle scattering properties using Fu et
  ! al. (1998) parameterization
  

  subroutine calc_ice_optics_fu_sw_GPU(nb, coeff, ice_wp, &
&  re, od, scat_od, g)


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp

real(jprb), intent(in) :: re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)
















end subroutine calc_ice_optics_fu_sw_GPU

  subroutine calc_ice_optics_fu_lw_GPU(nb, coeff, ice_wp, &
&  re, od, scat_od, g)


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp

real(jprb), intent(in) :: re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)
















end subroutine calc_ice_optics_fu_lw_GPU

  subroutine calc_ice_optics_fu_sw_CPU(nb, coeff, ice_wp, &
&  re, od, scat_od, g)


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp

real(jprb), intent(in) :: re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)















end subroutine calc_ice_optics_fu_sw_CPU

  subroutine calc_ice_optics_fu_lw_CPU(nb, coeff, ice_wp, &
&  re, od, scat_od, g)


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp

real(jprb), intent(in) :: re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)















end subroutine calc_ice_optics_fu_lw_CPU

end module radiation_ice_optics_fu

