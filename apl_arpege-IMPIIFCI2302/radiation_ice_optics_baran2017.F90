! radiation_ice_optics_baran2017.F90 - 2017 parameterization of Baran's ice optical properties
!
! (C) Copyright 2017- ECMWF.
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

module radiation_ice_optics_baran2017

  implicit none
  

  ! The number of ice coefficients depends on the parameterization
  integer, parameter :: NIceOpticsCoeffsBaran2017 = 9
  integer, parameter :: NIceOpticsGeneralCoeffsBaran2017 = 5

contains

  
  !---------------------------------------------------------------------
  ! Compute ice-particle scattering properties using a
  ! parameterization as a function of ice water mixing ratio and
  ! temperature
  

  subroutine calc_ice_optics_baran2017_GPU(nb, coeff_gen, coeff, ice_wp, &
&  qi, temperature, od, scat_od, g, lacc)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff_gen(:)

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp, qi

real(jprb), intent(in) :: temperature

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)


logical, intent (in) :: lacc










end subroutine calc_ice_optics_baran2017_GPU

  subroutine calc_ice_optics_baran2017_CPU(nb, coeff_gen, coeff, ice_wp, &
&  qi, temperature, od, scat_od, g)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff_gen(:)

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp, qi

real(jprb), intent(in) :: temperature

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)












end subroutine calc_ice_optics_baran2017_CPU

end module radiation_ice_optics_baran2017

