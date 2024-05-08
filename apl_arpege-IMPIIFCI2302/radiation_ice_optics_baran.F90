! radiation_ice_optics_fu.F90 - Scheme for ice optical properties adapted from Baran's data
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

module radiation_ice_optics_baran

  implicit none
  

  ! The number of ice coefficients depends on the parameterization
  integer, parameter :: NIceOpticsCoeffsBaran = 9
  integer, parameter :: NIceOpticsCoeffsBaran2016 = 5

contains

  
  !---------------------------------------------------------------------
  ! Compute ice-particle scattering properties using a
  ! parameterization as a function of ice water mixing ratio only
  


  !---------------------------------------------------------------------
  ! Compute ice-particle scattering properties using a
  ! parameterization as a function of ice water mixing ratio and
  ! temperature
  

  subroutine calc_ice_optics_baran_GPU(nb, coeff, ice_wp, &
&  qi, od, scat_od, g, lacc)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp, qi

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)
logical, intent (in) :: lacc










end subroutine calc_ice_optics_baran_GPU

  subroutine calc_ice_optics_baran2016_GPU(nb, coeff, ice_wp, &
&  qi, temperature, od, scat_od, g, lacc)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp, qi

real(jprb), intent(in) :: temperature

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)


logical, intent (in) :: lacc








end subroutine calc_ice_optics_baran2016_GPU

  subroutine calc_ice_optics_baran_CPU(nb, coeff, ice_wp, &
&  qi, od, scat_od, g)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp, qi

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)










end subroutine calc_ice_optics_baran_CPU

  subroutine calc_ice_optics_baran2016_CPU(nb, coeff, ice_wp, &
&  qi, temperature, od, scat_od, g)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp, qi

real(jprb), intent(in) :: temperature

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)










end subroutine calc_ice_optics_baran2016_CPU

end module radiation_ice_optics_baran

