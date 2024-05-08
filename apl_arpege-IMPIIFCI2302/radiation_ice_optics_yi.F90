! radiation_ice_optics_yi.F90 - Yi et al. (2013) ice optical properties
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
! Authors: Mark Fielding and Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! The reference for this ice optics parameterization is Yi, B.,
! P. Yang, B.A. Baum, T. L'Ecuyer, L. Oreopoulos, E.J. Mlawer,
! A.J. Heymsfield, and K. Liou, 2013: Influence of Ice Particle
! Surface Roughening on the Global Cloud Radiative
! Effect. J. Atmos. Sci., 70, 2794-2807,
! https://doi.org/10.1175/JAS-D-13-020.1

module radiation_ice_optics_yi

  implicit none
  

  ! The number of ice coefficients depends on the parameterization
  integer, parameter :: NIceOpticsCoeffsYiSW  = 69
  integer, parameter :: NIceOpticsCoeffsYiLW  = 69

  integer, parameter :: NSingleCoeffs = 23

contains

  !---------------------------------------------------------------------
  ! Compute shortwave ice-particle scattering properties using Yi et
  ! al. (2013) parameterization
  


  !---------------------------------------------------------------------
  ! Compute longwave ice-particle scattering properties using Yi et
  ! al. (2013) parameterization
  

  subroutine calc_ice_optics_yi_sw_GPU(nb, coeff, ice_wp, &
&  re, od, scat_od, g, lacc)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp

real(jprb), intent(in) :: re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)









logical, intent (in) :: lacc







 








end subroutine calc_ice_optics_yi_sw_GPU

  subroutine calc_ice_optics_yi_lw_GPU(nb, coeff, ice_wp, &
&  re, od, scat_od, g, lacc)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp

real(jprb), intent(in) :: re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)









logical, intent (in) :: lacc







 








end subroutine calc_ice_optics_yi_lw_GPU

  subroutine calc_ice_optics_yi_sw_CPU(nb, coeff, ice_wp, &
&  re, od, scat_od, g)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp

real(jprb), intent(in) :: re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)
















 








end subroutine calc_ice_optics_yi_sw_CPU

  subroutine calc_ice_optics_yi_lw_CPU(nb, coeff, ice_wp, &
&  re, od, scat_od, g)
use parkind1


integer, intent(in)  :: nb

real(jprb), intent(in) :: coeff(:,:)

real(jprb), intent(in) :: ice_wp

real(jprb), intent(in) :: re

real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)
















 








end subroutine calc_ice_optics_yi_lw_CPU

end module radiation_ice_optics_yi

