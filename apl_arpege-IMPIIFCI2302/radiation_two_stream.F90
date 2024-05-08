! radiation_two_stream.F90 - Compute two-stream coefficients
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
!   2017-05-04  P Dueben/R Hogan  Use JPRD where double precision essential
!   2017-07-12  R Hogan  Optimized LW coeffs in low optical depth case
!   2017-07-26  R Hogan  Added calc_frac_scattered_diffuse_sw routine
!   2017-10-23  R Hogan  Renamed single-character variables
!   2021-02-19  R Hogan  Security for shortwave singularity
!   2022-11-22  P Ukkonen/R Hogan  Single precision uses no double precision
!   2023-09-28  R Hogan  Increased security for single-precision SW "k"

#include "ecrad_config.h"

module radiation_two_stream

  use parkind1

  implicit none
  

  ! Elsasser's factor: the effective factor by which the zenith
  ! optical depth needs to be multiplied to account for longwave
  ! transmission at all angles through the atmosphere.  Alternatively
  ! think of acos(1/lw_diffusivity) to be the effective zenith angle
  ! of longwave radiation.
  real(jprd), parameter :: LwDiffusivity   = 1.66_jprd
  real(jprb), parameter :: LwDiffusivityWP = 1.66_jprb ! Working precision version

  ! Make minimum k value depend on precision, allowing to avoid JPRD altogether
  ! real(jprb), parameter :: KMin = 1.e4_jprb * epsilon(1._jprb)
  real(jprb), parameter :: KMinSw = merge (1.e-4_jprb, 1.e-12_jprb, jprb /= jprd)
  real(jprb), parameter :: KMinLw = merge (1.e-4_jprb, 1.e-12_jprb, jprb /= jprd)
  ! The routines in this module can be called millions of times, so
  ! calling Dr Hook for each one may be a significant overhead.
  ! Uncomment the following to turn Dr Hook on.
!#define DO_DR_HOOK_TWO_STREAM

  

contains

  !---------------------------------------------------------------------
  ! Calculate the two-stream coefficients gamma1 and gamma2 for the
  ! longwave
  


  !---------------------------------------------------------------------
  ! Calculate the two-stream coefficients gamma1-gamma4 in the
  ! shortwave
  


  !---------------------------------------------------------------------
  ! Compute the longwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver formulas, as well as the
  ! upward flux at the top and the downward flux at the base of the
  ! layer due to emission from within the layer assuming a linear
  ! variation of Planck function within the layer.
  


  !---------------------------------------------------------------------
  ! Compute the longwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver formulas, as well as the
  ! upward flux at the top and the downward flux at the base of the
  ! layer due to emission from within the layer assuming a linear
  ! variation of Planck function within the layer; this version
  ! computes gamma1 and gamma2 within the same routine.
  


  !---------------------------------------------------------------------
  ! Compute the longwave transmittance to diffuse radiation in the
  ! no-scattering case, as well as the upward flux at the top and the
  ! downward flux at the base of the layer due to emission from within
  ! the layer assuming a linear variation of Planck function within
  ! the layer.
  


  !---------------------------------------------------------------------
  ! Compute the shortwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver formulas, as well as the
  ! "direct" reflection and transmission, which really means the rate
  ! of transfer of direct solar radiation (into a plane perpendicular
  ! to the direct beam) into diffuse upward and downward streams at
  ! the top and bottom of the layer, respectively.  Finally,
  ! trans_dir_dir is the transmittance of the atmosphere to direct
  ! radiation with no scattering.
  


  !---------------------------------------------------------------------
  ! Compute the shortwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver formulas, as well as the
  ! "direct" reflection and transmission, which really means the rate
  ! of transfer of direct solar radiation (into a plane perpendicular
  ! to the direct beam) into diffuse upward and downward streams at
  ! the top and bottom of the layer, respectively.  Finally,
  ! trans_dir_dir is the transmittance of the atmosphere to direct
  ! radiation with no scattering. This version incorporates the
  ! calculation of the gamma terms.
  ! Faster version with all variables in working precision (jprb), inlined
  ! computations of gamma, and exponentials kept out of loops to facilitate vectorization
  ! Adapted from RTE (Radiative Transfer for Energetics) code (Robert Pincus)
  ! Optimizations by PU and RH
  


  !---------------------------------------------------------------------
  ! Compute the fraction of shortwave transmitted diffuse radiation
  ! that is scattered during its transmission, used to compute
  ! entrapment in SPARTACUS
  

  subroutine calc_two_stream_gammas_lw_GPU(ng, ssa, g, &
&                               gamma1, gamma2)
integer, intent(in) :: ng

real(jprb), intent(in),  dimension(ng) :: ssa, g
real(jprb), intent(out), dimension(ng) :: gamma1, gamma2







end subroutine calc_two_stream_gammas_lw_GPU

  subroutine calc_two_stream_gammas_sw_GPU(ng, mu0, ssa, g, &
&                               gamma1, gamma2, gamma3)
integer, intent(in) :: ng


real(jprb), intent(in)                :: mu0
real(jprb), intent(in),  dimension(ng) :: ssa, g
real(jprb), intent(out), dimension(ng) :: gamma1, gamma2, gamma3









end subroutine calc_two_stream_gammas_sw_GPU

  subroutine calc_reflectance_transmittance_lw_GPU(ng, &
&    od, gamma1, gamma2, planck_top, planck_bot, &
&    reflectance, transmittance, source_up, source_dn)
implicit none
integer, intent(in) :: ng

real(jprb), intent(in), dimension(ng) :: od



real(jprb), intent(in), dimension(ng) :: gamma1, gamma2


real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot



real(jprb), intent(out), dimension(ng) :: reflectance, transmittance


real(jprb), intent(out), dimension(ng) :: source_up, source_dn

  
 








end subroutine calc_reflectance_transmittance_lw_GPU

  subroutine calc_ref_trans_lw_GPU(ng, &
&    od, ssa, asymmetry, planck_top, planck_bot, &
&    reflectance, transmittance, source_up, source_dn, &
&    gamma1_out, gamma2_out)
integer, intent(in) :: ng

real(jprb), intent(in), dimension(ng) :: od

real(jprb), intent(in), dimension(ng) :: ssa, asymmetry


real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot



real(jprb), intent(out), dimension(ng) :: reflectance, transmittance


real(jprb), intent(out), dimension(ng) :: source_up, source_dn
real(jprb), intent(out), optional, target, dimension(ng) :: gamma1_out, gamma2_out






 








end subroutine calc_ref_trans_lw_GPU

  subroutine calc_no_scattering_transmittance_lw_GPU(ng, &
&    od, planck_top, planck_bot, transmittance, source_up, source_dn)
integer, intent(in) :: ng

real(jprb), intent(in), dimension(ng) :: od


real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot



real(jprb), intent(out), dimension(ng) :: transmittance


real(jprb), intent(out), dimension(ng) :: source_up, source_dn
 




end subroutine calc_no_scattering_transmittance_lw_GPU

  subroutine calc_reflectance_transmittance_sw_GPU(ng, mu0, od, ssa, &
&      gamma1, gamma2, gamma3, ref_diff, trans_diff, &
&      ref_dir, trans_dir_diff, trans_dir_dir)
integer, intent(in) :: ng

real(jprb), intent(in) :: mu0

real(jprb), intent(in), dimension(ng) :: od, ssa


real(jprb), intent(in), dimension(ng) :: gamma1, gamma2, gamma3




real(jprb), intent(out), dimension(ng) :: ref_dir, trans_dir_diff



real(jprb), intent(out), dimension(ng) :: ref_diff, trans_diff

real(jprb), intent(out), dimension(ng) :: trans_dir_dir

 
  
 














end subroutine calc_reflectance_transmittance_sw_GPU

  subroutine calc_ref_trans_sw_GPU(ng, mu0, od, ssa, &
&      asymmetry, ref_diff, trans_diff, &
&      ref_dir, trans_dir_diff, trans_dir_dir, &
&      gamma1_out, gamma2_out, gamma3_out)
implicit none
integer, intent(in) :: ng

real(jprb), intent(in) :: mu0

real(jprb), intent(in), dimension(ng) :: od, ssa, asymmetry




real(jprb), intent(out), dimension(ng) :: ref_dir, trans_dir_diff



real(jprb), intent(out), dimension(ng) :: ref_diff, trans_diff

real(jprb), intent(out), dimension(ng) :: trans_dir_dir


real(jprb), intent(out), optional, target, dimension(ng) :: gamma1_out, gamma2_out, gamma3_out


 

 










end subroutine calc_ref_trans_sw_GPU

  subroutine calc_frac_scattered_diffuse_sw_GPU(ng, od, &
&      gamma1, gamma2, frac_scat_diffuse, lacc)
integer, intent(in) :: ng

real(jprb), intent(in), dimension(ng) :: od


real(jprb), intent(in), dimension(ng) :: gamma1, gamma2


real(jprb), intent(out), dimension(ng) :: frac_scat_diffuse

  
 


logical, intent (in) :: lacc



end subroutine calc_frac_scattered_diffuse_sw_GPU

  subroutine calc_two_stream_gammas_lw_CPU(ng, ssa, g, &
&                               gamma1, gamma2)
integer, intent(in) :: ng

real(jprb), intent(in),  dimension(ng) :: ssa, g
real(jprb), intent(out), dimension(ng) :: gamma1, gamma2





end subroutine calc_two_stream_gammas_lw_CPU

  subroutine calc_two_stream_gammas_sw_CPU(ng, mu0, ssa, g, &
&                               gamma1, gamma2, gamma3)
integer, intent(in) :: ng


real(jprb), intent(in)                :: mu0
real(jprb), intent(in),  dimension(ng) :: ssa, g
real(jprb), intent(out), dimension(ng) :: gamma1, gamma2, gamma3







end subroutine calc_two_stream_gammas_sw_CPU

  subroutine calc_reflectance_transmittance_lw_CPU(ng, &
&    od, gamma1, gamma2, planck_top, planck_bot, &
&    reflectance, transmittance, source_up, source_dn)
implicit none
integer, intent(in) :: ng

real(jprb), intent(in), dimension(ng) :: od



real(jprb), intent(in), dimension(ng) :: gamma1, gamma2


real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot



real(jprb), intent(out), dimension(ng) :: reflectance, transmittance


real(jprb), intent(out), dimension(ng) :: source_up, source_dn

  
 





end subroutine calc_reflectance_transmittance_lw_CPU

  subroutine calc_ref_trans_lw_CPU(ng, &
&    od, ssa, asymmetry, planck_top, planck_bot, &
&    reflectance, transmittance, source_up, source_dn, &
&    gamma1_out, gamma2_out)
integer, intent(in) :: ng

real(jprb), intent(in), dimension(ng) :: od

real(jprb), intent(in), dimension(ng) :: ssa, asymmetry


real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot



real(jprb), intent(out), dimension(ng) :: reflectance, transmittance


real(jprb), intent(out), dimension(ng) :: source_up, source_dn
real(jprb), intent(out), optional, target, dimension(ng) :: gamma1_out, gamma2_out






 



end subroutine calc_ref_trans_lw_CPU

  subroutine calc_no_scattering_transmittance_lw_CPU(ng, &
&    od, planck_top, planck_bot, transmittance, source_up, source_dn)
integer, intent(in) :: ng

real(jprb), intent(in), dimension(ng) :: od


real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot



real(jprb), intent(out), dimension(ng) :: transmittance


real(jprb), intent(out), dimension(ng) :: source_up, source_dn
 



end subroutine calc_no_scattering_transmittance_lw_CPU

  subroutine calc_reflectance_transmittance_sw_CPU(ng, mu0, od, ssa, &
&      gamma1, gamma2, gamma3, ref_diff, trans_diff, &
&      ref_dir, trans_dir_diff, trans_dir_dir)
integer, intent(in) :: ng

real(jprb), intent(in) :: mu0

real(jprb), intent(in), dimension(ng) :: od, ssa


real(jprb), intent(in), dimension(ng) :: gamma1, gamma2, gamma3




real(jprb), intent(out), dimension(ng) :: ref_dir, trans_dir_diff



real(jprb), intent(out), dimension(ng) :: ref_diff, trans_diff

real(jprb), intent(out), dimension(ng) :: trans_dir_dir

 
  
 










end subroutine calc_reflectance_transmittance_sw_CPU

  subroutine calc_ref_trans_sw_CPU(ng, mu0, od, ssa, &
&      asymmetry, ref_diff, trans_diff, &
&      ref_dir, trans_dir_diff, trans_dir_dir, &
&      gamma1_out, gamma2_out, gamma3_out)
implicit none
integer, intent(in) :: ng

real(jprb), intent(in) :: mu0

real(jprb), intent(in), dimension(ng) :: od, ssa, asymmetry




real(jprb), intent(out), dimension(ng) :: ref_dir, trans_dir_diff



real(jprb), intent(out), dimension(ng) :: ref_diff, trans_diff

real(jprb), intent(out), dimension(ng) :: trans_dir_dir


real(jprb), intent(out), optional, target, dimension(ng) :: gamma1_out, gamma2_out, gamma3_out



 
  

 










end subroutine calc_ref_trans_sw_CPU

  subroutine calc_frac_scattered_diffuse_sw_CPU(ng, od, &
&      gamma1, gamma2, frac_scat_diffuse)
integer, intent(in) :: ng

real(jprb), intent(in), dimension(ng) :: od


real(jprb), intent(in), dimension(ng) :: gamma1, gamma2


real(jprb), intent(out), dimension(ng) :: frac_scat_diffuse

  
 





end subroutine calc_frac_scattered_diffuse_sw_CPU

end module radiation_two_stream

