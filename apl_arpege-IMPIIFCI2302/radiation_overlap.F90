! radiation_overlap.F90 - Module to compute cloud overlap quantities
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
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-10-05  R. Hogan  Generalized alpha overlap for non-equal regions
!   2018-10-08  R. Hogan  Removed calc_region_fractions

module radiation_overlap

  implicit none

  
  

  integer, parameter, private :: nreg = 3

contains


  ! This function now superceded by calc_region_properties in module
  ! radiation_regions
  ! !---------------------------------------------------------------------
  ! ! Return an array of length nreg containing the fraction of the
  ! ! gridbox occupied by each region for the specified cloud fraction.
  ! pure function calc_region_fractions(nreg, cloud_fraction)

  !   use parkind1, only : jprb

  !   integer,    intent(in)      :: nreg
  !   real(jprb), intent(in)      :: cloud_fraction

  !   real(jprb), dimension(nreg) :: calc_region_fractions
  !   integer :: jreg

  !   if (nreg == 1) then
  !     ! Only one region: must occupy all of gridbox
  !     calc_region_fractions(1) = 1.0_jprb
  !   else
  !     ! Two or more regions: the first is the cloud-free region
  !     calc_region_fractions(1) = 1.0_jprb - cloud_fraction

  !     do jreg = 2,nreg
  !       ! The cloudy regions are assumed to each have the same
  !       ! fraction - see Shonk and Hogan (2008) for justification
  !       calc_region_fractions(jreg) = cloud_fraction / (nreg - 1.0_jprb)
  !     end do
  !   end if

  ! end function calc_region_fractions

  !---------------------------------------------------------------------
  ! Calculate a matrix expressing the overlap of regions in adjacent
  ! layers, using the method of Shonk et al. (2010) in terms of their
  ! "beta" overlap parameter
  

  ! Double-precision version for SPARTACUS_SW
  
  !---------------------------------------------------------------------
  ! Calculate a matrix expressing the overlap of regions in adjacent
  ! layers, using the Hogan and Illingworth (2000) "alpha" overlap
  ! parameter, but allowing for the two cloudy regions in the
  ! Tripleclouds assumption to have different areas
  

  ! Double-precision version for SPARTACUS_SW
  

  !---------------------------------------------------------------------
  ! Calculate a matrix expressing the overlap of regions in adjacent
  ! layers, using the Hogan and Illingworth (2000) "alpha" overlap
  ! parameter, and assuming the two cloudy regions in the Tripleclouds
  ! assumption have the same area
  

  !---------------------------------------------------------------------
  ! Compute the upward and downward overlap matrices u_matrix and
  ! v_matrix, respectively, where u_matrix is defined such that
  ! y=u_matrix*x, where x is a vector of upwelling fluxes in each
  ! region just below an interface, and y is a vector of upwelling
  ! fluxes in each region just above that interface. For nlev model
  ! levels there are nlev+1 interfaces including the ground and
  ! top-of-atmosphere, and so that is one of the dimensions of
  ! u_matrix and v_matrix.
  

  ! Double-precision version for SPARTACUS_SW
  

  pure function calc_beta_overlap_matrix_GPU(nreg_in, op, frac_upper, frac_lower, &
&  frac_threshold, lacc) result(overlap_matrix)
use parkind1
integer, intent(in) :: nreg_in 


real(jprb), intent(in), dimension(nreg) :: op, frac_upper, frac_lower


real(jprb), intent(in) :: frac_threshold

real(jprb) :: overlap_matrix(nreg,nreg)







logical, intent (in) :: lacc









end function calc_beta_overlap_matrix_GPU

  pure function calc_beta_overlap_matrix_dp_GPU(nreg_in, op, frac_upper, frac_lower, &
&  frac_threshold, lacc) result(overlap_matrix)
use parkind1
integer, intent(in) :: nreg_in 


real(jprd), intent(in), dimension(nreg) :: op, frac_upper, frac_lower


real(jprd), intent(in) :: frac_threshold

real(jprd) :: overlap_matrix(nreg,nreg)







logical, intent (in) :: lacc









end function calc_beta_overlap_matrix_dp_GPU

  pure function calc_alpha_overlap_matrix_GPU(nreg_in, op, op_inhom, &
&  frac_upper, frac_lower, lacc) result(overlap_matrix)
use parkind1
integer, intent(in) :: nreg_in 


real(jprb), intent(in) :: op, op_inhom


real(jprb), intent(in), dimension(nreg) :: frac_upper, frac_lower

real(jprb) :: overlap_matrix(nreg,nreg)








logical, intent (in) :: lacc





end function calc_alpha_overlap_matrix_GPU

  pure function calc_alpha_overlap_matrix_dp_GPU(nreg_in, op, op_inhom, &
&  frac_upper, frac_lower, lacc) result(overlap_matrix)
use parkind1
integer, intent(in) :: nreg_in 


real(jprd), intent(in) :: op, op_inhom


real(jprd), intent(in), dimension(nreg) :: frac_upper, frac_lower

real(jprd) :: overlap_matrix(nreg,nreg)








logical, intent (in) :: lacc






end function calc_alpha_overlap_matrix_dp_GPU

  pure function calc_alpha_overlap_matrix_simple_GPU(nreg_in, op, op_inhom, &
&  cf_upper, cf_lower, lacc) result(overlap_matrix)
use parkind1
integer, intent(in) :: nreg_in 


real(jprb), intent(in) :: op, op_inhom

real(jprb), intent(in) :: cf_upper, cf_lower

real(jprb) :: overlap_matrix(nreg,nreg)



logical, intent (in) :: lacc




end function calc_alpha_overlap_matrix_simple_GPU

  subroutine calc_overlap_matrices_GPU(nlev, nreg_in, is_clear_sky_layer, &
&     region_fracs, overlap_param, v_matrix, u_matrix, decorrelation_scaling, &
&     cloud_fraction_threshold, cloud_cov, use_beta_overlap, lacc)
use parkind1
use yomhook

integer,  intent(in) :: nlev, nreg_in

logical, intent(in) :: is_clear_sky_layer(0:nlev+1)



real(jprb), intent(in), dimension(nreg,nlev)  :: region_fracs


real(jprb), intent(in), dimension(nlev-1)  :: overlap_param

real(jprb), intent(out), dimension(nreg,nreg,nlev+1) :: v_matrix

real(jprb), intent(out), optional, dimension(nreg,nreg,nlev+1) &
&  :: u_matrix





real(jprb), intent(in), optional :: decorrelation_scaling

real(jprb), intent(in), optional :: cloud_fraction_threshold

real(jprb), intent(out), optional :: cloud_cov

logical, intent(in), optional :: use_beta_overlap


















logical, intent (in) :: lacc














 



end subroutine calc_overlap_matrices_GPU

  subroutine calc_overlap_matrices_dp_GPU(nlev, nreg_in, is_clear_sky_layer, &
&     region_fracs, overlap_param, v_matrix, u_matrix, decorrelation_scaling, &
&     cloud_fraction_threshold, cloud_cov, use_beta_overlap, lacc)
use parkind1
use yomhook

integer,  intent(in) :: nlev, nreg_in

logical, intent(in) :: is_clear_sky_layer(0:nlev+1)



real(jprd), intent(in), dimension(nreg,nlev)  :: region_fracs


real(jprd), intent(in), dimension(nlev-1)  :: overlap_param

real(jprd), intent(out), dimension(nreg,nreg,nlev+1) :: v_matrix

real(jprd), intent(out), optional, dimension(nreg,nreg,nlev+1) &
&  :: u_matrix





real(jprd), intent(in), optional :: decorrelation_scaling

real(jprd), intent(in), optional :: cloud_fraction_threshold

real(jprd), intent(out), optional :: cloud_cov

logical, intent(in), optional :: use_beta_overlap


















logical, intent (in) :: lacc














 



end subroutine calc_overlap_matrices_dp_GPU

  pure function calc_beta_overlap_matrix_CPU(nreg_in, op, frac_upper, frac_lower, &
&  frac_threshold) result(overlap_matrix)
use parkind1
integer, intent(in) :: nreg_in 


real(jprb), intent(in), dimension(nreg) :: op, frac_upper, frac_lower


real(jprb), intent(in) :: frac_threshold

real(jprb) :: overlap_matrix(nreg,nreg)
















end function calc_beta_overlap_matrix_CPU

  pure function calc_beta_overlap_matrix_dp_CPU(nreg_in, op, frac_upper, frac_lower, &
&  frac_threshold) result(overlap_matrix)
use parkind1
integer, intent(in) :: nreg_in 


real(jprd), intent(in), dimension(nreg) :: op, frac_upper, frac_lower


real(jprd), intent(in) :: frac_threshold

real(jprd) :: overlap_matrix(nreg,nreg)
















end function calc_beta_overlap_matrix_dp_CPU

  pure function calc_alpha_overlap_matrix_CPU(nreg_in, op, op_inhom, &
&  frac_upper, frac_lower) result(overlap_matrix)
use parkind1
integer, intent(in) :: nreg_in 


real(jprb), intent(in) :: op, op_inhom


real(jprb), intent(in), dimension(nreg) :: frac_upper, frac_lower

real(jprb) :: overlap_matrix(nreg,nreg)













end function calc_alpha_overlap_matrix_CPU

  pure function calc_alpha_overlap_matrix_dp_CPU(nreg_in, op, op_inhom, &
&  frac_upper, frac_lower) result(overlap_matrix)
use parkind1
integer, intent(in) :: nreg_in 


real(jprd), intent(in) :: op, op_inhom


real(jprd), intent(in), dimension(nreg) :: frac_upper, frac_lower

real(jprd) :: overlap_matrix(nreg,nreg)














end function calc_alpha_overlap_matrix_dp_CPU

  pure function calc_alpha_overlap_matrix_simple_CPU(nreg_in, op, op_inhom, &
&  cf_upper, cf_lower) result(overlap_matrix)
use parkind1
integer, intent(in) :: nreg_in 


real(jprb), intent(in) :: op, op_inhom

real(jprb), intent(in) :: cf_upper, cf_lower

real(jprb) :: overlap_matrix(nreg,nreg)







end function calc_alpha_overlap_matrix_simple_CPU

  subroutine calc_overlap_matrices_CPU(nlev, nreg_in, is_clear_sky_layer, &
&     region_fracs, overlap_param, v_matrix, u_matrix, decorrelation_scaling, &
&     cloud_fraction_threshold, cloud_cov, use_beta_overlap)
use parkind1
use yomhook

integer,  intent(in) :: nlev, nreg_in

logical, intent(in) :: is_clear_sky_layer(0:nlev+1)



real(jprb), intent(in), dimension(nreg,nlev)  :: region_fracs


real(jprb), intent(in), dimension(nlev-1)  :: overlap_param

real(jprb), intent(out), dimension(nreg,nreg,nlev+1) :: v_matrix

real(jprb), intent(out), optional, dimension(nreg,nreg,nlev+1) &
&  :: u_matrix





real(jprb), intent(in), optional :: decorrelation_scaling

real(jprb), intent(in), optional :: cloud_fraction_threshold

real(jprb), intent(out), optional :: cloud_cov

logical, intent(in), optional :: use_beta_overlap
































 



end subroutine calc_overlap_matrices_CPU

  subroutine calc_overlap_matrices_dp_CPU(nlev, nreg_in, is_clear_sky_layer, &
&     region_fracs, overlap_param, v_matrix, u_matrix, decorrelation_scaling, &
&     cloud_fraction_threshold, cloud_cov, use_beta_overlap)
use parkind1
use yomhook

integer,  intent(in) :: nlev, nreg_in

logical, intent(in) :: is_clear_sky_layer(0:nlev+1)



real(jprd), intent(in), dimension(nreg,nlev)  :: region_fracs


real(jprd), intent(in), dimension(nlev-1)  :: overlap_param

real(jprd), intent(out), dimension(nreg,nreg,nlev+1) :: v_matrix

real(jprd), intent(out), optional, dimension(nreg,nreg,nlev+1) &
&  :: u_matrix





real(jprd), intent(in), optional :: decorrelation_scaling

real(jprd), intent(in), optional :: cloud_fraction_threshold

real(jprd), intent(out), optional :: cloud_cov

logical, intent(in), optional :: use_beta_overlap
































 



end subroutine calc_overlap_matrices_dp_CPU

end module radiation_overlap

