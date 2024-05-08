! radiation_cloud_generator_acc.F90 - Generate water-content or optical-depth scalings for McICA
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
! Generate clouds for McICA using a method modified from Raisanen et
! al. (2002)
! This is a copy of the original cloud_generator, that is better suited for OpenACC
!
! Modifications
!   2018-02-22  R. Hogan  Call masked version of PDF sampler for speed
!   2020-03-31  R. Hogan  More vectorizable version of Exp-Ran
!   2022-11-07  D. Hupp adaptation for ACC

module radiation_cloud_generator_acc

  

contains

  !---------------------------------------------------------------------
  ! Generate scaling factors for the cloud optical depth to represent
  ! cloud overlap, the overlap of internal cloud inhomogeneities and
  ! the fractional standard deviation of these inhomogeneities, for
  ! use in a Monte Carlo Independent Column Approximation radiation
  ! scheme. All returned profiles contain cloud, and the total cloud
  ! cover is also returned, so the calling function can then do a
  ! weighted average of clear and cloudy skies; this is a way to
  ! reduce the Monte Carlo noise in profiles with low cloud cover.
  

  subroutine cloud_generator_acc_GPU(ng, nlev, &
&  iseed, frac_threshold, &
&  frac, overlap_param, decorrelation_scaling, &
&  fractional_std, &
&  sample_ncdf, sample_nfsd, sample_fsd1, &
&  sample_inv_fsd_interval, sample_val, &
&  od_scaling, total_cloud_cover, &
&  ibegin, iend, &
&  cum_cloud_cover, &
&  pair_cloud_cover)
use parkind1
use radiation_random_numbers
use radiation_random_numbers
implicit none

integer, intent(in)     :: ng    
integer, intent(in)     :: nlev  
integer, intent(in)     :: iseed 


real(jprb), intent(in)  :: frac_threshold

real(jprb), intent(in)  :: frac(nlev)



real(jprb), intent(in)  :: overlap_param(nlev-1)

real(jprb), intent(in)  :: decorrelation_scaling

real(jprb), intent(in)  :: fractional_std(nlev)

integer, intent(in)  :: sample_ncdf, sample_nfsd
real(jprb), intent(in)  :: sample_fsd1, sample_inv_fsd_interval
real(jprb), intent(in), dimension(:,:)  :: sample_val


real(jprb), intent(out) :: od_scaling(ng,nlev)

real(jprb), intent(in) :: total_cloud_cover


real(jprb), intent(in) :: cum_cloud_cover(nlev)

integer, intent(in) :: ibegin, iend












real(jprb), intent(inout), dimension(nlev-1) :: pair_cloud_cover















 
end subroutine cloud_generator_acc_GPU

  subroutine cloud_generator_acc_CPU(ng, nlev, &
&  iseed, frac_threshold, &
&  frac, overlap_param, decorrelation_scaling, &
&  fractional_std, &
&  sample_ncdf, sample_nfsd, sample_fsd1, &
&  sample_inv_fsd_interval, sample_val, &
&  od_scaling, total_cloud_cover, &
&  ibegin, iend, &
&  cum_cloud_cover, &
&  pair_cloud_cover)
use parkind1
use radiation_random_numbers
use radiation_random_numbers
implicit none

integer, intent(in)     :: ng    
integer, intent(in)     :: nlev  
integer, intent(in)     :: iseed 


real(jprb), intent(in)  :: frac_threshold

real(jprb), intent(in)  :: frac(nlev)



real(jprb), intent(in)  :: overlap_param(nlev-1)

real(jprb), intent(in)  :: decorrelation_scaling

real(jprb), intent(in)  :: fractional_std(nlev)

integer, intent(in)  :: sample_ncdf, sample_nfsd
real(jprb), intent(in)  :: sample_fsd1, sample_inv_fsd_interval
real(jprb), intent(in), dimension(:,:)  :: sample_val


real(jprb), intent(out) :: od_scaling(ng,nlev)

real(jprb), intent(in) :: total_cloud_cover


real(jprb), intent(in) :: cum_cloud_cover(nlev)

integer, intent(in) :: ibegin, iend












real(jprb), intent(inout), dimension(nlev-1) :: pair_cloud_cover














 
end subroutine cloud_generator_acc_CPU

end module radiation_cloud_generator_acc

