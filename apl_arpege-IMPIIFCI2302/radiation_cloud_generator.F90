! radiation_cloud_generator.F90 - Generate water-content or optical-depth scalings for McICA
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
!
! Modifications
!   2018-02-22  R. Hogan  Call masked version of PDF sampler for speed
!   2020-03-31  R. Hogan  More vectorizable version of Exp-Ran

module radiation_cloud_generator

  

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
  


  !---------------------------------------------------------------------
  ! Generate a column of optical depth scalings using
  ! exponential-random overlap (which includes maximum-random overlap
  ! as a limiting case)
  


  !---------------------------------------------------------------------
  ! Generate a column of optical depth scalings using
  ! exponential-exponential overlap
  


  !---------------------------------------------------------------------
  ! Extract the value of a lognormal distribution with fractional
  ! standard deviation "fsd" corresponding to the cumulative
  ! distribution function value "cdf", and return it in x. Since this
  ! is an elemental subroutine, fsd, cdf and x may be arrays. SIMD version.
  


  !---------------------------------------------------------------------
  ! Generate columns of optical depth scalings using
  ! exponential-random overlap (which includes maximum-random overlap
  ! as a limiting case).  This version is intended to work better on
  ! hardware with long vector lengths.  As with all calculations in
  ! this file, we zoom into the fraction of the column with cloud at
  ! any height, so that all spectral intervals see a cloud somewhere.
  ! In the McICA solver, this is combined appropriately with the
  ! clear-sky calculation.
  

  subroutine cloud_generator_GPU(ng, nlev, i_overlap_scheme, &
&  iseed, frac_threshold, &
&  frac, overlap_param, decorrelation_scaling, &
&  fractional_std, pdf_sampler, &
&  od_scaling, total_cloud_cover, &
&  use_beta_overlap, use_vectorizable_generator, lacc)
use parkind1
use yomhook
use radiation_io
use random_numbers_mix
use radiation_pdf_sampler
use radiation_cloud_cover
implicit none

integer, intent(in)     :: ng    
integer, intent(in)     :: nlev  
integer, intent(in)     :: i_overlap_scheme
integer, intent(in)     :: iseed 


real(jprb), intent(in)  :: frac_threshold

real(jprb), intent(in)  :: frac(nlev)



real(jprb), intent(in)  :: overlap_param(nlev-1)

real(jprb), intent(in)  :: decorrelation_scaling

real(jprb), intent(in)  :: fractional_std(nlev)

type(pdf_sampler_type), intent(in) :: pdf_sampler





logical, intent(in), optional :: use_beta_overlap


logical, intent(in), optional :: use_vectorizable_generator


real(jprb), intent(out) :: od_scaling(ng,nlev)

real(jprb), intent(out) :: total_cloud_cover






















logical, intent (in) :: lacc


;



end subroutine cloud_generator_GPU

  subroutine generate_column_exp_ran_GPU(ng, nlev, ig, random_stream, pdf_sampler, &
&  frac, pair_cloud_cover, &
&  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
&  itrigger, iend, od_scaling, lacc)
use parkind1
use radiation_pdf_sampler
use random_numbers_mix
implicit none

integer, intent(in) :: ng, ig

integer, intent(in) :: nlev

type(randomnumberstream), intent(inout) :: random_stream

type(pdf_sampler_type), intent(in) :: pdf_sampler


real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std


real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom


integer, intent(in) :: itrigger, iend

real(jprb), intent(inout), dimension(ng,nlev) :: od_scaling










logical, intent (in) :: lacc









end subroutine generate_column_exp_ran_GPU

  subroutine generate_column_exp_exp_GPU(ng, nlev, ig, random_stream, pdf_sampler, &
&  frac, pair_cloud_cover, &
&  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
&  itrigger, iend, od_scaling, lacc)
use parkind1
use radiation_pdf_sampler
use random_numbers_mix
implicit none

integer, intent(in) :: ng, ig

integer, intent(in) :: nlev

type(randomnumberstream), intent(inout) :: random_stream

type(pdf_sampler_type), intent(in) :: pdf_sampler


real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std


real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom


integer, intent(in) :: itrigger, iend

real(jprb), intent(inout), dimension(ng,nlev) :: od_scaling











logical, intent (in) :: lacc


































end subroutine generate_column_exp_exp_GPU

  subroutine sample_from_pdf_simd_GPU(this, fsd, cdf, x, lacc)
use parkind1
use radiation_pdf_sampler
implicit none
type(pdf_sampler_type), intent(in)  :: this


real(jprb),              intent(in)  :: fsd, cdf

real(jprb),              intent(out) :: x




logical, intent (in) :: lacc








end subroutine sample_from_pdf_simd_GPU

  subroutine generate_columns_exp_ran_GPU(ng, nlev, iseed, pdf_sampler, &
&  total_cloud_cover, frac_threshold, frac, pair_cloud_cover, &
&  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
&  ibegin, iend, od_scaling, lacc)
use parkind1
use radiation_pdf_sampler
use radiation_random_numbers
implicit none

integer, intent(in) :: ng

integer, intent(in) :: nlev
integer, intent(in) :: iseed 




type(pdf_sampler_type), intent(in) :: pdf_sampler

real(jprb), intent(in) :: total_cloud_cover
real(jprb), intent(in) :: frac_threshold


real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std


real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom

integer, intent(inout) :: ibegin, iend

real(jprb), intent(inout), dimension(ng,nlev) :: od_scaling








    
  
 

logical, intent (in) :: lacc 

























end subroutine generate_columns_exp_ran_GPU

  subroutine cloud_generator_CPU(ng, nlev, i_overlap_scheme, &
&  iseed, frac_threshold, &
&  frac, overlap_param, decorrelation_scaling, &
&  fractional_std, pdf_sampler, &
&  od_scaling, total_cloud_cover, &
&  use_beta_overlap, use_vectorizable_generator)
use parkind1
use yomhook
use radiation_io
use random_numbers_mix
use radiation_pdf_sampler
use radiation_cloud_cover
implicit none

integer, intent(in)     :: ng    
integer, intent(in)     :: nlev  
integer, intent(in)     :: i_overlap_scheme
integer, intent(in)     :: iseed 


real(jprb), intent(in)  :: frac_threshold

real(jprb), intent(in)  :: frac(nlev)



real(jprb), intent(in)  :: overlap_param(nlev-1)

real(jprb), intent(in)  :: decorrelation_scaling

real(jprb), intent(in)  :: fractional_std(nlev)

type(pdf_sampler_type), intent(in) :: pdf_sampler





logical, intent(in), optional :: use_beta_overlap


logical, intent(in), optional :: use_vectorizable_generator


real(jprb), intent(out) :: od_scaling(ng,nlev)

real(jprb), intent(out) :: total_cloud_cover
























;



end subroutine cloud_generator_CPU

  subroutine generate_column_exp_ran_CPU(ng, nlev, ig, random_stream, pdf_sampler, &
&  frac, pair_cloud_cover, &
&  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
&  itrigger, iend, od_scaling)
use parkind1
use radiation_pdf_sampler
use random_numbers_mix
implicit none

integer, intent(in) :: ng, ig

integer, intent(in) :: nlev

type(randomnumberstream), intent(inout) :: random_stream

type(pdf_sampler_type), intent(in) :: pdf_sampler


real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std


real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom


integer, intent(in) :: itrigger, iend

real(jprb), intent(inout), dimension(ng,nlev) :: od_scaling



















end subroutine generate_column_exp_ran_CPU

  subroutine generate_column_exp_exp_CPU(ng, nlev, ig, random_stream, pdf_sampler, &
&  frac, pair_cloud_cover, &
&  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
&  itrigger, iend, od_scaling)
use parkind1
use radiation_pdf_sampler
use random_numbers_mix
implicit none

integer, intent(in) :: ng, ig

integer, intent(in) :: nlev

type(randomnumberstream), intent(inout) :: random_stream

type(pdf_sampler_type), intent(in) :: pdf_sampler


real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std


real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom


integer, intent(in) :: itrigger, iend

real(jprb), intent(inout), dimension(ng,nlev) :: od_scaling













































end subroutine generate_column_exp_exp_CPU

  subroutine sample_from_pdf_simd_CPU(this, fsd, cdf, x)
use parkind1
use radiation_pdf_sampler
implicit none
type(pdf_sampler_type), intent(in)  :: this


real(jprb),              intent(in)  :: fsd, cdf

real(jprb),              intent(out) :: x












end subroutine sample_from_pdf_simd_CPU

  subroutine generate_columns_exp_ran_CPU(ng, nlev, iseed, pdf_sampler, &
&  total_cloud_cover, frac_threshold, frac, pair_cloud_cover, &
&  cum_cloud_cover, overhang, fractional_std, overlap_param_inhom, &
&  ibegin, iend, od_scaling)
use parkind1
use radiation_pdf_sampler
use radiation_random_numbers
implicit none

integer, intent(in) :: ng

integer, intent(in) :: nlev
integer, intent(in) :: iseed 




type(pdf_sampler_type), intent(in) :: pdf_sampler

real(jprb), intent(in) :: total_cloud_cover
real(jprb), intent(in) :: frac_threshold


real(jprb), intent(in), dimension(nlev) :: frac, cum_cloud_cover, fractional_std


real(jprb), intent(in), dimension(nlev-1) :: pair_cloud_cover, overhang

real(jprb), intent(in), dimension(nlev-1) :: overlap_param_inhom

integer, intent(inout) :: ibegin, iend

real(jprb), intent(inout), dimension(ng,nlev) :: od_scaling








    
  
 
 

























end subroutine generate_columns_exp_ran_CPU

end module radiation_cloud_generator

