! radiation_cloud_cover.F90 - Compute cumulative cloud cover for McICA
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
! Generate profiles of the cumulative cloud cover as seen from TOA,
! used in the McICA cloud generator.
!
! Modifications
!   2020-10-07  R. Hogan  Ensure iobj1 initialized in case of alpha_obj==0

#include "ecrad_config.h"

module radiation_cloud_cover

  use parkind1

  

  ! Three overlap schemes.  Note that "Exponential" means that
  ! clear-sky regions have no special significance for computing the
  ! cumulative cloud cover: non-contiguous clouds are exponentially
  ! rather than randomly overlapped. This is the situaition in the
  ! McRad radiation scheme at ECMWF.
  enum, bind(c)
    enumerator IOverlapMaximumRandom, IOverlapExponentialRandom, &
         &     IOverlapExponential
  end enum
  character(len=*), parameter :: OverlapName(0:2) = (/ 'Max-Ran', &
       &                                               'Exp-Ran', &
       &                                               'Exp-Exp' /)

  ! Maximum cloud fraction distinguishable from 1
  real(jprb), parameter :: MaxCloudFrac = 1.0_jprb-epsilon(1.0_jprb)*10.0_jprb


contains

  !---------------------------------------------------------------------
  ! Convert "beta" overlap parameter of Shonk et al. (2010) to "alpha"
  ! overlap parameter of Hogan and Illingworth (2000)
  


  !---------------------------------------------------------------------
  ! Compute total cloud cover according to the specified overlap
  ! rule. This can be used to compute the high, mid and low cloud
  ! cover by passing in subsets of the cloud fraction array
  


  !---------------------------------------------------------------------
  ! Maximum-random overlap: Geleyn & Hollingsworth formula
  


  !---------------------------------------------------------------------
  ! Exponential-random overlap: exponential overlap for contiguous
  ! clouds, random overlap for non-contiguous clouds
  



  !---------------------------------------------------------------------
  ! Exponential-exponential overlap: exponential overlap for both
  ! contiguous and non-contiguous clouds. This is the result of the
  ! simple Raisanen cloud generator, but unfortunately it has no
  ! (known) analytic formula for the total cloud cover, or the
  ! cumulative cloud cover.  In partially cloudy columns, The McICA
  ! scheme needs this info in order to devote all the cloudy g-points
  ! to columns containing cloud, which reduces McICA noise. The
  ! following routine provides an approximate estimate of cumulative
  ! cloud cover consistent with the exponential-exponential scheme.
  

  elemental function beta2alpha_GPU(beta, frac1, frac2)
implicit none


real(jprb), intent(in) :: beta, frac1, frac2
real(jprb)             :: beta2alpha_GPU




end function beta2alpha_GPU

  function cloud_cover_GPU(nlev, i_overlap_scheme, frac, overlap_param, &
&               is_beta_overlap, lacc)
implicit none

integer, intent(in)    :: nlev, i_overlap_scheme


real(jprb), intent(in) :: frac(nlev), overlap_param(nlev-1)


logical, intent(in), optional :: is_beta_overlap

real(jprb)             :: cloud_cover_GPU




logical, intent (in) :: lacc


end function cloud_cover_GPU

  subroutine cum_cloud_cover_max_ran_GPU(nlev, frac, &
& cum_cloud_cover, pair_cloud_cover, lacc)
use yomhook
implicit none

integer, intent(in)     :: nlev  

real(jprb), intent(in)  :: frac(nlev)


real(jprb), intent(out) :: cum_cloud_cover(nlev)

real(jprb), intent(out) :: pair_cloud_cover(nlev-1)






logical, intent (in) :: lacc







end subroutine cum_cloud_cover_max_ran_GPU

  subroutine cum_cloud_cover_exp_ran_GPU(nlev, frac, overlap_param, &
& cum_cloud_cover, pair_cloud_cover, is_beta_overlap, lacc)
use yomhook
implicit none

integer, intent(in)     :: nlev  

real(jprb), intent(in)  :: frac(nlev)



real(jprb), intent(in)  :: overlap_param(nlev-1)





logical, intent(in), optional :: is_beta_overlap


real(jprb), intent(out) :: cum_cloud_cover(nlev)

real(jprb), intent(out) :: pair_cloud_cover(nlev-1)









logical, intent (in) :: lacc








end subroutine cum_cloud_cover_exp_ran_GPU

  subroutine cum_cloud_cover_exp_exp_GPU(nlev, frac, overlap_param, &
& cum_cloud_cover, pair_cloud_cover, is_beta_overlap, lacc)
use yomhook
implicit none

integer, intent(in)     :: nlev  

real(jprb), intent(in)  :: frac(nlev)



real(jprb), intent(in)  :: overlap_param(nlev-1)





logical, intent(in), optional :: is_beta_overlap


real(jprb), intent(out) :: cum_cloud_cover(nlev)

real(jprb), intent(out) :: pair_cloud_cover(nlev-1)






































logical, intent (in) :: lacc










 

end subroutine cum_cloud_cover_exp_exp_GPU

  elemental function beta2alpha_CPU(beta, frac1, frac2)
implicit none


real(jprb), intent(in) :: beta, frac1, frac2
real(jprb)             :: beta2alpha_CPU



end function beta2alpha_CPU

  function cloud_cover_CPU(nlev, i_overlap_scheme, frac, overlap_param, &
&               is_beta_overlap)
implicit none

integer, intent(in)    :: nlev, i_overlap_scheme


real(jprb), intent(in) :: frac(nlev), overlap_param(nlev-1)


logical, intent(in), optional :: is_beta_overlap

real(jprb)             :: cloud_cover_CPU






end function cloud_cover_CPU

  subroutine cum_cloud_cover_max_ran_CPU(nlev, frac, &
& cum_cloud_cover, pair_cloud_cover)
use yomhook
implicit none

integer, intent(in)     :: nlev  

real(jprb), intent(in)  :: frac(nlev)


real(jprb), intent(out) :: cum_cloud_cover(nlev)

real(jprb), intent(out) :: pair_cloud_cover(nlev-1)













end subroutine cum_cloud_cover_max_ran_CPU

  subroutine cum_cloud_cover_exp_ran_CPU(nlev, frac, overlap_param, &
& cum_cloud_cover, pair_cloud_cover, is_beta_overlap)
use yomhook
implicit none

integer, intent(in)     :: nlev  

real(jprb), intent(in)  :: frac(nlev)



real(jprb), intent(in)  :: overlap_param(nlev-1)





logical, intent(in), optional :: is_beta_overlap


real(jprb), intent(out) :: cum_cloud_cover(nlev)

real(jprb), intent(out) :: pair_cloud_cover(nlev-1)

















end subroutine cum_cloud_cover_exp_ran_CPU

  subroutine cum_cloud_cover_exp_exp_CPU(nlev, frac, overlap_param, &
& cum_cloud_cover, pair_cloud_cover, is_beta_overlap)
use yomhook
implicit none

integer, intent(in)     :: nlev  

real(jprb), intent(in)  :: frac(nlev)



real(jprb), intent(in)  :: overlap_param(nlev-1)





logical, intent(in), optional :: is_beta_overlap


real(jprb), intent(out) :: cum_cloud_cover(nlev)

real(jprb), intent(out) :: pair_cloud_cover(nlev-1)
















































 

end subroutine cum_cloud_cover_exp_exp_CPU

end module radiation_cloud_cover

