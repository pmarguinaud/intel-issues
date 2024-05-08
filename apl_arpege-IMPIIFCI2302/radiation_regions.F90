! radiation_regions.F90 -- Properties of horizontal regions in Tripleclouds & SPARTACUS
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
! Modifications
!   2017-07-14  R. Hogan  Incorporate gamma distribution option
!   2018-10-06  R. Hogan  Merged from radiation_optical_depth_scaling.h and radiation_overlap.F90

module radiation_regions

  implicit none

  

contains

  !---------------------------------------------------------------------
  ! Compute the optical depth scalings for the optically "thick" and
  ! "thin" regions of a Tripleclouds representation of a sub-grid PDF
  ! of cloud optical depth. Following Shonk and Hogan (2008), the 16th
  ! percentile is used for the thin region, and the formulas estimate
  ! this for both lognormal and gamma distributions. However, an
  ! adjustment is needed for the gamma distribution at large
  ! fractional standard deviations.
  

  subroutine calc_region_properties_GPU(nlev, nreg, istartcol, iendcol, do_gamma, &
&  cloud_fraction, frac_std, reg_fracs, od_scaling, cloud_fraction_threshold, lacc)
use parkind1
use yomhook
use radiation_io


















integer, intent(in) :: nlev, nreg

integer, intent(in) :: istartcol, iendcol

logical, intent(in) :: do_gamma


real(jprb), intent(in), dimension(:,:)  :: cloud_fraction 

real(jprb), intent(in), dimension(:,:)  :: frac_std       

real(jprb), intent(out) :: reg_fracs(1:nreg,nlev,istartcol:iendcol)

real(jprb), intent(out) :: od_scaling(2:nreg,nlev,istartcol:iendcol)

real(jprb), intent(in), optional :: cloud_fraction_threshold






logical, intent (in) :: lacc




end subroutine calc_region_properties_GPU

  subroutine calc_region_properties_CPU(nlev, nreg, istartcol, iendcol, do_gamma, &
&  cloud_fraction, frac_std, reg_fracs, od_scaling, cloud_fraction_threshold)
use parkind1
use yomhook
use radiation_io


















integer, intent(in) :: nlev, nreg

integer, intent(in) :: istartcol, iendcol

logical, intent(in) :: do_gamma


real(jprb), intent(in), dimension(:,:)  :: cloud_fraction 

real(jprb), intent(in), dimension(:,:)  :: frac_std       

real(jprb), intent(out) :: reg_fracs(1:nreg,nlev,istartcol:iendcol)

real(jprb), intent(out) :: od_scaling(2:nreg,nlev,istartcol:iendcol)

real(jprb), intent(in), optional :: cloud_fraction_threshold










end subroutine calc_region_properties_CPU

end module radiation_regions


