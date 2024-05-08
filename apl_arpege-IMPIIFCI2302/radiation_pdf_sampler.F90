! radiation_pdf_sampler.F90 - Get samples from a PDF for McICA
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

module radiation_pdf_sampler

  use parkind1

  implicit none
  

  !---------------------------------------------------------------------
  ! Derived type for sampling from a lognormal or gamma distribution,
  ! or other PDF, used to generate water content or optical depth
  ! scalings for use in the Monte Carlo Independent Column
  ! Approximation (McICA)
  type pdf_sampler_type
    ! Number of points in look-up table for cumulative distribution
    ! function (CDF) and fractional standard deviation (FSD)
    ! dimensions
    integer :: ncdf, nfsd

    ! First value of FSD and the reciprocal of the interval between
    ! FSD values (which are assumed to be uniformly distributed)
    real(jprb) :: fsd1, inv_fsd_interval

    ! Value of the distribution for each CDF and FSD bin
    real(jprb), allocatable, dimension(:,:) :: val

  contains

    
    
    
    
    
    

    
    
    
    

  procedure :: setup_GPU => setup_pdf_sampler_GPU

  procedure :: sample_GPU => sample_from_pdf_GPU

  procedure :: masked_sample_GPU => sample_from_pdf_masked_GPU

  procedure :: block_sample_GPU => sample_from_pdf_block_GPU

  procedure :: masked_block_sample_GPU => sample_from_pdf_masked_block_GPU

  procedure :: deallocate_GPU => deallocate_pdf_sampler_GPU

  procedure :: create_device_GPU

  procedure :: update_host_GPU

  procedure :: update_device_GPU

  procedure :: delete_device_GPU

  procedure :: setup_CPU => setup_pdf_sampler_CPU

  procedure :: sample_CPU => sample_from_pdf_CPU

  procedure :: masked_sample_CPU => sample_from_pdf_masked_CPU

  procedure :: block_sample_CPU => sample_from_pdf_block_CPU

  procedure :: masked_block_sample_CPU => sample_from_pdf_masked_block_CPU

  procedure :: deallocate_CPU => deallocate_pdf_sampler_CPU

  end type pdf_sampler_type

  private :: create_device_GPU, update_host_GPU, update_device_GPU, delete_device_GPU

contains

  !---------------------------------------------------------------------
  ! Load look-up table from a file
  

  !---------------------------------------------------------------------
  ! Deallocate data in pdf_sampler_type derived type
  


  !---------------------------------------------------------------------
  ! Extract the value from a PDF with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function value
  ! "cdf", and return it in val. Since this is an elemental
  ! subroutine, fsd, cdf and val may be arrays.
  


  !---------------------------------------------------------------------
  ! For true elements of mask, extract the values of a PDF with
  ! fractional standard deviation "fsd" corresponding to the
  ! cumulative distribution function values "cdf", and return in
  ! val. For false elements of mask, return zero in val.
  

  !---------------------------------------------------------------------
  ! Extract the values of a PDF with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function values
  ! "cdf", and return in val. This version works on 2D blocks of data.
  

  !---------------------------------------------------------------------
  ! Extract the values of a PDF with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function values
  ! "cdf", and return in val. This version works on 2D blocks of data.
  


  

  

  

  

  subroutine setup_pdf_sampler_GPU(this, file_name, iverbose, lacc)
use yomhook
use easy_netcdf
class(pdf_sampler_type), intent(inout) :: this
character(len=*),        intent(in)    :: file_name
integer, optional,       intent(in)    :: iverbose




logical, intent (in) :: lacc













end subroutine setup_pdf_sampler_GPU

  subroutine deallocate_pdf_sampler_GPU(this, lacc)
use yomhook
class(pdf_sampler_type), intent(inout) :: this

logical, intent (in) :: lacc



end subroutine deallocate_pdf_sampler_GPU

  elemental subroutine sample_from_pdf_GPU(this, fsd, cdf, val, lacc)
class(pdf_sampler_type), intent(in)  :: this


real(jprb),              intent(in)  :: fsd, cdf

real(jprb),              intent(out) :: val




logical, intent (in) :: lacc








end subroutine sample_from_pdf_GPU

  subroutine sample_from_pdf_masked_GPU(this, nsamp, fsd, cdf, val, mask, lacc)
class(pdf_sampler_type), intent(in)  :: this

integer,    intent(in) :: nsamp


real(jprb), intent(in)  :: fsd(nsamp), cdf(nsamp)

real(jprb), intent(out) :: val(:)

logical,    intent(in) :: mask(nsamp)






logical, intent (in) :: lacc

end subroutine sample_from_pdf_masked_GPU

  subroutine sample_from_pdf_block_GPU(this, nz, ng, fsd, cdf, val, lacc)
class(pdf_sampler_type), intent(in)  :: this

integer,    intent(in) :: nz, ng


real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

real(jprb), intent(out) :: val(:,:)






logical, intent (in) :: lacc

end subroutine sample_from_pdf_block_GPU

  subroutine sample_from_pdf_masked_block_GPU(this, nz, ng, fsd, cdf, val, mask, lacc)
class(pdf_sampler_type), intent(in)  :: this

integer,    intent(in) :: nz, ng


real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

real(jprb), intent(out) :: val(:,:)

logical,    intent(in), optional :: mask(nz)






logical, intent (in) :: lacc

end subroutine sample_from_pdf_masked_block_GPU

  subroutine create_device_GPU(this, lacc)
class(pdf_sampler_type), intent(inout) :: this
logical, intent (in) :: lacc

end subroutine create_device_GPU

  subroutine update_host_GPU(this, lacc)
class(pdf_sampler_type), intent(inout) :: this
logical, intent (in) :: lacc

end subroutine update_host_GPU

  subroutine update_device_GPU(this, lacc)
class(pdf_sampler_type), intent(inout) :: this
logical, intent (in) :: lacc

end subroutine update_device_GPU

  subroutine delete_device_GPU(this, lacc)
class(pdf_sampler_type), intent(inout) :: this
logical, intent (in) :: lacc

end subroutine delete_device_GPU

  subroutine setup_pdf_sampler_CPU(this, file_name, iverbose)
use yomhook
use easy_netcdf
class(pdf_sampler_type), intent(inout) :: this
character(len=*),        intent(in)    :: file_name
integer, optional,       intent(in)    :: iverbose

















end subroutine setup_pdf_sampler_CPU

  subroutine deallocate_pdf_sampler_CPU(this)
use yomhook
class(pdf_sampler_type), intent(inout) :: this




end subroutine deallocate_pdf_sampler_CPU

  elemental subroutine sample_from_pdf_CPU(this, fsd, cdf, val)
class(pdf_sampler_type), intent(in)  :: this


real(jprb),              intent(in)  :: fsd, cdf

real(jprb),              intent(out) :: val












end subroutine sample_from_pdf_CPU

  subroutine sample_from_pdf_masked_CPU(this, nsamp, fsd, cdf, val, mask)
class(pdf_sampler_type), intent(in)  :: this

integer,    intent(in) :: nsamp


real(jprb), intent(in)  :: fsd(nsamp), cdf(nsamp)

real(jprb), intent(out) :: val(:)

logical,    intent(in) :: mask(nsamp)







end subroutine sample_from_pdf_masked_CPU

  subroutine sample_from_pdf_block_CPU(this, nz, ng, fsd, cdf, val)
class(pdf_sampler_type), intent(in)  :: this

integer,    intent(in) :: nz, ng


real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

real(jprb), intent(out) :: val(:,:)







end subroutine sample_from_pdf_block_CPU

  subroutine sample_from_pdf_masked_block_CPU(this, nz, ng, fsd, cdf, val, mask)
class(pdf_sampler_type), intent(in)  :: this

integer,    intent(in) :: nz, ng


real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

real(jprb), intent(out) :: val(:,:)

logical,    intent(in), optional :: mask(nz)







end subroutine sample_from_pdf_masked_block_CPU

end module radiation_pdf_sampler

