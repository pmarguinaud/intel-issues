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

  use parkind1, only : jprb

  implicit none
  public

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

    procedure :: setup => setup_pdf_sampler
    procedure :: sample => sample_from_pdf
    procedure :: masked_sample => sample_from_pdf_masked
    procedure :: block_sample => sample_from_pdf_block
    procedure :: masked_block_sample => sample_from_pdf_masked_block
    procedure :: deallocate => deallocate_pdf_sampler

#ifdef _OPENACC
    procedure :: create_device
    procedure :: update_host
    procedure :: update_device
    procedure :: delete_device
#endif

  end type pdf_sampler_type

contains

  !---------------------------------------------------------------------
  ! Load look-up table from a file
  subroutine setup_pdf_sampler(this, file_name, iverbose)

    use yomhook,     only : lhook, dr_hook, jphook
    use easy_netcdf, only : netcdf_file

    class(pdf_sampler_type), intent(inout) :: this
    character(len=*),        intent(in)    :: file_name
    integer, optional,       intent(in)    :: iverbose

    type(netcdf_file)  :: file
    integer            :: iverb
    real(jprb), allocatable :: fsd(:)

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_pdf_sampler:setup',0,hook_handle)

    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    if (allocated(this%val)) then
      deallocate(this%val)
    end if

    call file%open(trim(file_name), iverbose=iverb)

    call file%get('fsd',fsd)
    call file%get('x',  this%val)

    call file%close()

    this%ncdf = size(this%val,1)
    this%nfsd = size(this%val,2)
    this%fsd1 = fsd(1)
    this%inv_fsd_interval = 1.0_jprb / (fsd(2)-fsd(1))

    deallocate(fsd)

    if (lhook) call dr_hook('radiation_pdf_sampler:setup',1,hook_handle)

  end subroutine setup_pdf_sampler

  !---------------------------------------------------------------------
  ! Deallocate data in pdf_sampler_type derived type
  subroutine deallocate_pdf_sampler(this)

    use yomhook,     only : lhook, dr_hook, jphook

    class(pdf_sampler_type), intent(inout) :: this
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_pdf_sampler:deallocate',0,hook_handle)

    if (allocated(this%val)) then
      deallocate(this%val)
    end if

    if (lhook) call dr_hook('radiation_pdf_sampler:deallocate',1,hook_handle)

  end subroutine deallocate_pdf_sampler


  !---------------------------------------------------------------------
  ! Extract the value from a PDF with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function value
  ! "cdf", and return it in val. Since this is an elemental
  ! subroutine, fsd, cdf and val may be arrays.
  elemental subroutine sample_from_pdf(this, fsd, cdf, val)

    class(pdf_sampler_type), intent(in)  :: this

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb),              intent(in)  :: fsd, cdf

    ! Sample from distribution
    real(jprb),              intent(out) :: val

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    ! Bilinear interpolation with bounds
    wcdf = cdf * (this%ncdf-1) + 1.0_jprb
    icdf = max(1, min(int(wcdf), this%ncdf-1))
    wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

    wfsd = (fsd-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
    ifsd = max(1, min(int(wfsd), this%nfsd-1))
    wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

    val =    (1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
         & + (1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
         & +           wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
         & +           wcdf *          wfsd  * this%val(icdf+1,ifsd+1)

  end subroutine sample_from_pdf


  !---------------------------------------------------------------------
  ! For true elements of mask, extract the values of a PDF with
  ! fractional standard deviation "fsd" corresponding to the
  ! cumulative distribution function values "cdf", and return in
  ! val. For false elements of mask, return zero in val.
  subroutine sample_from_pdf_masked(this, nsamp, fsd, cdf, val, mask)

    class(pdf_sampler_type), intent(in)  :: this

    ! Number of samples
    integer,    intent(in) :: nsamp

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nsamp), cdf(nsamp)

    ! Sample from distribution
    real(jprb), intent(out) :: val(:)

    ! Mask
    logical,    intent(in) :: mask(nsamp)

    ! Loop index
    integer    :: jsamp

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    do jsamp = 1,nsamp
      if (mask(jsamp)) then
        ! Bilinear interpolation with bounds
        wcdf = cdf(jsamp) * (this%ncdf-1) + 1.0_jprb
        icdf = max(1, min(int(wcdf), this%ncdf-1))
        wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

        wfsd = (fsd(jsamp)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
        ifsd = max(1, min(int(wfsd), this%nfsd-1))
        wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

        val(jsamp)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
             &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
             &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
             &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
      else
        val(jsamp) = 0.0_jprb
      end if
    end do
  end subroutine sample_from_pdf_masked

  !---------------------------------------------------------------------
  ! Extract the values of a PDF with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function values
  ! "cdf", and return in val. This version works on 2D blocks of data.
  subroutine sample_from_pdf_block(this, nz, ng, fsd, cdf, val)

    class(pdf_sampler_type), intent(in)  :: this

    ! Number of samples
    integer,    intent(in) :: nz, ng

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

    ! Sample from distribution
    real(jprb), intent(out) :: val(:,:)

    ! Loop index
    integer    :: jz, jg

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    do jz = 1,nz
      do jg = 1,ng
        if (cdf(jg, jz) > 0.0_jprb) then
          ! Bilinear interpolation with bounds
          wcdf = cdf(jg,jz) * (this%ncdf-1) + 1.0_jprb
          icdf = max(1, min(int(wcdf), this%ncdf-1))
          wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

          wfsd = (fsd(jz)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
          ifsd = max(1, min(int(wfsd), this%nfsd-1))
          wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

          val(jg,jz)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
               &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
               &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
               &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
        else
          val(jg,jz) = 0.0_jprb
        end if
      end do
    end do

  end subroutine sample_from_pdf_block

  !---------------------------------------------------------------------
  ! Extract the values of a PDF with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function values
  ! "cdf", and return in val. This version works on 2D blocks of data.
  subroutine sample_from_pdf_masked_block(this, nz, ng, fsd, cdf, val, mask)

    class(pdf_sampler_type), intent(in)  :: this

    ! Number of samples
    integer,    intent(in) :: nz, ng

    ! Fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

    ! Sample from distribution
    real(jprb), intent(out) :: val(:,:)

    ! Mask
    logical,    intent(in), optional :: mask(nz)

    ! Loop index
    integer    :: jz, jg

    ! Index to look-up table
    integer    :: ifsd, icdf

    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    do jz = 1,nz

      if (mask(jz)) then

        do jg = 1,ng
          if (cdf(jg, jz) > 0.0_jprb) then
            ! Bilinear interpolation with bounds
            wcdf = cdf(jg,jz) * (this%ncdf-1) + 1.0_jprb
            icdf = max(1, min(int(wcdf), this%ncdf-1))
            wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

            wfsd = (fsd(jz)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
            ifsd = max(1, min(int(wfsd), this%nfsd-1))
            wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

            val(jg,jz)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
                 &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
                 &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
                 &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
          else
            val(jg,jz) = 0.0_jprb
          end if
        end do

      end if

    end do

  end subroutine sample_from_pdf_masked_block

#ifdef _OPENACC

  subroutine create_device(this)
    class(pdf_sampler_type), intent(inout) :: this

    !$ACC ENTER DATA COPYIN(this%val) IF(allocated(this%val)) ASYNC(1)
  end subroutine create_device

  subroutine update_host(this)
    class(pdf_sampler_type), intent(inout) :: this

    !$ACC UPDATE HOST(this%val) IF(allocated(this%val)) ASYNC(1)
  end subroutine update_host

  subroutine update_device(this)
    class(pdf_sampler_type), intent(inout) :: this

    !$ACC UPDATE DEVICE(this%val) IF(allocated(this%val)) ASYNC(1)
  end subroutine update_device

  subroutine delete_device(this)
    class(pdf_sampler_type), intent(inout) :: this

    !$ACC EXIT DATA DELETE(this%val) IF(allocated(this%val)) ASYNC(1)
  end subroutine delete_device
#endif

end module radiation_pdf_sampler
