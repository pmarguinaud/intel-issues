! radiation_flux.F90 - Derived type to store the output fluxes
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
!   2017-09-08  R. Hogan  Store g-point fluxes
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2019-01-08  R. Hogan  Added "indexed_sum_profile"
!   2019-01-14  R. Hogan  out_of_physical_bounds calls routine in radiation_config
!   2021-01-20  R. Hogan  Added heating_rate_out_of_physical_bounds function
!   2022-12-07  R. Hogan  Added top-of-atmosphere spectral output

#include "ecrad_config.h"

module radiation_flux

  use parkind1

  implicit none
  

  !---------------------------------------------------------------------
  ! This derived type contains the output from the radiation
  ! calculation.  Currently this is solely flux profiles, but in
  ! future surface fluxes in each band may be stored in order that the
  ! calling program can compute surface-radiation such as
  ! photosynthetically active radiation and UV index.
  type flux_type
     ! All the following are broad-band fluxes in W m-2 with
     ! dimensions (ncol,nlev+1).  Note that only those fluxes that are
     ! requested will be used, so clear-sky and direct-beam arrays may
     ! not be allocated
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_up, lw_dn, &   ! Upwelling and downwelling longwave
          &  sw_up, sw_dn, &   ! Upwelling and downwelling shortwave
          &  sw_dn_direct, &   ! Direct-beam shortwave into a horizontal plane
          &  lw_up_clear, lw_dn_clear, & ! Clear-sky quantities...
          &  sw_up_clear, sw_dn_clear, &
          &  sw_dn_direct_clear
     ! As above but fluxes in each spectral band in W m-2 with
     ! dimensions (nband,ncol,nlev+1).  These are only allocated if
     ! config%do_save_spectral_flux==.true.
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  lw_up_band, lw_dn_band, &   ! Upwelling and downwelling longwave
          &  sw_up_band, sw_dn_band, &   ! Upwelling and downwelling shortwave
          &  sw_dn_direct_band, &        ! Direct-beam shortwave
          &  lw_up_clear_band, lw_dn_clear_band, & ! Clear-sky quantities...
          &  sw_up_clear_band, sw_dn_clear_band, &
          &  sw_dn_direct_clear_band
     ! Surface downwelling quantities at each g point, dimensioned
     ! (ng,ncol), that are always saved by the solver, except for the
     ! clear-sky ones that are only produced if
     ! config%do_clear==.true.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_dn_surf_g, lw_dn_surf_clear_g, &
          &  sw_dn_diffuse_surf_g, sw_dn_direct_surf_g, &
          &  sw_dn_diffuse_surf_clear_g, sw_dn_direct_surf_clear_g
     ! Top-of-atmosphere quantities at each g point, dimensioned
     ! (ng,ncol), that are always saved by the solver, except for the
     ! clear-sky ones that are only produced if
     ! config%do_clear==.true.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_up_toa_g, lw_up_toa_clear_g, &
          &  sw_dn_toa_g, sw_up_toa_g, sw_up_toa_clear_g
     ! Shortwave downwelling spectral fluxes in W m-2 at the surface,
     ! from which quantities such as photosynthetically active and UV
     ! radiation can be computed. Only allocated if
     ! config%do_surface_sw_spectral_flux==.true.  Note that the
     ! clear-sky quantities are only computed if
     ! config%do_clear==.true., but direct fluxes are computed whether
     ! or not do_direct==.true.. The dimensions are (nband,ncol).
     real(jprb), allocatable, dimension(:,:) :: &
          &  sw_dn_surf_band, sw_dn_direct_surf_band, &
          &  sw_dn_surf_clear_band, sw_dn_direct_surf_clear_band
     ! Top-of-atmosphere spectral fluxes in W m-2. Only allocated if
     ! config%do_toa_spectral_flux=.true.. Note that the clear-sky
     ! quantities are only computed if config%do_clear==.true.. The
     ! dimensions are (nband,ncol).
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_up_toa_band, lw_up_toa_clear_band, &
          &  sw_dn_toa_band, sw_up_toa_band, sw_up_toa_clear_band
     ! Surface downwelling fluxes in W m-2 at the spectral resolution
     ! needed by any subsequent canopy radiative transfer.  If
     ! config%use_canopy_full_spectrum_[sw|lw] then these will be at
     ! g-point resolution; otherwise they will be at
     ! config%n_albedo_bands and config%n_emiss_bands resolution.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_dn_surf_canopy, &
          &  sw_dn_diffuse_surf_canopy, sw_dn_direct_surf_canopy

     ! Diagnosed cloud cover from the short- and long-wave solvers
     real(jprb), allocatable, dimension(:) :: &
          &  cloud_cover_lw, cloud_cover_sw
     ! Longwave derivatives needed by Hogan and Bozzo (2015) method
     ! for approximate longwave updates in between the full radiation
     ! calls: rate of change of upwelling broad-band flux with respect
     ! to surface value, dimensioned (ncol,nlev+1)
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_derivatives

   contains
    
    
    
    
    
    
    
    
    
    
  procedure :: allocate_GPU   => allocate_flux_type_GPU

  procedure :: deallocate_GPU => deallocate_flux_type_GPU

  procedure :: calc_surface_spectral_GPU

  procedure :: calc_toa_spectral_GPU

  procedure :: out_of_physical_bounds_GPU

  procedure :: heating_rate_out_of_physical_bounds_GPU

  procedure :: create_device_GPU

  procedure :: update_host_GPU

  procedure :: update_device_GPU

  procedure :: delete_device_GPU

  procedure :: allocate_CPU   => allocate_flux_type_CPU

  procedure :: deallocate_CPU => deallocate_flux_type_CPU

  procedure :: calc_surface_spectral_CPU

  procedure :: calc_toa_spectral_CPU

  procedure :: out_of_physical_bounds_CPU

  procedure :: heating_rate_out_of_physical_bounds_CPU

  end type flux_type

  private :: create_device_GPU, update_host_GPU, update_device_GPU, delete_device_GPU

! Added for DWD (2020)
      logical, parameter :: use_indexed_sum_vec = .false.

contains

  !---------------------------------------------------------------------
  ! Allocate arrays for flux profiles, using config to define which
  ! fluxes are needed.  The arrays are dimensioned for columns between
  ! istartcol, iendcol and levels from 1 to nlev+1
  


  !---------------------------------------------------------------------
  ! Deallocate flux arrays
  


  !---------------------------------------------------------------------
  ! Calculate surface downwelling fluxes in each band using the
  ! downwelling surface fluxes at each g point
  


  !---------------------------------------------------------------------
  ! Calculate top-of-atmosphere fluxes in each band using the fluxes
  ! at each g point
  


  !---------------------------------------------------------------------
  ! Return .true. if the most important flux variables are out of a
  ! physically sensible range, optionally only considering columns
  ! between istartcol and iendcol
  

  !---------------------------------------------------------------------
  ! Return .true. if the heating rates are out of a physically
  ! sensible range, optionally only considering columns between
  ! istartcol and iendcol. This function allocates and deallocates
  ! memory due to the requirements for inputs of out_of_bounds_2d.
  


  !---------------------------------------------------------------------
  ! Sum elements of "source" into "dest" according to index "ind".
  ! "source" and "ind" should have the same size and bounds, and no
  ! element of "ind" should refer outside the bounds of "dest".  This
  ! version increments existing contents of "dest".
  


  !---------------------------------------------------------------------
  ! As "add_indexed_sum" but this version overwrites existing contents
  ! of "dest"
  

  !---------------------------------------------------------------------
  ! Vectorized version of "add_indexed_sum"
  

  !---------------------------------------------------------------------
  ! As "add_indexed_sum" but a whole vertical profiles
  


  !---------------------------------------------------------------------
  ! As "indexed_sum" but a whole vertical profiles
  

  !---------------------------------------------------------------------
  ! Creates fields on device
  
  !---------------------------------------------------------------------
  ! updates fields on host
  

  !---------------------------------------------------------------------
  ! updates fields on device
  

  !---------------------------------------------------------------------
  ! Deletes fields on device
  

  subroutine allocate_flux_type_GPU(this, config, istartcol, iendcol, nlev, lacc)
use yomhook
use radiation_io
use radiation_config
integer, intent(in)             :: istartcol, iendcol, nlev
class(flux_type), intent(inout) :: this
type(config_type), intent(in)   :: config


logical, intent (in) :: lacc












end subroutine allocate_flux_type_GPU

  subroutine deallocate_flux_type_GPU(this, lacc)
use yomhook
class(flux_type), intent(inout) :: this

logical, intent (in) :: lacc































end subroutine deallocate_flux_type_GPU

  subroutine calc_surface_spectral_GPU(this, config, istartcol, iendcol, lacc)
use yomhook
use radiation_io
use radiation_config
class(flux_type),  intent(inout) :: this
type(config_type), intent(in)    :: config
integer,           intent(in)    :: istartcol, iendcol





logical, intent (in) :: lacc


 

 


end subroutine calc_surface_spectral_GPU

  subroutine calc_toa_spectral_GPU(this, config, istartcol, iendcol, lacc)
use yomhook
use radiation_config
use radiation_io
class(flux_type),  intent(inout) :: this
type(config_type), intent(in)    :: config
integer,           intent(in)    :: istartcol, iendcol


logical, intent (in) :: lacc





end subroutine calc_toa_spectral_GPU

  function out_of_physical_bounds_GPU(this, istartcol, iendcol, lacc) result(is_bad)
use yomhook
use radiation_check
class(flux_type), intent(inout) :: this
integer, optional,intent(in) :: istartcol, iendcol
logical                      :: is_bad

logical, intent (in) :: lacc



end function out_of_physical_bounds_GPU

  function heating_rate_out_of_physical_bounds_GPU(this, nlev, istartcol, iendcol, pressure_hl, lacc) result(is_bad)
use radiation_check
use radiation_constants


class(flux_type), intent(inout) :: this
integer, intent(in) :: istartcol, iendcol, nlev
logical                      :: is_bad
real(jprb), intent(in) :: pressure_hl(:,:)


logical, intent (in) :: lacc









end function heating_rate_out_of_physical_bounds_GPU

  pure subroutine add_indexed_sum_GPU(source, ind, dest, lacc)
real(jprb), intent(in)    :: source(:)
integer,    intent(in)    :: ind(:)
real(jprb), intent(inout) :: dest(:)

logical, intent (in) :: lacc



end subroutine add_indexed_sum_GPU

  pure subroutine indexed_sum_GPU(source, ind, dest)
real(jprb), intent(in)  :: source(:)
integer,    intent(in)  :: ind(:)
real(jprb), intent(out) :: dest(:)







end subroutine indexed_sum_GPU

  subroutine indexed_sum_vec_GPU(source, ind, dest, ist, iend, lacc)
real(jprb), intent(in)  :: source(:,:)
integer,    intent(in)  :: ind(:)
real(jprb), intent(out) :: dest(:,:)
integer,    intent(in)  :: ist, iend

logical, intent (in) :: lacc


end subroutine indexed_sum_vec_GPU

  pure subroutine add_indexed_sum_profile_GPU(source, ind, dest, lacc)
real(jprb), intent(in)  :: source(:,:)
integer,    intent(in)  :: ind(:)
real(jprb), intent(out) :: dest(:,:)

logical, intent (in) :: lacc




end subroutine add_indexed_sum_profile_GPU

  pure subroutine indexed_sum_profile_GPU(source, ind, dest, lacc)
real(jprb), intent(in)  :: source(:,:)
integer,    intent(in)  :: ind(:)
real(jprb), intent(out) :: dest(:,:)

logical, intent (in) :: lacc





end subroutine indexed_sum_profile_GPU

  subroutine create_device_GPU(this, lacc)
class(flux_type), intent(inout) :: this
logical, intent (in) :: lacc














































end subroutine create_device_GPU

  subroutine update_host_GPU(this, lacc)
class(flux_type), intent(inout) :: this
logical, intent (in) :: lacc














































end subroutine update_host_GPU

  subroutine update_device_GPU(this, lacc)
class(flux_type), intent(inout) :: this
logical, intent (in) :: lacc














































end subroutine update_device_GPU

  subroutine delete_device_GPU(this, lacc)
class(flux_type), intent(inout) :: this
logical, intent (in) :: lacc














































end subroutine delete_device_GPU

  subroutine allocate_flux_type_CPU(this, config, istartcol, iendcol, nlev)
use yomhook
use radiation_io
use radiation_config
integer, intent(in)             :: istartcol, iendcol, nlev
class(flux_type), intent(inout) :: this
type(config_type), intent(in)   :: config














end subroutine allocate_flux_type_CPU

  subroutine deallocate_flux_type_CPU(this)
use yomhook
class(flux_type), intent(inout) :: this
































end subroutine deallocate_flux_type_CPU

  subroutine calc_surface_spectral_CPU(this, config, istartcol, iendcol)
use yomhook
use radiation_config
class(flux_type),  intent(inout) :: this
type(config_type), intent(in)    :: config
integer,           intent(in)    :: istartcol, iendcol






 

 


end subroutine calc_surface_spectral_CPU

  subroutine calc_toa_spectral_CPU(this, config, istartcol, iendcol)
use yomhook
use radiation_config
class(flux_type),  intent(inout) :: this
type(config_type), intent(in)    :: config
integer,           intent(in)    :: istartcol, iendcol






end subroutine calc_toa_spectral_CPU

  function out_of_physical_bounds_CPU(this, istartcol, iendcol) result(is_bad)
use yomhook
use radiation_check
class(flux_type), intent(inout) :: this
integer, optional,intent(in) :: istartcol, iendcol
logical                      :: is_bad




end function out_of_physical_bounds_CPU

  function heating_rate_out_of_physical_bounds_CPU(this, nlev, istartcol, iendcol, pressure_hl) result(is_bad)
use radiation_check
use radiation_constants


class(flux_type), intent(inout) :: this
integer, intent(in) :: istartcol, iendcol, nlev
logical                      :: is_bad
real(jprb), intent(in) :: pressure_hl(:,:)











end function heating_rate_out_of_physical_bounds_CPU

  pure subroutine add_indexed_sum_CPU(source, ind, dest)
real(jprb), intent(in)    :: source(:)
integer,    intent(in)    :: ind(:)
real(jprb), intent(inout) :: dest(:)




end subroutine add_indexed_sum_CPU

  pure subroutine indexed_sum_CPU(source, ind, dest)
real(jprb), intent(in)  :: source(:)
integer,    intent(in)  :: ind(:)
real(jprb), intent(out) :: dest(:)





end subroutine indexed_sum_CPU

  subroutine indexed_sum_vec_CPU(source, ind, dest, ist, iend)
real(jprb), intent(in)  :: source(:,:)
integer,    intent(in)  :: ind(:)
real(jprb), intent(out) :: dest(:,:)
integer,    intent(in)  :: ist, iend



end subroutine indexed_sum_vec_CPU

  pure subroutine add_indexed_sum_profile_CPU(source, ind, dest)
real(jprb), intent(in)  :: source(:,:)
integer,    intent(in)  :: ind(:)
real(jprb), intent(out) :: dest(:,:)





end subroutine add_indexed_sum_profile_CPU

  pure subroutine indexed_sum_profile_CPU(source, ind, dest)
real(jprb), intent(in)  :: source(:,:)
integer,    intent(in)  :: ind(:)
real(jprb), intent(out) :: dest(:,:)






end subroutine indexed_sum_profile_CPU

end module radiation_flux

