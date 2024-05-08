! radiation_aerosol.F90 - Derived type describing aerosol
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
!   2018-04-15  R. Hogan  Add "direct" option
!   2019-01-14  R. Hogan  Added out_of_physical_bounds routine

module radiation_aerosol

  use parkind1
  use radiation_io

  implicit none
  

  !---------------------------------------------------------------------
  ! Type describing the aerosol content in the atmosphere
  type aerosol_type
     ! The mass mixing ratio of config%n_aerosol_types different
     ! aerosol types dimensioned
     ! (ncol,istartlev:iendlev,config%n_aerosol_types), where ncol is
     ! the number of columns, istartlev:iendlev is the range of model
     ! levels where aerosols are present
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  mixing_ratio  ! mass mixing ratio (kg/kg)

     ! Alternatively, if is_direct=true, the optical properties are
     ! provided directly and are dimensioned
     ! (nband,istartlev:iendlev,ncol)
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  od_sw, ssa_sw, g_sw, & ! Shortwave optical properties
          &  od_lw, ssa_lw, g_lw    ! Longwave optical properties

     ! Range of levels in which the aerosol properties are provided
     integer :: istartlev, iendlev

     ! Are the optical properties going to be provided directly by the
     ! user?
     logical :: is_direct = .false.

   contains
      
      
      
      
      
      
      
      
  procedure :: allocate_GPU        => allocate_aerosol_arrays_GPU

  procedure :: allocate_direct_GPU => allocate_aerosol_arrays_direct_GPU

  procedure :: deallocate_GPU      => deallocate_aerosol_arrays_GPU

  procedure :: out_of_physical_bounds_GPU

  procedure :: create_device_GPU

  procedure :: update_host_GPU

  procedure :: update_device_GPU

  procedure :: delete_device_GPU

  procedure :: allocate_CPU        => allocate_aerosol_arrays_CPU

  procedure :: allocate_direct_CPU => allocate_aerosol_arrays_direct_CPU

  procedure :: deallocate_CPU      => deallocate_aerosol_arrays_CPU

  procedure :: out_of_physical_bounds_CPU

  end type aerosol_type

  private :: create_device_GPU, update_host_GPU, update_device_GPU, delete_device_GPU

contains

  !---------------------------------------------------------------------
  ! Allocate array for describing aerosols, although in the offline
  ! code these are allocated when they are read from the NetCDF file
  


  !---------------------------------------------------------------------
  ! Allocate arrays for describing aerosol optical properties
  


  !---------------------------------------------------------------------
  ! Deallocate arrays
  


  !---------------------------------------------------------------------
  ! Return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  

  !---------------------------------------------------------------------
  ! Creates fields on device
  

  !---------------------------------------------------------------------
  ! updates fields on host
  

  !---------------------------------------------------------------------
  ! updates fields on device
  

  !---------------------------------------------------------------------
  ! Deletes fields on device
  

  subroutine allocate_aerosol_arrays_GPU(this, ncol, istartlev, iendlev, ntype, lacc)
use yomhook
class(aerosol_type), intent(inout) :: this
integer, intent(in)                :: ncol  
integer, intent(in)                :: istartlev, iendlev 
integer, intent(in)                :: ntype 

logical, intent (in) :: lacc






end subroutine allocate_aerosol_arrays_GPU

  subroutine allocate_aerosol_arrays_direct_GPU(this, config, &
&                                    ncol, istartlev, iendlev, lacc)
use yomhook
use radiation_config
class(aerosol_type), intent(inout) :: this
type(config_type),   intent(in)    :: config
integer, intent(in)                :: ncol  
integer, intent(in)                :: istartlev, iendlev 


logical, intent (in) :: lacc







end subroutine allocate_aerosol_arrays_direct_GPU

  subroutine deallocate_aerosol_arrays_GPU(this, lacc)
use yomhook
class(aerosol_type), intent(inout) :: this

logical, intent (in) :: lacc









end subroutine deallocate_aerosol_arrays_GPU

  function out_of_physical_bounds_GPU(this, istartcol, iendcol, do_fix, lacc) result(is_bad)
use yomhook
use radiation_check
class(aerosol_type),   intent(inout) :: this
integer,      optional,intent(in) :: istartcol, iendcol
logical,      optional,intent(in) :: do_fix
logical                           :: is_bad


logical, intent (in) :: lacc




end function out_of_physical_bounds_GPU

  subroutine create_device_GPU(this, lacc)
class(aerosol_type), intent(inout) :: this
logical, intent (in) :: lacc







end subroutine create_device_GPU

  subroutine update_host_GPU(this, lacc)
class(aerosol_type), intent(inout) :: this
logical, intent (in) :: lacc







end subroutine update_host_GPU

  subroutine update_device_GPU(this, lacc)
class(aerosol_type), intent(inout) :: this
logical, intent (in) :: lacc







end subroutine update_device_GPU

  subroutine delete_device_GPU(this, lacc)
class(aerosol_type), intent(inout) :: this
logical, intent (in) :: lacc







end subroutine delete_device_GPU

  subroutine allocate_aerosol_arrays_CPU(this, ncol, istartlev, iendlev, ntype)
use yomhook
class(aerosol_type), intent(inout) :: this
integer, intent(in)                :: ncol  
integer, intent(in)                :: istartlev, iendlev 
integer, intent(in)                :: ntype 







end subroutine allocate_aerosol_arrays_CPU

  subroutine allocate_aerosol_arrays_direct_CPU(this, config, &
&                                    ncol, istartlev, iendlev)
use yomhook
use radiation_config
class(aerosol_type), intent(inout) :: this
type(config_type),   intent(in)    :: config
integer, intent(in)                :: ncol  
integer, intent(in)                :: istartlev, iendlev 









end subroutine allocate_aerosol_arrays_direct_CPU

  subroutine deallocate_aerosol_arrays_CPU(this)
use yomhook
class(aerosol_type), intent(inout) :: this










end subroutine deallocate_aerosol_arrays_CPU

  function out_of_physical_bounds_CPU(this, istartcol, iendcol, do_fix) result(is_bad)
use yomhook
use radiation_check
class(aerosol_type),   intent(inout) :: this
integer,      optional,intent(in) :: istartcol, iendcol
logical,      optional,intent(in) :: do_fix
logical                           :: is_bad






end function out_of_physical_bounds_CPU

end module radiation_aerosol

