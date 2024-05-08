! radiation_gas.F90 - Derived type to store the gas mixing ratios
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
!   2019-01-14  R. Hogan  Added out_of_physical_bounds routine

module radiation_gas

  use parkind1
  use radiation_gas_constants

  implicit none
  

  ! Available units
  enum, bind(c)
    enumerator IMassMixingRatio, IVolumeMixingRatio
  end enum

  !---------------------------------------------------------------------
  ! This derived type describes the gaseous composition of the
  ! atmosphere; gases may be stored as part of a 3D array (if their
  ! variation with height/column is to be represented) or one 1D array
  ! (if they are to be assumed globally well-mixed).
  type gas_type
    ! Units of each stored gas (or 0 if not present)
    integer :: iunits(NMaxGases) = 0

    ! Scaling factor that should be applied to each stored gas to get
    ! a dimensionless result, e.g. if iunits=IVolumeMixingRatio then
    ! 1.0e-6 is used to indicate the units are actually PPMV: need to
    ! multiply by 1e-6 to get mol/mol.
    real(jprb) :: scale_factor(NMaxGases) = 1.0_jprb

    ! Mixing ratios of variable gases, dimensioned (ncol, nlev,
    ! NMaxGases)
    real(jprb), allocatable, dimension(:,:,:) :: mixing_ratio

    ! Flag to indicate whether a gas is present
    logical :: is_present(NMaxGases) = .false.

    ! Flag to indicate whether a gas is well mixed
    logical :: is_well_mixed(NMaxGases) = .false.

    integer :: ntype          = 0 ! Number of gas types described

    integer :: ncol           = 0 ! Number of columns in mixing_ratio
    integer :: nlev           = 0 ! Number of levels  in mixing_ratio

    ! A list of length ntype of gases whose volume mixing ratios have
    ! been provided
    integer :: icode(NMaxGases) = 0

   contains
     
     
     
     
     
     
     
     
     
     
     
    
    
    
    

  procedure :: allocate_GPU   => allocate_gas_GPU

  procedure :: deallocate_GPU => deallocate_gas_GPU

  procedure :: put_GPU        => put_gas_GPU

  procedure :: put_well_mixed_GPU => put_well_mixed_gas_GPU

  procedure :: scale_GPU      => scale_gas_GPU

  procedure :: set_units_GPU  => set_units_gas_GPU

  procedure :: assert_units_GPU => assert_units_gas_GPU

  procedure :: get_GPU        => get_gas_GPU

  procedure :: get_scaling_GPU

  procedure :: reverse_GPU    => reverse_gas_GPU

  procedure :: out_of_physical_bounds_GPU

  procedure :: create_device_GPU

  procedure :: update_host_GPU

  procedure :: update_device_GPU

  procedure :: delete_device_GPU

  procedure :: allocate_CPU   => allocate_gas_CPU

  procedure :: deallocate_CPU => deallocate_gas_CPU

  procedure :: put_CPU        => put_gas_CPU

  procedure :: put_well_mixed_CPU => put_well_mixed_gas_CPU

  procedure :: scale_CPU      => scale_gas_CPU

  procedure :: set_units_CPU  => set_units_gas_CPU

  procedure :: assert_units_CPU => assert_units_gas_CPU

  procedure :: get_CPU        => get_gas_CPU

  procedure :: get_scaling_CPU

  procedure :: reverse_CPU    => reverse_gas_CPU

  procedure :: out_of_physical_bounds_CPU

  end type gas_type

  private :: create_device_GPU, update_host_GPU, update_device_GPU, delete_device_GPU

contains


  !---------------------------------------------------------------------
  ! Allocate a derived type for holding gas mixing ratios given the
  ! number of columns and levels
  


  !---------------------------------------------------------------------
  ! Deallocate memory and reset arrays
  


  !---------------------------------------------------------------------
  ! Put gas mixing ratio corresponding to gas ID "igas" with units
  ! "iunits"
  


  !---------------------------------------------------------------------
  ! Put well-mixed gas mixing ratio corresponding to gas ID "igas"
  ! with units "iunits"
  


  !---------------------------------------------------------------------
  ! Scale gas concentrations, e.g. igas=ICO2 and set scale_factor=2 to
  ! double CO2.  Note that this does not perform the scaling
  ! immediately, but changes the scale factor for the specified gas,
  ! ready to be used in set_units_gas.
  


  !---------------------------------------------------------------------
  ! Scale the gas concentrations so that they have the units "iunits"
  ! and are therefore ready to be used by the gas optics model within
  ! ecRad with no further scaling.  The existing scale_factor for each
  ! gas is applied.  If "igas" is present then apply only to gas with
  ! ID "igas", otherwise to all gases. Optional argument scale_factor
  ! specifies scaling that any subsequent access would need to apply
  ! to get a dimensionless result (consistent with definition of
  ! gas_type). So say that your gas optics model requires gas
  ! concentrations in PPMV, specify iunits=IVolumeMixingRatio and
  ! scale_factor=1.0e-6. If the gas concentrations were currently
  ! dimensionless volume mixing ratios, then the values would be
  ! internally divided by 1.0e-6.
  


  !---------------------------------------------------------------------
  ! Return a vector indicating the scaling that one would need to
  ! apply to each gas in order to obtain the dimension units in
  ! "iunits" (which can be IVolumeMixingRatio or IMassMixingRatio)
  


  !---------------------------------------------------------------------
  ! Assert that gas mixing ratio units are "iunits", applying to gas
  ! with ID "igas" if present, otherwise to all gases. Otherwise the
  ! program will exit, except if the optional argument "istatus" is
  ! provided in which case it will return true if the units are
  ! correct and false if they are not. Optional argument scale factor
  ! specifies any subsequent multiplication to apply; for PPMV one
  ! would use iunits=IVolumeMixingRatio and scale_factor=1.0e6.
  


  !---------------------------------------------------------------------
  ! Get gas mixing ratio corresponding to gas ID "igas" with units
  ! "iunits" and return as a 2D array of dimensions (ncol,nlev).  The
  ! array will contain zeros if the gas is not stored.
  


  !---------------------------------------------------------------------
  ! Copy data to "gas_rev", reversing the height ordering of the gas
  ! data
  

  !---------------------------------------------------------------------
  ! Return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  

  !---------------------------------------------------------------------
  ! creates fields on device
  

  !---------------------------------------------------------------------
  ! updates fields on host
  

  !---------------------------------------------------------------------
  ! updates fields on device
  

  !---------------------------------------------------------------------
  ! deletes fields on device
  

  subroutine allocate_gas_GPU(this, ncol, nlev, lacc)
use yomhook
use radiation_io
class(gas_type), intent(inout) :: this
integer,         intent(in)    :: ncol, nlev


logical, intent (in) :: lacc







end subroutine allocate_gas_GPU

  subroutine deallocate_gas_GPU(this, lacc)
use yomhook
class(gas_type), intent(inout) :: this

logical, intent (in) :: lacc











end subroutine deallocate_gas_GPU

  subroutine put_gas_GPU(this, igas, iunits, mixing_ratio, scale_factor, &
istartcol, lacc)
use yomhook
use radiation_io
class(gas_type),      intent(inout) :: this
integer,              intent(in)    :: igas
integer,              intent(in)    :: iunits
real(jprb),           intent(in)    :: mixing_ratio(:,:)
real(jprb), optional, intent(in)    :: scale_factor
integer,    optional, intent(in)    :: istartcol
logical,    optional, intent(in)    :: lacc

























end subroutine put_gas_GPU

  subroutine put_well_mixed_gas_GPU(this, igas, iunits, mixing_ratio, &
scale_factor, istartcol, iendcol, lacc)
use yomhook
use radiation_io
class(gas_type),      intent(inout) :: this
integer,              intent(in)    :: igas
integer,              intent(in)    :: iunits
real(jprb),           intent(in)    :: mixing_ratio
real(jprb), optional, intent(in)    :: scale_factor
integer,    optional, intent(in)    :: istartcol, iendcol
logical,    optional, intent(in)    :: lacc

























end subroutine put_well_mixed_gas_GPU

  subroutine scale_gas_GPU(this, igas, scale_factor, lverbose, lacc)
use radiation_io
class(gas_type),      intent(inout) :: this
integer,              intent(in)    :: igas
real(jprb),           intent(in)    :: scale_factor
logical,    optional, intent(in)    :: lverbose
logical, intent (in) :: lacc

end subroutine scale_gas_GPU

  recursive subroutine set_units_gas_GPU(this, iunits, igas, scale_factor, lacc)
class(gas_type),      intent(inout) :: this
integer,              intent(in)    :: iunits
integer,    optional, intent(in)    :: igas
real(jprb), optional, intent(in)    :: scale_factor
logical,    optional, intent(in)    :: lacc









end subroutine set_units_gas_GPU

  subroutine get_scaling_GPU(this, iunits, scaling, lacc)
class(gas_type), intent(in)  :: this
integer,         intent(in)  :: iunits
real(jprb),      intent(out) :: scaling(NMaxGases)

logical, intent (in) :: lacc


end subroutine get_scaling_GPU

  recursive subroutine assert_units_gas_GPU(this, iunits, igas, scale_factor, istatus, lacc)
use radiation_io
class(gas_type),      intent(in)  :: this
integer,              intent(in)  :: iunits
integer,    optional, intent(in)  :: igas
real(jprb), optional, intent(in)  :: scale_factor
logical,    optional, intent(out) :: istatus


logical, intent (in) :: lacc



end subroutine assert_units_gas_GPU

  subroutine get_gas_GPU(this, igas, iunits, mixing_ratio, scale_factor, &
&   istartcol, lacc)
use yomhook
use radiation_io
class(gas_type),      intent(in)  :: this
integer,              intent(in)  :: igas
integer,              intent(in)  :: iunits
real(jprb),           intent(out) :: mixing_ratio(:,:)
real(jprb), optional, intent(in)  :: scale_factor
integer,    optional, intent(in)  :: istartcol



logical, intent (in) :: lacc







end subroutine get_gas_GPU

  subroutine reverse_gas_GPU(this, istartcol, iendcol, gas_rev, lacc)
class(gas_type), intent(in) :: this
integer,        intent(in)  :: istartcol, iendcol
type(gas_type), intent(out) :: gas_rev
logical, intent (in) :: lacc










end subroutine reverse_gas_GPU

  function out_of_physical_bounds_GPU(this, istartcol, iendcol, do_fix, lacc) result(is_bad)
use yomhook
use radiation_check
class(gas_type),   intent(inout) :: this
integer,  optional,intent(in) :: istartcol, iendcol
logical,  optional,intent(in) :: do_fix
logical                       :: is_bad


logical, intent (in) :: lacc




end function out_of_physical_bounds_GPU

  subroutine create_device_GPU(this, lacc)
class(gas_type), intent(inout) :: this
logical, intent (in) :: lacc





end subroutine create_device_GPU

  subroutine update_host_GPU(this, lacc)
class(gas_type), intent(inout) :: this
logical, intent (in) :: lacc


end subroutine update_host_GPU

  subroutine update_device_GPU(this, lacc)
class(gas_type), intent(inout) :: this
logical, intent (in) :: lacc


end subroutine update_device_GPU

  subroutine delete_device_GPU(this, lacc)
class(gas_type), intent(inout) :: this
logical, intent (in) :: lacc


end subroutine delete_device_GPU

  subroutine allocate_gas_CPU(this, ncol, nlev)
use yomhook
use radiation_io
class(gas_type), intent(inout) :: this
integer,         intent(in)    :: ncol, nlev










end subroutine allocate_gas_CPU

  subroutine deallocate_gas_CPU(this)
use yomhook
class(gas_type), intent(inout) :: this












end subroutine deallocate_gas_CPU

  subroutine put_gas_CPU(this, igas, iunits, mixing_ratio, scale_factor, &
istartcol, lacc)
use yomhook
use radiation_io
class(gas_type),      intent(inout) :: this
integer,              intent(in)    :: igas
integer,              intent(in)    :: iunits
real(jprb),           intent(in)    :: mixing_ratio(:,:)
real(jprb), optional, intent(in)    :: scale_factor
integer,    optional, intent(in)    :: istartcol
logical,    optional, intent(in)    :: lacc




















end subroutine put_gas_CPU

  subroutine put_well_mixed_gas_CPU(this, igas, iunits, mixing_ratio, &
scale_factor, istartcol, iendcol, lacc)
use yomhook
use radiation_io
class(gas_type),      intent(inout) :: this
integer,              intent(in)    :: igas
integer,              intent(in)    :: iunits
real(jprb),           intent(in)    :: mixing_ratio
real(jprb), optional, intent(in)    :: scale_factor
integer,    optional, intent(in)    :: istartcol, iendcol
logical,    optional, intent(in)    :: lacc




















end subroutine put_well_mixed_gas_CPU

  subroutine scale_gas_CPU(this, igas, scale_factor, lverbose)
use radiation_io
class(gas_type),      intent(inout) :: this
integer,              intent(in)    :: igas
real(jprb),           intent(in)    :: scale_factor
logical,    optional, intent(in)    :: lverbose

end subroutine scale_gas_CPU

  recursive subroutine set_units_gas_CPU(this, iunits, igas, scale_factor, lacc)
class(gas_type),      intent(inout) :: this
integer,              intent(in)    :: iunits
integer,    optional, intent(in)    :: igas
real(jprb), optional, intent(in)    :: scale_factor
logical,    optional, intent(in)    :: lacc









end subroutine set_units_gas_CPU

  subroutine get_scaling_CPU(this, iunits, scaling)
class(gas_type), intent(in)  :: this
integer,         intent(in)  :: iunits
real(jprb),      intent(out) :: scaling(NMaxGases)



end subroutine get_scaling_CPU

  recursive subroutine assert_units_gas_CPU(this, iunits, igas, scale_factor, istatus)
use radiation_io
class(gas_type),      intent(in)  :: this
integer,              intent(in)  :: iunits
integer,    optional, intent(in)  :: igas
real(jprb), optional, intent(in)  :: scale_factor
logical,    optional, intent(out) :: istatus





end subroutine assert_units_gas_CPU

  subroutine get_gas_CPU(this, igas, iunits, mixing_ratio, scale_factor, &
&   istartcol)
use yomhook
use radiation_io
class(gas_type),      intent(in)  :: this
integer,              intent(in)  :: igas
integer,              intent(in)  :: iunits
real(jprb),           intent(out) :: mixing_ratio(:,:)
real(jprb), optional, intent(in)  :: scale_factor
integer,    optional, intent(in)  :: istartcol













end subroutine get_gas_CPU

  subroutine reverse_gas_CPU(this, istartcol, iendcol, gas_rev)
class(gas_type), intent(in) :: this
integer,        intent(in)  :: istartcol, iendcol
type(gas_type), intent(out) :: gas_rev










end subroutine reverse_gas_CPU

  function out_of_physical_bounds_CPU(this, istartcol, iendcol, do_fix) result(is_bad)
use yomhook
use radiation_check
class(gas_type),   intent(inout) :: this
integer,  optional,intent(in) :: istartcol, iendcol
logical,  optional,intent(in) :: do_fix
logical                       :: is_bad






end function out_of_physical_bounds_CPU

end module radiation_gas

