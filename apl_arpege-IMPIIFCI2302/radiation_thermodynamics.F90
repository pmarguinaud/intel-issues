! radiation_thermodynamics.F90 - Derived type for pressure & temperature
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
!   2017-05-11  R. Hogan  Fix startcol/endcol for get_layer_mass
!   2019-01-14  R. Hogan  Added out_of_physical_bounds routine
!   2019-01-14  R. Hogan  Capped h2o_sat_liq at 1

module radiation_thermodynamics

  use parkind1

  implicit none
  

  !---------------------------------------------------------------------
  ! Derived type for storing pressure and temperature at half and full levels
  type thermodynamics_type
     real(jprb), allocatable, dimension(:,:) :: &
          &  pressure_hl,    & ! (ncol,nlev+1) pressure (Pa)
          &  temperature_hl, & ! (ncol,nlev+1) temperature (K)
          &  pressure_fl,    & ! (ncol,nlev) pressure (Pa)
          &  temperature_fl    ! (ncol,nlev) temperature (K)

     ! The following is a function of pressure and temperature: you
     ! can calculate it according to your favourite formula, or the
     ! calc_saturation_wrt_liquid subroutine can be used to do this
     ! approximately
     real(jprb), allocatable, dimension(:,:) :: &
          &  h2o_sat_liq ! (ncol,nlev) specific humidity at liquid
                         ! saturation (kg/kg)

     ! Using the interpolation method for temperature and pressure from half levels
     ! to full levels that is used in the subroutine gas_optics in radiation_ifs_rrtm
     ! can result in values for ind1 that exceed the bounds of absa within
     ! srtm_taumol16 for ecRad inside ICON. This can be avoided by directly
     ! passing pressure_fl and temperature_fl from ICON to ecRad. With
     ! rrtm_pass_temppres_fl = .TRUE., the fields pressure_fl and temperature_fl
     ! are allocated and used within gas_optics in radiation_ifs_rrtm.
     logical :: &
          &  rrtm_pass_temppres_fl
   contains
     
     
     
     
     
     
    
    
    
    

  procedure :: allocate_GPU   => allocate_thermodynamics_arrays_GPU

  procedure :: deallocate_GPU => deallocate_thermodynamics_arrays_GPU

  procedure :: calc_saturation_wrt_liquid_GPU

  procedure :: get_layer_mass_GPU

  procedure :: get_layer_mass_column_GPU

  procedure :: out_of_physical_bounds_GPU

  procedure :: create_device_GPU

  procedure :: update_host_GPU

  procedure :: update_device_GPU

  procedure :: delete_device_GPU

  procedure :: allocate_CPU   => allocate_thermodynamics_arrays_CPU

  procedure :: deallocate_CPU => deallocate_thermodynamics_arrays_CPU

  procedure :: calc_saturation_wrt_liquid_CPU

  procedure :: get_layer_mass_CPU

  procedure :: get_layer_mass_column_CPU

  procedure :: out_of_physical_bounds_CPU

  end type thermodynamics_type

  private :: create_device_GPU, update_host_GPU, update_device_GPU, delete_device_GPU

contains


  !---------------------------------------------------------------------
  ! Allocate variables with specified dimensions
  


  !---------------------------------------------------------------------
  ! Deallocate variables
  


  !---------------------------------------------------------------------
  ! Calculate approximate saturation with respect to liquid
  


  !---------------------------------------------------------------------
  ! Calculate the dry mass of each layer, neglecting humidity effects.
  ! The first version is for all columns.
  

  !---------------------------------------------------------------------
  ! Calculate the dry mass of each layer, neglecting humidity effects.
  ! The second version is for one column, the one numbered "icol".
  


  !---------------------------------------------------------------------
  ! Estimate the separation between the mid-points of model layers
  ! given the half-level pressure and temperature.  This is not in
  ! terms of the "thermodynamics" type as it is useful for computing
  ! overlap decorrelation lengths and hence cloud cover outside the
  ! radiation scheme.
  


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
  

  subroutine allocate_thermodynamics_arrays_GPU(this, ncol, nlev, &
&                                    use_h2o_sat, rrtm_pass_temppres_fl, lacc)
use yomhook
class(thermodynamics_type), intent(inout) :: this
integer, intent(in)           :: ncol  
integer, intent(in)           :: nlev  
logical, intent(in), optional :: use_h2o_sat 
logical, intent(in), optional :: rrtm_pass_temppres_fl 



logical, intent (in) :: lacc










end subroutine allocate_thermodynamics_arrays_GPU

  subroutine deallocate_thermodynamics_arrays_GPU(this, lacc)
use yomhook
class(thermodynamics_type), intent(inout) :: this

logical, intent (in) :: lacc







end subroutine deallocate_thermodynamics_arrays_GPU

  subroutine calc_saturation_wrt_liquid_GPU(this,istartcol,iendcol, lacc)
use yomhook
class(thermodynamics_type), intent(inout) :: this
integer, intent(in)                       :: istartcol, iendcol
logical, optional, intent(in)             :: lacc





 
 











end subroutine calc_saturation_wrt_liquid_GPU

  subroutine get_layer_mass_GPU(this,istartcol,iendcol,layer_mass, lacc)
use yomhook
use radiation_constants
class(thermodynamics_type), intent(in)  :: this
integer,                    intent(in)  :: istartcol, iendcol
real(jprb),                 intent(out) :: layer_mass(:,:)



logical, intent (in) :: lacc








end subroutine get_layer_mass_GPU

  subroutine get_layer_mass_column_GPU(this, icol, layer_mass, lacc)
use yomhook
use radiation_constants
class(thermodynamics_type), intent(in)  :: this
integer,                    intent(in)  :: icol
real(jprb),                 intent(out) :: layer_mass(:)



logical, intent (in) :: lacc





end subroutine get_layer_mass_column_GPU

  subroutine get_layer_separation_GPU(pressure_hl, temperature_hl, layer_separation, lacc)
use yomhook
use radiation_constants



real(jprb), dimension(:,:), intent(in)  :: pressure_hl, temperature_hl

real(jprb), dimension(:,:), intent(out) :: layer_separation






logical, intent (in) :: lacc






end subroutine get_layer_separation_GPU

  function out_of_physical_bounds_GPU(this, istartcol, iendcol, do_fix, lacc) result(is_bad)
use yomhook
use radiation_check
class(thermodynamics_type), intent(inout) :: this
integer,           optional,intent(in) :: istartcol, iendcol
logical,           optional,intent(in) :: do_fix
logical                                :: is_bad


logical, intent (in) :: lacc





end function out_of_physical_bounds_GPU

  subroutine create_device_GPU(this, lacc)
class(thermodynamics_type), intent(inout) :: this
logical, intent (in) :: lacc










end subroutine create_device_GPU

  subroutine update_host_GPU(this, lacc)
class(thermodynamics_type), intent(inout) :: this
logical, intent (in) :: lacc










end subroutine update_host_GPU

  subroutine update_device_GPU(this, lacc)
class(thermodynamics_type), intent(inout) :: this
logical, intent (in) :: lacc










end subroutine update_device_GPU

  subroutine delete_device_GPU(this, lacc)
class(thermodynamics_type), intent(inout) :: this
logical, intent (in) :: lacc










end subroutine delete_device_GPU

  subroutine allocate_thermodynamics_arrays_CPU(this, ncol, nlev, &
&                                    use_h2o_sat, rrtm_pass_temppres_fl)
use yomhook
class(thermodynamics_type), intent(inout) :: this
integer, intent(in)           :: ncol  
integer, intent(in)           :: nlev  
logical, intent(in), optional :: use_h2o_sat 
logical, intent(in), optional :: rrtm_pass_temppres_fl 













end subroutine allocate_thermodynamics_arrays_CPU

  subroutine deallocate_thermodynamics_arrays_CPU(this)
use yomhook
class(thermodynamics_type), intent(inout) :: this








end subroutine deallocate_thermodynamics_arrays_CPU

  subroutine calc_saturation_wrt_liquid_CPU(this,istartcol,iendcol, lacc)
use yomhook
class(thermodynamics_type), intent(inout) :: this
integer, intent(in)                       :: istartcol, iendcol
logical, optional, intent(in)             :: lacc





 
 








end subroutine calc_saturation_wrt_liquid_CPU

  subroutine get_layer_mass_CPU(this,istartcol,iendcol,layer_mass)
use yomhook
use radiation_constants
class(thermodynamics_type), intent(in)  :: this
integer,                    intent(in)  :: istartcol, iendcol
real(jprb),                 intent(out) :: layer_mass(:,:)








end subroutine get_layer_mass_CPU

  subroutine get_layer_mass_column_CPU(this, icol, layer_mass)
use yomhook
use radiation_constants
class(thermodynamics_type), intent(in)  :: this
integer,                    intent(in)  :: icol
real(jprb),                 intent(out) :: layer_mass(:)








end subroutine get_layer_mass_column_CPU

  subroutine get_layer_separation_CPU(pressure_hl, temperature_hl, layer_separation)
use yomhook
use radiation_constants



real(jprb), dimension(:,:), intent(in)  :: pressure_hl, temperature_hl

real(jprb), dimension(:,:), intent(out) :: layer_separation












end subroutine get_layer_separation_CPU

  function out_of_physical_bounds_CPU(this, istartcol, iendcol, do_fix) result(is_bad)
use yomhook
use radiation_check
class(thermodynamics_type), intent(inout) :: this
integer,           optional,intent(in) :: istartcol, iendcol
logical,           optional,intent(in) :: do_fix
logical                                :: is_bad







end function out_of_physical_bounds_CPU

end module radiation_thermodynamics

