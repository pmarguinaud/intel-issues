! radiation_cloud.F90 - Derived type to store cloud/precip properties
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
!   2019-01-14  R. Hogan  Added inv_inhom_effective_size variable
!   2019-01-14  R. Hogan  Added out_of_physical_bounds routine
!   2019-06-14  R. Hogan  Added capability to store any number of cloud/precip types

module radiation_cloud

  use parkind1

  implicit none
  

  !---------------------------------------------------------------------
  ! The intention is that all variables describing clouds and
  ! radiatively-active precipitation are contained in this derived
  ! type, and if cloud variables are to be added in future, they can
  ! be added to this type without requiring extra variables to be
  ! passed between subroutines elsewhere in the program.
  type cloud_type
    ! For maximum flexibility, an arbitrary number "ntype" of
    ! hydrometeor types can be stored, dimensioned (ncol,nlev,ntype)
    integer                                   :: ntype = 0
    real(jprb), allocatable, dimension(:,:,:) :: &
         &  mixing_ratio, &  ! mass mixing ratio (kg/kg)
         &  effective_radius ! (m)

    ! For backwards compatibility, we also allow for the two
    ! traditional cloud types, liquid cloud droplets and ice cloud
    ! particles, dimensioned (ncol,nlev)
    real(jprb), pointer, dimension(:,:) :: &
         &  q_liq,  q_ice,  & ! mass mixing ratio (kg/kg)
         &  re_liq, re_ice    ! effective radius (m)

    ! For the moment, the different types of hydrometeor are assumed
    ! to be mixed with each other, so there is just one cloud fraction
    ! variable varying from 0 to 1
    real(jprb), allocatable, dimension(:,:) :: fraction

    ! The fractional standard deviation of cloud optical depth in the
    ! cloudy part of the gridbox.  In the Tripleclouds representation
    ! of cloud inhomogeneity, this is implemented by splitting the
    ! cloudy part of the gridbox into two equal-area regions, one
    ! with the cloud optical depth scaled by 1+fractional_std and the
    ! other scaled by 1-fractional_std. This variable is dimensioned
    ! (ncol,nlev)
    real(jprb), allocatable, dimension(:,:) :: fractional_std

    ! The inverse of the effective horizontal size of the clouds in
    ! the gridbox, used to compute the cloud edge length per unit
    ! gridbox area for use in representing 3D effects. This variable
    ! is dimensioned (ncol,nlev).
    real(jprb), allocatable, dimension(:,:) :: inv_cloud_effective_size ! m-1

    ! Similarly for the in-cloud heterogeneities, used to compute the
    ! edge length between the optically thin and thick cloudy regions
    ! of the gridbox.
    real(jprb), allocatable, dimension(:,:) :: inv_inhom_effective_size ! m-1

    ! The following variable describes the overlap of cloud boundaries
    ! in adjacent layers, with dimensions (ncol,nlev-1): 1 corresponds
    ! to maximum overlap and 0 to random overlap. Depending on the
    ! ecRad configuration, it may be the "alpha" overlap parameter of
    ! Hogan and Illingworth (2000) or the "beta" overlap parameter of
    ! Shonk et al. (2010).
    real(jprb), allocatable, dimension(:,:) :: overlap_param

  contains
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

  procedure :: allocate_GPU   => allocate_cloud_arrays_GPU

  procedure :: deallocate_GPU => deallocate_cloud_arrays_GPU

  procedure :: set_overlap_param_fix_GPU

  procedure :: set_overlap_param_var_GPU

  procedure :: set_overlap_param_approx_GPU

  procedure :: create_fractional_std_GPU

  procedure :: create_inv_cloud_effective_size_GPU

  procedure :: create_inv_cloud_effective_size_eta_GPU

  procedure :: param_cloud_effective_separation_eta_GPU

  procedure :: crop_cloud_fraction_GPU

  procedure :: out_of_physical_bounds_GPU

  procedure :: create_device_GPU

  procedure :: update_host_GPU

  procedure :: update_device_GPU

  procedure :: delete_device_GPU

  procedure :: allocate_CPU   => allocate_cloud_arrays_CPU

  procedure :: deallocate_CPU => deallocate_cloud_arrays_CPU

  procedure :: set_overlap_param_fix_CPU

  procedure :: set_overlap_param_var_CPU

  procedure :: set_overlap_param_approx_CPU

  procedure :: create_fractional_std_CPU

  procedure :: create_inv_cloud_effective_size_CPU

  procedure :: create_inv_cloud_effective_size_eta_CPU

  procedure :: param_cloud_effective_separation_eta_CPU

  procedure :: crop_cloud_fraction_CPU

  procedure :: out_of_physical_bounds_CPU

  generic   :: set_overlap_param_GPU => set_overlap_param_fix_GPU, set_overlap_param_var_GPU

  generic   :: set_overlap_param_CPU => set_overlap_param_fix_CPU, set_overlap_param_var_CPU

  end type cloud_type

  private :: create_device_GPU, update_host_GPU, update_device_GPU, delete_device_GPU

contains

  !---------------------------------------------------------------------
  ! Allocate arrays for describing clouds and precipitation, although
  ! in the offline code these are allocated when they are read from
  ! the NetCDF file
  


  !---------------------------------------------------------------------
  ! Deallocate arrays
  


  !---------------------------------------------------------------------
  ! Compute and store the overlap parameter from the provided overlap
  ! decorrelation length (in metres).  If istartcol and/or iendcol are
  ! provided then only columns in this range are computed.  If the
  ! overlap_param array has not been allocated then it will be
  ! allocated to be of the correct size relative to the pressure
  ! field. This version assumes a fixed decorrelation_length for all
  ! columns.
  


  !---------------------------------------------------------------------
  ! Compute and store the overlap parameter from the provided overlap
  ! decorrelation length (in metres), which may vary with column. Only
  ! columns from istartcol to iendcol are computed.  If the
  ! overlap_param array has not been allocated then it will be
  ! allocated to be of the correct size relative to the pressure
  ! field.
  


  !---------------------------------------------------------------------
  ! Compute and store the overlap parameter from the provided overlap
  ! decorrelation length (in metres).  If istartcol and/or iendcol are
  ! provided then only columns in this range are computed.  If the
  ! overlap_param array has not been allocated then it will be
  ! allocated to be of the correct size relative to the pressure
  ! field. This is the APPROXIMATE method as it assumes a fixed
  ! atmospheric scale height, which leads to differences particularly
  ! in low cloud.
  


  !---------------------------------------------------------------------
  ! Create a matrix of constant fractional standard deviations
  ! (dimensionless)
  


  !---------------------------------------------------------------------
  ! Create a matrix of constant inverse cloud effective size (m-1)
  


  !---------------------------------------------------------------------
  ! Create a matrix of inverse cloud effective size (m-1) according to
  ! the value of eta (=pressure divided by surface pressure)
  


  !---------------------------------------------------------------------
  ! Create a matrix of inverse cloud and inhomogeneity effective size
  ! (m-1) parameterized according to the value of eta (=pressure
  ! divided by surface pressure): effective_separation =
  ! coeff_a + coeff_b*exp(-(eta**power)).
  


  !---------------------------------------------------------------------
  ! Remove "ghost" clouds: those with a cloud fraction that is too
  ! small to treat sensibly (e.g. because it implies that the
  ! "in-cloud" water content is too high), or with a cloud water
  ! content that is too small.  We do this in one place to ensure that
  ! all subsequent subroutines can assume that if cloud_fraction > 0.0
  ! then cloud is really present and should be treated.
  


  !---------------------------------------------------------------------
  ! Return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  

  
  !---------------------------------------------------------------------
  ! updates fields on host
  

  !---------------------------------------------------------------------
  ! updates fields on device
  

  

  subroutine allocate_cloud_arrays_GPU(this, ncol, nlev, ntype, use_inhom_effective_size, lacc)
use yomhook



class(cloud_type), intent(inout), target :: this
integer, intent(in)              :: ncol   
integer, intent(in)              :: nlev   




integer, intent(in), optional    :: ntype
logical, intent(in), optional    :: use_inhom_effective_size

logical, intent (in) :: lacc















end subroutine allocate_cloud_arrays_GPU

  subroutine deallocate_cloud_arrays_GPU(this, lacc)
use yomhook
class(cloud_type), intent(inout) :: this

logical, intent (in) :: lacc













end subroutine deallocate_cloud_arrays_GPU

  subroutine set_overlap_param_fix_GPU(this, thermodynamics, decorrelation_length, &
&  istartcol, iendcol, lacc)
use yomhook
use radiation_thermodynamics
use radiation_constants
class(cloud_type),         intent(inout) :: this
type(thermodynamics_type), intent(in)    :: thermodynamics
real(jprb),                intent(in)    :: decorrelation_length 
integer,         optional, intent(in)    :: istartcol, iendcol
logical, optional, intent(in) :: lacc























end subroutine set_overlap_param_fix_GPU

  subroutine set_overlap_param_var_GPU(this, thermodynamics, decorrelation_length, &
&                           istartcol, iendcol, lacc)
use yomhook
use radiation_thermodynamics
use radiation_constants
use radiation_io
class(cloud_type),         intent(inout) :: this
type(thermodynamics_type), intent(in)    :: thermodynamics
integer,                   intent(in)    :: istartcol, iendcol
real(jprb),                intent(in)    :: decorrelation_length(istartcol:iendcol) 
logical, optional, intent(in) :: lacc























end subroutine set_overlap_param_var_GPU

  subroutine set_overlap_param_approx_GPU(this, thermodynamics, decorrelation_length, &
&  istartcol, iendcol, lacc)
use yomhook
use radiation_thermodynamics
class(cloud_type),         intent(inout) :: this
type(thermodynamics_type), intent(in)    :: thermodynamics
real(jprb),                intent(in)    :: decorrelation_length 
integer,         optional, intent(in)    :: istartcol, iendcol











logical, intent (in) :: lacc










end subroutine set_overlap_param_approx_GPU

  subroutine create_fractional_std_GPU(this, ncol, nlev, frac_std, lacc)
use yomhook
class(cloud_type), intent(inout) :: this
integer,           intent(in)    :: ncol, nlev
real(jprb),        intent(in)    :: frac_std
logical, optional, intent(in) :: lacc
















end subroutine create_fractional_std_GPU

  subroutine create_inv_cloud_effective_size_GPU(this, ncol, nlev, inv_eff_size, lacc)
use yomhook
class(cloud_type), intent(inout) :: this
integer,           intent(in)    :: ncol, nlev
real(jprb),        intent(in)    :: inv_eff_size

logical, intent (in) :: lacc





end subroutine create_inv_cloud_effective_size_GPU

  subroutine create_inv_cloud_effective_size_eta_GPU(this, ncol, nlev, &
&  pressure_hl, inv_eff_size_low, inv_eff_size_mid, inv_eff_size_high, &
&  eta_low_mid, eta_mid_high, istartcol, iendcol, lacc)
use yomhook
class(cloud_type), intent(inout) :: this
integer,           intent(in)    :: ncol, nlev

real(jprb),        intent(in)    :: pressure_hl(:,:)

real(jprb),        intent(in)    :: inv_eff_size_low
real(jprb),        intent(in)    :: inv_eff_size_mid
real(jprb),        intent(in)    :: inv_eff_size_high

real(jprb),        intent(in)    :: eta_low_mid, eta_mid_high
integer, optional, intent(in)    :: istartcol, iendcol







logical, intent (in) :: lacc









end subroutine create_inv_cloud_effective_size_eta_GPU

  subroutine param_cloud_effective_separation_eta_GPU(this, ncol, nlev, &
&  pressure_hl, separation_surf, separation_toa, power, &
&  inhom_separation_factor, istartcol, iendcol, lacc)
use yomhook
class(cloud_type), intent(inout) :: this
integer,           intent(in)    :: ncol, nlev

real(jprb),        intent(in)    :: pressure_hl(:,:)


real(jprb),           intent(in) :: separation_surf 
real(jprb),           intent(in) :: separation_toa 
real(jprb),           intent(in) :: power
real(jprb), optional, intent(in) :: inhom_separation_factor
integer,    optional, intent(in) :: istartcol, iendcol











logical, intent (in) :: lacc















end subroutine param_cloud_effective_separation_eta_GPU

  subroutine crop_cloud_fraction_GPU(this, istartcol, iendcol, &
&    cloud_fraction_threshold, cloud_mixing_ratio_threshold, lacc)
use yomhook
class(cloud_type), intent(inout) :: this
integer,           intent(in)    :: istartcol, iendcol


real(jprb) :: cloud_fraction_threshold, cloud_mixing_ratio_threshold


logical, intent (in) :: lacc








end subroutine crop_cloud_fraction_GPU

  function out_of_physical_bounds_GPU(this, istartcol, iendcol, do_fix, lacc) result(is_bad)
use yomhook
use radiation_check
class(cloud_type), intent(inout) :: this
integer,  optional,intent(in) :: istartcol, iendcol
logical,  optional,intent(in) :: do_fix
logical                       :: is_bad


logical, intent (in) :: lacc




end function out_of_physical_bounds_GPU

  subroutine create_device_GPU(this, lacc)
class(cloud_type), intent(inout) :: this
logical, intent (in) :: lacc











end subroutine create_device_GPU

  subroutine update_host_GPU(this, lacc)
class(cloud_type), intent(inout) :: this
logical, intent (in) :: lacc







end subroutine update_host_GPU

  subroutine update_device_GPU(this, lacc)
#ifdef _OPENACC
use openacc
#endif
class(cloud_type), intent(inout) :: this
logical, intent (in) :: lacc


#ifdef __PGI

#endif
#ifdef __PGI

#endif
#ifdef __PGI

#endif
#ifdef __PGI

#endif





end subroutine update_device_GPU

  subroutine delete_device_GPU(this, lacc)
class(cloud_type), intent(inout) :: this
logical, intent (in) :: lacc











end subroutine delete_device_GPU

  subroutine allocate_cloud_arrays_CPU(this, ncol, nlev, ntype, use_inhom_effective_size)
use yomhook



class(cloud_type), intent(inout), target :: this
integer, intent(in)              :: ncol   
integer, intent(in)              :: nlev   




integer, intent(in), optional    :: ntype
logical, intent(in), optional    :: use_inhom_effective_size
















end subroutine allocate_cloud_arrays_CPU

  subroutine deallocate_cloud_arrays_CPU(this)
use yomhook
class(cloud_type), intent(inout) :: this














end subroutine deallocate_cloud_arrays_CPU

  subroutine set_overlap_param_fix_CPU(this, thermodynamics, decorrelation_length, &
&  istartcol, iendcol, lacc)
use yomhook
use radiation_thermodynamics
use radiation_constants
class(cloud_type),         intent(inout) :: this
type(thermodynamics_type), intent(in)    :: thermodynamics
real(jprb),                intent(in)    :: decorrelation_length 
integer,         optional, intent(in)    :: istartcol, iendcol
logical, optional, intent(in) :: lacc




















end subroutine set_overlap_param_fix_CPU

  subroutine set_overlap_param_var_CPU(this, thermodynamics, decorrelation_length, &
&                           istartcol, iendcol, lacc)
use yomhook
use radiation_thermodynamics
use radiation_constants
class(cloud_type),         intent(inout) :: this
type(thermodynamics_type), intent(in)    :: thermodynamics
integer,                   intent(in)    :: istartcol, iendcol
real(jprb),                intent(in)    :: decorrelation_length(istartcol:iendcol) 
logical, optional, intent(in) :: lacc






















end subroutine set_overlap_param_var_CPU

  subroutine set_overlap_param_approx_CPU(this, thermodynamics, decorrelation_length, &
&  istartcol, iendcol)
use yomhook
use radiation_thermodynamics
class(cloud_type),         intent(inout) :: this
type(thermodynamics_type), intent(in)    :: thermodynamics
real(jprb),                intent(in)    :: decorrelation_length 
integer,         optional, intent(in)    :: istartcol, iendcol





















end subroutine set_overlap_param_approx_CPU

  subroutine create_fractional_std_CPU(this, ncol, nlev, frac_std, lacc)
use yomhook
class(cloud_type), intent(inout) :: this
integer,           intent(in)    :: ncol, nlev
real(jprb),        intent(in)    :: frac_std
logical, optional, intent(in) :: lacc













end subroutine create_fractional_std_CPU

  subroutine create_inv_cloud_effective_size_CPU(this, ncol, nlev, inv_eff_size)
use yomhook
class(cloud_type), intent(inout) :: this
integer,           intent(in)    :: ncol, nlev
real(jprb),        intent(in)    :: inv_eff_size






end subroutine create_inv_cloud_effective_size_CPU

  subroutine create_inv_cloud_effective_size_eta_CPU(this, ncol, nlev, &
&  pressure_hl, inv_eff_size_low, inv_eff_size_mid, inv_eff_size_high, &
&  eta_low_mid, eta_mid_high, istartcol, iendcol)
use yomhook
class(cloud_type), intent(inout) :: this
integer,           intent(in)    :: ncol, nlev

real(jprb),        intent(in)    :: pressure_hl(:,:)

real(jprb),        intent(in)    :: inv_eff_size_low
real(jprb),        intent(in)    :: inv_eff_size_mid
real(jprb),        intent(in)    :: inv_eff_size_high

real(jprb),        intent(in)    :: eta_low_mid, eta_mid_high
integer, optional, intent(in)    :: istartcol, iendcol
















end subroutine create_inv_cloud_effective_size_eta_CPU

  subroutine param_cloud_effective_separation_eta_CPU(this, ncol, nlev, &
&  pressure_hl, separation_surf, separation_toa, power, &
&  inhom_separation_factor, istartcol, iendcol)
use yomhook
class(cloud_type), intent(inout) :: this
integer,           intent(in)    :: ncol, nlev

real(jprb),        intent(in)    :: pressure_hl(:,:)


real(jprb),           intent(in) :: separation_surf 
real(jprb),           intent(in) :: separation_toa 
real(jprb),           intent(in) :: power
real(jprb), optional, intent(in) :: inhom_separation_factor
integer,    optional, intent(in) :: istartcol, iendcol


























end subroutine param_cloud_effective_separation_eta_CPU

  subroutine crop_cloud_fraction_CPU(this, istartcol, iendcol, &
&    cloud_fraction_threshold, cloud_mixing_ratio_threshold)
use yomhook
class(cloud_type), intent(inout) :: this
integer,           intent(in)    :: istartcol, iendcol


real(jprb) :: cloud_fraction_threshold, cloud_mixing_ratio_threshold







end subroutine crop_cloud_fraction_CPU

  function out_of_physical_bounds_CPU(this, istartcol, iendcol, do_fix) result(is_bad)
use yomhook
use radiation_check
class(cloud_type), intent(inout) :: this
integer,  optional,intent(in) :: istartcol, iendcol
logical,  optional,intent(in) :: do_fix
logical                       :: is_bad






end function out_of_physical_bounds_CPU

end module radiation_cloud

