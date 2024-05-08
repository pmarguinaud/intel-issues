! radiation_spectral_definition.F90 - Derived type to describe a spectral definition
!
! (C) Copyright 2020- ECMWF.
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
! License: see the COPYING file for details
!

#include "ecrad_config.h"

module radiation_spectral_definition

  use parkind1

  implicit none

  

  real(jprb), parameter :: SolarReferenceTemperature       = 5777.0_jprb ! K
  real(jprb), parameter :: TerrestrialReferenceTemperature = 273.15_jprb ! K

  !---------------------------------------------------------------------
  ! A derived type describing the contribution of the g points of a
  ! correlated k-distribution gas-optics model from each part of the
  ! spectrum. This is used primarily to map the cloud and aerosol
  ! optical properties on to the gas g points.
  type spectral_definition_type
    
    ! Spectral mapping of g points

    ! Number of wavenumber intervals
    integer :: nwav = 0
    ! Number of k terms / g points
    integer :: ng   = 0
    ! Start and end wavenumber (cm-1), dimensioned (nwav)
    real(jprb), allocatable :: wavenumber1(:)
    real(jprb), allocatable :: wavenumber2(:)
    ! Fraction of each g point in each wavenumber interval,
    ! dimensioned (nwav, ng)
    real(jprb), allocatable :: gpoint_fraction(:,:)

    ! Spectral weighting information for generating mappings to/from
    ! different spectral grids: this can be in terms of a reference
    ! temperature (K) to generate a Planck function, or the
    ! solar_spectral_irradiance (W m-2) if available in the gas-optics
    ! file.
    real(jprb) :: reference_temperature = -1.0_jprb
    real(jprb), allocatable :: solar_spectral_irradiance(:)
    
    ! Band information

    ! Number of bands
    integer :: nband = 0
    ! Lower and upper bounds of wavenumber bands (cm-1), dimensioned
    ! (nband)
    real(jprb), allocatable :: wavenumber1_band(:)
    real(jprb), allocatable :: wavenumber2_band(:)
    ! Band (one based) to which each g point belongs
    integer,    allocatable :: i_band_number(:)

  contains
    
    
    
    
    
    
    
    
    

  procedure :: read_GPU => read_spectral_definition_GPU

  procedure :: allocate_bands_only_GPU

  procedure :: deallocate_GPU

  procedure :: find_GPU => find_wavenumber_GPU

  procedure :: calc_mapping_GPU

  procedure :: calc_mapping_from_bands_GPU

  procedure :: calc_mapping_from_wavenumber_bands_GPU

  procedure :: print_mapping_from_bands_GPU

  procedure :: min_wavenumber_GPU, max_wavenumber_GPU

  procedure :: read_CPU => read_spectral_definition_CPU

  procedure :: allocate_bands_only_CPU

  procedure :: deallocate_CPU

  procedure :: find_CPU => find_wavenumber_CPU

  procedure :: calc_mapping_CPU

  procedure :: calc_mapping_from_bands_CPU

  procedure :: calc_mapping_from_wavenumber_bands_CPU

  procedure :: print_mapping_from_bands_CPU

  procedure :: min_wavenumber_CPU, max_wavenumber_CPU

  end type spectral_definition_type

contains

  !---------------------------------------------------------------------
  ! Read the description of a spectral definition from a NetCDF
  ! file of the type used to describe an ecCKD model
  


  !---------------------------------------------------------------------
  ! Store a simple band description by copying over the reference
  ! temperature and the lower and upper wavenumbers of each band
  


  !---------------------------------------------------------------------
  ! Deallocate memory inside a spectral definition object
  


  !---------------------------------------------------------------------
  ! Find the index to the highest wavenumber in the spectral
  ! definition that is lower than or equal to "wavenumber", used for
  ! implementing look-up tables
  


  !---------------------------------------------------------------------
  ! Compute a mapping matrix "mapping" that can be used in an
  ! expression y=matmul(mapping,x) where x is a variable containing
  ! optical properties at each input "wavenumber", and y is this
  ! variable mapped on to the spectral intervals in the spectral
  ! definition "this". 
  


  !---------------------------------------------------------------------
  ! Under normal operation (if use_fluxes is .false. or not present),
  ! compute a mapping matrix "mapping" that can be used in an
  ! expression y=matmul(mapping^T,x) where x is a variable containing
  ! optical properties in input bands (e.g. albedo in shortwave albedo
  ! bands), and y is this variable mapped on to the spectral intervals
  ! in the spectral definition "this". Note that "mapping" is here
  ! transposed from the convention in the calc_mapping routine.  Under
  ! the alternative operation (if use_fluxes is present and .true.),
  ! the mapping works in the reverse sense: if y contains fluxes in
  ! each ecRad band or g-point, then x=matmul(mapping,y) would return
  ! fluxes in x averaged to user-supplied "input" bands. In this
  ! version, the bands are described by their wavelength bounds
  ! (wavelength_bound, which must be increasing and exclude the end
  ! points) and the index of the mapping matrix that each band
  ! corresponds to (i_intervals, which has one more element than
  ! wavelength_bound and can have duplicated values if an
  ! albedo/emissivity value is to be associated with more than one
  ! discontinuous ranges of the spectrum).
  


  !---------------------------------------------------------------------
  ! As calc_mapping_from_bands but in terms of wavenumber bounds from
  ! wavenumber1 to wavenumber2
  


  !---------------------------------------------------------------------
  ! Print out the mapping computed by calc_mapping_from_bands
  


  !---------------------------------------------------------------------
  ! Return the minimum wavenumber of this object in cm-1
  


  !---------------------------------------------------------------------
  ! Return the maximum wavenumber of this object in cm-1
  


  !---------------------------------------------------------------------
  ! Return the Planck function (in W m-2 (cm-1)-1) for a given
  ! wavenumber (cm-1) and temperature (K), ensuring double precision
  ! for internal calculation.  If temperature is 0 or less then unity
  ! is returned; since this function is primarily used to weight an
  ! integral by the Planck function, a temperature of 0 or less means
  ! no weighting is to be applied.
  

  subroutine read_spectral_definition_GPU(this, file, lacc)
use easy_netcdf
use yomhook
class(spectral_definition_type), intent(inout) :: this
type(netcdf_file),               intent(inout) :: file

logical, intent (in) :: lacc















;


end subroutine read_spectral_definition_GPU

  subroutine allocate_bands_only_GPU(this, reference_temperature, wavenumber1, wavenumber2, lacc)
use yomhook
class(spectral_definition_type), intent(inout) :: this
real(jprb),                      intent(in)    :: reference_temperature    
real(jprb),        dimension(:), intent(in)    :: wavenumber1, wavenumber2 

logical, intent (in) :: lacc









end subroutine allocate_bands_only_GPU

  subroutine deallocate_GPU(this, lacc)
class(spectral_definition_type), intent(inout) :: this
logical, intent (in) :: lacc










end subroutine deallocate_GPU

  pure function find_wavenumber_GPU(this, wavenumber, lacc)
class(spectral_definition_type), intent(in) :: this
real(jprb),                      intent(in) :: wavenumber 
integer                                     :: find_wavenumber_GPU
logical, intent (in) :: lacc

end function find_wavenumber_GPU

  subroutine calc_mapping_GPU(this, wavenumber, mapping, weighting_temperature, use_bands, lacc)
use yomhook
use radiation_io
class(spectral_definition_type), intent(in)    :: this
real(jprb),                      intent(in)    :: wavenumber(:) 
real(jprb), allocatable,         intent(inout) :: mapping(:,:)
real(jprb), optional,            intent(in)    :: weighting_temperature 
logical,    optional,            intent(in)    :: use_bands





 








logical, intent (in) :: lacc







end subroutine calc_mapping_GPU

  subroutine calc_mapping_from_bands_GPU(this, &
&  wavelength_bound, i_intervals, mapping, use_bands, use_fluxes, lacc)
use yomhook
use radiation_io
class(spectral_definition_type), intent(in)    :: this



real(jprb),                      intent(in)    :: wavelength_bound(:)

integer,                         intent(in)    :: i_intervals(:)
real(jprb), allocatable,         intent(inout) :: mapping(:,:)
logical,    optional,            intent(in)    :: use_bands
logical,    optional,            intent(in)    :: use_fluxes


         
 

























logical, intent (in) :: lacc













end subroutine calc_mapping_from_bands_GPU

  subroutine calc_mapping_from_wavenumber_bands_GPU(this, &
&  wavenumber1, wavenumber2, mapping, use_bands, use_fluxes, lacc)
use yomhook
class(spectral_definition_type), intent(in)    :: this
real(jprb), intent(in)    :: wavenumber1(:), wavenumber2(:)
real(jprb), allocatable,         intent(inout) :: mapping(:,:)
logical,    optional,            intent(in)    :: use_bands
logical,    optional,            intent(in)    :: use_fluxes





















logical, intent (in) :: lacc







end subroutine calc_mapping_from_wavenumber_bands_GPU

  subroutine print_mapping_from_bands_GPU(this, mapping, use_bands, lacc)
use radiation_io
class(spectral_definition_type), intent(in) :: this
real(jprb), allocatable,         intent(in) :: mapping(:,:) 
logical,    optional,            intent(in) :: use_bands



logical, intent (in) :: lacc




end subroutine print_mapping_from_bands_GPU

  pure function min_wavenumber_GPU(this, lacc)
class(spectral_definition_type), intent(in)    :: this
real(jprb) :: min_wavenumber_GPU
logical, intent (in) :: lacc

end function min_wavenumber_GPU

  pure function max_wavenumber_GPU(this, lacc)
class(spectral_definition_type), intent(in)    :: this
real(jprb) :: max_wavenumber_GPU
logical, intent (in) :: lacc

end function max_wavenumber_GPU

  elemental function calc_planck_function_wavenumber_GPU(wavenumber, temperature, lacc)
use parkind1
use radiation_constants
real(jprb), intent(in) :: wavenumber  
real(jprb), intent(in) :: temperature 
real(jprb) :: calc_planck_function_wavenumber_GPU
 

logical, intent (in) :: lacc 

end function calc_planck_function_wavenumber_GPU

  subroutine read_spectral_definition_CPU(this, file)
use easy_netcdf
use yomhook
class(spectral_definition_type), intent(inout) :: this
type(netcdf_file),               intent(inout) :: file
















;


end subroutine read_spectral_definition_CPU

  subroutine allocate_bands_only_CPU(this, reference_temperature, wavenumber1, wavenumber2)
use yomhook
class(spectral_definition_type), intent(inout) :: this
real(jprb),                      intent(in)    :: reference_temperature    
real(jprb),        dimension(:), intent(in)    :: wavenumber1, wavenumber2 










end subroutine allocate_bands_only_CPU

  subroutine deallocate_CPU(this)
class(spectral_definition_type), intent(inout) :: this










end subroutine deallocate_CPU

  pure function find_wavenumber_CPU(this, wavenumber)
class(spectral_definition_type), intent(in) :: this
real(jprb),                      intent(in) :: wavenumber 
integer                                     :: find_wavenumber_CPU

end function find_wavenumber_CPU

  subroutine calc_mapping_CPU(this, wavenumber, mapping, weighting_temperature, use_bands)
use yomhook
use radiation_io
class(spectral_definition_type), intent(in)    :: this
real(jprb),                      intent(in)    :: wavenumber(:) 
real(jprb), allocatable,         intent(inout) :: mapping(:,:)
real(jprb), optional,            intent(in)    :: weighting_temperature 
logical,    optional,            intent(in)    :: use_bands





 















end subroutine calc_mapping_CPU

  subroutine calc_mapping_from_bands_CPU(this, &
&  wavelength_bound, i_intervals, mapping, use_bands, use_fluxes)
use yomhook
use radiation_io
class(spectral_definition_type), intent(in)    :: this



real(jprb),                      intent(in)    :: wavelength_bound(:)

integer,                         intent(in)    :: i_intervals(:)
real(jprb), allocatable,         intent(inout) :: mapping(:,:)
logical,    optional,            intent(in)    :: use_bands
logical,    optional,            intent(in)    :: use_fluxes


         
 






































end subroutine calc_mapping_from_bands_CPU

  subroutine calc_mapping_from_wavenumber_bands_CPU(this, &
&  wavenumber1, wavenumber2, mapping, use_bands, use_fluxes)
use yomhook
class(spectral_definition_type), intent(in)    :: this
real(jprb), intent(in)    :: wavenumber1(:), wavenumber2(:)
real(jprb), allocatable,         intent(inout) :: mapping(:,:)
logical,    optional,            intent(in)    :: use_bands
logical,    optional,            intent(in)    :: use_fluxes




























end subroutine calc_mapping_from_wavenumber_bands_CPU

  subroutine print_mapping_from_bands_CPU(this, mapping, use_bands)
use radiation_io
class(spectral_definition_type), intent(in) :: this
real(jprb), allocatable,         intent(in) :: mapping(:,:) 
logical,    optional,            intent(in) :: use_bands







end subroutine print_mapping_from_bands_CPU

  pure function min_wavenumber_CPU(this)
class(spectral_definition_type), intent(in)    :: this
real(jprb) :: min_wavenumber_CPU

end function min_wavenumber_CPU

  pure function max_wavenumber_CPU(this)
class(spectral_definition_type), intent(in)    :: this
real(jprb) :: max_wavenumber_CPU

end function max_wavenumber_CPU

  elemental function calc_planck_function_wavenumber_CPU(wavenumber, temperature)
use parkind1
use radiation_constants
real(jprb), intent(in) :: wavenumber  
real(jprb), intent(in) :: temperature 
real(jprb) :: calc_planck_function_wavenumber_CPU
 
 

end function calc_planck_function_wavenumber_CPU

end module radiation_spectral_definition

