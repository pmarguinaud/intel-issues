! radiation_general_cloud_optics_data.F90 - Type to store generalized cloud optical properties
!
! (C) Copyright 2019- ECMWF.
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

module radiation_general_cloud_optics_data

  use parkind1

  implicit none

  

  !---------------------------------------------------------------------
  ! This type holds the configuration information to compute optical
  ! properties for a particular type of cloud or hydrometeor in one of
  ! the shortwave or longwave
  type general_cloud_optics_type
    ! Band-specific (or g-point-specific) values as a look-up table
    ! versus effective radius dimensioned (nband,n_effective_radius)

    ! Extinction coefficient per unit mass (m2 kg-1)
    real(jprb), allocatable, dimension(:,:) :: &
         &  mass_ext

    ! Single-scattering albedo and asymmetry factor (dimensionless)
    real(jprb), allocatable, dimension(:,:) :: &
         &  ssa, asymmetry

    ! Number of effective radius coefficients, start value and
    ! interval in look-up table
    integer    :: n_effective_radius = 0
    real(jprb) :: effective_radius_0, d_effective_radius

    ! Name of cloud/precip type and scattering model
    ! (e.g. "mie_droplet", "fu-muskatel_ice"). These are used to
    ! generate the name of the data file from which the coefficients
    ! are read.
    character(len=511) :: type_name

    ! Do we use bands or g-points?
    logical :: use_bands = .false.

   contains
     
     
     

  procedure :: setup_GPU => setup_general_cloud_optics_data_GPU

  procedure :: add_optical_properties_GPU

  procedure :: save_GPU => save_general_cloud_optics_data_GPU

  procedure :: setup_CPU => setup_general_cloud_optics_data_CPU

  procedure :: add_optical_properties_CPU

  procedure :: save_CPU => save_general_cloud_optics_data_CPU

  end type general_cloud_optics_type

contains

  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Setup cloud optics coefficients by reading them from a file
  


  !---------------------------------------------------------------------
  ! Add the optical properties of a particular cloud type to the
  ! accumulated optical properties of all cloud types
  


  !---------------------------------------------------------------------
  ! Return the Planck function (in W m-2 (cm-1)-1) for a given
  ! wavenumber (cm-1) and temperature (K), ensuring double precision
  ! for internal calculation
  

  !---------------------------------------------------------------------
  ! Save cloud optical properties in the named file
  

  subroutine setup_general_cloud_optics_data_GPU(this, file_name, specdef, &
&                                use_bands, use_thick_averaging, &
&                                weighting_temperature, &
&                                iverbose, lacc)
use yomhook
use easy_netcdf
use radiation_spectral_definition
use radiation_io
class(general_cloud_optics_type), intent(inout)    :: this
character(len=*), intent(in)               :: file_name
type(spectral_definition_type), intent(in) :: specdef
logical, intent(in), optional              :: use_bands, use_thick_averaging
real(jprb), intent(in), optional           :: weighting_temperature 
integer, intent(in), optional              :: iverbose







       
 









  
 


logical, intent (in) :: lacc











































end subroutine setup_general_cloud_optics_data_GPU

  subroutine add_optical_properties_GPU(this, ng_in, nlev, ncol, &
&                            cloud_fraction, &
&                            water_path, effective_radius, &
&                            od, scat_od, scat_asymmetry, lacc)
use yomhook
class(general_cloud_optics_type), intent(in) :: this



integer, intent(in) :: ng_in, nlev, ncol

real(jprb), intent(in) :: cloud_fraction(ncol,nlev)
real(jprb), intent(in) :: water_path(ncol,nlev)       
real(jprb), intent(in) :: effective_radius(ncol,nlev) 


real(jprb), intent(inout), dimension(ng_in,nlev,ncol) &
&  :: od             
real(jprb), intent(inout), dimension(ng_in,nlev,ncol), optional &
&  :: scat_od, &     
&     scat_asymmetry 





logical, intent (in) :: lacc



end subroutine add_optical_properties_GPU

  elemental function calc_planck_function_wavenumber_GPU(wavenumber, temperature, lacc)
use parkind1
use radiation_constants
real(jprb), intent(in) :: wavenumber  
real(jprb), intent(in) :: temperature 
real(jprb) :: calc_planck_function_wavenumber_GPU
 

logical, intent (in) :: lacc 



end function calc_planck_function_wavenumber_GPU

  subroutine save_general_cloud_optics_data_GPU(this, file_name, iverbose, lacc)
use yomhook
use easy_netcdf
class(general_cloud_optics_type), intent(in) :: this
character(len=*),                 intent(in) :: file_name
integer,                optional, intent(in) :: iverbose





logical, intent (in) :: lacc






















end subroutine save_general_cloud_optics_data_GPU

  subroutine setup_general_cloud_optics_data_CPU(this, file_name, specdef, &
&                                use_bands, use_thick_averaging, &
&                                weighting_temperature, &
&                                iverbose)
use yomhook
use easy_netcdf
use radiation_spectral_definition
use radiation_io
class(general_cloud_optics_type), intent(inout)    :: this
character(len=*), intent(in)               :: file_name
type(spectral_definition_type), intent(in) :: specdef
logical, intent(in), optional              :: use_bands, use_thick_averaging
real(jprb), intent(in), optional           :: weighting_temperature 
integer, intent(in), optional              :: iverbose







       
 









  
 













































end subroutine setup_general_cloud_optics_data_CPU

  subroutine add_optical_properties_CPU(this, ng_in, nlev, ncol, &
&                            cloud_fraction, &
&                            water_path, effective_radius, &
&                            od, scat_od, scat_asymmetry)
use yomhook
class(general_cloud_optics_type), intent(in) :: this



integer, intent(in) :: ng_in, nlev, ncol

real(jprb), intent(in) :: cloud_fraction(ncol,nlev)
real(jprb), intent(in) :: water_path(ncol,nlev)       
real(jprb), intent(in) :: effective_radius(ncol,nlev) 


real(jprb), intent(inout), dimension(ng_in,nlev,ncol) &
&  :: od             
real(jprb), intent(inout), dimension(ng_in,nlev,ncol), optional &
&  :: scat_od, &     
&     scat_asymmetry 








end subroutine add_optical_properties_CPU

  elemental function calc_planck_function_wavenumber_CPU(wavenumber, temperature)
use parkind1
use radiation_constants
real(jprb), intent(in) :: wavenumber  
real(jprb), intent(in) :: temperature 
real(jprb) :: calc_planck_function_wavenumber_CPU
 
 



end function calc_planck_function_wavenumber_CPU

  subroutine save_general_cloud_optics_data_CPU(this, file_name, iverbose)
use yomhook
use easy_netcdf
class(general_cloud_optics_type), intent(in) :: this
character(len=*),                 intent(in) :: file_name
integer,                optional, intent(in) :: iverbose



























end subroutine save_general_cloud_optics_data_CPU

end module radiation_general_cloud_optics_data

