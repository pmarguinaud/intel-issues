! radiation_ecckd.F90 - ecCKD generalized gas optics model
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

module radiation_ecckd

  use parkind1
  use radiation_gas_constants
  use radiation_ecckd_gas
  use radiation_spectral_definition

  implicit none

  

  !---------------------------------------------------------------------
  ! This derived type contains all the data needed to describe a
  ! correlated k-distribution gas optics model created using the ecCKD
  ! tool
  type ckd_model_type

    ! Gas information

    ! Number of gases
    integer :: ngas = 0
    ! Array of individual-gas data objects
    type(ckd_gas_type), allocatable :: single_gas(:)
    ! Mapping from the "gas codes" in the radiation_gas_constants
    ! module to an index to the single_gas array, where zero means gas
    ! not present (or part of a composite gas)
    integer :: i_gas_mapping(0:NMaxGases)

    ! Coordinates of main look-up table for absorption coeffts

    ! Number of pressure and temperature points
    integer :: npress = 0
    integer :: ntemp  = 0
    ! Natural logarithm of first (lowest) pressure (Pa) and increment
    real(jprb) :: log_pressure1, d_log_pressure
    ! First temperature profile (K), dimensioned (npress)
    real(jprb), allocatable :: temperature1(:)
    ! Temperature increment (K)
    real(jprb) :: d_temperature

    ! Look-up table for Planck function

    ! Number of entries
    integer :: nplanck = 0
    ! Temperature of first element of look-up table and increment (K)
    real(jprb), allocatable :: temperature1_planck
    real(jprb), allocatable :: d_temperature_planck
    ! Planck function (black body flux into a horizontal plane) in W
    ! m-2, dimensioned (ng,nplanck)
    real(jprb), allocatable :: planck_function(:,:)

    ! Normalized solar irradiance in each g point, dimensioned (ng)
    real(jprb), allocatable :: norm_solar_irradiance(:)

    ! Normalized amplitude of variations in the solar irradiance
    ! through the solar cycle in each g point, dimensioned (ng).
    ! Since the user always provides the solar irradiance SI
    ! integrated across the spectrum, this variable must sum to zero:
    ! this ensures that the solar irradiance in each g-point is
    ! SSI=SI*(norm_solar_irradiance +
    ! A*norm_amplitude_solar_irradiance) for any A, where A denotes
    ! the amplitude of deviations from the mean solar spectrum,
    ! typically between -1.0 and 1.0 and provided by
    ! single_level%solar_spectral_multiplier.
    real(jprb), allocatable :: norm_amplitude_solar_irradiance(:)

    ! Rayleigh molar scattering coefficient in m2 mol-1 in each g
    ! point
    real(jprb), allocatable :: rayleigh_molar_scat(:)

    ! ! Spectral mapping of g points

    ! ! Number of wavenumber intervals
    ! integer :: nwav = 0

    ! Number of k terms / g points
    integer :: ng   = 0

    ! Spectral definition describing bands and g points
    type(spectral_definition_type) :: spectral_def

    ! Shortwave: true, longwave: false
    logical :: is_sw

  contains

    
    
! Vectorized version of the optical depth look-up performs better on
! NEC, but slower on x86
    
    
    
    
!    procedure :: deallocate => deallocate_ckd_model

  procedure :: read_GPU => read_ckd_model_GPU

  procedure :: read_spectral_solar_cycle_GPU

  procedure :: calc_optical_depth_GPU => calc_optical_depth_ckd_model_GPU

  procedure :: print_GPU => print_ckd_model_GPU

  procedure :: calc_planck_function_GPU

  procedure :: calc_incoming_sw_GPU

  procedure :: read_CPU => read_ckd_model_CPU

  procedure :: read_spectral_solar_cycle_CPU

  procedure :: calc_optical_depth_CPU => calc_optical_depth_ckd_model_CPU

  procedure :: print_CPU => print_ckd_model_CPU

  procedure :: calc_planck_function_CPU

  procedure :: calc_incoming_sw_CPU

  end type ckd_model_type


contains

  !---------------------------------------------------------------------
  ! Read a complete ecCKD gas optics model from a NetCDF file
  ! "filename"
  

  !---------------------------------------------------------------------
  ! Print a description of the correlated k-distribution model to the
  ! "nulout" unit
  


  !---------------------------------------------------------------------
  ! Read the amplitude of the spectral variations associated with the
  ! solar cycle and map to g-points
  


  !---------------------------------------------------------------------
  ! Compute layerwise optical depth for each g point for ncol columns
  ! at nlev layers
  


  !---------------------------------------------------------------------
  ! Vectorized variant of above routine
  


  !---------------------------------------------------------------------
  ! Calculate the Planck function integrated across each of the g
  ! points of this correlated k-distribution model, for a given
  ! temperature, where Planck function is defined as the flux emitted
  ! by a black body (rather than radiance)
  


  !---------------------------------------------------------------------
  ! Return the spectral solar irradiance integrated over each g point
  ! of a solar correlated k-distribution model, given the
  ! total_solar_irradiance
  

  subroutine read_ckd_model_GPU(this, filename, iverbose, lacc)
use easy_netcdf

use yomhook
class(ckd_model_type), intent(inout) :: this
character(len=*),      intent(in)    :: filename
integer, optional,     intent(in)    :: iverbose










logical, intent (in) :: lacc









;









;











end subroutine read_ckd_model_GPU

  subroutine print_ckd_model_GPU(this, lacc)
use radiation_io
use radiation_gas_constants
class(ckd_model_type), intent(in)  :: this

logical, intent (in) :: lacc





end subroutine print_ckd_model_GPU

  subroutine read_spectral_solar_cycle_GPU(this, filename, iverbose, use_updated_solar_spectrum, lacc)
use easy_netcdf
use radiation_io
use yomhook


class(ckd_model_type), intent(inout) :: this
character(len=*),      intent(in)    :: filename
integer, optional,     intent(in)    :: iverbose


logical, optional,     intent(in)    :: use_updated_solar_spectrum



 
 
 
















logical, intent (in) :: lacc
































end subroutine read_spectral_solar_cycle_GPU

  subroutine calc_optical_depth_ckd_model_GPU(this, ncol, nlev, istartcol, iendcol, nmaxgas, &
&  pressure_hl, temperature_fl, mole_fraction_fl, &
&  optical_depth_fl, rayleigh_od_fl, concentration_scaling, lacc)
use yomhook
use radiation_constants

class(ckd_model_type), intent(in), target  :: this

integer,               intent(in)  :: ncol, nlev, nmaxgas, istartcol, iendcol

real(jprb),            intent(in)  :: pressure_hl(ncol,nlev+1)

real(jprb),            intent(in)  :: temperature_fl(istartcol:iendcol,nlev)

real(jprb),            intent(in)  :: mole_fraction_fl(ncol,nlev,nmaxgas)

real(jprb), optional,  intent(in)  :: concentration_scaling(nmaxgas)


real(jprb),            intent(out) :: optical_depth_fl(this%ng,nlev,istartcol:iendcol)

real(jprb),  optional, intent(out) :: rayleigh_od_fl(this%ng,nlev,istartcol:iendcol)



















logical, intent (in) :: lacc




end subroutine calc_optical_depth_ckd_model_GPU

  subroutine calc_optical_depth_ckd_model_vec_GPU(this, ncol, nlev, istartcol, iendcol, nmaxgas, &
&  pressure_hl, temperature_fl, mole_fraction_fl, &
&  optical_depth_fl, rayleigh_od_fl, lacc)
use yomhook
use radiation_constants

class(ckd_model_type), intent(in), target  :: this

integer,               intent(in)  :: ncol, nlev, nmaxgas, istartcol, iendcol

real(jprb),            intent(in)  :: pressure_hl(ncol,nlev+1)

real(jprb),            intent(in)  :: temperature_fl(istartcol:iendcol,nlev)

real(jprb),            intent(in)  :: mole_fraction_fl(ncol,nlev,nmaxgas)


real(jprb),            intent(out) :: optical_depth_fl(this%ng,nlev,istartcol:iendcol)

real(jprb),  optional, intent(out) :: rayleigh_od_fl(this%ng,nlev,istartcol:iendcol)




















logical, intent (in) :: lacc






end subroutine calc_optical_depth_ckd_model_vec_GPU

  subroutine calc_planck_function_GPU(this, nt, temperature, planck, lacc)
class(ckd_model_type), intent(in)  :: this
integer,    intent(in)  :: nt
real(jprb), intent(in)  :: temperature(:) 
real(jprb), intent(out) :: planck(this%ng,nt) 


logical, intent (in) :: lacc

end subroutine calc_planck_function_GPU

  subroutine calc_incoming_sw_GPU(this, total_solar_irradiance, &
&                      spectral_solar_irradiance, &
&                      solar_spectral_multiplier, lacc)
use radiation_io
class(ckd_model_type), intent(in)    :: this
real(jprb),            intent(in)    :: total_solar_irradiance 
real(jprb),            intent(inout) :: spectral_solar_irradiance(:,:) 
real(jprb), optional,  intent(in)    :: solar_spectral_multiplier
logical, intent (in) :: lacc

end subroutine calc_incoming_sw_GPU

  subroutine read_ckd_model_CPU(this, filename, iverbose)
use easy_netcdf

use yomhook
class(ckd_model_type), intent(inout) :: this
character(len=*),      intent(in)    :: filename
integer, optional,     intent(in)    :: iverbose



















;









;











end subroutine read_ckd_model_CPU

  subroutine print_ckd_model_CPU(this)
use radiation_io
use radiation_gas_constants
class(ckd_model_type), intent(in)  :: this






end subroutine print_ckd_model_CPU

  subroutine read_spectral_solar_cycle_CPU(this, filename, iverbose, use_updated_solar_spectrum)
use easy_netcdf
use radiation_io
use yomhook


class(ckd_model_type), intent(inout) :: this
character(len=*),      intent(in)    :: filename
integer, optional,     intent(in)    :: iverbose


logical, optional,     intent(in)    :: use_updated_solar_spectrum



 
 
 
















































end subroutine read_spectral_solar_cycle_CPU

  subroutine calc_optical_depth_ckd_model_CPU(this, ncol, nlev, istartcol, iendcol, nmaxgas, &
&  pressure_hl, temperature_fl, mole_fraction_fl, &
&  optical_depth_fl, rayleigh_od_fl, concentration_scaling)
use yomhook
use radiation_constants

class(ckd_model_type), intent(in), target  :: this

integer,               intent(in)  :: ncol, nlev, nmaxgas, istartcol, iendcol

real(jprb),            intent(in)  :: pressure_hl(ncol,nlev+1)

real(jprb),            intent(in)  :: temperature_fl(istartcol:iendcol,nlev)

real(jprb),            intent(in)  :: mole_fraction_fl(ncol,nlev,nmaxgas)

real(jprb), optional,  intent(in)  :: concentration_scaling(nmaxgas)


real(jprb),            intent(out) :: optical_depth_fl(this%ng,nlev,istartcol:iendcol)

real(jprb),  optional, intent(out) :: rayleigh_od_fl(this%ng,nlev,istartcol:iendcol)























end subroutine calc_optical_depth_ckd_model_CPU

  subroutine calc_optical_depth_ckd_model_vec_CPU(this, ncol, nlev, istartcol, iendcol, nmaxgas, &
&  pressure_hl, temperature_fl, mole_fraction_fl, &
&  optical_depth_fl, rayleigh_od_fl)
use yomhook
use radiation_constants

class(ckd_model_type), intent(in), target  :: this

integer,               intent(in)  :: ncol, nlev, nmaxgas, istartcol, iendcol

real(jprb),            intent(in)  :: pressure_hl(ncol,nlev+1)

real(jprb),            intent(in)  :: temperature_fl(istartcol:iendcol,nlev)

real(jprb),            intent(in)  :: mole_fraction_fl(ncol,nlev,nmaxgas)


real(jprb),            intent(out) :: optical_depth_fl(this%ng,nlev,istartcol:iendcol)

real(jprb),  optional, intent(out) :: rayleigh_od_fl(this%ng,nlev,istartcol:iendcol)


























end subroutine calc_optical_depth_ckd_model_vec_CPU

  subroutine calc_planck_function_CPU(this, nt, temperature, planck)
class(ckd_model_type), intent(in)  :: this
integer,    intent(in)  :: nt
real(jprb), intent(in)  :: temperature(:) 
real(jprb), intent(out) :: planck(this%ng,nt) 



end subroutine calc_planck_function_CPU

  subroutine calc_incoming_sw_CPU(this, total_solar_irradiance, &
&                      spectral_solar_irradiance, &
&                      solar_spectral_multiplier)
use radiation_io
class(ckd_model_type), intent(in)    :: this
real(jprb),            intent(in)    :: total_solar_irradiance 
real(jprb),            intent(inout) :: spectral_solar_irradiance(:,:) 
real(jprb), optional,  intent(in)    :: solar_spectral_multiplier

end subroutine calc_incoming_sw_CPU

end module radiation_ecckd

