! radiation_aerosol_optics_data.F90 - Type to store aerosol optical properties
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
! Modifications
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-04-20  A. Bozzo  Read optical properties at selected wavelengths

#include "ecrad_config.h"

module radiation_aerosol_optics_data

  use parkind1
  use radiation_io

  implicit none
  

  

  ! The following provide possible values for
  ! aerosol_optics_type%iclass, which is used to map the user's type
  ! index to the hydrophobic or hydrophilic types loaded from the
  ! aerosol optics file. Initially iclass is equal to
  ! AerosolClassUndefined, which will throw an error if ever the user
  ! tries to use this aerosol type. The user may specify that an
  ! aerosol type is to be ignored in the radiation calculation, in
  ! which case iclass will be set equal to AerosolClassIgnored.
  enum, bind(c)
     enumerator IAerosolClassUndefined,   IAerosolClassIgnored, &
          &     IAerosolClassHydrophobic, IAerosolClassHydrophilic
  end enum

  integer, parameter :: NMaxStringLength = 2000
  integer, parameter :: NMaxLineLength   = 200

  !---------------------------------------------------------------------
  ! This type holds the configuration information to compute
  ! aerosol optical properties
  type aerosol_optics_type
     ! A vector of length ntype, iclass maps user-defined types on to
     ! the hydrophilic or hydrophobic aerosol classes using the
     ! enumerators above
     integer, allocatable, dimension(:) :: iclass

     ! Also a vector of length ntype, itype maps user-defined types on
     ! to specific hydrophilic or hydrophobic aerosol types
     integer, allocatable, dimension(:) :: itype

     ! Wavenumber (cm-1) upper and lower bounds of each spectral
     ! interval, which if used in the RRTMG gas optics scheme should
     ! match its band bounds
     real(jprb), allocatable, dimension(:) :: wavenumber1_sw, wavenumber2_sw
     real(jprb), allocatable, dimension(:) :: wavenumber1_lw, wavenumber2_lw

     ! Scattering properties are provided separately in the shortwave
     ! and longwave for hydrophobic and hydrophilic aerosols.
     ! Hydrophobic aerosols are dimensioned (nband,n_type_phobic):
     real(jprb), allocatable, dimension(:,:) :: &
          &  mass_ext_sw_phobic, & ! Mass-extinction coefficient (m2 kg-1)
          &  ssa_sw_phobic,      & ! Single scattering albedo
          &  g_sw_phobic,        & ! Asymmetry factor
!          &  ssa_g_sw_phobic,    & ! ssa*g
          &  mass_ext_lw_phobic, & ! Mass-extinction coefficient (m2 kg-1)
!          &  mass_abs_lw_phobic, & ! Mass-absorption coefficient (m2 kg-1)
          &  ssa_lw_phobic,      & ! Single scattering albedo
          &  g_lw_phobic           ! Asymmetry factor

     ! Hydrophilic aerosols are dimensioned (nband,nrh,n_type_philic):
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  mass_ext_sw_philic, & ! Mass-extinction coefficient (m2 kg-1)
          &  ssa_sw_philic,      & ! Single scattering albedo
          &  g_sw_philic,        & ! Asymmetry factor
 !         &  ssa_g_sw_philic,    & ! ssa*g
          &  mass_ext_lw_philic, & ! Mass-extinction coefficient (m2 kg-1)
 !         &  mass_abs_lw_philic, & ! Mass-absorption coefficient (m2 kg-1)
          &  ssa_lw_philic,      & ! Single scattering albedo
          &  g_lw_philic           ! Asymmetry factor

     ! Wavelengths at which monochromatic properties are stored,
     ! dimensioned (n_mono_wl), units metres
     real(jprb), allocatable :: wavelength_mono(:)

     ! Scattering properties at selected monochromatic wavelengths
     ! (n_mono_wl,n_type_phobic)
     real(jprb), allocatable, dimension(:,:) :: &
          &  mass_ext_mono_phobic, & ! Mass-extinction coefficient (m2 kg-1)
          &  ssa_mono_phobic,      & ! Single scattering albedo
          &  g_mono_phobic,        & ! Asymmetry factor
          &  lidar_ratio_mono_phobic ! Lidar Ratio
     ! ...hydrophilic aerosols dimensioned (n_mono_wl,nrh,n_type_philic):
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  mass_ext_mono_philic, & ! Mass-extinction coefficient (m2 kg-1)
          &  ssa_mono_philic,      & ! Single scattering albedo
          &  g_mono_philic,        & ! Asymmetry factor
          &  lidar_ratio_mono_philic ! Lidar Ratio

     ! For hydrophilic aerosols, the lower bounds of the relative
     ! humidity bins is a vector of length nrh:
     real(jprb), allocatable, dimension(:) :: &
          &  rh_lower    ! Dimensionless (1.0 = 100% humidity)

     ! Strings describing the aerosol types
     character(len=NMaxStringLength) :: description_phobic_str = ' '
     character(len=NMaxStringLength) :: description_philic_str = ' '

     ! The number of user-defined aerosol types
     integer :: ntype

     ! The number of hydrophobic and hydrophilic types read from the
     ! aerosol optics file
     integer :: n_type_phobic = 0
     integer :: n_type_philic = 0

     ! Number of relative humidity bins
     integer :: nrh = 0

     ! Number of longwave and shortwave bands of the data in the file,
     ! and monochromatic wavelengths
     integer :: n_bands_lw = 0, n_bands_sw = 0, n_mono_wl = 0

     ! Do we have any hydrophilic types?
     logical :: use_hydrophilic = .true.

     ! Do we have monochromatic optical properties
     logical :: use_monochromatic = .false.

   contains
     
     
     
     
     
     
     
     
     
     

     
     
     
     

  procedure :: setup_GPU => setup_aerosol_optics_data_GPU

  procedure :: save_GPU  => save_aerosol_optics_GPU

  procedure :: allocate_GPU

  procedure :: initialize_types_GPU

  procedure :: set_hydrophobic_type_GPU

  procedure :: set_hydrophilic_type_GPU

  procedure :: set_empty_type_GPU

  procedure :: set_types_GPU

  procedure :: calc_rh_index_GPU

  procedure :: print_description_GPU

  procedure :: create_device_GPU

  procedure :: update_host_GPU

  procedure :: update_device_GPU

  procedure :: delete_device_GPU

  procedure :: setup_CPU => setup_aerosol_optics_data_CPU

  procedure :: save_CPU  => save_aerosol_optics_CPU

  procedure :: allocate_CPU

  procedure :: initialize_types_CPU

  procedure :: set_hydrophobic_type_CPU

  procedure :: set_hydrophilic_type_CPU

  procedure :: set_empty_type_CPU

  procedure :: set_types_CPU

  procedure :: calc_rh_index_CPU

  procedure :: print_description_CPU

  end type aerosol_optics_type

  private :: create_device_GPU, update_host_GPU, update_device_GPU, delete_device_GPU

contains


  !---------------------------------------------------------------------
  ! Setup aerosol optics coefficients by reading them from a file
  


  !---------------------------------------------------------------------
  ! Initialize the arrays describing the user's aerosol types
  

  !---------------------------------------------------------------------
  ! Allocate arrays for aerosol optics data type
  


  !---------------------------------------------------------------------
  ! Save aerosol optical properties in the named file
  


  !---------------------------------------------------------------------
  ! Map user type "itype" onto stored hydrophobic type "i_type_phobic"
  


  !---------------------------------------------------------------------
  ! Map user type "itype" onto stored hydrophilic type "i_type_philic"
  


  !---------------------------------------------------------------------
  ! Set a user type "itype" to be ignored in the radiation scheme
  


  !---------------------------------------------------------------------
  ! Set user types "itypes" to map onto the stored hydrophobic and
  ! hydrophilic types according to its sign and value, with a value of
  ! 0 indicating that this type is to be ignored.  Thus if itypes=(/
  ! 3, 4, -6, 0 /) then user types 1 and 2 map on to hydrophobic types
  ! 3 and 4, user type 3 maps on to hydrophilic type 6 and user type 4
  ! is ignored.
  


  !---------------------------------------------------------------------
  ! Return an index to the relative-humdity array, or zero if no
  ! hydrophilic types are present. This function does so little that
  ! it is best to remove the Dr Hook call.
  


  !---------------------------------------------------------------------
  ! Print a description of the aerosol types to nulout
  


  !---------------------------------------------------------------------
  ! Private helper function for print_description
  


  

  

  

  


  subroutine setup_aerosol_optics_data_GPU(this, file_name, iverbose, lacc)
use yomhook
use easy_netcdf
use radiation_io
class(aerosol_optics_type), intent(inout) :: this
character(len=*), intent(in)              :: file_name
integer, intent(in), optional             :: iverbose





logical, intent (in) :: lacc



































end subroutine setup_aerosol_optics_data_GPU

  subroutine initialize_types_GPU(this, ntype, lacc)
class(aerosol_optics_type), intent(inout) :: this
integer,                    intent(in)    :: ntype
logical, intent (in) :: lacc






end subroutine initialize_types_GPU

  subroutine allocate_GPU(this, n_type_phobic, n_type_philic, nrh, &
&              n_bands_lw, n_bands_sw, n_mono_wl, lacc)
use yomhook
class(aerosol_optics_type), intent(inout) :: this
integer, intent(in) :: n_type_phobic, n_type_philic, nrh
integer, intent(in) :: n_bands_lw, n_bands_sw, n_mono_wl

logical, intent (in) :: lacc













end subroutine allocate_GPU

  subroutine save_aerosol_optics_GPU(this, file_name, iverbose, lacc)
use yomhook
use easy_netcdf
class(aerosol_optics_type), intent(inout) :: this
character(len=*),           intent(in)    :: file_name
integer,          optional, intent(in)    :: iverbose



logical, intent (in) :: lacc












































end subroutine save_aerosol_optics_GPU

  subroutine set_hydrophobic_type_GPU(this, itype, i_type_phobic, lacc)
use yomhook
class(aerosol_optics_type), intent(inout) :: this
integer, intent(in)                       :: itype, i_type_phobic

logical, intent (in) :: lacc






end subroutine set_hydrophobic_type_GPU

  subroutine set_hydrophilic_type_GPU(this, itype, i_type_philic, lacc)
use yomhook
class(aerosol_optics_type), intent(inout) :: this
integer, intent(in)                       :: itype, i_type_philic

logical, intent (in) :: lacc







end subroutine set_hydrophilic_type_GPU

  subroutine set_empty_type_GPU(this, itype, lacc)
use yomhook
class(aerosol_optics_type), intent(inout) :: this
integer, intent(in)                       :: itype

logical, intent (in) :: lacc




end subroutine set_empty_type_GPU

  subroutine set_types_GPU(this, itypes, lacc)
use yomhook
class(aerosol_optics_type), intent(inout) :: this
integer, dimension(:), intent(in)         :: itypes



logical, intent (in) :: lacc





end subroutine set_types_GPU

  function calc_rh_index_GPU(this, rh)


class(aerosol_optics_type), intent(in)    :: this
real(jprb),                 intent(in)    :: rh
integer                                   :: calc_rh_index_GPU




end function calc_rh_index_GPU

  subroutine print_description_GPU(this, i_type_map, lacc)
use radiation_io
class(aerosol_optics_type), intent(in) :: this
integer,                    intent(in) :: i_type_map(:)

logical, intent (in) :: lacc


end subroutine print_description_GPU

  pure function get_line_GPU(str,iline, lacc) result(line_str)
character(len=*), intent(in)  :: str
integer,          intent(in)  :: iline
character(len=NMaxLineLength) :: line_str


logical, intent (in) :: lacc








end function get_line_GPU

  subroutine create_device_GPU(this, lacc)
class(aerosol_optics_type), intent(inout) :: this
logical, intent (in) :: lacc




























end subroutine create_device_GPU

  subroutine update_host_GPU(this, lacc)
class(aerosol_optics_type), intent(inout) :: this
logical, intent (in) :: lacc




























end subroutine update_host_GPU

  subroutine update_device_GPU(this, lacc)
class(aerosol_optics_type), intent(inout) :: this
logical, intent (in) :: lacc




























end subroutine update_device_GPU

  subroutine delete_device_GPU(this, lacc)
class(aerosol_optics_type), intent(inout) :: this
logical, intent (in) :: lacc




























end subroutine delete_device_GPU

  subroutine setup_aerosol_optics_data_CPU(this, file_name, iverbose)
use yomhook
use easy_netcdf
use radiation_io
class(aerosol_optics_type), intent(inout) :: this
character(len=*), intent(in)              :: file_name
integer, intent(in), optional             :: iverbose








































end subroutine setup_aerosol_optics_data_CPU

  subroutine initialize_types_CPU(this, ntype)
class(aerosol_optics_type), intent(inout) :: this
integer,                    intent(in)    :: ntype






end subroutine initialize_types_CPU

  subroutine allocate_CPU(this, n_type_phobic, n_type_philic, nrh, &
&              n_bands_lw, n_bands_sw, n_mono_wl)
use yomhook
class(aerosol_optics_type), intent(inout) :: this
integer, intent(in) :: n_type_phobic, n_type_philic, nrh
integer, intent(in) :: n_bands_lw, n_bands_sw, n_mono_wl














end subroutine allocate_CPU

  subroutine save_aerosol_optics_CPU(this, file_name, iverbose)
use yomhook
use easy_netcdf
class(aerosol_optics_type), intent(inout) :: this
character(len=*),           intent(in)    :: file_name
integer,          optional, intent(in)    :: iverbose















































end subroutine save_aerosol_optics_CPU

  subroutine set_hydrophobic_type_CPU(this, itype, i_type_phobic)
use yomhook
class(aerosol_optics_type), intent(inout) :: this
integer, intent(in)                       :: itype, i_type_phobic







end subroutine set_hydrophobic_type_CPU

  subroutine set_hydrophilic_type_CPU(this, itype, i_type_philic)
use yomhook
class(aerosol_optics_type), intent(inout) :: this
integer, intent(in)                       :: itype, i_type_philic








end subroutine set_hydrophilic_type_CPU

  subroutine set_empty_type_CPU(this, itype)
use yomhook
class(aerosol_optics_type), intent(inout) :: this
integer, intent(in)                       :: itype





end subroutine set_empty_type_CPU

  subroutine set_types_CPU(this, itypes)
use yomhook
class(aerosol_optics_type), intent(inout) :: this
integer, dimension(:), intent(in)         :: itypes








end subroutine set_types_CPU

  function calc_rh_index_CPU(this, rh)

class(aerosol_optics_type), intent(in)    :: this
real(jprb),                 intent(in)    :: rh
integer                                   :: calc_rh_index_CPU




end function calc_rh_index_CPU

  subroutine print_description_CPU(this, i_type_map)
use radiation_io
class(aerosol_optics_type), intent(in) :: this
integer,                    intent(in) :: i_type_map(:)



end subroutine print_description_CPU

  pure function get_line_CPU(str,iline) result(line_str)
character(len=*), intent(in)  :: str
integer,          intent(in)  :: iline
character(len=NMaxLineLength) :: line_str










end function get_line_CPU

end module radiation_aerosol_optics_data

