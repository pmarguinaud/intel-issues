! radiation_aerosol_optics.F90 - Computing aerosol optical properties
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
!   2018-04-15  R. Hogan  Add "direct" option
!   2020-11-14  R. Hogan  Add setup_general_aerosol_optics for ecCKD compatibility
!   2022-03-27  R. Hogan  Add setup_general_aerosol_optics_legacy to use RRTM aerosol files with ecCKD
!   2022-11-22  P. Ukkonen / R. Hogan  Optimizations to enhance vectorization

#include "ecrad_config.h"

module radiation_aerosol_optics

  implicit none
  

contains

  ! Provides the elemental function "delta_eddington_extensive"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Load aerosol scattering data; this subroutine delegates to one
  ! in radiation_aerosol_optics_data.F90
  


  !---------------------------------------------------------------------
  ! Read file containing high spectral resolution optical properties
  ! and average to the spectral intervals of the current gas-optics
  ! scheme
  


  !---------------------------------------------------------------------
  ! Read file containing legacy-style band-wise aerosol optical
  ! properties and average to the spectral intervals of the current
  ! gas-optics scheme
  


  !---------------------------------------------------------------------
  ! Compute aerosol optical properties and add to existing gas optical
  ! depth and scattering properties
  


  !---------------------------------------------------------------------
  ! Add precomputed optical properties to gas optical depth and
  ! scattering properties
  


  !---------------------------------------------------------------------
  ! Sometimes it is useful to specify aerosol in terms of its optical
  ! depth at a particular wavelength.  This function returns the dry
  ! mass-extinction coefficient, i.e. the extinction cross section per
  ! unit mass, for aerosol of type "itype" at the specified wavelength
  ! (m). For hydrophilic types, the value at the first relative
  ! humidity bin is taken.
  


  !---------------------------------------------------------------------
  ! Compute aerosol extinction coefficient at a particular wavelength
  ! and a single height - this is useful for visibility diagnostics
  

  subroutine setup_aerosol_optics_GPU(config, lacc)
use parkind1
use yomhook
use radiation_config
use radiation_aerosol_optics_data
use radiation_io
type(config_type), intent(inout) :: config

logical, intent (in) :: lacc




end subroutine setup_aerosol_optics_GPU

  subroutine setup_general_aerosol_optics_GPU(config, lacc)
use parkind1
use yomhook
use easy_netcdf
use radiation_config
use radiation_aerosol_optics_data
use radiation_spectral_definition
use radiation_io
type(config_type), intent(inout), target :: config



 

    
         
           
 

      
           
             
   















logical, intent (in) :: lacc



end subroutine setup_general_aerosol_optics_GPU

  subroutine setup_general_aerosol_optics_legacy_GPU(config, file_name, lacc)
use parkind1
use yomhook
use easy_netcdf
use radiation_config
use radiation_aerosol_optics_data
use radiation_spectral_definition
type(config_type), intent(inout), target :: config

character(len=*), intent(in) :: file_name











logical, intent (in) :: lacc




















end subroutine setup_general_aerosol_optics_legacy_GPU

  subroutine add_aerosol_optics_GPU(ng_sw_in, ng_lw_in, nbnd_sw_in, nbnd_lw_in, &
&  nlev, istartcol, iendcol, &
&  config, thermodynamics, gas, aerosol, &
&  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw, lacc)
use parkind1
use radiation_io
use yomhook
use radiation_config
use radiation_thermodynamics
use radiation_gas
use radiation_aerosol
use radiation_constants
use radiation_aerosol_optics_data
use radiation_aerosol_optics_data, only : calc_rh_index_GPU


integer, intent(in) :: ng_sw_in, ng_lw_in 
integer, intent(in) :: nbnd_sw_in, nbnd_lw_in 
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in)  :: config
type(thermodynamics_type),intent(in)  :: thermodynamics
type(gas_type),           intent(in)  :: gas
type(aerosol_type),       intent(in)  :: aerosol





real(jprb), dimension(ng_lw_in,nlev,istartcol:iendcol), &
&   intent(inout) :: od_lw
real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
&   intent(out)   :: ssa_lw, g_lw
real(jprb), dimension(ng_sw_in,nlev,istartcol:iendcol), &
&   intent(inout) :: od_sw, ssa_sw
real(jprb), dimension(ng_sw_in,nlev,istartcol:iendcol), &
&   intent(inout)   :: g_sw









 

















logical, intent (in) :: lacc



end subroutine add_aerosol_optics_GPU

  subroutine add_aerosol_optics_direct_GPU(nlev,istartcol,iendcol, &
&  config, aerosol, &
&  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw, lacc)
use parkind1
use radiation_io
use yomhook
use radiation_config
use radiation_aerosol
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in), target :: config
type(aerosol_type),       intent(in)  :: aerosol





real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), &
&   intent(inout) :: od_lw
real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
&   intent(out)   :: ssa_lw, g_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
&   intent(inout) :: od_sw, ssa_sw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
&   intent(out)   :: g_sw

















logical, intent (in) :: lacc






end subroutine add_aerosol_optics_direct_GPU

  function dry_aerosol_mass_extinction_GPU(config, itype, wavelength, lacc)
use parkind1
use radiation_io
use radiation_config
use radiation_aerosol_optics_data
type(config_type), intent(in), target :: config

integer, intent(in) :: itype

real(jprb), intent(in) :: wavelength
real(jprb) :: dry_aerosol_mass_extinction_GPU




logical, intent (in) :: lacc




end function dry_aerosol_mass_extinction_GPU

  subroutine aerosol_extinction_GPU(ncol,istartcol,iendcol, &
&  config, wavelength, mixing_ratio, relative_humidity, extinction, lacc)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_aerosol_optics_data
integer, intent(in) :: ncol               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in), target :: config
real(jprb), intent(in)  :: wavelength 
real(jprb), intent(in)  :: mixing_ratio(ncol,config%n_aerosol_types)
real(jprb), intent(in)  :: relative_humidity(ncol)
real(jprb), intent(out) :: extinction(ncol)











logical, intent (in) :: lacc








end subroutine aerosol_extinction_GPU

  subroutine setup_aerosol_optics_CPU(config)
use parkind1
use yomhook
use radiation_config
use radiation_aerosol_optics_data
use radiation_io
type(config_type), intent(inout) :: config





end subroutine setup_aerosol_optics_CPU

  subroutine setup_general_aerosol_optics_CPU(config)
use parkind1
use yomhook
use easy_netcdf
use radiation_config
use radiation_aerosol_optics_data
use radiation_spectral_definition
use radiation_io
type(config_type), intent(inout), target :: config



 

    
         
           
 

      
           
             
   


















end subroutine setup_general_aerosol_optics_CPU

  subroutine setup_general_aerosol_optics_legacy_CPU(config, file_name)
use parkind1
use yomhook
use easy_netcdf
use radiation_config
use radiation_aerosol_optics_data
use radiation_spectral_definition
type(config_type), intent(inout), target :: config

character(len=*), intent(in) :: file_name































end subroutine setup_general_aerosol_optics_legacy_CPU

  subroutine add_aerosol_optics_CPU(ng_sw_in, ng_lw_in, nbnd_sw_in, nbnd_lw_in, &
&  nlev, istartcol, iendcol, &
&  config, thermodynamics, gas, aerosol, &
&  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)
use parkind1
use radiation_io
use yomhook
use radiation_config
use radiation_thermodynamics
use radiation_gas
use radiation_aerosol
use radiation_constants
use radiation_aerosol_optics_data
use radiation_aerosol_optics_data, only : calc_rh_index_CPU


integer, intent(in) :: ng_sw_in, ng_lw_in 
integer, intent(in) :: nbnd_sw_in, nbnd_lw_in 
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type),        intent(in)  :: config
type(thermodynamics_type),intent(in)  :: thermodynamics
type(gas_type),           intent(in)  :: gas
type(aerosol_type),       intent(in)  :: aerosol





real(jprb), dimension(ng_lw_in,nlev,istartcol:iendcol), &
&   intent(inout) :: od_lw
real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
&   intent(out)   :: ssa_lw, g_lw
real(jprb), dimension(ng_sw_in,nlev,istartcol:iendcol), &
&   intent(inout) :: od_sw, ssa_sw
real(jprb), dimension(ng_sw_in,nlev,istartcol:iendcol), &
&   intent(inout)   :: g_sw









 




















end subroutine add_aerosol_optics_CPU

  subroutine add_aerosol_optics_direct_CPU(nlev,istartcol,iendcol, &
&  config, aerosol, &
&  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)
use parkind1
use radiation_io
use yomhook
use radiation_config
use radiation_aerosol
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in), target :: config
type(aerosol_type),       intent(in)  :: aerosol





real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), &
&   intent(inout) :: od_lw
real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
&   intent(out)   :: ssa_lw, g_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
&   intent(inout) :: od_sw, ssa_sw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), &
&   intent(out)   :: g_sw





















end subroutine add_aerosol_optics_direct_CPU

  function dry_aerosol_mass_extinction_CPU(config, itype, wavelength)
use parkind1
use radiation_io
use radiation_config
use radiation_aerosol_optics_data
type(config_type), intent(in), target :: config

integer, intent(in) :: itype

real(jprb), intent(in) :: wavelength
real(jprb) :: dry_aerosol_mass_extinction_CPU








end function dry_aerosol_mass_extinction_CPU

  subroutine aerosol_extinction_CPU(ncol,istartcol,iendcol, &
&  config, wavelength, mixing_ratio, relative_humidity, extinction)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_aerosol_optics_data
integer, intent(in) :: ncol               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in), target :: config
real(jprb), intent(in)  :: wavelength 
real(jprb), intent(in)  :: mixing_ratio(ncol,config%n_aerosol_types)
real(jprb), intent(in)  :: relative_humidity(ncol)
real(jprb), intent(out) :: extinction(ncol)



















end subroutine aerosol_extinction_CPU

end module radiation_aerosol_optics

