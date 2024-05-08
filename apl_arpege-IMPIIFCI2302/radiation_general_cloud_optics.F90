! radiation_general_cloud_optics.F90 - Computing generalized cloud optical properties
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

module radiation_general_cloud_optics

  implicit none

  
  
contains

  ! Provides elemental function "delta_eddington_scat_od"
#include "radiation_delta_eddington.h"


  !---------------------------------------------------------------------
  ! Load cloud scattering data; this subroutine delegates to one
  ! in radiation_general_cloud_optics_data.F90
  

  !---------------------------------------------------------------------
  ! Compute cloud optical properties
  


  !---------------------------------------------------------------------
  ! Save all the cloud optics look-up tables for sw/lw and for each
  ! hydrometeor type
  

  subroutine setup_general_cloud_optics_GPU(config, lacc)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_spectral_definition
type(config_type), intent(inout) :: config

 


logical, intent (in) :: lacc














end subroutine setup_general_cloud_optics_GPU

  subroutine general_cloud_optics_GPU(nlev,istartcol,iendcol, &
&  config, thermodynamics, cloud, &
&  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
&  od_sw_cloud, ssa_sw_cloud, g_sw_cloud, lacc)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_thermodynamics
use radiation_cloud
use radiation_constants

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in), target :: config
type(thermodynamics_type),intent(in)  :: thermodynamics
type(cloud_type),   intent(in)        :: cloud




real(jprb), dimension(config%n_bands_lw,nlev,istartcol:iendcol), intent(out) :: &
&   od_lw_cloud
real(jprb), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol), &
&   intent(out) :: ssa_lw_cloud, g_lw_cloud


real(jprb), dimension(config%n_bands_sw,nlev,istartcol:iendcol), intent(out) :: &
&   od_sw_cloud, ssa_sw_cloud, g_sw_cloud




logical, intent (in) :: lacc















end subroutine general_cloud_optics_GPU

  subroutine save_general_cloud_optics_GPU(config, file_prefix, iverbose, lacc)
use yomhook
use easy_netcdf
use radiation_config
type(config_type),  intent(in) :: config
character(len=*),   intent(in) :: file_prefix
integer,  optional, intent(in) :: iverbose


logical, intent (in) :: lacc



end subroutine save_general_cloud_optics_GPU

  subroutine setup_general_cloud_optics_CPU(config)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_spectral_definition
type(config_type), intent(inout) :: config

 
















end subroutine setup_general_cloud_optics_CPU

  subroutine general_cloud_optics_CPU(nlev,istartcol,iendcol, &
&  config, thermodynamics, cloud, &
&  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
&  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_thermodynamics
use radiation_cloud
use radiation_constants

integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in), target :: config
type(thermodynamics_type),intent(in)  :: thermodynamics
type(cloud_type),   intent(in)        :: cloud




real(jprb), dimension(config%n_bands_lw,nlev,istartcol:iendcol), intent(out) :: &
&   od_lw_cloud
real(jprb), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol), &
&   intent(out) :: ssa_lw_cloud, g_lw_cloud


real(jprb), dimension(config%n_bands_sw,nlev,istartcol:iendcol), intent(out) :: &
&   od_sw_cloud, ssa_sw_cloud, g_sw_cloud



















end subroutine general_cloud_optics_CPU

  subroutine save_general_cloud_optics_CPU(config, file_prefix, iverbose)
use yomhook
use easy_netcdf
use radiation_config
type(config_type),  intent(in) :: config
character(len=*),   intent(in) :: file_prefix
integer,  optional, intent(in) :: iverbose





end subroutine save_general_cloud_optics_CPU

end module radiation_general_cloud_optics

