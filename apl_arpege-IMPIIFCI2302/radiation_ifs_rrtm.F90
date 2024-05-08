! radiation_ifs_rrtm.F90 - Interface to IFS implementation of RRTM-G
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
!   2017-04-11  R. Hogan  Receive "surface" dummy argument
!   2017-09-08  R. Hogan  Reverted some changes
!   2017-10-18  R. Hogan  Added planck_function public function
!   2018-01-11  R. Hogan  Added optional spectral scaling of incoming solar radiation
!   2018-02-22  R. Hogan  Optimized reverse indexing of heights
!   2018-05-05  R. Hogan  gas_optics can be called for reduced number of levels
!   2019-01-02  R. Hogan  Initialize shortwave props to zero in case sun below horizon

module radiation_ifs_rrtm

  implicit none

  

contains

  !---------------------------------------------------------------------
  ! Setup the IFS implementation of RRTM-G gas absorption model
  


  !---------------------------------------------------------------------
  ! Scale gas mixing ratios according to required units
  


  !---------------------------------------------------------------------
  ! Compute gas optical depths, shortwave scattering, Planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  


  !---------------------------------------------------------------------
  ! Compute Planck function of the atmosphere
  


  !---------------------------------------------------------------------
  ! Compute Planck function of the surface
  


  !---------------------------------------------------------------------
  ! Externally facing function for computing the Planck function
  ! without reference to any gas profile; typically this would be used
  ! for computing the emission by facets of a complex surface.  Note
  ! that this uses fixed "PFRAC" values, obtained by averaging over
  ! those derived from RRTM-G for near-surface conditions over a line
  ! of meridian from the ECMWF model.
  

  subroutine setup_gas_optics_rrtm_GPU(config, directory, lacc)
use yoerrtm,   only : jpglw
use yoesrtm,   only : jpgsw
use yoerrtftr, only : ngb_lw => ngb
use yoesrtm,   only : ngb_sw => ngbsw
use yomhook
use radiation_config
use radiation_spectral_definition
type(config_type), intent(inout), target :: config
character(len=*), intent(in)     :: directory
logical, intent (in) :: lacc
end subroutine setup_gas_optics_rrtm_GPU

  subroutine set_gas_units_rrtm_GPU(gas, lacc)
use radiation_gas
type(gas_type),    intent(inout) :: gas
logical, optional, intent(in)    :: lacc
end subroutine set_gas_units_rrtm_GPU

  subroutine gas_optics_rrtm_GPU(ncol,nlev,istartcol,iendcol, &
&  config, single_level, thermodynamics, gas, &
&  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
&  incoming_sw, lacc)
use parkind1
use radiation_io
USE PARRRTM
USE YOERRTM  , ONLY : JPGPT_LW => JPGPT
USE YOESRTM  , ONLY : JPGPT_SW => JPGPT
use yomhook
use radiation_config
use radiation_thermodynamics
use radiation_single_level
use radiation_gas
integer, intent(in) :: ncol
integer, intent(in) :: nlev
integer, intent(in) :: istartcol, iendcol
type(config_type), intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(gas_type),           intent(in) :: gas
real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&  intent(in), optional :: lw_albedo
real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
&   od_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(inout) :: &
&   od_sw, ssa_sw
real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), &
&   intent(out), optional :: planck_hl
real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&   intent(out), optional :: lw_emission
real(jprb), dimension(config%n_g_sw,istartcol:iendcol), &
&   intent(out), optional :: incoming_sw
logical, intent (in) :: lacc
end subroutine gas_optics_rrtm_GPU

  subroutine planck_function_atmos_GPU(nlev,istartcol,iendcol, &
config, thermodynamics, PFRAC, &
planck_hl, lacc)
use parkind1
USE YOERRTM  , ONLY : JPGPT_LW => JPGPT
use yoerrtwn, only : totplnk, delwave
use yomhook
use radiation_config
use radiation_thermodynamics
integer, intent(in) :: nlev
integer, intent(in) :: istartcol, iendcol
type(config_type), intent(in) :: config
type(thermodynamics_type),intent(in) :: thermodynamics
real(jprb), intent(in) :: PFRAC(istartcol:iendcol,JPGPT_LW,nlev)
real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out) :: &
&   planck_hl
logical, intent (in) :: lacc
end subroutine planck_function_atmos_GPU

  subroutine planck_function_surf_GPU(istartcol, iendcol, config, temperature, PFRAC, &
&  planck_surf, lacc)
use parkind1
USE YOERRTM  , ONLY : JPGPT_LW => JPGPT
use yoerrtwn, only : totplnk, delwave
use yomhook
use radiation_config
integer, intent(in) :: istartcol, iendcol
type(config_type), intent(in) :: config
real(jprb), intent(in) :: temperature(:)
real(jprb), intent(in) :: PFRAC(istartcol:iendcol,JPGPT_LW)
real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&  intent(out) :: planck_surf
logical, intent (in) :: lacc
end subroutine planck_function_surf_GPU

  subroutine planck_function_GPU(config, temperature, planck_surf, lacc)
use parkind1
use radiation_config
type(config_type), intent(in) :: config
real(jprb), intent(in) :: temperature
real(jprb), dimension(config%n_g_lw), &
&  intent(out) :: planck_surf
logical, intent (in) :: lacc
end subroutine planck_function_GPU

  subroutine setup_gas_optics_rrtm_CPU(config, directory)
use yoerrtm,   only : jpglw
use yoesrtm,   only : jpgsw
use yoerrtftr, only : ngb_lw => ngb
use yoesrtm,   only : ngb_sw => ngbsw
use yomhook
use radiation_config
use radiation_spectral_definition
type(config_type), intent(inout), target :: config
character(len=*), intent(in)     :: directory
end subroutine setup_gas_optics_rrtm_CPU

  subroutine set_gas_units_rrtm_CPU(gas, lacc)
use radiation_gas
type(gas_type),    intent(inout) :: gas
logical, optional, intent(in)    :: lacc
end subroutine set_gas_units_rrtm_CPU

  subroutine gas_optics_rrtm_CPU(ncol,nlev,istartcol,iendcol, &
&  config, single_level, thermodynamics, gas, &
&  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
&  incoming_sw)
use parkind1
USE PARRRTM
USE YOERRTM  , ONLY : JPGPT_LW => JPGPT
USE YOESRTM  , ONLY : JPGPT_SW => JPGPT
use yomhook
use radiation_config
use radiation_thermodynamics
use radiation_single_level
use radiation_gas
integer, intent(in) :: ncol
integer, intent(in) :: nlev
integer, intent(in) :: istartcol, iendcol
type(config_type), intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(gas_type),           intent(in) :: gas
real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&  intent(in), optional :: lw_albedo
real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
&   od_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(inout) :: &
&   od_sw, ssa_sw
real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), &
&   intent(out), optional :: planck_hl
real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&   intent(out), optional :: lw_emission
real(jprb), dimension(config%n_g_sw,istartcol:iendcol), &
&   intent(out), optional :: incoming_sw
end subroutine gas_optics_rrtm_CPU

  subroutine planck_function_atmos_CPU(nlev,istartcol,iendcol, &
config, thermodynamics, PFRAC, &
planck_hl)
use parkind1
USE YOERRTM  , ONLY : JPGPT_LW => JPGPT
use yoerrtwn, only : totplnk, delwave
use yomhook
use radiation_config
use radiation_thermodynamics
integer, intent(in) :: nlev
integer, intent(in) :: istartcol, iendcol
type(config_type), intent(in) :: config
type(thermodynamics_type),intent(in) :: thermodynamics
real(jprb), intent(in) :: PFRAC(istartcol:iendcol,JPGPT_LW,nlev)
real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out) :: &
&   planck_hl
end subroutine planck_function_atmos_CPU

  subroutine planck_function_surf_CPU(istartcol, iendcol, config, temperature, PFRAC, &
&  planck_surf)
use parkind1
USE YOERRTM  , ONLY : JPGPT_LW => JPGPT
use yoerrtwn, only : totplnk, delwave
use yomhook
use radiation_config
integer, intent(in) :: istartcol, iendcol
type(config_type), intent(in) :: config
real(jprb), intent(in) :: temperature(:)
real(jprb), intent(in) :: PFRAC(istartcol:iendcol,JPGPT_LW)
real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
&  intent(out) :: planck_surf
end subroutine planck_function_surf_CPU

  subroutine planck_function_CPU(config, temperature, planck_surf)
use parkind1
use radiation_config
type(config_type), intent(in) :: config
real(jprb), intent(in) :: temperature
real(jprb), dimension(config%n_g_lw), &
&  intent(out) :: planck_surf
end subroutine planck_function_CPU

end module radiation_ifs_rrtm

