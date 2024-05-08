! radiation_ifs_rrtmgp.F90 - Interface to IFS implementation of RRTMGP-NN,
! which can use either RRTMGP k-distributions loaded from a netCDF file,
! or faster NN emulators of these distributions
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
! Author:  Peter U kkonen
! Email:   peterukk@gmail.com
!
! Modifications
!   2020-01-30  P. Ukkonen  Initial version
!   2022-02-01  P. Ukkonen  Reworked to incur minimal changes to existing ecRAD code

module radiation_ifs_rrtmgp

  implicit none
  
  

  ! List of character for case conversion, need this for changing an entry
  ! in a string array
  character(len=26), parameter :: LOWER_CASE_CHARS = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter :: UPPER_CASE_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

contains

  !---------------------------------------------------------------------
  ! Setup the IFS implementation of RRTMGP gas absorption model
  

  !---------------------------------------------------------------------
  ! Scale gas mixing ratios according to required units
  


  !---------------------------------------------------------------------
  ! Compute gas optical depths, shortwave scattering, Planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  


  

  


  


    subroutine setup_gas_optics_ifs_rrtmgp_GPU(config, lacc)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_spectral_definition
use mo_load_coefficients
use mo_gas_concentrations
use mo_gas_optics_rrtmgp
type(config_type),      intent(inout), target    :: config
 





logical, intent (in) :: lacc




















end subroutine setup_gas_optics_ifs_rrtmgp_GPU

  subroutine set_gas_units_rrtmgp_GPU(gas, lacc)
use radiation_gas
type(gas_type),    intent(inout) :: gas
logical, intent (in) :: lacc

end subroutine set_gas_units_rrtmgp_GPU

  subroutine gas_optics_ifs_rrtmgp_GPU(ncol,nlev,istartcol,iendcol, &
&  config, single_level, thermodynamics, gas, &
&  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
&  incoming_sw, lacc)
use parkind1
use yomhook
use radiation_config
use radiation_thermodynamics
use radiation_single_level
use radiation_gas
use radiation_io
use mo_gas_optics_rrtmgp
use mo_gas_concentrations
use radiation_gas_constants
use mo_rte_util_array
integer, intent(in) :: ncol               
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(gas_type),           intent(in) :: gas

real(jprb), dimension(config%n_g_lw,istartcol:iendcol), intent(in)  :: lw_albedo



real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: od_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: od_sw, ssa_sw


real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out)   :: planck_hl

real(jprb), dimension(config%n_g_lw,istartcol:iendcol),      intent(out)     :: lw_emission



real(jprb), dimension(config%n_g_sw,istartcol:iendcol),  intent(out)         :: incoming_sw



 
 





logical, intent (in) :: lacc







































end subroutine gas_optics_ifs_rrtmgp_GPU

  pure function string_loc_in_array_GPU(string, array, lacc)
character(len=*),               intent(in) :: string
character(len=*), dimension(:), intent(in) :: array
integer                                    :: string_loc_in_array_GPU


logical, intent (in) :: lacc



end function string_loc_in_array_GPU

  pure function lower_case_GPU( input_string, lacc ) result( output_string )
character(len=*), intent(in) :: input_string
character(len=len(input_string)) :: output_string

logical, intent (in) :: lacc




end function

  subroutine stop_on_err_GPU(error_msg, lacc)
use iso_c_binding
use radiation_io
character(len=*), intent(in) :: error_msg
logical, intent (in) :: lacc

end subroutine stop_on_err_GPU

  subroutine setup_gas_optics_ifs_rrtmgp_CPU(config)
use parkind1
use yomhook
use radiation_io
use radiation_config
use radiation_spectral_definition
use mo_load_coefficients
use mo_gas_concentrations
use mo_gas_optics_rrtmgp
type(config_type),      intent(inout), target    :: config
 

























end subroutine setup_gas_optics_ifs_rrtmgp_CPU

  subroutine set_gas_units_rrtmgp_CPU(gas)
use radiation_gas
type(gas_type),    intent(inout) :: gas

end subroutine set_gas_units_rrtmgp_CPU

  subroutine gas_optics_ifs_rrtmgp_CPU(ncol,nlev,istartcol,iendcol, &
&  config, single_level, thermodynamics, gas, &
&  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
&  incoming_sw)
use parkind1
use yomhook
use radiation_config
use radiation_thermodynamics
use radiation_single_level
use radiation_gas
use radiation_io
use mo_gas_optics_rrtmgp
use mo_gas_concentrations
use radiation_gas_constants
use mo_rte_util_array
integer, intent(in) :: ncol               
integer, intent(in) :: nlev               
integer, intent(in) :: istartcol, iendcol 
type(config_type), intent(in) :: config
type(single_level_type),  intent(in) :: single_level
type(thermodynamics_type),intent(in) :: thermodynamics
type(gas_type),           intent(in) :: gas

real(jprb), dimension(config%n_g_lw,istartcol:iendcol), intent(in)  :: lw_albedo



real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: od_lw
real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: od_sw, ssa_sw


real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out)   :: planck_hl

real(jprb), dimension(config%n_g_lw,istartcol:iendcol),      intent(out)     :: lw_emission



real(jprb), dimension(config%n_g_sw,istartcol:iendcol),  intent(out)         :: incoming_sw



 
 












































end subroutine gas_optics_ifs_rrtmgp_CPU

  pure function string_loc_in_array_CPU(string, array)
character(len=*),               intent(in) :: string
character(len=*), dimension(:), intent(in) :: array
integer                                    :: string_loc_in_array_CPU





end function string_loc_in_array_CPU

  pure function lower_case_CPU( input_string ) result( output_string )
character(len=*), intent(in) :: input_string
character(len=len(input_string)) :: output_string





end function

  subroutine stop_on_err_CPU(error_msg)
use iso_c_binding
use radiation_io
character(len=*), intent(in) :: error_msg

end subroutine stop_on_err_CPU

end module radiation_ifs_rrtmgp

