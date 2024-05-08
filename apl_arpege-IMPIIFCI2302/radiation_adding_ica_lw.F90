! radiation_adding_ica_lw.F90 - Longwave adding method in independent column approximation
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
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-07-12  R. Hogan  Fast adding method for if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_adding_ica_lw

  

contains

  !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering, by successively adding the contribution of
  ! layers starting from the surface to compute the total albedo and
  ! total upward emission of the increasingly larger block of
  ! atmospheric layers.
  


  !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering in cloudy layers only.
  


  !---------------------------------------------------------------------
  ! If there is no scattering then fluxes may be computed simply by
  ! passing down through the atmosphere computing the downwelling
  ! fluxes from the transmission and emission of each layer, and then
  ! passing back up through the atmosphere to compute the upwelling
  ! fluxes in the same way.
  

  subroutine adding_ica_lw_GPU(ncol, nlev, &
&  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
&  flux_up, flux_dn, lacc)
use parkind1
use yomhook
implicit none

integer, intent(in) :: ncol 
integer, intent(in) :: nlev 

real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn


real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn











logical, intent (in) :: lacc


















end subroutine adding_ica_lw_GPU

  subroutine fast_adding_ica_lw_GPU(ncol, nlev, &
&  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
&  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
&  flux_up, flux_dn, albedo, source, inv_denominator)
use parkind1
use yomhook
implicit none

integer, intent(in) :: ncol 
integer, intent(in) :: nlev 

real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

logical, intent(in) :: is_clear_sky_layer(nlev)

integer, intent(in) :: i_cloud_top

real(jprb), intent(in), dimension(ncol, nlev+1)  :: flux_dn_clear


real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn


real(jprb), intent(inout), dimension(ncol, nlev+1) :: albedo


real(jprb), intent(inout), dimension(ncol, nlev+1) :: source

real(jprb), intent(inout), dimension(ncol, nlev)   :: inv_denominator
























end subroutine fast_adding_ica_lw_GPU

  subroutine calc_fluxes_no_scattering_lw_GPU(ncol, nlev, &
&  transmittance, source_up, source_dn, emission_surf, albedo_surf, flux_up, flux_dn)
use parkind1
implicit none

integer, intent(in) :: ncol 
integer, intent(in) :: nlev 

real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

real(jprb), intent(in),  dimension(ncol, nlev)   :: transmittance

real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn


real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn





















end subroutine calc_fluxes_no_scattering_lw_GPU

  subroutine adding_ica_lw_CPU(ncol, nlev, &
&  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
&  flux_up, flux_dn)
use parkind1
use yomhook
implicit none

integer, intent(in) :: ncol 
integer, intent(in) :: nlev 

real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn


real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn





























end subroutine adding_ica_lw_CPU

  subroutine fast_adding_ica_lw_CPU(ncol, nlev, &
&  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
&  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
&  flux_up, flux_dn, albedo, source, inv_denominator)
use parkind1
use yomhook
implicit none

integer, intent(in) :: ncol 
integer, intent(in) :: nlev 

real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

logical, intent(in) :: is_clear_sky_layer(nlev)

integer, intent(in) :: i_cloud_top

real(jprb), intent(in), dimension(ncol, nlev+1)  :: flux_dn_clear


real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn


real(jprb), intent(inout), dimension(ncol, nlev+1) :: albedo


real(jprb), intent(inout), dimension(ncol, nlev+1) :: source

real(jprb), intent(inout), dimension(ncol, nlev)   :: inv_denominator





















end subroutine fast_adding_ica_lw_CPU

  subroutine calc_fluxes_no_scattering_lw_CPU(ncol, nlev, &
&  transmittance, source_up, source_dn, emission_surf, albedo_surf, flux_up, flux_dn)
use parkind1
use yomhook
implicit none

integer, intent(in) :: ncol 
integer, intent(in) :: nlev 

real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

real(jprb), intent(in),  dimension(ncol, nlev)   :: transmittance

real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn


real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn



















end subroutine calc_fluxes_no_scattering_lw_CPU

end module radiation_adding_ica_lw

