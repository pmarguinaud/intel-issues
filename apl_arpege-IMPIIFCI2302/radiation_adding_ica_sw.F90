! radiation_adding_ica_sw.F90 - Shortwave adding method in independent column approximation
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

module radiation_adding_ica_sw

  

contains

  

  subroutine adding_ica_sw_GPU(ncol, nlev, incoming_toa, &
&  albedo_surf_diffuse, albedo_surf_direct, cos_sza, &
&  reflectance, transmittance, ref_dir, trans_dir_diff, trans_dir_dir, &
&  flux_up, flux_dn_diffuse, flux_dn_direct, &
&  albedo, source, inv_denominator)
use parkind1
use yomhook
implicit none

integer, intent(in) :: ncol 
integer, intent(in) :: nlev 

real(jprb), intent(in),  dimension(ncol)         :: incoming_toa

real(jprb), intent(in),  dimension(ncol)         :: albedo_surf_diffuse, &
&                                              albedo_surf_direct

real(jprb), intent(in)                           :: cos_sza

real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance



real(jprb), intent(in),  dimension(ncol, nlev)   :: ref_dir, trans_dir_diff


real(jprb), intent(in),  dimension(ncol, nlev)   :: trans_dir_dir


real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn_diffuse, &
&                                              flux_dn_direct


real(jprb), intent(inout), dimension(ncol, nlev+1) :: albedo


real(jprb), intent(inout), dimension(ncol, nlev+1) :: source

real(jprb), intent(inout), dimension(ncol, nlev)   :: inv_denominator

































end subroutine adding_ica_sw_GPU

  subroutine adding_ica_sw_CPU(ncol, nlev, incoming_toa, &
&  albedo_surf_diffuse, albedo_surf_direct, cos_sza, &
&  reflectance, transmittance, ref_dir, trans_dir_diff, trans_dir_dir, &
&  flux_up, flux_dn_diffuse, flux_dn_direct, &
&  albedo, source, inv_denominator)
use parkind1
use yomhook
implicit none

integer, intent(in) :: ncol 
integer, intent(in) :: nlev 

real(jprb), intent(in),  dimension(ncol)         :: incoming_toa

real(jprb), intent(in),  dimension(ncol)         :: albedo_surf_diffuse, &
&                                              albedo_surf_direct

real(jprb), intent(in)                           :: cos_sza

real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance



real(jprb), intent(in),  dimension(ncol, nlev)   :: ref_dir, trans_dir_diff


real(jprb), intent(in),  dimension(ncol, nlev)   :: trans_dir_dir


real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn_diffuse, &
&                                              flux_dn_direct


real(jprb), intent(inout), dimension(ncol, nlev+1) :: albedo


real(jprb), intent(inout), dimension(ncol, nlev+1) :: source

real(jprb), intent(inout), dimension(ncol, nlev)   :: inv_denominator































end subroutine adding_ica_sw_CPU

end module radiation_adding_ica_sw

