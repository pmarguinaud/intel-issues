! radiation_lw_derivatives.F90 - Compute longwave derivatives for Hogan and Bozzo (2015) method
!
! (C) Copyright 2016- ECMWF.
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
! This module provides routines to compute the rate of change of
! broadband upwelling longwave flux at each half level with respect to
! the surface broadband upwelling flux.  This is done from the surface
! spectral fluxes and the spectral transmittance of each atmospheric
! layer, assuming no longwave scattering. The result may be used to
! perform approximate updates to the longwave flux profile in between
! calls to the full radiation scheme, accounting for the change in
! skin temperature, following the method of Hogan and Bozzo (JAMES
! 2015).  Separate routines are provided for each solver.
!
! Note that currently a more approximate calculation is performed from
! the exact one in Hogan and Bozzo (2015); here we assume that a
! change in temperature increases the spectral fluxes in proportion,
! when in reality there is a change in shape of the Planck function in
! addition to an overall increase in the total emission.
!
! Modifications
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2022-11-22  P. Ukkonen / R. Hogan  Optimized calc_lw_derivatives_region

module radiation_lw_derivatives

  

! Allow size of inner dimension (number of g-points) to be known at compile time if NG_LW is defined

  integer, parameter, private :: nreg = 3


contains

  !---------------------------------------------------------------------
  ! Calculation for the Independent Column Approximation
  


  !---------------------------------------------------------------------
  ! Calculation for the Independent Column Approximation
  



  !---------------------------------------------------------------------
  ! Calculation for solvers involving multiple regions and matrices
  


  !---------------------------------------------------------------------
  ! Calculation for solvers involving multiple regions but no 3D
  ! effects: the difference from calc_lw_derivatives_matrix is that transmittance
  ! has one fewer dimensions
  


  subroutine calc_lw_derivatives_ica_GPU(ng_lw_in, nlev, icol, transmittance, flux_up_surf, lw_derivatives)
use parkind1
implicit none

integer,    intent(in) :: ng_lw_in   
integer,    intent(in) :: nlev 
integer,    intent(in) :: icol 
real(jprb), intent(in) :: transmittance(ng_lw_in,nlev)
real(jprb), intent(in) :: flux_up_surf(ng_lw_in) 

real(jprb), intent(out) :: lw_derivatives(:,:) 

















end subroutine calc_lw_derivatives_ica_GPU

  subroutine modify_lw_derivatives_ica_GPU(ng_lw_in, nlev, icol, transmittance, &
&                               flux_up_surf, weight, lw_derivatives)
use parkind1
use yomhook
implicit none

integer,    intent(in) :: ng_lw_in   
integer,    intent(in) :: nlev 
integer,    intent(in) :: icol 
real(jprb), intent(in) :: transmittance(ng_lw_in,nlev)
real(jprb), intent(in) :: flux_up_surf(ng_lw_in) 
real(jprb), intent(in) :: weight 

real(jprb), intent(inout) :: lw_derivatives(:,:) 


















end subroutine modify_lw_derivatives_ica_GPU

  subroutine calc_lw_derivatives_matrix_GPU(ng_lw_in, nlev, nreg_in, icol, transmittance, &
&                                u_matrix, flux_up_surf, lw_derivatives, lacc)
use parkind1
use yomhook
use radiation_matrix
implicit none

integer,    intent(in) :: ng_lw_in   
integer,    intent(in) :: nlev 
integer,    intent(in) :: nreg_in 
integer,    intent(in) :: icol 
real(jprb), intent(in) :: transmittance(ng_lw_in,nreg,nreg,nlev)
real(jprb), intent(in) :: u_matrix(nreg,nreg,nlev+1) 
real(jprb), intent(in) :: flux_up_surf(ng_lw_in) 

real(jprb), intent(out) :: lw_derivatives(:,:) 







logical, intent (in) :: lacc











end subroutine calc_lw_derivatives_matrix_GPU

  subroutine calc_lw_derivatives_region_GPU(ng_lw_in, nlev, nreg_in, icol, transmittance, &
&                                u_matrix, flux_up_surf, lw_derivatives, lacc)
use parkind1
use yomhook
use radiation_matrix
implicit none

integer,    intent(in) :: ng_lw_in   
integer,    intent(in) :: nlev 
integer,    intent(in) :: nreg_in 
integer,    intent(in) :: icol 
real(jprb), intent(in) :: transmittance(ng_lw_in,nreg,nlev)
real(jprb), intent(in) :: u_matrix(nreg,nreg,nlev+1) 
real(jprb), intent(in) :: flux_up_surf(ng_lw_in) 

real(jprb), intent(out) :: lw_derivatives(:,:) 






logical, intent (in) :: lacc









end subroutine calc_lw_derivatives_region_GPU

  subroutine calc_lw_derivatives_ica_CPU(ng_lw_in, nlev, icol, transmittance, flux_up_surf, lw_derivatives)
use parkind1
use yomhook
implicit none

integer,    intent(in) :: ng_lw_in   
integer,    intent(in) :: nlev 
integer,    intent(in) :: icol 
real(jprb), intent(in) :: transmittance(ng_lw_in,nlev)
real(jprb), intent(in) :: flux_up_surf(ng_lw_in) 

real(jprb), intent(out) :: lw_derivatives(:,:) 
















end subroutine calc_lw_derivatives_ica_CPU

  subroutine modify_lw_derivatives_ica_CPU(ng_lw_in, nlev, icol, transmittance, &
&                               flux_up_surf, weight, lw_derivatives)
use parkind1
use yomhook
implicit none

integer,    intent(in) :: ng_lw_in   
integer,    intent(in) :: nlev 
integer,    intent(in) :: icol 
real(jprb), intent(in) :: transmittance(ng_lw_in,nlev)
real(jprb), intent(in) :: flux_up_surf(ng_lw_in) 
real(jprb), intent(in) :: weight 

real(jprb), intent(inout) :: lw_derivatives(:,:) 

















end subroutine modify_lw_derivatives_ica_CPU

  subroutine calc_lw_derivatives_matrix_CPU(ng_lw_in, nlev, nreg_in, icol, transmittance, &
&                                u_matrix, flux_up_surf, lw_derivatives)
use parkind1
use yomhook
use radiation_matrix
implicit none

integer,    intent(in) :: ng_lw_in   
integer,    intent(in) :: nlev 
integer,    intent(in) :: nreg_in 
integer,    intent(in) :: icol 
real(jprb), intent(in) :: transmittance(ng_lw_in,nreg,nreg,nlev)
real(jprb), intent(in) :: u_matrix(nreg,nreg,nlev+1) 
real(jprb), intent(in) :: flux_up_surf(ng_lw_in) 

real(jprb), intent(out) :: lw_derivatives(:,:) 


















end subroutine calc_lw_derivatives_matrix_CPU

  subroutine calc_lw_derivatives_region_CPU(ng_lw_in, nlev, nreg_in, icol, transmittance, &
&                                u_matrix, flux_up_surf, lw_derivatives)
use parkind1
use yomhook
use radiation_matrix
implicit none

integer,    intent(in) :: ng_lw_in   
integer,    intent(in) :: nlev 
integer,    intent(in) :: nreg_in 
integer,    intent(in) :: icol 
real(jprb), intent(in) :: transmittance(ng_lw_in,nreg,nlev)
real(jprb), intent(in) :: u_matrix(nreg,nreg,nlev+1) 
real(jprb), intent(in) :: flux_up_surf(ng_lw_in) 

real(jprb), intent(out) :: lw_derivatives(:,:) 















end subroutine calc_lw_derivatives_region_CPU

end module radiation_lw_derivatives

