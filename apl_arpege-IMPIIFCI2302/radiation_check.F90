! radiation_check.F90 - Checking routines
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

module radiation_check

  use parkind1, only : jprb

  implicit none
  public

contains

  !---------------------------------------------------------------------
  ! Return .true. if 1D allocatable array "var" is out of physical
  ! range specified by boundmin and boundmax, and issue a warning.
  ! "do_fix" determines whether erroneous values are fixed to lie
  ! within the physical range. To check only a subset of the array,
  ! specify i1 and i2 for the range.
  function out_of_bounds_1d(var, var_name, boundmin, boundmax, do_fix, i1, i2) result (is_bad)
use radiation_io,     only : nulout
real(jprb), allocatable, intent(inout) :: var(:)
character(len=*),        intent(in) :: var_name
real(jprb),              intent(in) :: boundmin, boundmax
logical,                 intent(in) :: do_fix
integer,       optional, intent(in) :: i1, i2
logical                       :: is_bad



end function out_of_bounds_1d


  !---------------------------------------------------------------------
  ! Return .true. if 2D allocatable array "var" is out of physical
  ! range specified by boundmin and boundmax, and issue a warning.  To
  ! check only a subset of the array, specify i1 and i2 for the range
  ! of the first dimension and j1 and j2 for the range of the second.
  function out_of_bounds_2d(var, var_name, boundmin, boundmax, do_fix, &
&                    i1, i2, j1, j2) result (is_bad)
use radiation_io,     only : nulout
real(jprb), allocatable, intent(inout) :: var(:,:)
character(len=*),        intent(in) :: var_name
real(jprb),              intent(in) :: boundmin, boundmax
logical,                 intent(in) :: do_fix
integer,       optional, intent(in) :: i1, i2, j1, j2


logical                       :: is_bad



end function out_of_bounds_2d


  !---------------------------------------------------------------------
  ! Return .true. if 3D allocatable array "var" is out of physical
  ! range specified by boundmin and boundmax, and issue a warning.  To
  ! check only a subset of the array, specify i1 and i2 for the range
  ! of the first dimension, j1 and j2 for the second and k1 and k2 for
  ! the third.
  function out_of_bounds_3d(var, var_name, boundmin, boundmax, do_fix, &
&                    i1, i2, j1, j2, k1, k2) result (is_bad)
use radiation_io,     only : nulout
real(jprb), allocatable, intent(inout) :: var(:,:,:)
character(len=*),        intent(in) :: var_name
real(jprb),              intent(in) :: boundmin, boundmax
logical,                 intent(in) :: do_fix
integer,       optional, intent(in) :: i1, i2, j1, j2, k1, k2


logical                       :: is_bad



end function out_of_bounds_3d

end module radiation_check
