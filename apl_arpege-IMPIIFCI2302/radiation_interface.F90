! radiation_interface.F90 - Public interface to radiation scheme
!
! (C) Copyright 2014- ECMWF.
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
!   2017-04-11  R. Hogan  Changes to enable generalized surface description
!   2017-09-08  R. Hogan  Reverted some changes
!   2022-01-18  P. Ukkonen Added support for RRTMGP gas optics
!
! To use the radiation scheme, create a configuration_type object,
! call "setup_radiation" on it once to load the look-up-tables and
! data describing how gas and hydrometeor absorption/scattering are to
! be represented, and call "radiation" multiple times on different
! input profiles.

module radiation_interface

  use radiation_interface_GPU
  use radiation_interface_CPU

  implicit none

end module radiation_interface

