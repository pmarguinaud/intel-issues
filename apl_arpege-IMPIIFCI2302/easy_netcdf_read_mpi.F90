! easy_netcdf_read_mpi.f90 - Read netcdf file on one task and share with other tasks
!
! Copyright (C) 2017- ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module easy_netcdf_read_mpi

  use easy_netcdf,   only : netcdf_file_raw => netcdf_file
  use parkind1,      only : jpim, jprb
  use radiation_io,  only : nulout, nulerr, my_abort => radiation_abort

  implicit none

  ! MPI tag for radiation and physics communication
  integer(kind=jpim), parameter :: mtagrad = 2800

  !---------------------------------------------------------------------
  ! An object of this type provides convenient read or write access to
  ! a NetCDF file
  type netcdf_file
    type(netcdf_file_raw) :: file
    logical :: is_master_task = .true.
    logical :: multi_open = .false.
    logical :: mpi_enabled    = .false.
  contains
    procedure :: open => open_netcdf_file
    procedure :: open_active => open_netcdf_file_active
    procedure :: close => close_netcdf_file
    procedure :: get_rank
    procedure :: get_real_scalar
    procedure :: get_int_scalar
    procedure :: get_real_vector
    procedure :: get_real_vector_active
    procedure :: get_int_vector
    procedure :: get_real_matrix
    procedure :: get_real_matrix_active
    procedure :: get_real_array3
    procedure :: get_real_array3_active
    procedure :: get_real_array3_indexed
    procedure :: get_real_array3_indexed2
    procedure :: get_real_array4
    procedure :: get_real_array4_active
    generic   :: get => get_real_scalar, get_int_scalar, &
         &              get_real_vector, get_int_vector, &
         &              get_real_matrix, get_real_array3, &
         &              get_real_array4, get_real_array3_indexed, &
         &              get_real_array3_indexed2
    generic   :: get_active => get_real_vector_active, get_real_matrix_active, &
         &                     get_real_array3_active, get_real_array4_active
    procedure :: get_global_attribute

    procedure :: set_verbose
    procedure :: transpose_matrices
    procedure :: exists
  end type netcdf_file

contains

  ! --- GENERIC SUBROUTINES ---

  !---------------------------------------------------------------------
  ! Open a NetCDF file with name "file_name", optionally specifying the
  ! verbosity level (0-5)
  subroutine open_netcdf_file(this, file_name, iverbose)


    class(netcdf_file)            :: this
    character(len=*), intent(in)  :: file_name
    integer, intent(in), optional :: iverbose

    integer                       :: istatus

    ! Store verbosity level in object
    if (present(iverbose)) then
      this%file%iverbose = iverbose
    else
      ! By default announce files being opened and closed, but not
      ! variables read/written
      this%file%iverbose = 2
    end if

    ! determine whether we are in an MPI run (standard) or not (Single Column Model)
    call mpi_initialized(this%mpi_enabled, istatus)
    if (istatus /= 0) call abor1 ('easynetcdf MPI_initialized() check failed')

    ! By default we don't transpose 2D arrays on read
    this%file%do_transpose_2d = .false.

    if (this%mpi_enabled) then 
    else !! not an mpi job so is_master_task is set to true
      this%is_master_task = .true.
      call this%file%open(file_name, iverbose)
    end if

  end subroutine open_netcdf_file

  !---------------------------------------------------------------------
  ! Open a NetCDF file with name "file_name", optionally specifying the
  ! verbosity level (0-5), all tasks open file to enable multiple mpi
  ! ranks to collaborate on file reading (without using pNetCDF)
  subroutine open_netcdf_file_active(this, file_name, iverbose)


    class(netcdf_file)            :: this
    character(len=*), intent(in)  :: file_name
    integer, intent(in), optional :: iverbose

    integer                       :: istatus

    ! Store verbosity level in object
    if (present(iverbose)) then
      this%file%iverbose = iverbose
    else
      ! By default announce files being opened and closed, but not
      ! variables read/written
      this%file%iverbose = 2
    end if

    ! By default we don't transpose 2D arrays on read
    this%file%do_transpose_2d = .false.

    call this%file%open(file_name, iverbose)
    this%multi_open = .true.

  end subroutine open_netcdf_file_active


  !---------------------------------------------------------------------
  ! Close the NetCDF file
  subroutine close_netcdf_file(this)
    class(netcdf_file) :: this
    integer            :: istatus

    if (this%is_master_task .or. this%multi_open) then
      call this%file%close()
    end if

  end subroutine close_netcdf_file


  !---------------------------------------------------------------------
  ! Set the verbosity level from 0 to 5, where the codes have the
  ! following meaning: 0=errors only, 1=warning, 2=info, 3=progress,
  ! 4=detailed, 5=debug
  subroutine set_verbose(this, ival)
    class(netcdf_file) :: this
    integer, optional  :: ival

    if (present(ival)) then
      this%file%iverbose = ival
    else
      this%file%iverbose = 2
    end if

  end subroutine set_verbose



  !---------------------------------------------------------------------
  ! Specify whether 2D arrays should be transposed on read
  subroutine transpose_matrices(this, do_transpose)
    class(netcdf_file) :: this
    logical, optional  :: do_transpose

    if (present(do_transpose)) then
      this%file%do_transpose_2d = do_transpose
    else
      this%file%do_transpose_2d = .true.
    end if

  end subroutine transpose_matrices



  ! --- READING SUBROUTINES ---

  !---------------------------------------------------------------------
  ! Return true if the variable is present, false otherwise
  function exists(this, var_name) result(is_present)


    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name

    logical :: is_present

    if (this%is_master_task) then
      is_present = this%file%exists(var_name)
    end if

    if (this%mpi_enabled) then
    end if

  end function exists


  !---------------------------------------------------------------------
  ! Return the number of dimensions of variable with name var_name, or
  ! -1 if the variable is not found
  function get_rank(this, var_name) result(ndims)


    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    integer :: ndims

    if (this%is_master_task) then
      ndims = this%file%get_rank(var_name)
    end if

  end function get_rank


  !---------------------------------------------------------------------
  ! The method "get" will read either a scalar, vector or matrix
  ! depending on the rank of the output argument. This version reads a
  ! scalar.
  subroutine get_real_scalar(this, var_name, scalar)


    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    real(jprb), intent(out)      :: scalar

    if (this%is_master_task) then
      call this%file%get(var_name, scalar)
    end if

    if (this%mpi_enabled) then
    end if

  end subroutine get_real_scalar


  !---------------------------------------------------------------------
  ! Read an integer scalar
  subroutine get_int_scalar(this, var_name, scalar)


    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    integer,          intent(out):: scalar

    if (this%is_master_task) then
      call this%file%get(var_name, scalar)
    end if


  end subroutine get_int_scalar


  !---------------------------------------------------------------------
  ! Read a 1D array into "vector", which must be allocatable and will
  ! be reallocated if necessary
  subroutine get_real_vector(this, var_name, vector)


    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    real(jprb), allocatable, intent(out) :: vector(:)

    integer                      :: n  ! Length of vector

    n = 0

    if (this%is_master_task) then
      call this%file%get(var_name, vector)
      n = size(vector)
    end if

    if (this%mpi_enabled) then
    end if

  end subroutine get_real_vector

  ! version with active rank specified
  ! irequest : for non_blocking broadcasts, the message handle that will have to be waited on
  ! imp_type : choose between JP_(BLOCKING/NON-BLOCKING)_(STANDARD/BUFFERED)
  subroutine get_real_vector_active(this, var_name, vector, iactive_rank, irequest, imp_type)


    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: vector(:)
    integer(jpim), intent(in)            :: iactive_rank
    integer, optional, intent(inout)     :: irequest
    integer, optional, intent(in)        :: imp_type

    integer :: n  ! Length of vector

    n = 0

      call this%file%get(var_name, vector)
      n = size(vector)


  end subroutine get_real_vector_active


  !---------------------------------------------------------------------
  ! Read a 1D integer array into "vector", which must be allocatable
  ! and will be reallocated if necessary
  subroutine get_int_vector(this, var_name, vector)


    class(netcdf_file)                :: this
    character(len=*),     intent(in)  :: var_name
    integer, allocatable, intent(out) :: vector(:)

    integer                      :: n  ! Length of vector

    n = 0

    if (this%is_master_task) then
      call this%file%get(var_name, vector)
      n = size(vector)
    end if

  end subroutine get_int_vector


  !---------------------------------------------------------------------
  ! Read 2D array into "matrix", which must be allocatable and will be
  ! reallocated if necessary.  Whether to transpose is specifed by the
  ! final optional argument, but can also be specified by the
  ! do_transpose_2d class data member.
  subroutine get_real_matrix(this, var_name, matrix, do_transp)


    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    real(jprb), allocatable, intent(out) :: matrix(:,:)
    logical, optional, intent(in):: do_transp ! Transpose data?

    integer                      :: n(2)

    n = 0

    if (this%is_master_task) then
      call this%file%get(var_name, matrix, do_transp)
      n = shape(matrix)
    end if

    if (this%mpi_enabled) then
    end if

  end subroutine get_real_matrix

  ! version with active rank specified
  ! irequest : for non_blocking broadcasts, the message handle that will have to be waited on
  ! imp_type : choose between JP_(BLOCKING/NON-BLOCKING)_(STANDARD/BUFFERED)
  subroutine get_real_matrix_active(this, var_name, matrix, do_transp, iactive_rank, irequest, imp_type)


    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: matrix(:,:)
    logical, optional, intent(in)        :: do_transp ! Transpose data?
    integer(jpim), intent(in)            :: iactive_rank
    integer, optional, intent(inout)     :: irequest
    integer, optional, intent(in)        :: imp_type

    integer                      :: n(2)

    n = 0

      call this%file%get(var_name, matrix, do_transp)
      n = shape(matrix)



  end subroutine get_real_matrix_active


  !---------------------------------------------------------------------
  ! Read 3D array into "var", which must be allocatable and will be
  ! reallocated if necessary.  Whether to pemute is specifed by the
  ! final optional argument
  subroutine get_real_array3(this, var_name, var, ipermute)


    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: var(:,:,:)
    integer, optional, intent(in)        :: ipermute(3)

    integer                              :: n(3)

    n = 0

    if (this%is_master_task) then
      call this%file%get(var_name, var, ipermute)
      n = shape(var)
    end if

    if (this%mpi_enabled) then
    end if

  end subroutine get_real_array3

  ! version with active rank specified
  ! irequest : for non_blocking broadcasts, the message handle that will have to be waited on
  ! imp_type : choose between JP_(BLOCKING/NON-BLOCKING)_(STANDARD/BUFFERED)
  subroutine get_real_array3_active(this, var_name, var, iactive_rank, ipermute, irequest, imp_type)


    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: var(:,:,:)
    integer(jpim), intent(in)            :: iactive_rank
    integer, optional, intent(in)        :: ipermute(3)
    integer, optional, intent(inout)     :: irequest
    integer, optional, intent(in)        :: imp_type

    integer                              :: n(3)

    n = 0

      call this%file%get(var_name, var, ipermute)
      n = shape(var)


  end subroutine get_real_array3_active


  !---------------------------------------------------------------------
  ! Read 3D array into "var", which must be allocatable and will be
  ! reallocated if necessary.  Whether to pemute is specifed by the
  ! final optional argument
  subroutine get_real_array3_indexed(this, var_name, var, index, ipermute)


    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: var(:,:,:)
    integer, intent(in)                  :: index
    integer, optional, intent(in)        :: ipermute(3)

    integer                              :: n(3)

    n = 0

    if (this%is_master_task) then
      call this%file%get(var_name, var, index, ipermute)
      n = shape(var)
    end if

  end subroutine get_real_array3_indexed


  !---------------------------------------------------------------------
  ! Read 3D array into "var", which must be allocatable and will be
  ! reallocated if necessary.  Whether to pemute is specifed by the
  ! final optional argument
  subroutine get_real_array3_indexed2(this, var_name, var, index4, index5, ipermute)


    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: var(:,:,:)
    integer, intent(in)                  :: index4, index5
    integer, optional, intent(in)        :: ipermute(3)

    integer                              :: n(3)

    n = 0

    if (this%is_master_task) then
      call this%file%get(var_name, var, index4, index5, ipermute)
      n = shape(var)
    end if

  end subroutine get_real_array3_indexed2


  !---------------------------------------------------------------------
  ! Read 4D array into "var", which must be allocatable and will be
  ! reallocated if necessary.  Whether to pemute is specifed by the
  ! final optional argument
  subroutine get_real_array4(this, var_name, var, ipermute)


    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: var(:,:,:,:)
    integer, optional, intent(in)        :: ipermute(4)

    integer                              :: n(4)

    n = 0

    if (this%is_master_task) then
      call this%file%get(var_name, var, ipermute)
      n = shape(var)
    end if

    if (this%mpi_enabled) then
    end if

  end subroutine get_real_array4

  ! version with active rank specified
  ! irequest : for non_blocking broadcasts, the message handle that will have to be waited on
  ! imp_type : choose between JP_(BLOCKING/NON-BLOCKING)_(STANDARD/BUFFERED)
  subroutine get_real_array4_active(this, var_name, var, iactive_rank, ipermute, irequest, imp_type)


    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: var(:,:,:,:)
    integer(jpim), intent(in)            :: iactive_rank
    integer, optional, intent(in)        :: ipermute(4)
    integer, optional, intent(inout)     :: irequest
    integer, optional, intent(in)        :: imp_type

    integer                              :: n(4)

    n = 0

      call this%file%get(var_name, var, ipermute)
      n = shape(var)

  end subroutine get_real_array4_active


  !---------------------------------------------------------------------
  ! Get a global attribute as a character string
  subroutine get_global_attribute(this, attr_name, attr_str)


    class(netcdf_file) :: this

    character(len=*), intent(in)    :: attr_name
    character(len=*), intent(inout) :: attr_str

    if (this%is_master_task) then
      call this%file%get_global_attribute(attr_name, attr_str)
    end if

    if (this%mpi_enabled) then
    end if

  end subroutine get_global_attribute

end module easy_netcdf_read_mpi
