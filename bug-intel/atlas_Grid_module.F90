module atlas_Grid_module

implicit none

public :: atlas_Grid
public :: atlas_StructuredGrid

private

TYPE :: atlas_Grid
  character*32 :: id
END TYPE atlas_Grid

interface atlas_Grid
  module procedure atlas_Grid__ctor_id
end interface

TYPE, extends(atlas_Grid) :: atlas_StructuredGrid
END TYPE atlas_StructuredGrid

interface atlas_StructuredGrid
  module procedure atlas_StructuredGrid__ctor_id
end interface

contains

function atlas_Grid__ctor_id(identifier) result(this)
  type(atlas_Grid) :: this
  character(len=*), intent(in) :: identifier
  this%id = identifier
end function

function atlas_StructuredGrid__ctor_id(identifier) result(this)
  type(atlas_StructuredGrid) :: this
  character(len=*), intent(in) :: identifier
  this%id = identifier
end function

end module atlas_Grid_module
