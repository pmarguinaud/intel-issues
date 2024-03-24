MODULE tt_owned_mod

type :: tt_owned
  integer :: pp
contains 
  final :: tt_owned_object__final_auto
end type

CONTAINS

subroutine tt_owned_object__final_auto (this)
type (tt_owned), intent (inout) :: this
end subroutine

END

