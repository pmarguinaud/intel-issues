

IMPLICIT NONE

INTEGER, PARAMETER :: N = 10
REAL, POINTER :: X (:)

X => NULL ()

CALL TOTO (N, X)

ALLOCATE (X (N))

CALL TOTO (N, X)

CONTAINS

SUBROUTINE TOTO (N, X)

INTEGER :: N
REAL :: X (N)

CALL TITI (N, X)

END SUBROUTINE

SUBROUTINE TITI (N, X)

INTEGER :: N
REAL, OPTIONAL :: X (N)

PRINT *, PRESENT (X)


END SUBROUTINE

END
