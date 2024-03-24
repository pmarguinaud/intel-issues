USE OMP_LIB
IMPLICIT NONE

! OpenMP 3.1

! For a firstprivate clause on a parallel or task construct, the initial value of
! the new list item is the value of the original list item that exists immediately
! prior to the construct in the task region where the construct is encountered.

! If the original list item does not have the POINTER attribute, initialization of
! the new list items occurs as if by intrinsic assignment, unless the original list item
! has the allocation status of not currently allocated, in which case the new list items
! will have the same status.
! If the original list item has the POINTER attribute, the new list items receive
! the same association status of the original list item as if by pointer assignment.


TYPE T1
  INTEGER :: K
END TYPE T1

TYPE T2
  TYPE (T1), ALLOCATABLE :: YT1 (:)
END TYPE

INTEGER, PARAMETER :: N = 40, M = 20

INTEGER :: J
TYPE (T2) :: YLT2

ALLOCATE (YLT2%YT1 (M))


!$OMP PARALLEL FIRSTPRIVATE (YLT2)

!$OMP DO PRIVATE (J)
  DO J = 1, N
    YLT2%YT1 (:)%K = J
    PRINT *, OMP_GET_THREAD_NUM (), J, LOC (YLT2%YT1 (1)), YLT2%YT1 (3)%K == J
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

PRINT *, '-------'
CALL PAR1 (YLT2)
PRINT *, '-------'
CALL PAR2 (YLT2)
PRINT *, '-------'

CONTAINS

SUBROUTINE PAR1 (YDT2)

TYPE (T2) :: YDT2
TYPE (T2) :: YLT2
INTEGER :: J

!$OMP PARALLEL PRIVATE (YLT2)

YLT2 = YDT2

!$OMP DO PRIVATE (J)
  DO J = 1, N
    YLT2%YT1 (:)%K = J
    PRINT *, OMP_GET_THREAD_NUM (), J, LOC (YLT2%YT1 (1)), YLT2%YT1 (3)%K == J
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE PAR2 (YDT2)

TYPE (T2) :: YDT2

INTEGER :: J

!$OMP PARALLEL FIRSTPRIVATE (YDT2)

!$OMP DO PRIVATE (J)
  DO J = 1, N
    YDT2%YT1 (:)%K = J
    PRINT *, OMP_GET_THREAD_NUM (), J, LOC (YDT2%YT1 (1)), YDT2%YT1 (3)%K == J
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

END
