
MODULE FIELD_5RM_FACTORY_MODULE


USE FIELD_MODULE
USE FIELD_GANG_MODULE
USE PARKIND1, ONLY : JPRM, JPRB, JPRD, JPIM, JPLM

IMPLICIT NONE

PRIVATE

INTERFACE FIELD_NEW
  MODULE PROCEDURE FIELD_5RM_NEW_OWNER
  MODULE PROCEDURE FIELD_5RM_NEW_WRAPPER
  MODULE PROCEDURE FIELD_5RM_NEW_GANG_WRAPPER
  MODULE PROCEDURE FIELD_5RM_NEW_GANG_OWNER
END INTERFACE

PUBLIC :: FIELD_NEW

INTERFACE FIELD_DELETE
  MODULE PROCEDURE FIELD_5RM_DELETE
END INTERFACE FIELD_DELETE

PUBLIC :: FIELD_DELETE

INTERFACE FIELD_RESIZE
  MODULE PROCEDURE FIELD_5RM_RESIZE
END INTERFACE FIELD_RESIZE

PUBLIC :: FIELD_RESIZE

CONTAINS

SUBROUTINE FIELD_5RM_NEW_OWNER (FIELD_PTR, UBOUNDS, LBOUNDS, PERSISTENT, DELAYED, INIT_VALUE, PINNED, MAP_DEVPTR)
CLASS(FIELD_5RM), POINTER :: FIELD_PTR

INTEGER(KIND=JPIM), INTENT(IN) :: UBOUNDS (5)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: LBOUNDS (5)
LOGICAL, OPTIONAL, INTENT(IN) :: PERSISTENT
LOGICAL, OPTIONAL,  INTENT(IN) :: DELAYED
REAL(KIND=JPRM), OPTIONAL, INTENT(IN) :: INIT_VALUE
LOGICAL, OPTIONAL,  INTENT(IN) :: PINNED
LOGICAL, OPTIONAL,  INTENT(IN) :: MAP_DEVPTR



END SUBROUTINE

SUBROUTINE FIELD_5RM_NEW_WRAPPER (FIELD_PTR, LBOUNDS, PERSISTENT, DATA, MAP_DEVPTR, SYNC_ON_FINAL)
CLASS(FIELD_5RM), POINTER :: FIELD_PTR
REAL(KIND=JPRM), TARGET, INTENT (IN) :: DATA (:,:,:,:,:)

INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: LBOUNDS (5)
LOGICAL, OPTIONAL, INTENT(IN) :: PERSISTENT
LOGICAL, OPTIONAL,  INTENT(IN) :: MAP_DEVPTR
LOGICAL, OPTIONAL,  INTENT(IN) :: SYNC_ON_FINAL



END SUBROUTINE

SUBROUTINE FIELD_5RM_NEW_GANG_WRAPPER (FIELD_PTR, CHILDREN, LBOUNDS, PERSISTENT, DATA, SYNC_ON_FINAL)
CLASS(FIELD_5RM), POINTER :: FIELD_PTR
TYPE(FIELD_4RM_PTR), ALLOCATABLE :: CHILDREN (:)
REAL(KIND=JPRM), TARGET, INTENT (IN) :: DATA (:,:,:,:,:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: LBOUNDS (5)
LOGICAL, OPTIONAL, INTENT(IN) :: PERSISTENT
LOGICAL, OPTIONAL,  INTENT(IN) :: SYNC_ON_FINAL







END SUBROUTINE

SUBROUTINE FIELD_5RM_NEW_GANG_OWNER (FIELD_PTR, CHILDREN, UBOUNDS, LBOUNDS, PERSISTENT, DELAYED, INIT_VALUE)
CLASS(FIELD_5RM), POINTER :: FIELD_PTR
TYPE(FIELD_4RM_PTR), ALLOCATABLE :: CHILDREN (:)
INTEGER(KIND=JPIM), INTENT(IN) :: UBOUNDS (5)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: LBOUNDS (5)
LOGICAL, OPTIONAL, INTENT(IN) :: PERSISTENT
LOGICAL, OPTIONAL,  INTENT(IN) :: DELAYED
REAL(KIND=JPRM), OPTIONAL, INTENT(IN) :: INIT_VALUE







END SUBROUTINE


SUBROUTINE FIELD_5RM_DELETE (FIELD_PTR)
CLASS(FIELD_5RM), POINTER :: FIELD_PTR



END SUBROUTINE

SUBROUTINE FIELD_5RM_RESIZE (FIELD_PTR, UBOUNDS, LBOUNDS, PERSISTENT)
CLASS(FIELD_5RM), POINTER :: FIELD_PTR
INTEGER(KIND=JPIM), INTENT(IN) :: UBOUNDS (5)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: LBOUNDS (5)
LOGICAL, OPTIONAL, INTENT(IN) :: PERSISTENT

END SUBROUTINE FIELD_5RM_RESIZE


END MODULE
