
MODULE FIELD_4RB_UTIL_MODULE


USE FIELD_MODULE
USE FIELD_ACCESS_MODULE

IMPLICIT NONE

INTERFACE LOAD
  MODULE PROCEDURE LOAD_FIELD_4RB
  MODULE PROCEDURE LOAD_FIELD_4RB_PTR
  MODULE PROCEDURE LOAD_FIELD_4RB_VIEW
END INTERFACE

INTERFACE SAVE
  MODULE PROCEDURE SAVE_FIELD_4RB
  MODULE PROCEDURE SAVE_FIELD_4RB_PTR
  MODULE PROCEDURE SAVE_FIELD_4RB_VIEW
END INTERFACE

INTERFACE DIFF
  MODULE PROCEDURE DIFF_FIELD_4RB
END INTERFACE

INTERFACE COPY
  MODULE PROCEDURE COPY_FIELD_4RB
  MODULE PROCEDURE COPY_FIELD_4RB_PTR
  MODULE PROCEDURE COPY_FIELD_4RB_VIEW
END INTERFACE

INTERFACE WIPE
  MODULE PROCEDURE WIPE_FIELD_4RB
  MODULE PROCEDURE WIPE_FIELD_4RB_PTR
  MODULE PROCEDURE WIPE_FIELD_4RB_VIEW
END INTERFACE

INTERFACE HOST
  MODULE PROCEDURE HOST_FIELD_4RB
  MODULE PROCEDURE HOST_FIELD_4RB_PTR
  MODULE PROCEDURE HOST_FIELD_4RB_VIEW
END INTERFACE

INTERFACE CRC64
  MODULE PROCEDURE CRC64_FIELD_4RB
END INTERFACE


CONTAINS

INTEGER*8 FUNCTION CRC64_FIELD_4RB (YD)
CLASS (FIELD_4RB), POINTER :: YD






END FUNCTION

SUBROUTINE LOAD_FIELD_4RB (KLUN, YD)
USE FIELD_ABORT_MODULE
INTEGER (KIND=JPIM), INTENT (IN) :: KLUN
CLASS (FIELD_4RB), POINTER :: YD

END SUBROUTINE

SUBROUTINE SAVE_FIELD_4RB (KLUN, YD)
USE FIELD_ABORT_MODULE
INTEGER (KIND=JPIM), INTENT (IN) :: KLUN
CLASS (FIELD_4RB), POINTER :: YD

END SUBROUTINE

SUBROUTINE DIFF_FIELD_4RB (CDMESS, YD, YO)
USE FIELD_ABORT_MODULE
CHARACTER (LEN=*), INTENT(IN) :: CDMESS
CLASS (FIELD_4RB), POINTER :: YD, YO

END SUBROUTINE

SUBROUTINE COPY_FIELD_4RB (SELF, LDCREATED)
USE FIELD_ABORT_MODULE
CLASS (FIELD_4RB), POINTER :: SELF
LOGICAL, INTENT (IN), OPTIONAL :: LDCREATED

END SUBROUTINE 

SUBROUTINE WIPE_FIELD_4RB (SELF, LDDELETED)
USE FIELD_ABORT_MODULE
CLASS (FIELD_4RB) :: SELF
LOGICAL, INTENT (IN), OPTIONAL :: LDDELETED

END SUBROUTINE 

SUBROUTINE HOST_FIELD_4RB (SELF)
CLASS (FIELD_4RB), POINTER :: SELF

END SUBROUTINE 


SUBROUTINE LOAD_FIELD_4RB_VIEW (KLUN, YD)
INTEGER (KIND=JPIM), INTENT (IN) :: KLUN
CLASS (FIELD_4RB_VIEW) :: YD

END SUBROUTINE

SUBROUTINE SAVE_FIELD_4RB_VIEW (KLUN, YD)
INTEGER (KIND=JPIM), INTENT (IN) :: KLUN
CLASS (FIELD_4RB_VIEW) :: YD

END SUBROUTINE

SUBROUTINE COPY_FIELD_4RB_VIEW (SELF, LDCREATED)
CLASS (FIELD_4RB_VIEW) :: SELF
LOGICAL, INTENT (IN), OPTIONAL :: LDCREATED

END SUBROUTINE 

SUBROUTINE WIPE_FIELD_4RB_VIEW (SELF, LDDELETED)
CLASS (FIELD_4RB_VIEW) :: SELF
LOGICAL, INTENT (IN), OPTIONAL :: LDDELETED

END SUBROUTINE 

SUBROUTINE HOST_FIELD_4RB_VIEW (SELF)
CLASS (FIELD_4RB_VIEW) :: SELF

END SUBROUTINE 


SUBROUTINE LOAD_FIELD_4RB_PTR (KLUN, YD)
USE FIELD_ABORT_MODULE
INTEGER (KIND=JPIM), INTENT (IN) :: KLUN
CLASS (FIELD_4RB_PTR) :: YD

END SUBROUTINE

SUBROUTINE SAVE_FIELD_4RB_PTR (KLUN, YD)
USE FIELD_ABORT_MODULE
INTEGER (KIND=JPIM), INTENT (IN) :: KLUN
CLASS (FIELD_4RB_PTR) :: YD

END SUBROUTINE

SUBROUTINE COPY_FIELD_4RB_PTR (SELF, LDCREATED)
USE FIELD_ABORT_MODULE
CLASS (FIELD_4RB_PTR) :: SELF
LOGICAL, INTENT (IN), OPTIONAL :: LDCREATED

END SUBROUTINE 

SUBROUTINE WIPE_FIELD_4RB_PTR (SELF, LDDELETED)
USE FIELD_ABORT_MODULE
CLASS (FIELD_4RB_PTR) :: SELF
LOGICAL, INTENT (IN), OPTIONAL :: LDDELETED

END SUBROUTINE 

SUBROUTINE HOST_FIELD_4RB_PTR (SELF)
CLASS (FIELD_4RB_PTR) :: SELF

END SUBROUTINE 


END MODULE
