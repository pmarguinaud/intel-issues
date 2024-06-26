
MODULE FIELD_2RM_ACCESS_MODULE


USE FIELD_MODULE
USE PARKIND1, ONLY : JPRM, JPRB, JPRD, JPIM, JPLM

IMPLICIT NONE

PRIVATE


INTERFACE GET_DEVICE_DATA_RDONLY
  MODULE PROCEDURE :: GET_DEVICE_DATA_RDONLY_FIELD_2RM
END INTERFACE GET_DEVICE_DATA_RDONLY

PUBLIC :: GET_DEVICE_DATA_RDONLY


INTERFACE GET_HOST_DATA_RDONLY
  MODULE PROCEDURE :: GET_HOST_DATA_RDONLY_FIELD_2RM
END INTERFACE GET_HOST_DATA_RDONLY

PUBLIC :: GET_HOST_DATA_RDONLY


INTERFACE GET_DEVICE_DATA_RDWR
  MODULE PROCEDURE :: GET_DEVICE_DATA_RDWR_FIELD_2RM
END INTERFACE GET_DEVICE_DATA_RDWR

PUBLIC :: GET_DEVICE_DATA_RDWR


INTERFACE GET_HOST_DATA_RDWR
  MODULE PROCEDURE :: GET_HOST_DATA_RDWR_FIELD_2RM
END INTERFACE GET_HOST_DATA_RDWR

PUBLIC :: GET_HOST_DATA_RDWR


REAL(KIND=JPRM), TARGET, SAVE :: DUMMY_FIELD_2RM (1, 1)
!$acc declare create (DUMMY_FIELD_2RM)


CONTAINS


  FUNCTION GET_DEVICE_DATA_RDONLY_FIELD_2RM (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_2RM), POINTER :: FIELD_PTR
REAL(KIND=JPRM), POINTER :: PTR(:,:)

END FUNCTION GET_DEVICE_DATA_RDONLY_FIELD_2RM



  FUNCTION GET_HOST_DATA_RDONLY_FIELD_2RM (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_2RM), POINTER :: FIELD_PTR
REAL(KIND=JPRM), POINTER :: PTR(:,:)

END FUNCTION GET_HOST_DATA_RDONLY_FIELD_2RM



  FUNCTION GET_DEVICE_DATA_RDWR_FIELD_2RM (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_2RM), POINTER :: FIELD_PTR
REAL(KIND=JPRM), POINTER :: PTR(:,:)

END FUNCTION GET_DEVICE_DATA_RDWR_FIELD_2RM



  FUNCTION GET_HOST_DATA_RDWR_FIELD_2RM (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_2RM), POINTER :: FIELD_PTR
REAL(KIND=JPRM), POINTER :: PTR(:,:)

END FUNCTION GET_HOST_DATA_RDWR_FIELD_2RM



END MODULE FIELD_2RM_ACCESS_MODULE
