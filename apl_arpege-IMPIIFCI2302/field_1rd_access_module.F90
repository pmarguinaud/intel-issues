
MODULE FIELD_1RD_ACCESS_MODULE


USE FIELD_MODULE
USE PARKIND1, ONLY : JPRM, JPRB, JPRD, JPIM, JPLM

IMPLICIT NONE

PRIVATE


INTERFACE GET_DEVICE_DATA_RDONLY
  MODULE PROCEDURE :: GET_DEVICE_DATA_RDONLY_FIELD_1RD
END INTERFACE GET_DEVICE_DATA_RDONLY

PUBLIC :: GET_DEVICE_DATA_RDONLY


INTERFACE GET_HOST_DATA_RDONLY
  MODULE PROCEDURE :: GET_HOST_DATA_RDONLY_FIELD_1RD
END INTERFACE GET_HOST_DATA_RDONLY

PUBLIC :: GET_HOST_DATA_RDONLY


INTERFACE GET_DEVICE_DATA_RDWR
  MODULE PROCEDURE :: GET_DEVICE_DATA_RDWR_FIELD_1RD
END INTERFACE GET_DEVICE_DATA_RDWR

PUBLIC :: GET_DEVICE_DATA_RDWR


INTERFACE GET_HOST_DATA_RDWR
  MODULE PROCEDURE :: GET_HOST_DATA_RDWR_FIELD_1RD
END INTERFACE GET_HOST_DATA_RDWR

PUBLIC :: GET_HOST_DATA_RDWR


REAL(KIND=JPRD), TARGET, SAVE :: DUMMY_FIELD_1RD (1)
!$acc declare create (DUMMY_FIELD_1RD)


CONTAINS


  FUNCTION GET_DEVICE_DATA_RDONLY_FIELD_1RD (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_1RD), POINTER :: FIELD_PTR
REAL(KIND=JPRD), POINTER :: PTR(:)

END FUNCTION GET_DEVICE_DATA_RDONLY_FIELD_1RD



  FUNCTION GET_HOST_DATA_RDONLY_FIELD_1RD (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_1RD), POINTER :: FIELD_PTR
REAL(KIND=JPRD), POINTER :: PTR(:)

END FUNCTION GET_HOST_DATA_RDONLY_FIELD_1RD



  FUNCTION GET_DEVICE_DATA_RDWR_FIELD_1RD (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_1RD), POINTER :: FIELD_PTR
REAL(KIND=JPRD), POINTER :: PTR(:)

END FUNCTION GET_DEVICE_DATA_RDWR_FIELD_1RD



  FUNCTION GET_HOST_DATA_RDWR_FIELD_1RD (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_1RD), POINTER :: FIELD_PTR
REAL(KIND=JPRD), POINTER :: PTR(:)

END FUNCTION GET_HOST_DATA_RDWR_FIELD_1RD



END MODULE FIELD_1RD_ACCESS_MODULE
