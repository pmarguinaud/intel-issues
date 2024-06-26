
MODULE FIELD_5LM_ACCESS_MODULE


USE FIELD_MODULE
USE PARKIND1, ONLY : JPRM, JPRB, JPRD, JPIM, JPLM

IMPLICIT NONE

PRIVATE


INTERFACE GET_DEVICE_DATA_RDONLY
  MODULE PROCEDURE :: GET_DEVICE_DATA_RDONLY_FIELD_5LM
END INTERFACE GET_DEVICE_DATA_RDONLY

PUBLIC :: GET_DEVICE_DATA_RDONLY


INTERFACE GET_HOST_DATA_RDONLY
  MODULE PROCEDURE :: GET_HOST_DATA_RDONLY_FIELD_5LM
END INTERFACE GET_HOST_DATA_RDONLY

PUBLIC :: GET_HOST_DATA_RDONLY


INTERFACE GET_DEVICE_DATA_RDWR
  MODULE PROCEDURE :: GET_DEVICE_DATA_RDWR_FIELD_5LM
END INTERFACE GET_DEVICE_DATA_RDWR

PUBLIC :: GET_DEVICE_DATA_RDWR


INTERFACE GET_HOST_DATA_RDWR
  MODULE PROCEDURE :: GET_HOST_DATA_RDWR_FIELD_5LM
END INTERFACE GET_HOST_DATA_RDWR

PUBLIC :: GET_HOST_DATA_RDWR


LOGICAL(KIND=JPLM), TARGET, SAVE :: DUMMY_FIELD_5LM (1, 1, 1, 1, 1)
!$acc declare create (DUMMY_FIELD_5LM)


CONTAINS


  FUNCTION GET_DEVICE_DATA_RDONLY_FIELD_5LM (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_5LM), POINTER :: FIELD_PTR
LOGICAL(KIND=JPLM), POINTER :: PTR(:,:,:,:,:)

END FUNCTION GET_DEVICE_DATA_RDONLY_FIELD_5LM



  FUNCTION GET_HOST_DATA_RDONLY_FIELD_5LM (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_5LM), POINTER :: FIELD_PTR
LOGICAL(KIND=JPLM), POINTER :: PTR(:,:,:,:,:)

END FUNCTION GET_HOST_DATA_RDONLY_FIELD_5LM



  FUNCTION GET_DEVICE_DATA_RDWR_FIELD_5LM (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_5LM), POINTER :: FIELD_PTR
LOGICAL(KIND=JPLM), POINTER :: PTR(:,:,:,:,:)

END FUNCTION GET_DEVICE_DATA_RDWR_FIELD_5LM



  FUNCTION GET_HOST_DATA_RDWR_FIELD_5LM (FIELD_PTR) RESULT (PTR)
CLASS (FIELD_5LM), POINTER :: FIELD_PTR
LOGICAL(KIND=JPLM), POINTER :: PTR(:,:,:,:,:)

END FUNCTION GET_HOST_DATA_RDWR_FIELD_5LM



END MODULE FIELD_5LM_ACCESS_MODULE
