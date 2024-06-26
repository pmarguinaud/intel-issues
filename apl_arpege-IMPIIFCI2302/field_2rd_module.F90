
MODULE FIELD_2RD_MODULE


USE IEEE_ARITHMETIC, ONLY: IEEE_SIGNALING_NAN
USE DEV_ALLOC_MODULE
USE HOST_ALLOC_MODULE
USE FIELD_BASIC_MODULE
USE FIELD_CONSTANTS_MODULE
USE FIELD_DEFAULTS_MODULE
#ifdef _OPENACC
USE OPENACC
#endif
USE PARKIND1, ONLY : JPRM, JPRB, JPRD, JPIM, JPLM

USE FIELD_2RD_DATA_MODULE, ONLY : FIELD_2RD_COPY_INTF

IMPLICIT NONE

PRIVATE

TYPE, ABSTRACT, EXTENDS (FIELD_BASIC) :: FIELD_2RD
  REAL(KIND=JPRD), POINTER :: PTR(:,:) => NULL()
  REAL(KIND=JPRD), POINTER, CONTIGUOUS :: DEVPTR(:,:) => NULL()
  PROCEDURE (FIELD_2RD_COPY_INTF), POINTER, NOPASS :: COPY_FUNC => NULL ()
CONTAINS

  PROCEDURE :: FINAL => FIELD_2RD_FINAL
  PROCEDURE :: FIELD_2RD_FINAL
  PROCEDURE :: DELETE_DEVICE_DATA => FIELD_2RD_DELETE_DEVICE_DATA
  PROCEDURE :: GET_VIEW => FIELD_2RD_GET_VIEW
  PROCEDURE :: GET_DEVICE_DATA_RDONLY => FIELD_2RD_GET_DEVICE_DATA_RDONLY
  PROCEDURE :: GET_DEVICE_DATA_RDWR => FIELD_2RD_GET_DEVICE_DATA_RDWR
  PROCEDURE :: GET_HOST_DATA_RDONLY => FIELD_2RD_GET_HOST_DATA_RDONLY
  PROCEDURE :: GET_HOST_DATA_RDWR => FIELD_2RD_GET_HOST_DATA_RDWR
  PROCEDURE :: SYNC_HOST_RDWR => FIELD_2RD_SYNC_HOST_RDWR
  PROCEDURE :: SYNC_HOST_RDONLY => FIELD_2RD_SYNC_HOST_RDONLY
  PROCEDURE :: SYNC_DEVICE_RDWR => FIELD_2RD_SYNC_DEVICE_RDWR
  PROCEDURE :: SYNC_DEVICE_RDONLY => FIELD_2RD_SYNC_DEVICE_RDONLY
  PROCEDURE :: COPY_OBJECT => FIELD_2RD_COPY_OBJECT
  PROCEDURE :: WIPE_OBJECT => FIELD_2RD_WIPE_OBJECT
  PROCEDURE(GET_DIMS), DEFERRED :: GET_DIMS
  PROCEDURE(RESIZE), DEFERRED :: RESIZE

  PROCEDURE :: GET_DEVICE_DATA => FIELD_2RD_GET_DEVICE_DATA
  PROCEDURE :: GET_HOST_DATA => FIELD_2RD_GET_HOST_DATA
  PROCEDURE, PRIVATE :: FIELD_2RD_GET_HOST_DATA
  PROCEDURE, PRIVATE :: FIELD_2RD_GET_DEVICE_DATA
  PROCEDURE, PRIVATE :: COPY_DATA =>  FIELD_2RD_COPY_DATA
  PROCEDURE :: CREATE_DEVICE_DATA => FIELD_2RD_CREATE_DEVICE_DATA
#ifdef __PGI
  PROCEDURE :: SET_STATUS => FIELD_2RD_SET_STATUS
#endif
END TYPE FIELD_2RD

ABSTRACT INTERFACE
  SUBROUTINE GET_DIMS(SELF, LBOUNDS, UBOUNDS)
    USE PARKIND1, ONLY : JPRM, JPRB, JPRD, JPIM, JPLM
    IMPORT ::  FIELD_2RD
    CLASS(FIELD_2RD),               INTENT(IN) :: SELF
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: LBOUNDS(2)
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: UBOUNDS(2)
  END SUBROUTINE GET_DIMS
  SUBROUTINE RESIZE (SELF, UBOUNDS, LBOUNDS, PERSISTENT)
    USE PARKIND1, ONLY : JPRM, JPRB, JPRD, JPIM, JPLM
    IMPORT ::  FIELD_2RD
    CLASS(FIELD_2RD),               INTENT(IN) :: SELF
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: UBOUNDS(2)
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: LBOUNDS(2)
    LOGICAL, OPTIONAL,            INTENT(IN) :: PERSISTENT
  END SUBROUTINE RESIZE
END INTERFACE

PUBLIC :: FIELD_2RD

TYPE, EXTENDS(FIELD_2RD) :: FIELD_2RD_WRAPPER
  LOGICAL :: SYNC_ON_FINAL = .TRUE.
CONTAINS
  PROCEDURE :: INIT => FIELD_2RD_WRAPPER_INIT
  PROCEDURE :: FINAL => FIELD_2RD_WRAPPER_FINAL
  PROCEDURE :: GET_DIMS => FIELD_2RD_WRAPPER_GET_DIMS
  PROCEDURE :: RESIZE => FIELD_2RD_WRAPPER_RESIZE
END TYPE FIELD_2RD_WRAPPER

PUBLIC :: FIELD_2RD_WRAPPER

TYPE, EXTENDS(FIELD_2RD) :: FIELD_2RD_OWNER
  INTEGER(KIND=JPIM) :: LBOUNDS(2), UBOUNDS(2)
  LOGICAL :: HAS_INIT_VALUE = .FALSE.
  LOGICAL :: PINNED = .FALSE.
  REAL(KIND=JPRD) :: INIT_VALUE
CONTAINS
  PROCEDURE :: INIT => FIELD_2RD_OWNER_INIT
  PROCEDURE :: FINAL => FIELD_2RD_OWNER_FINAL
  PROCEDURE, PRIVATE :: CREATE_HOST_DATA => FIELD_2RD_CREATE_HOST_DATA
  PROCEDURE :: GET_HOST_DATA => FIELD_2RD_OWNER_GET_HOST_DATA
  PROCEDURE :: GET_DEVICE_DATA => FIELD_2RD_OWNER_GET_DEVICE_DATA
  PROCEDURE :: GET_DIMS => FIELD_2RD_OWNER_GET_DIMS
  PROCEDURE :: RESIZE => FIELD_2RD_OWNER_RESIZE
END TYPE FIELD_2RD_OWNER

PUBLIC :: FIELD_2RD_OWNER

TYPE FIELD_2RD_PTR
  CLASS(FIELD_2RD), POINTER :: PTR => NULL()
END TYPE FIELD_2RD_PTR

PUBLIC :: FIELD_2RD_PTR

TYPE FIELD_2RD_VIEW
  REAL(KIND=JPRD), POINTER :: P(:) => NULL()
END TYPE FIELD_2RD_VIEW

PUBLIC :: FIELD_2RD_VIEW


CONTAINS

  SUBROUTINE FIELD_2RD_WRAPPER_INIT(SELF, DATA, PERSISTENT, LBOUNDS, MAP_DEVPTR, SYNC_ON_FINAL)
USE FIELD_ABORT_MODULE
USE FIELD_DEFAULTS_MODULE
USE FIELD_2RD_DATA_MODULE, ONLY : FIELD_2RD_COPY_FUNC

CLASS(FIELD_2RD_WRAPPER) :: SELF
REAL(KIND=JPRD), TARGET, INTENT(IN) :: DATA(:,:)
LOGICAL, INTENT(IN), OPTIONAL :: PERSISTENT
LOGICAL, INTENT(IN), OPTIONAL :: MAP_DEVPTR
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: LBOUNDS(2)
LOGICAL, INTENT(IN), OPTIONAL :: SYNC_ON_FINAL














END SUBROUTINE FIELD_2RD_WRAPPER_INIT

  SUBROUTINE FIELD_2RD_OWNER_INIT (SELF, LBOUNDS, UBOUNDS, PERSISTENT, DELAYED, INIT_VALUE, PINNED, MAP_DEVPTR)
USE FIELD_ABORT_MODULE
USE FIELD_2RD_DATA_MODULE, ONLY : FIELD_2RD_COPY_FUNC
CLASS(FIELD_2RD_OWNER) :: SELF
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: LBOUNDS(2)
INTEGER(KIND=JPIM), INTENT(IN) :: UBOUNDS(2)
LOGICAL, OPTIONAL,  INTENT(IN) :: PERSISTENT
LOGICAL, OPTIONAL,  INTENT(IN) :: DELAYED
LOGICAL, OPTIONAL,  INTENT(IN) :: PINNED
LOGICAL, OPTIONAL,  INTENT(IN) :: MAP_DEVPTR
REAL(KIND=JPRD), OPTIONAL, INTENT(IN) :: INIT_VALUE


















END SUBROUTINE FIELD_2RD_OWNER_INIT

  SUBROUTINE FIELD_2RD_CREATE_HOST_DATA (SELF)

CLASS(FIELD_2RD_OWNER) :: SELF


END SUBROUTINE FIELD_2RD_CREATE_HOST_DATA

  FUNCTION FIELD_2RD_GET_VIEW(SELF, BLOCK_INDEX, ZERO) RESULT(VIEW_PTR)
USE FIELD_ABORT_MODULE
CLASS(FIELD_2RD) :: SELF
REAL(KIND=JPRD), POINTER :: VIEW_PTR(:)
INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
LOGICAL, OPTIONAL,  INTENT(IN) :: ZERO









END FUNCTION FIELD_2RD_GET_VIEW

  SUBROUTINE FIELD_2RD_DELETE_DEVICE_DATA(SELF)

CLASS(FIELD_2RD) :: SELF

END SUBROUTINE FIELD_2RD_DELETE_DEVICE_DATA

  SUBROUTINE FIELD_2RD_FINAL(SELF)

CLASS(FIELD_2RD) :: SELF


END SUBROUTINE FIELD_2RD_FINAL

  SUBROUTINE FIELD_2RD_WRAPPER_FINAL(SELF)

CLASS(FIELD_2RD_WRAPPER) :: SELF



END SUBROUTINE FIELD_2RD_WRAPPER_FINAL

  SUBROUTINE FIELD_2RD_OWNER_FINAL(SELF)

CLASS(FIELD_2RD_OWNER) :: SELF


END SUBROUTINE FIELD_2RD_OWNER_FINAL

  SUBROUTINE FIELD_2RD_COPY_OBJECT (SELF, LDCREATED)
USE FIELD_ABORT_MODULE
CLASS(FIELD_2RD) :: SELF
LOGICAL, INTENT (IN), OPTIONAL :: LDCREATED





END SUBROUTINE FIELD_2RD_COPY_OBJECT

  SUBROUTINE FIELD_2RD_WIPE_OBJECT (SELF, LDDELETED)
USE FIELD_ABORT_MODULE
CLASS(FIELD_2RD) :: SELF
LOGICAL, INTENT (IN), OPTIONAL :: LDDELETED





END SUBROUTINE FIELD_2RD_WIPE_OBJECT

  SUBROUTINE FIELD_2RD_COPY_DATA (SELF, KDIR, QUEUE)
CLASS(FIELD_2RD) :: SELF
INTEGER (KIND=JPIM),           INTENT(IN) :: KDIR
INTEGER (KIND=JPIM), OPTIONAL, INTENT(IN) :: QUEUE





END SUBROUTINE FIELD_2RD_COPY_DATA

  SUBROUTINE FIELD_2RD_GET_HOST_DATA (SELF, MODE, PTR, QUEUE)
CLASS(FIELD_2RD) :: SELF
INTEGER (KIND=JPIM),           INTENT(IN)    :: MODE
REAL(KIND=JPRD), POINTER,          INTENT(INOUT) :: PTR(:,:)
INTEGER (KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE





END SUBROUTINE FIELD_2RD_GET_HOST_DATA

  SUBROUTINE FIELD_2RD_OWNER_GET_HOST_DATA (SELF, MODE, PTR, QUEUE)
CLASS(FIELD_2RD_OWNER) :: SELF
INTEGER (KIND=JPIM),           INTENT(IN)    :: MODE
REAL(KIND=JPRD), POINTER,          INTENT(INOUT) :: PTR(:,:)
INTEGER (KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE


END SUBROUTINE FIELD_2RD_OWNER_GET_HOST_DATA

  SUBROUTINE FIELD_2RD_GET_HOST_DATA_RDONLY (SELF, PPTR, QUEUE)
CLASS(FIELD_2RD) :: SELF
REAL(KIND=JPRD), POINTER,         INTENT(INOUT) :: PPTR(:,:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE

END SUBROUTINE FIELD_2RD_GET_HOST_DATA_RDONLY

  SUBROUTINE FIELD_2RD_SYNC_HOST_RDONLY (SELF, QUEUE)
CLASS(FIELD_2RD) :: SELF
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE


END SUBROUTINE FIELD_2RD_SYNC_HOST_RDONLY

  SUBROUTINE FIELD_2RD_GET_HOST_DATA_RDWR (SELF, PPTR, QUEUE)
CLASS(FIELD_2RD) :: SELF
REAL(KIND=JPRD), POINTER,         INTENT(INOUT) :: PPTR(:,:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE

END SUBROUTINE FIELD_2RD_GET_HOST_DATA_RDWR

  SUBROUTINE FIELD_2RD_SYNC_HOST_RDWR (SELF, QUEUE)
CLASS(FIELD_2RD) :: SELF
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE


END SUBROUTINE FIELD_2RD_SYNC_HOST_RDWR

  SUBROUTINE FIELD_2RD_CREATE_DEVICE_DATA (SELF)
CLASS(FIELD_2RD) :: SELF

END SUBROUTINE

  SUBROUTINE FIELD_2RD_GET_DEVICE_DATA (SELF, MODE, PTR, QUEUE)
CLASS(FIELD_2RD) :: SELF
INTEGER (KIND=JPIM),           INTENT(IN)    :: MODE
REAL(KIND=JPRD), POINTER,          INTENT(INOUT) :: PTR(:,:)
INTEGER (KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE






END SUBROUTINE FIELD_2RD_GET_DEVICE_DATA

  SUBROUTINE FIELD_2RD_OWNER_GET_DEVICE_DATA (SELF, MODE, PTR, QUEUE)
CLASS(FIELD_2RD_OWNER) :: SELF
INTEGER (KIND=JPIM),           INTENT(IN)    :: MODE
REAL(KIND=JPRD), POINTER,          INTENT(INOUT) :: PTR(:,:)
INTEGER (KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE


END SUBROUTINE FIELD_2RD_OWNER_GET_DEVICE_DATA

  SUBROUTINE FIELD_2RD_GET_DEVICE_DATA_RDONLY (SELF, PPTR, QUEUE)
CLASS(FIELD_2RD) :: SELF
REAL(KIND=JPRD), POINTER,         INTENT(INOUT) :: PPTR(:,:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE

END SUBROUTINE FIELD_2RD_GET_DEVICE_DATA_RDONLY

  SUBROUTINE FIELD_2RD_SYNC_DEVICE_RDONLY (SELF, QUEUE)
CLASS(FIELD_2RD) :: SELF
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE


END SUBROUTINE FIELD_2RD_SYNC_DEVICE_RDONLY

  SUBROUTINE FIELD_2RD_GET_DEVICE_DATA_RDWR (SELF, PPTR, QUEUE)
CLASS(FIELD_2RD) :: SELF
REAL(KIND=JPRD), POINTER,         INTENT(INOUT) :: PPTR(:,:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE

END SUBROUTINE FIELD_2RD_GET_DEVICE_DATA_RDWR

  SUBROUTINE FIELD_2RD_SYNC_DEVICE_RDWR (SELF, QUEUE)
CLASS(FIELD_2RD) :: SELF
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: QUEUE


END SUBROUTINE FIELD_2RD_SYNC_DEVICE_RDWR

  SUBROUTINE FIELD_2RD_WRAPPER_GET_DIMS (SELF, LBOUNDS, UBOUNDS)
CLASS(FIELD_2RD_WRAPPER),       INTENT(IN) :: SELF
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: LBOUNDS(2)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: UBOUNDS(2)


END SUBROUTINE FIELD_2RD_WRAPPER_GET_DIMS

  SUBROUTINE FIELD_2RD_OWNER_GET_DIMS (SELF, LBOUNDS, UBOUNDS)
CLASS(FIELD_2RD_OWNER),         INTENT(IN) :: SELF
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: LBOUNDS(2)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: UBOUNDS(2)


END SUBROUTINE FIELD_2RD_OWNER_GET_DIMS

  SUBROUTINE FIELD_2RD_WRAPPER_RESIZE (SELF, UBOUNDS, LBOUNDS, PERSISTENT)
USE FIELD_ABORT_MODULE
CLASS(FIELD_2RD_WRAPPER),       INTENT(IN) :: SELF
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: UBOUNDS(2)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: LBOUNDS(2)
LOGICAL, OPTIONAL,            INTENT(IN) :: PERSISTENT

END SUBROUTINE FIELD_2RD_WRAPPER_RESIZE

  SUBROUTINE FIELD_2RD_OWNER_RESIZE (SELF, UBOUNDS, LBOUNDS, PERSISTENT)
CLASS(FIELD_2RD_OWNER),         INTENT(IN) :: SELF
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: UBOUNDS(2)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: LBOUNDS(2)
LOGICAL, OPTIONAL,            INTENT(IN) :: PERSISTENT






END SUBROUTINE FIELD_2RD_OWNER_RESIZE

#ifdef __PGI
  SUBROUTINE FIELD_2RD_SET_STATUS (SELF, KSTATUS)
CLASS (FIELD_2RD) :: SELF
INTEGER (KIND=JPIM), INTENT (IN) :: KSTATUS

END SUBROUTINE FIELD_2RD_SET_STATUS
#endif


END MODULE FIELD_2RD_MODULE
