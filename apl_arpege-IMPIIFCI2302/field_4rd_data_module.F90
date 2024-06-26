
MODULE FIELD_4RD_DATA_MODULE


USE FIELD_CONSTANTS_MODULE
USE PARKIND1, ONLY : JPRM, JPRB, JPRD, JPIM, JPLM

IMPLICIT NONE

PRIVATE


PUBLIC :: FIELD_4RD_COPY
PUBLIC :: FIELD_4RD_COPY_FUNC
PUBLIC :: FIELD_4RD_COPY_INTF

ABSTRACT INTERFACE
  SUBROUTINE FIELD_4RD_COPY_INTF (HST, DEV, MAP_DEVPTR, KDIR, QUEUE)
    IMPORT :: JPIM, JPRD
    REAL(KIND=JPRD), POINTER :: HST (:,:,:,:), DEV (:,:,:,:)
    LOGICAL,                       INTENT (IN) :: MAP_DEVPTR
    INTEGER (KIND=JPIM),           INTENT (IN) :: KDIR
    INTEGER (KIND=JPIM), OPTIONAL, INTENT (IN) :: QUEUE
  END SUBROUTINE
END INTERFACE

CONTAINS


  FUNCTION FIELD_4RD_COPY_FUNC (HST, DEV) RESULT (FUNC)
USE FIELD_ABORT_MODULE
PROCEDURE (FIELD_4RD_COPY_INTF), POINTER :: FUNC
REAL(KIND=JPRD), POINTER, OPTIONAL :: HST (:,:,:,:), DEV (:,:,:,:)




END FUNCTION

  SUBROUTINE FIELD_4RD_COPY (HST, DEV, MAP_DEVPTR, KDIR, QUEUE)
USE FIELD_ABORT_MODULE
REAL(KIND=JPRD), POINTER :: HST (:,:,:,:), DEV (:,:,:,:)
LOGICAL,                       INTENT (IN) :: MAP_DEVPTR
INTEGER (KIND=JPIM),           INTENT (IN) :: KDIR
INTEGER (KIND=JPIM), OPTIONAL, INTENT (IN) :: QUEUE



END SUBROUTINE

  SUBROUTINE FIELD_4RD_COPY_DIM0_CONTIGUOUS (HST, DEV, MAP_DEVPTR, KDIR, QUEUE)
#ifdef _OPENACC
USE OPENACC
#endif
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : INT64
REAL(KIND=JPRD), POINTER :: HST (:,:,:,:), DEV (:,:,:,:)
LOGICAL,                       INTENT (IN) :: MAP_DEVPTR
INTEGER (KIND=JPIM),           INTENT (IN) :: KDIR
INTEGER (KIND=JPIM), OPTIONAL, INTENT (IN) :: QUEUE


#ifdef _OPENACC

#endif

END SUBROUTINE

  SUBROUTINE FIELD_4RD_COPY_DIM1_CONTIGUOUS (HST, DEV, MAP_DEVPTR, KDIR, QUEUE)
#ifdef _OPENACC
USE OPENACC
#endif
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : INT64
REAL(KIND=JPRD), POINTER :: HST (:,:,:,:), DEV (:,:,:,:)
LOGICAL,                       INTENT (IN) :: MAP_DEVPTR
INTEGER (KIND=JPIM),           INTENT (IN) :: KDIR
INTEGER (KIND=JPIM), OPTIONAL, INTENT (IN) :: QUEUE


#ifdef _OPENACC

#endif

END SUBROUTINE

  SUBROUTINE FIELD_4RD_COPY_DIM2_CONTIGUOUS (HST, DEV, MAP_DEVPTR, KDIR, QUEUE)
#ifdef _OPENACC
USE OPENACC
#endif
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : INT64
REAL(KIND=JPRD), POINTER :: HST (:,:,:,:), DEV (:,:,:,:)
LOGICAL,                       INTENT (IN) :: MAP_DEVPTR
INTEGER (KIND=JPIM),           INTENT (IN) :: KDIR
INTEGER (KIND=JPIM), OPTIONAL, INTENT (IN) :: QUEUE


#ifdef _OPENACC

#endif

END SUBROUTINE

  SUBROUTINE FIELD_4RD_COPY_DIM3_CONTIGUOUS (HST, DEV, MAP_DEVPTR, KDIR, QUEUE)
#ifdef _OPENACC
USE OPENACC
#endif
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : INT64
REAL(KIND=JPRD), POINTER :: HST (:,:,:,:), DEV (:,:,:,:)
LOGICAL,                       INTENT (IN) :: MAP_DEVPTR
INTEGER (KIND=JPIM),           INTENT (IN) :: KDIR
INTEGER (KIND=JPIM), OPTIONAL, INTENT (IN) :: QUEUE


#ifdef _OPENACC

#endif

END SUBROUTINE

  SUBROUTINE FIELD_4RD_COPY_DIM4_CONTIGUOUS (HST, DEV, MAP_DEVPTR, KDIR, QUEUE)
#ifdef _OPENACC
USE OPENACC
#endif
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : INT64
REAL(KIND=JPRD), POINTER :: HST (:,:,:,:), DEV (:,:,:,:)
LOGICAL,                       INTENT (IN) :: MAP_DEVPTR
INTEGER (KIND=JPIM),           INTENT (IN) :: KDIR
INTEGER (KIND=JPIM), OPTIONAL, INTENT (IN) :: QUEUE


#ifdef _OPENACC

#endif
#ifdef _OPENACC

#endif


END SUBROUTINE





  INTEGER (KIND=JPIM) FUNCTION FIELD_4RD_GET_LAST_CONTIGUOUS_DIMENSION (PTR, AFTER) RESULT (JDIM)
REAL(KIND=JPRD), POINTER :: PTR (:,:,:,:)
INTEGER (KIND=JPIM) :: AFTER










END FUNCTION FIELD_4RD_GET_LAST_CONTIGUOUS_DIMENSION


END MODULE FIELD_4RD_DATA_MODULE
