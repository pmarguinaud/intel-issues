MODULE FIELD_CONTAINER_SP_MOD
USE PARKIND1                , ONLY : JPIM, JPRB
USE YOMHOOK                 , ONLY : LHOOK, DR_HOOK, JPHOOK
USE FIELD_CONTAINER_BASE_MOD
USE FIELD_DEFINITIONS_BASE  , ONLY : FIELD_ACCESS_BASE, FIELD_METADATA_BASE

IMPLICIT NONE

PRIVATE
PUBLIC :: TYPE_ITERATOR
TYPE TYPE_FIELD_INDEX
  LOGICAL :: L_ON = .FALSE.
END TYPE TYPE_FIELD_INDEX
!!$TYPE TYPE_STORAGE_ARRAYS
!!$  REAL(KIND=JPRB),POINTER :: STORE1D(:) => NULL()
!!$  REAL(KIND=JPRB),POINTER :: STORE2D(:,:) => NULL()
!!$  REAL(KIND=JPRB),POINTER :: STORE3D(:,:,:) => NULL()
!!$END type TYPE_STORAGE_ARRAYS
!!$TYPE, PUBLIC :: TYPE_ITERATOR
!!$INTEGER :: ICURRENT=0
!!$END type TYPE_ITERATOR
TYPE,PUBLIC,EXTENDS(FIELD_CONTAINER_BASE) :: FIELD_CONTAINER_SP
  CONTAINS
  PROCEDURE :: FIELD_CREATE
  PROCEDURE :: FIELD_DESTROY
  PROCEDURE :: FIELD_OPEN_FAC
  PROCEDURE :: FIELD_CLOSE_FAC
  PROCEDURE :: FIELD_ITERATED
  PROCEDURE :: FIELD_OPEN_ONE_1D
  PROCEDURE :: FIELD_OPEN_ONE_2D
  GENERIC,PUBLIC :: FIELD_OPEN_ONE => FIELD_OPEN_ONE_1D
  GENERIC,PUBLIC :: FIELD_OPEN_ONE => FIELD_OPEN_ONE_2D
  PROCEDURE :: FIELD_CLOSE_ONE_1D
  PROCEDURE :: FIELD_CLOSE_ONE_2D
  GENERIC,PUBLIC :: FIELD_CLOSE_ONE => FIELD_CLOSE_ONE_1D
  GENERIC,PUBLIC :: FIELD_CLOSE_ONE => FIELD_CLOSE_ONE_2D
  PROCEDURE :: ZERO
END TYPE FIELD_CONTAINER_SP

CONTAINS
!=======================================================================================
SUBROUTINE FIELD_CREATE(SELF,NAMESPACE,KPOINTS,KLEVELS,KFIDS,KATTACH)
CLASS(FIELD_CONTAINER_SP), INTENT(INOUT) :: SELF
CLASS(FIELD_METADATA_BASE),   INTENT(IN)    :: NAMESPACE  ! Metadata defining what namespace we're using (Main, radiation, etc.)
INTEGER(KIND=JPIM),           INTENT(IN)    :: KPOINTS    ! size per level of field on this PE
INTEGER(KIND=JPIM),           INTENT(IN)    :: KLEVELS(:) ! vertical dimensions
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KFIDS(:)    ! list of field IDs (otherwise make all fields)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KATTACH     ! Attachment mode (GFL backwards compatibility)

INTEGER(KIND=JPIM) :: ID
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_CREATE',0,ZHOOK_HANDLE)
CALL SELF%FIELD_CREATE_BASE(NAMESPACE,KPOINTS,KLEVELS,KFIDS,KATTACH)
CALL SELF%SETUP_ACTIVE_FIELDS(KFIDS)
IF(SELF%GET_STORE_TYPE() == JP_STORAGE_ARR) THEN
  CALL SELF%ALLOCATE_STORAGE_ARRAYS()
ELSEIF(SELF%GET_STORE_TYPE() == JP_STORAGEOLD) THEN
ENDIF
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_CREATE',1,ZHOOK_HANDLE)
END SUBROUTINE FIELD_CREATE
!===============================================================================
SUBROUTINE FIELD_DESTROY(SELF)
CLASS(FIELD_CONTAINER_SP), INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM) :: ID
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_DESTROY',0,ZHOOK_HANDLE)
IF(SELF%GET_STORE_TYPE() == JP_STORAGE_ARR) THEN
  DO ID=1,SELF%NFIELDS
    IF (SELF%FIELD_INDEX(ID)%L_ON) THEN
      IF (SELF%METADATA(ID)%NDIMS == 0) THEN
        DEALLOCATE(SELF%STORAGE_ARRAYS(ID)%STORE1D)
      ELSEIF (SELF%METADATA(ID)%NDIMS == 1) THEN
        DEALLOCATE(SELF%STORAGE_ARRAYS(ID)%STORE1D)
      ELSEIF (SELF%METADATA(ID)%NDIMS == 2) THEN
        DEALLOCATE(SELF%STORAGE_ARRAYS(ID)%STORE2D)
      ELSEIF (SELF%METADATA(ID)%NDIMS == 3) THEN
        DEALLOCATE(SELF%STORAGE_ARRAYS(ID)%STORE3D)
      ENDIF
    ENDIF
  ENDDO
  DEALLOCATE(SELF%STORAGE_ARRAYS)
ENDIF
DEALLOCATE(SELF%NLEVELS)
DEALLOCATE(SELF%FIELD_INDEX)
DEALLOCATE(SELF%MALL_FIELDS)
DEALLOCATE(SELF%METADATA)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_DESTROY',1,ZHOOK_HANDLE)
END SUBROUTINE FIELD_DESTROY
!===============================================================================
SUBROUTINE FIELD_OPEN_FAC(SELF, KFIDS, FAC)
CLASS(FIELD_CONTAINER_SP), INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM), TARGET, OPTIONAL,INTENT(IN) :: KFIDS(:)
CLASS(FIELD_ACCESS_BASE),           INTENT(OUT) :: FAC

INTEGER(KIND=JPIM) :: IFIELDS_USED,ID,J
INTEGER(KIND=JPIM),ALLOCATABLE :: IUSED_FIDS(:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!LOGICAL :: LLON
!!$IFIELDS_USED = 0
!!$DO ID=1,SELF%NFIELDS
!!$  LLON = .FALSE.
!!$  IF (SELF%FIELD_INDEX(ID)%L_ON) THEN
!!$    IF (PRESENT(KFIDS)) THEN
!!$      IF (ANY(ID==KFIDS)) THEN
!!$        LLON = .TRUE.
!!$      ENDIF
!!$    ELSE
!!$      LLON = .TRUE.
!!$    ENDIF
!!$  ENDIF
!!$
!!$  IF(LLON) THEN
!!$    IFIELDS_USED = IFIELDS_USED + 1
!!$    IUSED_FIDS(IFIELDS_USED) = ID
!!$  ENDIF
!!$
!!$ENDDO
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_OPEN_FAC',0,ZHOOK_HANDLE)
CALL SELF%GET_USED_FIELDS(IUSED_FIDS,KFIDS)
IFIELDS_USED = SIZE(IUSED_FIDS)
IF(PRESENT(KFIDS)) THEN
  IF(IFIELDS_USED < SIZE(KFIDS)) THEN
    WRITE(0,*) 'FIELD_OPEN - KFIDS MISSING '
    DO J=1,SIZE(KFIDS)
      IF(.NOT. ANY(IUSED_FIDS(1:IFIELDS_USED) == KFIDS(J))) WRITE(0,*) 'FID ',J,KFIDS(J),' MISSING'
    ENDDO
  ENDIF
ENDIF


ALLOCATE(FAC%FIELDS(IFIELDS_USED))
FAC%FIELDS(:) = IUSED_FIDS(1:IFIELDS_USED)

CALL MAP_NAMES_TO_STORAGE(SELF, FAC%FIELDS, MFAC=FAC)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_OPEN_FAC',1,ZHOOK_HANDLE)
END SUBROUTINE FIELD_OPEN_FAC
!===============================================================================

!===============================================================================
SUBROUTINE MAP_NAMES_TO_STORAGE(SELF, KFIDS, MFAC, LD_RELEASE)
CLASS(FIELD_CONTAINER_SP),TARGET,INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM),            INTENT(IN)    :: KFIDS(:)
CLASS(FIELD_ACCESS_BASE), OPTIONAL, INTENT(INOUT) :: MFAC
LOGICAL, OPTIONAL,             INTENT(IN)    :: LD_RELEASE ! Indicates when access is being closed
REAL(KIND=JPRB), POINTER :: FIELD_1D(:)
REAL(KIND=JPRB), POINTER :: FIELD_2D(:,:)
REAL(KIND=JPRB), POINTER :: FIELD_3D(:,:,:)
INTEGER(KIND=JPIM)       :: ID,INFIDS,J
LOGICAL :: LL_NULLIFY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:MAP_NAMES_TO_STORAGE',0,ZHOOK_HANDLE)
LL_NULLIFY = .FALSE.
IF(PRESENT(LD_RELEASE)) LL_NULLIFY=LD_RELEASE
INFIDS = SIZE(KFIDS)
DO J=1,INFIDS
  ID = KFIDS(J)
  IF (SELF%FIELD_INDEX(ID)%L_ON) THEN
    IF(SELF%GET_STORE_TYPE() == JP_STORAGE_ARR) THEN
      IF (SELF%METADATA(ID)%NDIMS == 0) THEN
        FIELD_1D => SELF%STORAGE_ARRAYS(ID)%STORE1D
        CALL MFAC%FIELD_MAP_STORAGE( ID, STORAGE_1D=FIELD_1D)
      ELSEIF (SELF%METADATA(ID)%NDIMS == 1) THEN
        FIELD_1D => SELF%STORAGE_ARRAYS(ID)%STORE1D
        CALL MFAC%FIELD_MAP_STORAGE( ID, STORAGE_1D=FIELD_1D)
      ELSEIF (SELF%METADATA(ID)%NDIMS == 2) THEN
        FIELD_2D => SELF%STORAGE_ARRAYS(ID)%STORE2D
        CALL MFAC%FIELD_MAP_STORAGE( ID, STORAGE_2D=FIELD_2D)
      ELSEIF (SELF%METADATA(ID)%NDIMS == 3) THEN
        FIELD_3D =>  SELF%STORAGE_ARRAYS(ID)%STORE3D
        CALL MFAC%FIELD_MAP_STORAGE( ID, STORAGE_3D=FIELD_3D)
      ENDIF
    ENDIF
  ENDIF
ENDDO
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:MAP_NAMES_TO_STORAGE',1,ZHOOK_HANDLE)

END SUBROUTINE MAP_NAMES_TO_STORAGE
!===============================================================================
SUBROUTINE FIELD_CLOSE_FAC(SELF, FAC)
CLASS(FIELD_CONTAINER_SP), INTENT(INOUT) :: SELF
CLASS(FIELD_ACCESS_BASE),    INTENT(INOUT) :: FAC

INTEGER(KIND=JPIM) :: J,ID,INFIDS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_CLOSE_FAC',0,ZHOOK_HANDLE)
INFIDS=SIZE(FAC%FIELDS)
DO J=1,INFIDS
  ID = FAC%FIELDS(J)
  CALL FAC%FIELD_MAP_STORAGE( ID,ld_nullify=.TRUE.)
ENDDO
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_CLOSE_FAC',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_CLOSE_FAC
!===============================================================================
FUNCTION FIELD_ITERATED(SELF, ITERATOR,KFIDS, FIELD_1D, FIELD_2D,KOTID)
INTEGER(KIND=JPIM) :: FIELD_ITERATED
CLASS(FIELD_CONTAINER_SP),TARGET, INTENT(IN) :: SELF
TYPE(TYPE_ITERATOR), INTENT(INOUT) :: ITERATOR
INTEGER(KIND=JPIM), OPTIONAL, TARGET,INTENT(IN)    :: KFIDS(:)    ! list of field IDs (otherwise make all fields)
REAL(KIND=JPRB), POINTER, OPTIONAL, INTENT(OUT)     :: FIELD_1D(:)   ! 1D field
REAL(KIND=JPRB), POINTER, OPTIONAL, INTENT(OUT)     :: FIELD_2D(:,:) ! 2D field
INTEGER(KIND=JPIM), OPTIONAL,         INTENT(OUT)   :: KOTID         ! ID of the field being returned

INTEGER(KIND=JPIM), POINTER :: IUSE_FIDS(:)
INTEGER(KIND=JPIM) :: INO_FLDS_ITER,ID
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_ITERATED',0,ZHOOK_HANDLE)
IF(PRESENT(KFIDS)) THEN
  IUSE_FIDS=>KFIDS ! AJGDB - is this guaranteed not to disapper? No!
ELSE
  IUSE_FIDS=>SELF%MALL_FIELDS
ENDIF
INO_FLDS_ITER = SIZE(IUSE_FIDS)

IF(ITERATOR%ICURRENT < INO_FLDS_ITER) THEN
  ITERATOR%ICURRENT = ITERATOR%ICURRENT+1
  ID = IUSE_FIDS(ITERATOR%ICURRENT)
  IF(PRESENT(KOTID)) KOTID = ID
  IF(SELF%FIELD_INDEX(ID)%L_ON) THEN
    IF(SELF%METADATA(ID)%NDIMS==1 .AND. PRESENT(FIELD_1D)) THEN
      IF(SELF%GET_STORE_TYPE() == JP_STORAGE_ARR) THEN
        FIELD_1D => SELF%STORAGE_ARRAYS(ID)%STORE1D
      ENDIF
      IF( PRESENT(FIELD_2D) ) NULLIFY(FIELD_2D)
    ELSEIF(SELF%METADATA(ID)%NDIMS==2 .AND. PRESENT(FIELD_2D)) THEN
      IF(SELF%GET_STORE_TYPE() == JP_STORAGE_ARR) THEN
        FIELD_2D => SELF%STORAGE_ARRAYS(ID)%STORE2D
      ENDIF
      IF( PRESENT(FIELD_1D) ) NULLIFY(FIELD_1D)
    ENDIF
    FIELD_ITERATED = 1
  ELSE
    IF(PRESENT(FIELD_1D)) NULLIFY(FIELD_1D)
    IF(PRESENT(FIELD_2D)) NULLIFY(FIELD_2D)
    FIELD_ITERATED = 0
  ENDIF
ELSE
  IF(PRESENT(FIELD_1D)) NULLIFY(FIELD_1D)
  IF(PRESENT(FIELD_2D)) NULLIFY(FIELD_2D)
  FIELD_ITERATED = 0
  ITERATOR%ICURRENT = 0
ENDIF
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_ITERATED',1,ZHOOK_HANDLE)

END FUNCTION FIELD_ITERATED



SUBROUTINE FIELD_OPEN_ONE_1D(SELF, KFID, PFIELD)
! -------------------------------------------------
! Get access to one model field via a pointer
! There are versions for 1D and 2D fields
! -------------------------------------------------

CLASS(FIELD_CONTAINER_SP),      INTENT(IN) :: SELF
INTEGER(KIND=JPIM),         INTENT(IN)    :: KFID      ! field ID
REAL(KIND=JPRB), POINTER, INTENT(OUT)     :: PFIELD(:)

INTEGER(KIND=JPIM) :: IFIDS(1)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_OPEN_ONE_1D',0,ZHOOK_HANDLE)
IFIDS(1) = KFID

! Track the use of this storage via an "access token"

CALL FMAP_1D(SELF, KFID, PFIELD)

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_OPEN_ONE_1D',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_OPEN_ONE_1D
!=======================================================================================
SUBROUTINE FIELD_OPEN_ONE_2D(SELF, KFID, PFIELD)

CLASS(FIELD_CONTAINER_SP),        INTENT(IN) :: SELF
INTEGER(KIND=JPIM),           INTENT(IN)    :: KFID ! field ID
REAL(KIND=JPRB), POINTER, INTENT(OUT)       :: PFIELD(:,:)

INTEGER(KIND=JPIM) :: IFIDS(1)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_OPEN_ONE_2D',0,ZHOOK_HANDLE)
IFIDS(1) = KFID

CALL FMAP_2D(SELF, KFID, PFIELD)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_OPEN_ONE_2D',1,ZHOOK_HANDLE)


END SUBROUTINE FIELD_OPEN_ONE_2D
!=======================================================================================
SUBROUTINE FIELD_CLOSE_ONE_1D(SELF, KFID, PFIELD)
! -------------------------------------------------
! Finish accessing model fields (one field)
! -------------------------------------------------

CLASS(FIELD_CONTAINER_SP), INTENT(IN)    :: SELF
INTEGER(KIND=JPIM),        INTENT(IN)    :: KFID ! field ID
REAL(KIND=JPRB), POINTER,  INTENT(INOUT) :: PFIELD(:)

INTEGER(KIND=JPIM) :: IFIDS(1)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_CLOSE_ONE_1D',0,ZHOOK_HANDLE)
IFIDS(1) = KFID

! Release the "lock" on these fields and nullify the pointer
NULLIFY(PFIELD)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_CLOSE_ONE_1D',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_CLOSE_ONE_1D
!=======================================================================================
SUBROUTINE FIELD_CLOSE_ONE_2D(SELF, KFID, PFIELD)

CLASS(FIELD_CONTAINER_SP), INTENT(IN)    :: SELF
INTEGER(KIND=JPIM),        INTENT(IN)    :: KFID ! field ID
REAL(KIND=JPRB), POINTER,  INTENT(INOUT) :: PFIELD(:,:)

INTEGER(KIND=JPIM) :: IFIDS(1)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_CLOSE_ONE_2D',0,ZHOOK_HANDLE)
IFIDS(1) = KFID

NULLIFY(PFIELD)
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FIELD_CLOSE_ONE_2D',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_CLOSE_ONE_2D
!=======================================================================================
SUBROUTINE FMAP_1D(SELF, KFID, PFIELD)
CLASS(FIELD_CONTAINER_SP),TARGET,INTENT(IN) :: SELF
INTEGER(KIND=JPIM),           INTENT(IN)    :: KFID ! field ID
REAL(KIND=JPRB), POINTER, INTENT(OUT)       :: PFIELD(:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FMAP_1D',0,ZHOOK_HANDLE)
IF(SELF%METADATA(KFID)%NDIMS==1) THEN
  IF(SELF%GET_STORE_TYPE() == JP_STORAGE_ARR) THEN
   PFIELD => SELF%STORAGE_ARRAYS(KFID)%STORE1D
  ENDIF
ENDIF
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FMAP_1D',1,ZHOOK_HANDLE)
END SUBROUTINE FMAP_1D
!=======================================================================================
SUBROUTINE FMAP_2D(SELF, KFID, PFIELD)
CLASS(FIELD_CONTAINER_SP),TARGET,INTENT(IN) :: SELF
INTEGER(KIND=JPIM),           INTENT(IN)    :: KFID ! field ID
REAL(KIND=JPRB), POINTER, INTENT(OUT)       :: PFIELD(:,:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FMAP_2D',0,ZHOOK_HANDLE)
IF(SELF%METADATA(KFID)%NDIMS==2) THEN
  IF(SELF%GET_STORE_TYPE() == JP_STORAGE_ARR) THEN
   PFIELD => SELF%STORAGE_ARRAYS(KFID)%STORE2D
  ENDIF
ENDIF
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:FMAP_2D',1,ZHOOK_HANDLE)
END SUBROUTINE FMAP_2D
!=======================================================================================
SUBROUTINE ZERO(SELF)
CLASS(FIELD_CONTAINER_SP),        INTENT(INOUT) :: SELF
INTEGER(KIND=JPIM) :: ID
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:ZERO',0,ZHOOK_HANDLE)
IF(SELF%GET_STORE_TYPE() == JP_STORAGE_ARR) THEN
  DO ID=1,SELF%NFIELDS
    IF (SELF%FIELD_INDEX(ID)%L_ON) THEN
      IF (SELF%METADATA(ID)%NDIMS == 0) THEN
        SELF%STORAGE_ARRAYS(ID)%STORE1D(:) = 0.0_JPRB
      ELSEIF (SELF%METADATA(ID)%NDIMS == 1) THEN
        SELF%STORAGE_ARRAYS(ID)%STORE1D(:) = 0.0_JPRB
      ELSEIF (SELF%METADATA(ID)%NDIMS == 2) THEN
        SELF%STORAGE_ARRAYS(ID)%STORE2D(:,:) = 0.0_JPRB
      ELSEIF (SELF%METADATA(ID)%NDIMS == 3) THEN
        SELF%STORAGE_ARRAYS(ID)%STORE3D(:,:,:) = 0.0_JPRB
      ENDIF
    ENDIF
  ENDDO
ELSEIF(SELF%GET_STORE_TYPE() == JP_STORAGEOLD) THEN
ENDIF
IF (LHOOK) CALL DR_HOOK('FIELD_CONTAINER_SP_MOD:ZERO',1,ZHOOK_HANDLE)
END SUBROUTINE ZERO
!=======================================================================================
END MODULE FIELD_CONTAINER_SP_MOD
