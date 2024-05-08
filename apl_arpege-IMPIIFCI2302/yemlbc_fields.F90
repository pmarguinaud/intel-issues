MODULE YEMLBC_FIELDS 

! Purpose :
! -------
!    Forcing a LAM model by another model: part 0C
!    - forcing by lateral boundary conditions
!    - pressure tendency coupling
!    - spectral nudging

! Interface :
! ---------
!    Empty.

! External :
! --------
!    None.

! Method :
! ------
!    See Documentation.

! Reference :
! ---------

! Author :
! ------
!    Bogdan Bochenek - part from ELBC0B
! Original : Sep 2014

! Modifications :
! -------------
! H. Dhouioui (Sep 2017) Renamed from elbc0c_mod.F90
!-----------------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0    , ONLY : LALLOPR, NCONF
USE YOMLUN    , ONLY : NULOUT
USE YOMGMV    , ONLY : TGMV
USE YEMLBC_MODEL , ONLY : TELBC_MODEL
IMPLICIT NONE
SAVE

!=============================================================================

!      1.    DECLARATIONS
!            ------------


TYPE :: TELBC_FIELDS 


!      1.2   Other variables for grid-point coupling.

! GMVCPL, GMVSCPL, GFLCPL: buffers containing the LBC for GMV, GMVS, GFL
! GMVUCPL, GMVSUCPL, GFLUCPL: buffers containing the UBC for GMV, GMVS, GFL
! GMVSTENC       : cf. GMVSCPL but for LTENC.
! MGMV0          : index for "GMV" memory transfers between GMV and coupling buffers.
! MGMV1          : index for "GMV" memory transfers between GMVT1 and coupling buffers.
! MGMVS0         : index for "GMVS" memory transfers between GMVS and coupling buffers.
! MGMVS1         : index for "GMVS" memory transfers between GMVT1S and coupling buffers.
! CCFIELD_GMV, CCFIELD_GMVS, CCFIELD_GFL: fields names for EWRLSGRAD.

REAL(KIND=JPRB),ALLOCATABLE:: GMVCPL(:,:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: GMVSCPL(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: GFLCPL(:,:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: GMVUCPL(:,:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: GMVSUCPL(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: GFLUCPL(:,:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: GMVSTENC(:,:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: MGMV0(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: MGMV1(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: MGMVS0(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: MGMVS1(:)
CHARACTER(LEN=12),ALLOCATABLE:: CCFIELD_GMV(:)
CHARACTER(LEN=12),ALLOCATABLE:: CCFIELD_GMVS(:)
CHARACTER(LEN=12),ALLOCATABLE:: CCFIELD_GFL(:)

!      1.3   Other variables for spectral nudging.

! GT3SPBUF   : buffer for spectral boundary fields
! NOFFGT3BSP : Supposed to be half the size of GT3SPBUF (which contains 2 sets of fields)

REAL(KIND=JPRB),ALLOCATABLE :: GT3SPBUF(:)
INTEGER(KIND=JPIM) :: NOFFGT3BSP

END TYPE TELBC_FIELDS 


!=============================================================================

CONTAINS

!      2.    SET-UP

SUBROUTINE SUELBC_FIELDS(YDELBC_FIELDS,YDML_LBC,YDGEOMETRY,YDDYNA,YDGMV,YGFL,KNFD2D,KNS3D)




USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOM_YGFL     , ONLY : TYPE_GFLD
USE YOMDYNA      , ONLY : TDYNA
TYPE(TELBC_FIELDS)     , INTENT(INOUT), POINTER :: YDELBC_FIELDS
TYPE(TELBC_MODEL)     , INTENT(INOUT) :: YDML_LBC
TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
TYPE(TDYNA)       , INTENT(IN)    :: YDDYNA
TYPE(TGMV)        , INTENT(INOUT) :: YDGMV
TYPE(TYPE_GFLD)   , INTENT(INOUT) :: YGFL
INTEGER(KIND=JPIM), INTENT(IN)    :: KNFD2D
INTEGER(KIND=JPIM), INTENT(IN)    :: KNS3D








END SUBROUTINE SUELBC_FIELDS 

SUBROUTINE DEALLOCATE_ELBC0C(YDELBC_FIELDS)



IMPLICIT NONE
TYPE(TELBC_FIELDS) , INTENT(INOUT) :: YDELBC_FIELDS





















END SUBROUTINE DEALLOCATE_ELBC0C

!=============================================================================

END MODULE YEMLBC_FIELDS 
