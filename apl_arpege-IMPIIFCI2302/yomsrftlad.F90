MODULE YOMSRFTLAD

USE PARKIND1  ,ONLY : JPRB,JPIM

IMPLICIT NONE

SAVE

! ------ SKIN TEMPERATURE FOR LINEARIZED PHYSICS

! GPSURF  : BUFFER FOR PERTURBATION OF TOP LAYER SURF. FIELDS 
!           (SKIN TEMPERATURE,...)
! NGSKIN  : NUMBER OF TOP SOIL FIELDS FOR PERTURBATION

! LREGSF  : .TRUE. if the regularization for SURF computation is used

TYPE :: TSRFTLAD
REAL(KIND=JPRB),ALLOCATABLE :: GPTSKIN0(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: GPTSKIN9(:,:,:)

INTEGER(KIND=JPIM) :: NGSKIN

LOGICAL :: LREGSF
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TSRFTLAD
!============================================================================

!!TYPE(TSRFTLAD), POINTER :: YRSRFTLAD => NULL()

!     ------------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TSRFTLAD), INTENT(IN) :: SELF
  INTEGER        , INTENT(IN) :: KDEPTH
  INTEGER        , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH+2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_slin%yrephli : '
  IF (ALLOCATED(SELF%GPTSKIN0)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GPTSKIN0 allocated of shape ', &
 &        SHAPE(SELF%GPTSKIN0),' and sum ',SUM(SELF%GPTSKIN0)
  IF (ALLOCATED(SELF%GPTSKIN9)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GPTSKIN9 allocated of shape ', &
 &        SHAPE(SELF%GPTSKIN9),' and sum ',SUM(SELF%GPTSKIN9)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGSKIN = ', SELF%NGSKIN
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LREGSF = ', SELF%LREGSF
 
END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOMSRFTLAD
