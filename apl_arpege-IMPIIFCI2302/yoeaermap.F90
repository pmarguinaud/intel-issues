MODULE YOEAERMAP

USE PARKIND1  ,ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAERMAP* - GEOGRAPHICALLY DEPENDENT PARAMETERS
!     ------------------------------------------------------------------

TYPE :: TEAERMAP
INTEGER(KIND=JPIM) :: NDUSRCP(50)

REAL(KIND=JPRB)    :: RDDUAER(50)
REAL(KIND=JPRB)    :: RDUSRCP(10,2)
!---------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TEAERMAP
!=====================================================================

!!TYPE(TEAERMAP), POINTER :: YREAERMAP => NULL()

!     ------------------------------------------------------------------
! NDUSRCP ! INTEGER ! INDEX FOR DUST SOURCE TYPE
! RDUSRCP : REAL    ! 1: DUST SOURCE-RELATED REFERENCE THRESHOLD SPEED
! RDUSRCP : REAL    ! 2: DUST SOURCE-RELATED REFERENCE PARTICLE RADIUS
! RDDUAER ! REAL    ! DUST EMISSION POTENTIAL SCALING FACTOR
!     ------------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TEAERMAP), INTENT(IN) :: SELF
INTEGER     , INTENT(IN) :: KDEPTH
INTEGER     , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_aer%yreaermap : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDUSRCP sum = ', SUM(SELF%NDUSRCP)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDDUAER sum = ', SUM(SELF%RDDUAER)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RDUSRCP sum = ', SUM(SELF%RDUSRCP)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOEAERMAP



