MODULE YOMRADF

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE

! EMTD    - longwave net flux
! TRSW    - shortwave net transmissivity (multiply by incoming SW to get flux)
! EMTC    - clear-sky net longwave flux
! TRSC    - clear-sky net shortwave transmissivity

! TAUAER  - prognostic aerosol variable for radiation and clouds

! SRSWD   - downward SW radiation at the surface 
! SRLWD   - downward LW radiation at the surface
! SRLWDC  - clear-sky downward LW radiation at the surface
! SRSWDC  - clear-sky downward SW radiation at the surface
! SRSWDCS - clear-sky NET SW radiation at the surface
! SRLWDCS - clear-sky NET LW radiation at the surface
! SRSWDV  - downward SW visible radiation at the surface
! SRSWDUV - downward SW ultraviolet/visible radiation at the surface
! SRSWPAR - downward SW PAR radiation at the surface
! SRSWUVB - downward UV-B radiation at the surface
! SRSWPARC- downward clear-sky SW PAR radiation at the surface
! SRSWTINC- TOA incident solar radiation 
! RMOON   - M-F military application
! SRSWTINC- TOA incident solar radiation
! SRFDIR  - total sky direct downward SW radiation
! SRCDIR  - clear-sky direct downward SW radiation
! DerivativeLw - derivative to update LW radiation between calls to full radiation scheme

! This type is used to pass radiative fluxes from the radiation scheme
! (see RADDRV) to the heating-rate calculation (see RADHEATN), and
! the arrays are allocated in SUECRAD
TYPE :: TRADF
! Dimesioned NPROMA,NFLEVG+1,NGPBLKS:
REAL(KIND=JPRB),ALLOCATABLE :: EMTD(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: TRSW(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: EMTC(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: TRSC(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: EMTU(:,:,:)

! Dimesioned NPROMA,NFLEVG,6,NGPBLKS:
REAL(KIND=JPRB),ALLOCATABLE :: TAUAER(:,:,:,:)

! Dimesioned NPROMA,NGPBLKS:
REAL(KIND=JPRB),ALLOCATABLE :: SRSWD(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRLWDC(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRLWD(:,:,:) ! NPROMA,NLWOUT,NGPBLKS
REAL(KIND=JPRB),ALLOCATABLE :: SRSWDC(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWDCS(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRLWDCS(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWDV(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWDUV(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: EDRO(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWPAR(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWUVB(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWPARC(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRSWTINC(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRFDIR(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SRCDIR(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RMOON(:,:)
! The following is a pointer rather than an allocatable in order that
! pointers to it may be passed to subroutines (specifically to
! RADHEATN) even if it is not associated with any data, which occurs
! if LApproxLwUpdate is set to .FALSE.
REAL(KIND=JPRB),POINTER      :: DERIVATIVELW(:,:,:)=>NULL()
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TRADF
!============================================================================

!!TYPE(TRADF), POINTER :: YRRADF => NULL()

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TRADF), INTENT(IN) :: SELF
INTEGER     , INTENT(IN) :: KDEPTH
INTEGER     , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMRADF:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)
IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_rad%yrradf : '
IF (ALLOCATED(SELF%EMTD)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'EMTD ALLOCATED OF SHAPE ', SHAPE(SELF%EMTD), ' SUM ',SUM(SELF%EMTD)
IF (ALLOCATED(SELF%TRSW)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TRSW ALLOCATED OF SHAPE ', SHAPE(SELF%TRSW), ' SUM ',SUM(SELF%TRSW)
IF (ALLOCATED(SELF%EMTC)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'EMTC ALLOCATED OF SHAPE ', SHAPE(SELF%EMTC), ' SUM ',SUM(SELF%EMTC)
IF (ALLOCATED(SELF%TRSC)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TRSC ALLOCATED OF SHAPE ', SHAPE(SELF%TRSC), ' SUM ',SUM(SELF%TRSC)
IF (ALLOCATED(SELF%EMTU)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'EMTU ALLOCATED OF SHAPE ', SHAPE(SELF%EMTU), ' SUM ',SUM(SELF%EMTU)
IF (ALLOCATED(SELF%TAUAER)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TAUAER ALLOCATED OF SHAPE ', SHAPE(SELF%TAUAER), ' SUM ',SUM(SELF%TAUAER)
IF (ALLOCATED(SELF%SRSWD)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRSWD ALLOCATED OF SHAPE ', SHAPE(SELF%SRSWD), ' SUM ',SUM(SELF%SRSWD)
IF (ALLOCATED(SELF%SRLWDC)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRLWDC ALLOCATED OF SHAPE ', SHAPE(SELF%SRLWDC), ' SUM ',SUM(SELF%SRLWDC)
IF (ALLOCATED(SELF%SRLWD)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRLWD ALLOCATED OF SHAPE ', SHAPE(SELF%SRLWD), ' SUM ',SUM(SELF%SRLWD)
IF (ALLOCATED(SELF%SRSWDC)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRSWDC ALLOCATED OF SHAPE ', SHAPE(SELF%SRSWDC), ' SUM ',SUM(SELF%SRSWDC)
IF (ALLOCATED(SELF%SRSWDCS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRSWDCS ALLOCATED OF SHAPE ', SHAPE(SELF%SRSWDCS), ' SUM ',SUM(SELF%SRSWDCS)
IF (ALLOCATED(SELF%SRLWDCS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRLWDCS ALLOCATED OF SHAPE ', SHAPE(SELF%SRLWDCS), ' SUM ',SUM(SELF%SRLWDCS)
IF (ALLOCATED(SELF%SRSWDV)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRSWDV ALLOCATED OF SHAPE ', SHAPE(SELF%SRSWDV), ' SUM ',SUM(SELF%SRSWDV)
IF (ALLOCATED(SELF%SRSWDUV)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRSWDUV ALLOCATED OF SHAPE ', SHAPE(SELF%SRSWDUV), ' SUM ',SUM(SELF%SRSWDUV)
IF (ALLOCATED(SELF%EDRO)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'EDRO ALLOCATED OF SHAPE ', SHAPE(SELF%EDRO), ' SUM ',SUM(SELF%EDRO)
IF (ALLOCATED(SELF%SRSWPAR)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRSWPAR ALLOCATED OF SHAPE ', SHAPE(SELF%SRSWPAR), ' SUM ',SUM(SELF%SRSWPAR)
IF (ALLOCATED(SELF%SRSWUVB)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRSWUVB ALLOCATED OF SHAPE ', SHAPE(SELF%SRSWUVB), ' SUM ',SUM(SELF%SRSWUVB)
IF (ALLOCATED(SELF%SRSWPARC)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRSWPARC ALLOCATED OF SHAPE ', SHAPE(SELF%SRSWPARC), ' SUM ',SUM(SELF%SRSWPARC)
IF (ALLOCATED(SELF%SRSWTINC)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRSWTINC ALLOCATED OF SHAPE ', SHAPE(SELF%SRSWTINC), ' SUM ',SUM(SELF%SRSWTINC)
IF (ALLOCATED(SELF%SRFDIR)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRFDIR ALLOCATED OF SHAPE ', SHAPE(SELF%SRFDIR), ' SUM ',SUM(SELF%SRFDIR)
IF (ALLOCATED(SELF%SRCDIR)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SRCDIR ALLOCATED OF SHAPE ', SHAPE(SELF%SRCDIR), ' SUM ',SUM(SELF%SRCDIR)
IF (ALLOCATED(SELF%RMOON)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RMOON ALLOCATED OF SHAPE ', SHAPE(SELF%RMOON), ' SUM ',SUM(SELF%RMOON)
IF (LHOOK) CALL DR_HOOK('YOMRADF:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOMRADF
