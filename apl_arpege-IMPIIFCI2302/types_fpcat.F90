MODULE TYPES_FPCAT

! Purpose :
! -------
!    To define the user type "TYPES_FPCAT" for the post-processing of individual CAT/MTW indices
!    and the final combinaison of indices to compute diagnostic EDR.

!     - %CLNAME        :  the  indice name,
!     - %ZDEL          :  Delta for EDR transpose
!     - %ZRAT          :  Range for EDR transpose
!     - %LCONVEDR      :  Indicate conversion into EDR unit when indice is post-processed
!     - %LCOMB         :  Used for final EDR diag combinaison for each classes or layer

!    Also define a super-type 'ALL_FPCAT_TYPES' with all TYPE_FPCAT objets

!     - %XXX           :  Individual indice of CAT
!      
!     - %LCONVEDRALL   : Global key to control conversion of individual indices into EDR unit when they are post-processed.
!                        When LCONVEDR =.T., if at least one individual indice LCONVEDR is .T.
         
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
!    Fullpos technical & users guide.

! Author :
! ------
!    Olivier Jaron *METEO-FRANCE*

! Modifications :
! -------------
! Original : 2019-04-02

!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE
SAVE

TYPE TYPE_FPCAT

LOGICAL  :: LCONVEDR = .FALSE.
CHARACTER(LEN=16) :: CLNAME   = 'XXXXXXXX        '
REAL(KIND=JPRB)   :: ZDEL(1:3)  = 0._JPRB ! Value of delta in log_log transpose for HIGH/MEDIUM/LOW levels
REAL(KIND=JPRB)   :: ZRAT(1:3)  = 1._JPRB ! Value of ratio in log_log transpose for HIGH/MEDIUM/LOW levels
LOGICAL           :: LCOMB(1:3) = .FALSE. ! Combine the index for HIGH/MEDIUM/LOW levels

END TYPE TYPE_FPCAT

TYPE ALL_FPCAT_TYPES

LOGICAL :: LCONVEDRALL    = .FALSE. 
REAL(KIND=JPRB):: PPRESH  =  7170._JPRB  ! Pression de sommet de la couche HIGH  # FL600
REAL(KIND=JPRB):: PPRESM  = 46569._JPRB  ! Pression de transition MEDIUM -> HIGH # FL200
REAL(KIND=JPRB):: PPRESL  = 69686._JPRB  ! Pression de transition LOW -> MEDIUM  # FL100
REAL(KIND=JPRB):: PPDELTA =  5000._JPRB  ! Epaisseur (Pa) de la transition
REAL(KIND=JPRB):: PPTHLDBVF =10.E-5_JPRB ! Threshold applied to Brunt-Vaissala Frequ.
REAL(KIND=JPRB):: PPTHLDRI  =1.E-3_JPRB  ! Threshold applied to Ri
LOGICAL :: LPSMORI         = .TRUE.      ! Apply a 1-2-1 points vertical mean to Ri
LOGICAL :: LPVMEAN         = .TRUE.      ! Apply vertical means ti 1/Ri and DTheta_virtual/Dz
LOGICAL :: LPFINDIF        = .TRUE.      ! Use vertical finite difference on a non regular grid
REAL(KIND=JPRB):: PHLLX    = 3010._JPRB  ! Maximum of height to compute low levels EDRDC when postprocess height levels

! Individual Diagnostics
TYPE(TYPE_FPCAT) :: BR1        ! Brown1       : Brown1 
TYPE(TYPE_FPCAT) :: BR2        ! Brown2       : Brown 2
TYPE(TYPE_FPCAT) :: DUT        ! Dutton       : Dutton
TYPE(TYPE_FPCAT) :: LAZ        ! LAZ          : Laikthman and Alter-Zalik TKE
TYPE(TYPE_FPCAT) :: TI1        ! Ellrod1      : Ellrod 1
TYPE(TYPE_FPCAT) :: TI2        ! Ellrod2      : Ellrod 2
TYPE(TYPE_FPCAT) :: CP         ! CP           : Colson-Panofsky TKE
TYPE(TYPE_FPCAT) :: RI         ! Ri           : Richardson Number
TYPE(TYPE_FPCAT) :: INVRI      ! 1/Ri         : Inverse Ri 
TYPE(TYPE_FPCAT) :: INVRIM     ! 1/Ri_m       : Inverse Ri moist
TYPE(TYPE_FPCAT) :: HS         ! HS           : | Horizontal shear of horizontal wind |      
TYPE(TYPE_FPCAT) :: LRT        ! LRT          : Lapse rate of temperature
TYPE(TYPE_FPCAT) :: DF         ! DEFSQ        : | Total Deformations | **2 
TYPE(TYPE_FPCAT) :: INVRITW    ! 1/Ri_tw      : Inverse Ri with vertical shear from thermal wind relation
TYPE(TYPE_FPCAT) :: AGI        ! AGI          : Anomalous gradient instability
TYPE(TYPE_FPCAT) :: F2D        ! F2D          : 2D frontogenesis function on z constant 
TYPE(TYPE_FPCAT) :: NG1        ! NGM1         : Wind Speed x | Deformation |
TYPE(TYPE_FPCAT) :: NG2        ! NGM2         : | Vertical gradient temperature| x | Deformation |
TYPE(TYPE_FPCAT) :: TG         ! TEMPG        : Horizontal temperature gradient
TYPE(TYPE_FPCAT) :: IAW        ! IAWIND       : Inertial advective wind
TYPE(TYPE_FPCAT) :: SV         ! VWS          :  Vertical wind shear
TYPE(TYPE_FPCAT) :: BV         ! BV           : Brunt Vaissala frequency
TYPE(TYPE_FPCAT) :: BVM        ! BVM          :Brunt Vaissala frequency (moist)
TYPE(TYPE_FPCAT) :: III        ! Stone        : Inertial instability index (-Stone)
TYPE(TYPE_FPCAT) :: VORSQ      ! VORTSQ       : | Vertical vorticity | **2
TYPE(TYPE_FPCAT) :: VVSQ       ! VSQ          : Vertical Velocity (Omega)  ** 2 
TYPE(TYPE_FPCAT) :: WND        ! Wind         : FF
TYPE(TYPE_FPCAT) :: EDR        ! EDR          : Subgrid EDR
TYPE(TYPE_FPCAT) :: TKE        ! TKE          : Subgrid TKE
!TYPE(TYPE_FPCAT) :: LUN        ! EDRLUN       : Epsilon^1/3 from simplified DRi/Dt

! Diagnostics divided by the Richardson Number
TYPE(TYPE_FPCAT) :: BR1R       ! Brown1/Ri    : Brown1 / Ri 
TYPE(TYPE_FPCAT) :: HSR        ! HS / Ri      : | Horizontal shear of horizontal wind | / Ri 
TYPE(TYPE_FPCAT) :: EDL        ! ENDLICH      : Wind Speed x tuning 
TYPE(TYPE_FPCAT) :: DFR        ! DEFSQ/Ri     : | Total Deformations | **2 / Ri
TYPE(TYPE_FPCAT) :: DVR        ! |DIV|/Ri     : | Horizontal divergence | / Ri
TYPE(TYPE_FPCAT) :: F2DR       ! F2D/Ri       : 2D frontogenesis function on z constant / Ri
TYPE(TYPE_FPCAT) :: NG1R       ! NGM1/Ri      : Wind Speed x | Deformation | / Ri
TYPE(TYPE_FPCAT) :: NG2R       ! NGM2/Ri      : | Vertical gradient temperature| x | Deformation | / Ri
TYPE(TYPE_FPCAT) :: TKER       ! TKE/Ri       : Subgrid TKE / Ri
TYPE(TYPE_FPCAT) :: TGR        ! TEMPG/Ri     : Horizontal temperature gradient / Ri
TYPE(TYPE_FPCAT) :: VORSQR     ! VORTSQ/Ri    : | Vertical vorticity | **2 / Ri
TYPE(TYPE_FPCAT) :: VVSQRI     ! VSQ/Ri       : Vertical Velocity (Omega)  ** 2 / Ri
TYPE(TYPE_FPCAT) :: IAWR       ! IAWIND/Ri    : Inertial advective wind / Ri

! Moutain Waves diagnotics
TYPE(TYPE_FPCAT) :: VVSQMW   
TYPE(TYPE_FPCAT) :: F2DMW
TYPE(TYPE_FPCAT) :: WNDMW
TYPE(TYPE_FPCAT) :: DVRMW
TYPE(TYPE_FPCAT) :: NG1MW
TYPE(TYPE_FPCAT) :: IAWMW
TYPE(TYPE_FPCAT) :: DFMW

! Final CAT/MTW/MAX combinaison
TYPE(TYPE_FPCAT) :: EDRDC

! 2D Diagnostics and max over the vertical
TYPE(TYPE_FPCAT) :: DS
TYPE(TYPE_FPCAT) :: TI1H
TYPE(TYPE_FPCAT) :: TI1M
TYPE(TYPE_FPCAT) :: EDRDCH
TYPE(TYPE_FPCAT) :: EDRDCM
TYPE(TYPE_FPCAT) :: EDRDCL

END TYPE ALL_FPCAT_TYPES
!-----------------------

CONTAINS

SUBROUTINE SUCATIND(YDCAT_X,CDNAME,LDEDR,PDELTA,PRATIO,LDCOMB)

! Purpose :
! -------
!    To set default values to the type TYPE_FPCAT

!     PDELRA  : Delta for the three classes of level
!     LDCOMB  : Indicates if indece should be combined, for the three classes of level
!     LDEDR   : Indicates if indice shloud be converted into EDR units for outputs, following clim. log-normal law

TYPE(TYPE_FPCAT) , INTENT(INOUT)          :: YDCAT_X
CHARACTER(LEN=*), INTENT(IN)              :: CDNAME
LOGICAL , INTENT(IN)           , OPTIONAL :: LDEDR 
REAL(KIND=JPRB) , INTENT(IN)   , OPTIONAL :: PDELTA(3)
REAL(KIND=JPRB) , INTENT(IN)   , OPTIONAL :: PRATIO(3)
LOGICAL , INTENT(IN)           , OPTIONAL :: LDCOMB(3)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TYPES_FPCAT:SUCATIND',0,ZHOOK_HANDLE)

YDCAT_X%CLNAME=CDNAME
IF (PRESENT(LDEDR))  YDCAT_X%LCONVEDR=LDEDR
IF (PRESENT(PDELTA)) YDCAT_X%ZDEL=PDELTA
IF (PRESENT(PRATIO)) YDCAT_X%ZRAT=PRATIO
IF (PRESENT(LDCOMB)) YDCAT_X%LCOMB=LDCOMB

IF (LHOOK) CALL DR_HOOK('TYPES_FPCAT:SUCATIND',1,ZHOOK_HANDLE)

END SUBROUTINE SUCATIND

SUBROUTINE PRINT_CONFIG_EDRD(YDCAT_X, CDNAME)
USE YOMLUN       , ONLY : NULOUT
! Purpose :
! -------
!   Print configuration of EDRD combinaison (coefficients, canditate for final combinaison, 
!   convert or not in EDRD unit)
CHARACTER(LEN=12), INTENT(IN) :: CDNAME
TYPE(TYPE_FPCAT) , INTENT(IN) :: YDCAT_X
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TYPES_FPCAT:PRINT_CONFIG_EDRD',0,ZHOOK_HANDLE)

WRITE(NULOUT, FMT ='(1X,A7,'' : '',A16,1X,6(F9.6,1X), 2X, 3(L1,1X),3X,L1)') &  
  & CDNAME,YDCAT_X%CLNAME,YDCAT_X%ZDEL(1),YDCAT_X%ZDEL(2),YDCAT_X%ZDEL(3),  &
  & YDCAT_X%ZRAT(1) ,YDCAT_X%ZRAT(2) ,YDCAT_X%ZRAT(3),&
  & YDCAT_X%LCOMB(1),YDCAT_X%LCOMB(2),YDCAT_X%LCOMB(3),YDCAT_X%LCONVEDR

IF (LHOOK) CALL DR_HOOK('TYPES_FPCAT:PRINT_CONFIG_EDRD',1,ZHOOK_HANDLE)

ENDSUBROUTINE PRINT_CONFIG_EDRD


END MODULE TYPES_FPCAT
