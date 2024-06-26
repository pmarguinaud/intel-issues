MODULE YOMTOPH

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE TTOPH
!*
!     ------------------------------------------------------------------
!     Top limits of parametrization call
!     we have one by parameterization,
!              ETXXXX : Top pressure
!              NTXXXX : Corresponding level in standard atmosphere

!       ETQSAT,ETDIFU,ETCOEF,ETDRAG,ETCVIM,ET850,ET950,ETPLUI,ETRADI,ETNEBU
!       ETOZON,ETDRME,ETCOET,ETAJUC,NTAJUC
!       NTQSAT,NTDIFU,NTCOEF,NTDRAG,NTCVIM,NT850,NT950,NTPLUI,NTRADI,NTNEBU
!       NTOZON,NTDRME,NTCOET
!
!       1D model MUSC Same bottom pressure for the NUDGING level 
!               for T, Qv and Wind :
!       ETRELAX, NTRELAX

!     Constants for Mesospheric drag parameterization

!       XDRMUK  : Mesospheric drag coefficient for U and V
!       XDRMUX  : Maximum mesospheric drag coefficient for U and V
!       XDRMUP  : Bottom pressure for U and V
!       XDRMTK  : Mesospheric drag coefficient for T
!       XDRMTX  : Maximum mesospheric drag coefficient for T
!       XDRMTP  : Bottom pressure for T
!       XDRMQK  : Mesospheric drag coefficient for Q
!       XDRMQP  : Bottom pressure for Q
!       XDRMQX  : Maximum mesospheric drag coefficient for Q

!       RMESOU(NFLEVG) : Vertical mesospheric drag profil for U and V
!       RMESOT(NFLEVG) : Vertical mesospheric drag profil for T
!       RMESOQ(NFLEV)  : Vertical mesospheric drag profil for Q

!       RFMESOQ   : Mesospheric water reference

!       RCLX      : Multiplicator for Cl term in ozone photochemistry

!       TPSCLIM   : Temperature threshold for activation of heterogeneous
!                   chemistry (polar stratospheric clouds temperature formation)

!       ETELEM    : Molecule conservation is ensured at pressures below ETELEM (arpclim_repro and arpclim_relacs chemistry)
!       ETRELA    : Concentrations are relaxed at pressures higher than ETRELA (arpclim_repro and arpclim_relacs chemistry)
!       ETPSCT    : Heterogenous chemistry is computed between the ETPSCT and ETPSCB pressures (arpclim_repro and arpclim_relacs chemistry)
!       ETPSCB    :
!       ETHOXF    : HOx as a family at pressures below ETHOXF (arpclim_repro and arpclim_relacs chemistry)
!       ETCHEMNO  : Vertical profiles imposed to chemical species at pressures below ETCHEMNO (arpclim_repro and arpclim_relacs chemistry)

REAL(KIND=JPRB),ALLOCATABLE:: RMESOU(:)
REAL(KIND=JPRB),ALLOCATABLE:: RMESOT(:)
REAL(KIND=JPRB),ALLOCATABLE:: RMESOQ(:)

!       RUREL    : PROFIL VERTICAL DE RAPPEL DE U
!                   RELAXATION PROFILE FOR U
!       RVREL    : PROFIL VERTICAL DE RAPPEL DE V
!                   RELAXATION PROFILE FOR V
!       RTREL    : PROFIL VERTICAL DE RAPPEL DE T
!                   RELAXATION PROFILE FOR T
!       RQREL    : PROFIL VERTICAL DE RAPPEL DE Q
!                   RELAXATION PROFILE FOR Q

REAL(KIND=JPRB), ALLOCATABLE :: RUREL(:)
REAL(KIND=JPRB), ALLOCATABLE :: RVREL(:)
REAL(KIND=JPRB), ALLOCATABLE :: RTREL(:)
REAL(KIND=JPRB), ALLOCATABLE :: RQREL(:)

REAL(KIND=JPRB) :: RFMESOQ

INTEGER(KIND=JPIM) :: NTQSAT
INTEGER(KIND=JPIM) :: NTDIFU
INTEGER(KIND=JPIM) :: NTCOEF
INTEGER(KIND=JPIM) :: NTDRAG
INTEGER(KIND=JPIM) :: NTCVIM
INTEGER(KIND=JPIM) :: NT850
INTEGER(KIND=JPIM) :: NT950
INTEGER(KIND=JPIM) :: NTPLUI
INTEGER(KIND=JPIM) :: NTRADI
INTEGER(KIND=JPIM) :: NTNEBU
INTEGER(KIND=JPIM) :: NTOZON
INTEGER(KIND=JPIM) :: NTDRME
INTEGER(KIND=JPIM) :: NTCOET
INTEGER(KIND=JPIM) :: NTAJUC
INTEGER(KIND=JPIM) :: NTRELAXT
INTEGER(KIND=JPIM) :: NTRELAXQ
INTEGER(KIND=JPIM) :: NTRELAXU
INTEGER(KIND=JPIM) :: CHEM_HPSC_NTOP  ! arpclim_repro and arpclim_relacs heterogeneous chemistry
INTEGER(KIND=JPIM) :: CHEM_HPSC_NBOT  ! arpclim_repro and arpclim_relacs heterogeneous chemistry
INTEGER(KIND=JPIM) :: CHEM_HOXF_NBOT  ! arpclim_repro and arpclim_relacs chemistry
INTEGER(KIND=JPIM) :: CHEM_ELEM_NTOP  ! arpclim_repro and arpclim_relacs chemistry
INTEGER(KIND=JPIM) :: CHEM_RELA_NIV   ! arpclim_repro and arpclim_relacs chemistry
INTEGER(KIND=JPIM) :: CHEM_NO_NIV     ! arpclim_repro and arpclim_relacs chemistry  level above which the chemistry is not computed (in the mesophere)
REAL(KIND=JPRB) :: ETQSAT
REAL(KIND=JPRB) :: ETDIFU
REAL(KIND=JPRB) :: ETCOEF
REAL(KIND=JPRB) :: ETDRAG
REAL(KIND=JPRB) :: ETCVIM
REAL(KIND=JPRB) :: ET850
REAL(KIND=JPRB) :: ET950
REAL(KIND=JPRB) :: ETPLUI
REAL(KIND=JPRB) :: ETRADI
REAL(KIND=JPRB) :: ETNEBU
REAL(KIND=JPRB) :: ETOZON
REAL(KIND=JPRB) :: ETDRME
REAL(KIND=JPRB) :: ETCOET
REAL(KIND=JPRB) :: ETAJUC
REAL(KIND=JPRB) :: ETRELAXT
REAL(KIND=JPRB) :: ETRELAXQ
REAL(KIND=JPRB) :: ETRELAXU
REAL(KIND=JPRB) :: XDRMUK
REAL(KIND=JPRB) :: XDRMUX
REAL(KIND=JPRB) :: XDRMUP
REAL(KIND=JPRB) :: XDRMTK
REAL(KIND=JPRB) :: XDRMTX
REAL(KIND=JPRB) :: XDRMTP
REAL(KIND=JPRB) :: XDRMQK
REAL(KIND=JPRB) :: XDRMQP
REAL(KIND=JPRB) :: XDRMQX

REAL(KIND=JPRB) :: RCLX

REAL(KIND=JPRB) :: TPSCLIM

REAL(KIND=JPRB) :: ETELEM
REAL(KIND=JPRB) :: ETRELA
REAL(KIND=JPRB) :: ETPSCT 
REAL(KIND=JPRB) :: ETPSCB
REAL(KIND=JPRB) :: ETHOXF
REAL(KIND=JPRB) :: ETCHEMNO

END TYPE TTOPH

!! TYPE(TTOPH), POINTER :: YRTOPH => NULL()


!     ------------------------------------------------------------------
END MODULE YOMTOPH
