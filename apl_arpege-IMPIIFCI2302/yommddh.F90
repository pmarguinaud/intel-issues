MODULE YOMMDDH

USE PARKIND1 , ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!     DIAGNOSTIQUES DOMAINES HORIZONTAUX
!     ----------------------------------

!     DIMENSIONS DES DIAGNOSTIQUES DDH
!     ---------------------------------
! === basic dimensions for reading NAMDDH ===
INTEGER(KIND=JPIM), PARAMETER :: JPDHNOX=10000
INTEGER(KIND=JPIM), PARAMETER :: JPDHXPU=5

!     JPDHNOX: Maximum number of user area
!     JPDHXPU: Maximum number of user or mask memory plans

TYPE :: TMDDH
INTEGER(KIND=JPIM) :: NDHKD
INTEGER(KIND=JPIM) :: NDHNPU
INTEGER(KIND=JPIM) :: NDHBPU(JPDHXPU)
INTEGER(KIND=JPIM) :: NDHBPX
INTEGER(KIND=JPIM) :: NDHNOM
INTEGER(KIND=JPIM) :: NDHDDX
INTEGER(KIND=JPIM) :: NDHIDH
INTEGER(KIND=JPIM) :: NDHCS
INTEGER(KIND=JPIM) :: NDHCV
INTEGER(KIND=JPIM) :: NDHCVSU
INTEGER(KIND=JPIM) :: NDHCSSU
INTEGER(KIND=JPIM) :: NDHCVSUN
INTEGER(KIND=JPIM) :: NDHCVSUL
!     NDHKD : NOMBRE DE BANDES DE LATITUDES
!     NDHNPU : NOMBRE DE MASQUES UTILISES NDHNPU .LE. JPDHXPU
!     NDHBPU(JPDHXPU) : NOMBRE DE DOMAINES CONTENU DANS CHAQUE MASQUE
!     NDHBPX : NOMBRE DE DOMAINE MAXIMAL PAR MASQUE
!                                         NDHBPX = MAX( NDHBPU )
!     NDHNOM : NOMBRE TOTALE DE DOMAINES LIMITES EN NOMENCLATURE
!                                         NDHNOM .LE. JPDHNOX
!     NDHDDX : NOMBRE MAXIMAL DE DOMAINES INTERNES POSSIBLES
!     NDHIDH : NOMBRE DE DOMAINES INTERNES ( OR DOMAINE 0 EVENTUEL )
!              INITIALISE PAR CALCUL DANS SUMDDH
!     NDHCS  : NOMBRE TOTAL DE CHAMPS AU SOL
!     NDHCV  : NOMBRE TOTAL DE CHAMPS EN PROFILS VERTICAUX
!     NDHCVSU: DIMENSION TABLEAUX LOCAUX DE CHAMPS VERTICAUX
!              NDHCVSU EST IMPAIR ET VAUT 1 AU MOINS
!     NDHCSSU: DIMENSION TABLEAUX LOCAUX DE CHAMPS AU SOL

!     INDICATEURS D ORGANISATION DE LA REPARTITION DES CHAMPS
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: NDHVV
INTEGER(KIND=JPIM) :: NDHFVD
INTEGER(KIND=JPIM) :: NDHFVP
INTEGER(KIND=JPIM) :: NDHVS
INTEGER(KIND=JPIM) :: NDHFSD
INTEGER(KIND=JPIM) :: NDHFSP
INTEGER(KIND=JPIM) :: NDHFFS
INTEGER(KIND=JPIM) :: NDHVFS
INTEGER(KIND=JPIM) :: NFSVAR_AERO
INTEGER(KIND=JPIM) :: NFSFLX_AERO

INTEGER(KIND=JPIM) :: NDHVTLS
INTEGER(KIND=JPIM) :: NDHFTLS
INTEGER(KIND=JPIM) :: NDHVTSS
INTEGER(KIND=JPIM) :: NDHFTSS
INTEGER(KIND=JPIM) :: NDHVTTS
INTEGER(KIND=JPIM) :: NDHFTTS
INTEGER(KIND=JPIM) :: NDHVTIS
INTEGER(KIND=JPIM) :: NDHFTIS
INTEGER(KIND=JPIM) :: NDHVSSS
INTEGER(KIND=JPIM) :: NDHFSSS
INTEGER(KIND=JPIM) :: NDHVIIS
INTEGER(KIND=JPIM) :: NDHFIIS
INTEGER(KIND=JPIM) :: NDHVWLS
INTEGER(KIND=JPIM) :: NDHFWLS

INTEGER(KIND=JPIM) :: NDHTHK
INTEGER(KIND=JPIM) :: NDHVHK
INTEGER(KIND=JPIM) :: NDHFHKD
INTEGER(KIND=JPIM) :: NDHFHKP
INTEGER(KIND=JPIM) :: NDHTMC
INTEGER(KIND=JPIM) :: NDHVMC
INTEGER(KIND=JPIM) :: NDHFMCD
INTEGER(KIND=JPIM) :: NDHFMCP
INTEGER(KIND=JPIM) :: NDHTEN
INTEGER(KIND=JPIM) :: NDHVEN
INTEGER(KIND=JPIM) :: NDHFEND
INTEGER(KIND=JPIM) :: NDHFENP
INTEGER(KIND=JPIM) :: NDHAVD
INTEGER(KIND=JPIM) :: NDHBVD
INTEGER(KIND=JPIM) :: NDHAVP
INTEGER(KIND=JPIM) :: NDHBVP
INTEGER(KIND=JPIM) :: NDHAHKD
INTEGER(KIND=JPIM) :: NDHBHKD
INTEGER(KIND=JPIM) :: NDHAHKP
INTEGER(KIND=JPIM) :: NDHBHKP
INTEGER(KIND=JPIM) :: NDHAMCD
INTEGER(KIND=JPIM) :: NDHBMCD
INTEGER(KIND=JPIM) :: NDHAMCP
INTEGER(KIND=JPIM) :: NDHBMCP
INTEGER(KIND=JPIM) :: NDHAEND
INTEGER(KIND=JPIM) :: NDHBEND
INTEGER(KIND=JPIM) :: NDHAENP
INTEGER(KIND=JPIM) :: NDHBENP

!   * LES INDICATEURS PRECEDES DE * SONT DES NOMBRES ATTACHES
!     AU CODE FORTRAN DE CALCUL SCIENTIFIQUE, ET NON DES
!     PARAMETRES AJUSTABLES AUX CIRCONSTANCES

!     NDHVV  : NOMBRE DE VARIABLES SUR LA VERTICALE
!     NDHFVD : NOMBRE DE FLUX+TENDANCES DYNAMIQUES SUR LA VERTICALE
!           DETAIL CI DESSOUS
!     NDHFVP : NOMBRE DE FLUX+TENDANCES PHYSIQUES SUR LA VERTICALE
!           DETAIL CI DESSOUS

!   * NDHVS  : NOMBRE DE VARIABLES AU SOL
!   * NDHFSD : NOMBRE DE FLUX DYNAMIQUES AU SOL
!   * NDHFSP : NOMBRE DE FLUX PHYSIQUES AU SOL
!   * NDHFFS : Number of free-style fluxes
!   * NDHVFS : Number of free-style variables
!   * NFSVAR_AERO : Index of first aerosol variable in the list of free-style vars.
!   * NFSFLX_AERO : Index of first aerosol flux in the list of free-style fluxes

!   * NDHVTLS: Number of variables for individual tiles
!   * NDHFTLS: Number of fluxes for individual tiles
!   * NDHVTSS: Number of variables for snow energy budget
!   * NDHFTSS: Number of fluxes for snow energy budget
!   * NDHVTTS: Number of variables for soil energy budget
!   * NDHFTTS: Number of fluxes for soil energy budget
!   * NDHVTIS: Number of variables for sea ice energy budget
!   * NDHFTIS: Number of fluxes for sea ice energy budget
!   * NDHVSSS: Number of variables for snow water budget
!   * NDHFSSS: Number of fluxes for snow water budget
!   * NDHVIIS: Number of variables for interception layer water budget
!   * NDHFIIS: Number of fluxes for interception layer water budget
!   * NDHVWLS: Number of variables for soil water budget
!   * NDHFWLS: Number of fluxes for soil water budget

!     NDHTHK : NOMBRE TOTAL DE CHAMPS VERTICAUX SOUS LHDHKS
!   * NDHVHK : NOMBRE DE VARIABLES SOUS LHDHKS
!   * NDHFHKD: NOMBRE DE FLUX+TENDANCES DYNAMIQUES SOUS LHDHKS
!   * NDHFHKP: NOMBRE DE FLUX+TENDANCES PHYSIQUES SOUS LHDHKS

!     NDHTMC : NOMBRE TOTAL DE CHAMPS VERTICAUX SOUS LHDMCI
!   * NDHVMC : NOMBRE DE VARIABLES SOUS LHDMCI
!   * NDHFMCD: NOMBRE DE FLUX+TENDANCES DYNAMIQUES SOUS LHDMCI
!   * NDHFMCP: NOMBRE DE FLUX+TENDANCES PHYSIQUES SOUS LHDMCI

!     NDHTEN : NOMBRE TOTAL DE CHAMPS SOUS LHDENT
!   * NDHVEN : NOMBRE DE VARIABLES SOUS LHDENT
!   * NDHFEND: NOMBRE DE FLUX+TENDANCES DYNAMIQUES SOUS LHDENT
!   * NDHFENP: NOMBRE DE FLUX+TENDANCES PHYSIQUES SOUS LHDENT

!     NDHAVD : NOMBRE TOTAL DE TENDANCES DYNAMIQUES
!     NDHBVD : NOMBRE TOTAL DE FLUX DYNAMIQUES
!     NDHAVP : NOMBRE TOTAL DE FLUX PHYSIQUES
!     NDHBVP : NOMBRE TOTAL DE TENDANCES PHYSIQUES
!   * NDHAHKD : NOMBRE TOTAL DE TENDANCES DYNAMIQUES, OPTION LHDHKS
!   * NDHBHKD : NOMBRE TOTAL DE FLUX DYNAMIQUES
!   * NDHAHKP : NOMBRE TOTAL DE FLUX PHYSIQUES
!   * NDHBHKP : NOMBRE TOTAL DE TENDANCES PHYSIQUES
!   * NDHAMCD : NOMBRE TOTAL DE TENDANCES DYNAMIQUES, OPTION LHDMCI
!   * NDHBMCD : NOMBRE TOTAL DE FLUX DYNAMIQUES
!   * NDHAMCP : NOMBRE TOTAL DE FLUX PHYSIQUES
!   * NDHBMCP : NOMBRE TOTAL DE TENDANCES PHYSIQUES
!   * NDHAEND : NOMBRE TOTAL DE TENDANCES DYNAMIQUES, OPTION LHDENT
!   * NDHBEND : NOMBRE TOTAL DE FLUX DYNAMIQUES
!   * NDHAENP : NOMBRE TOTAL DE FLUX PHYSIQUES
!   * NDHBENP : NOMBRE TOTAL DE TENDANCES PHYSIQUES

!     IDENTIFICATION DES DOMAINES DE L UTILISATEUR
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: NDHZPR
REAL(KIND=JPRB) :: FNODDH(11,JPDHNOX)
REAL(KIND=JPRB) :: BDEDDH(10,JPDHNOX)
REAL(KIND=JPRB) :: HDSFGL

!     FNODDH : DESCRIPTEURS DES DOMAINES LIMITES UTILISATEUR ACCEPTES
!           (1,-) : PLAN MEMOIRE (MASQUE)
!           (2,-) : NUMERO D ORDRE DANS LE PLAN
!          CES DEUX NOMBRES FORMENT LES COORDONNEES INTERNES DU DOMAINE
!           (3,-) A (10,-) : COORDONNEES (EN GAL LON, SIN(LAT) EN RD)
!                            SUIVANT LE TYPE (VOIR DOCUMENTATION)
!           (11,-) : TYPE DE DOMAINE

!     BDEDDH : DESCRIPTEUR BRUT LU SUR NAMELIST
!              (CERTAINS DOMAINES PEUVENT AVOIR ETE REJETES)
!           (1,-) : TYPE DE DOMAINE
!           (2,-) : PLAN MEMOIRE
!           (3,-) A (10,-) : COORDONNEES (EN GAL LON, LAT EN DEG)
!                            SUIVANT LE TYPE (VOIR DOCUMENTATION)
!     HDSFGL : VOIR CI DESSOUS
!     NDHZPR : INDICE DE LA BANDE ZONALE EVENTUELLEMENT IMPRIMEE

!     MASQUES DES DIAGNOSTIQUES DDH, POIDS POUR LES MOYENNES
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM),ALLOCATABLE:: NDDHLA(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NDDHPU(:,:)

INTEGER(KIND=JPIM),ALLOCATABLE:: NDDHI(:)
!     -----
INTEGER(KIND=JPIM),ALLOCATABLE:: NLRDDH(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NURDDH(:,:,:)
!     -----
INTEGER(KIND=JPIM),ALLOCATABLE:: NLXDDH(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NUXDDH(:,:)
!     -----
REAL(KIND=JPRB),ALLOCATABLE:: HDSFLA(:)
REAL(KIND=JPRB),ALLOCATABLE:: HDSFDU(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: HDSF(:)

!     NDDHLA : MASQUE BANDES DE LATITUDES
!     NDDHPU : MASQUES POUR LES DOMAINES LIMITES ET LES POINTS

!     NDDHI : MASQUE DES DOMAINES INTERNES

!     NLRDDH : DISTRIBUTION DE CHAQUE BANDE DE LATITUDE EN DOMAINES INTERNES
!     NURDDH : DISTRIBUTION DE CHAQUE DOMAINE DE CHAQUE PLAN EN DOMAINES INTERNE

!     NLXDDH : BORNE DE LECTURE DE NLRDDH POUR CHAQUE BANDE
!     NUXDDH : BORNE DE LECTURE DE NURDDH POUR CHAQUE DOMAINE DE CHAQUE PLAN

!     HDSFLA : POIDS DE CHAQUE BANDE DE LATITUDE
!     HDSFDU : POIDS DE CHAQUE DOMAINE UTILISATEUR
!     HDSFGL : POIDS POUR LES MOYENNES GLOBALES
!     HDSF : POIDS DE CHAQUE POINT DE LA GRILLE PHYSIQUE

!     POINTEURS

! ---- DDH

!     ------------------------------------------------------------------

! ----DIMENSIONS OF THE GMV AND GFL ARRAYS AND POINTERS
!     FOR DYNAMICAL SI AND HD DDH TENDENCIES

INTEGER(KIND=JPIM) :: NDIMHDGFL
INTEGER(KIND=JPIM) :: NDIMSIGMV

INTEGER(KIND=JPIM) :: MSIDDH_U1
INTEGER(KIND=JPIM) :: MSIDDH_V1
INTEGER(KIND=JPIM) :: MSIDDH_T1
INTEGER(KIND=JPIM) :: MSIDDH_PD1
INTEGER(KIND=JPIM) :: MSIDDH_VD1

INTEGER(KIND=JPIM) :: MSIDDH_U0
INTEGER(KIND=JPIM) :: MSIDDH_V0
INTEGER(KIND=JPIM) :: MSIDDH_T0
INTEGER(KIND=JPIM) :: MSIDDH_PD0
INTEGER(KIND=JPIM) :: MSIDDH_VD0

INTEGER(KIND=JPIM) :: MSIDDH_U9
INTEGER(KIND=JPIM) :: MSIDDH_V9
INTEGER(KIND=JPIM) :: MSIDDH_T9
INTEGER(KIND=JPIM) :: MSIDDH_PD9
INTEGER(KIND=JPIM) :: MSIDDH_VD9

INTEGER(KIND=JPIM) :: MHDDDH_U
INTEGER(KIND=JPIM) :: MHDDDH_V
INTEGER(KIND=JPIM) :: MHDDDH_T
INTEGER(KIND=JPIM) :: MHDDDH_Q
INTEGER(KIND=JPIM) :: MHDDDH_PD
INTEGER(KIND=JPIM) :: MHDDDH_VD
INTEGER(KIND=JPIM) :: MHDDDH_NHX

! pointers for DDH arrays used in the semi-Lagrangian scheme
INTEGER(KIND=JPIM) :: MSLDDH_U
INTEGER(KIND=JPIM) :: MSLDDH_V
INTEGER(KIND=JPIM) :: MSLDDH_T
INTEGER(KIND=JPIM) :: MSLDDH_PD
INTEGER(KIND=JPIM) :: MSLDDH_VD
INTEGER(KIND=JPIM) :: MSLDDH_NHX

!     MTNDSI_DDH_U : POINTEUR SUR LA COMPOSANTE-U DU VENT DANS LE TABLEAU 
!                    CONTENANT LES TENDANCES SEMI-IMPLICITES 
!     MTNDSI_DDH_V : POINTEUR SUR LA COMPOSANTE-V DU VENT DANS LE TABLEAU 
!                    CONTENANT LES TENDANCES SEMI-IMPLICITES
!     MTNDSI_DDH_T : POINTEUR SUR LA TEMPERATURE DANS LE TABLEAU 
!                    CONTENANT LES TENDANCES SEMI-IMPLICITES
!     MTNDHD_DDH_U : POINTEUR SUR LA COMPOSANTE-U DU VENT DANS LE TABLEAU GMV
!                    CONTENANT LES TENDANCES LIEES A LA DIFFUSION HORIZONTALE
!     MTNDHD_DDH_V : POINTEUR SUR LA COMPOSANTE-V DU VENT DANS LE TABLEAU GMV
!                    CONTENANT LES TENDANCES LIEES A LA DIFFUSION HORIZONTALE
!     MTNDHD_DDH_T : POINTEUR SUR LA TEMPERATURE DANS LE TABLEAU GMV
!                    CONTENANT LES TENDANCES LIEES A LA DIFFUSION HORIZONTALE
!     MTNDHD_DDH_Q : POINTEUR SUR L'HUMIDITE SPECIFIQUE DANS LE TABLEAU GFL
!                    CONTENANT LES TENDANCES LIEES A LA DIFFUSION HORIZONTALE
!     NTNDHD_GFL   : NOMBRE TOTAL DE  GFL TRAITES EN SPECTRAL ET DIFFUSSES

! CFPATHDDH   : directory name (ending with a '/') of DDH output files

CHARACTER(LEN=200) :: CFPATHDDH

INTEGER(KIND=JPIM) :: NFIELDS3D_AUTO
INTEGER(KIND=JPIM) :: NFIELDS3D_OFFSET
INTEGER(KIND=JPIM) :: NFIELDSMAX
INTEGER(KIND=JPIM) :: NFIELDS2D_AUTO
INTEGER(KIND=JPIM) :: NFIELDS2D_OFFSET


END TYPE TMDDH

!!TYPE(TMDDH), POINTER :: YRMDDH => NULL()

!     ------------------------------------------------------------------
END MODULE YOMMDDH
