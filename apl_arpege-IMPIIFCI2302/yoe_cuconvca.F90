MODULE YOE_CUCONVCA

! Modified : 
!    R. El Khatib 09-Mar-2012 Allocate RCUCONVCA/RNLCONVCA for safe bound checkings later
!    L. Bengtsson 05-Aug-2014 Correction if using with GOL
!    L. Gerard    31-Mar-2016 Add RCADELX
USE PARKIND1 , ONLY : JPIM, JPRB
USE RANDOM_NUMBERS_MIX, ONLY: RANDOMNUMBERSTREAM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------

!*    control parameters for cellular automaton for convection

!      yd_random_stream_CA: random number stream for CA update
!      NIJH         : Multiplicative for number of points in high resolution grid
!      NBmax        : maximum possible number of neighbours for every gridcell
!      RCADELX      : Length of individual CA cell (meters) to compute NIJH
!      RCUCONVCA    : Array for interaction with the physics
!      RNLCONVCA    : Array for interaction with the physics
!      RCAPECONVCA  : Array to store "good" CAPE between its calculation [MOD(KSTEP*TSPHY,3600._JPRB)]
!      NCELLCU      : Lat-lon grid for cellular automaton
!      NFERTCU      : Lat-lon grid for cellular automaton (fertile cell indicator)
!      NBLS         : neighbours for large scale grid (= model red gaussian grid)
!      NBSS         : neighbours for small scale grid (= NIJH*NIJH cells in every LS cell)
!      NRNBLS       : number of neighbours for large scale cells
!      NRNBSS       : number of neighbours for small scale cells
!      NBNS         : neighbour north/south or on same lat (1/-1/0)
!      NBEW         : neighbour east/west or on same long (1/-1/0)
!      RWASALIVE    : Array indicating grid-points where CA was alive before
!      RWGHTCU      : Weighted and smoothed CA pattern
!      LCUCONV_CA   : TRUE IF CELLULAR AUTOMATON USED IN CONVECTION
!      LCA_ADVECT   : switch for "kind of" semi-Lagrangian advection
!      LCA_GLOBAL   : switch for global CA instead ov "convective" CA
!      LCA_SMOOTH   : swith for smoothing CA on large scale grid
!      LCA_RANTROP  : switch for random coupling of CA to deep convection (concerns initialization)
!      LCA_TEST     : switch for initialize CA at single point
!      LCA_ADVTEST  : switch for advection test (CA only evolved first 5 steps)
!      NTESTPROC    : set on which processor the single point should lie
!      NTESTGP      : set which gridpoint the single point is on
!      NSPINUP      : set number of spin-up cycles for global pattern
!      CA_FORC      : switch to choose convective forcing
!      CA_WIND      : switch to choose CA-wind (real/idealized)
!      CA_PROB      : switch to choose probabilities
!      NLIVES       : switch to choose number of lives (scaled by CAPE in callpar)
!      NFERTYRS     : switch to choose max number of steps a cell can be fertile
!      NFRCASEED    : Frequency of seeding CA where physics diagnoses deep convection
!      RCA_SEEDPROB : Probability of random seeding for global CA

!      RPROB_SURVIVE: probabilities for cell survival
!      RPROB_BIRTH: probabilities for cell birth
!      RPROB_FUN_FERT_S: function for specifying probabilities for survival as function of upwind(uw)+downwind(dw) neighbours
!      RPROB_FUN_UW_DW_S: function for specifying probabilities for survival as function of upwind(uw)-downwind(dw) neighbours
!      RPROB_FUN_FERT_B: function for specifying probabilities for birth as function of upwind(uw)+downwind(dw) neighbours
!      RPROB_FUN_UW_DW_B: function for specifying probabilities for birth as function of upwind(uw)-downwind(dw) neighbours

!      SLLONG(YRSL%NASLB1)       : Array containing the longitude for each gridpoint in a semi/lagrangian buffer array
!      SLLAT(YRSL%NASLB1)        : Array containing the latitude for each gridpoint in a semi/lagrangian buffer array
!      SLDLONG(YRSL%NASLB1)      : Array containing the width of each gridbox in longitudinal direction 
!                             in a semi/lagrangian buffer array
!      SLDLAT(YRSL%NASLB1)       : Array containing the width of each gridbox in latitudinal direction 
!                             in a semi/lagrangian buffer array
!      SLDDLAT(YRSL%NASLB1)      : Array containing the difference of the gridbox's northern boundary latitude 
!                             and the latitude of the gp (necessary because latitudes are not uniformly distributed)

INTEGER(KIND=JPIM)            :: NIJH=4
INTEGER(KIND=JPIM), PARAMETER :: NBMAX=8

TYPE :: TECUCONVCA
TYPE (RANDOMNUMBERSTREAM)     :: YD_RANDOM_STREAM_CA
INTEGER(KIND=JPIM)            :: NLIVES
INTEGER(KIND=JPIM)            :: NFERTYRS
INTEGER(KIND=JPIM)            :: NSPINUP

INTEGER(KIND=JPIM),ALLOCATABLE:: NCELLCU(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NFERTCU(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NBLS(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NBSS(:,:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRNBLS(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NRNBSS(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NBEW(:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NBNS(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RCUCONVCA(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RNLCONVCA(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RCAPECONVCA(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RWASALIVE(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RWGHTCU(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RLATLONNBLS(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RLATLONNBSS(:,:,:,:)

REAL(KIND=JPRB)               :: RPROB_SURVIVE(0:8,0:8)
REAL(KIND=JPRB)               :: RPROB_BIRTH(0:8,0:8)
REAL(KIND=JPRB)               :: RPROB_FUN_FERT_S(0:16)
REAL(KIND=JPRB)               :: RPROB_FUN_UW_DW_S(0:16)
REAL(KIND=JPRB)               :: RPROB_FUN_FERT_B(0:16)
REAL(KIND=JPRB)               :: RPROB_FUN_UW_DW_B(0:16)

REAL(KIND=JPRB)   ,ALLOCATABLE:: RLONDEP(:,:,:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: RLATDEP(:,:,:)

REAL(KIND=JPRB)               :: RCA_SEEDPROB
REAL(KIND=JPRB)               :: RCADELX

LOGICAL                       :: LCUCONV_CA
LOGICAL                       :: LCA_ADVECT
LOGICAL                       :: LCA_GLOBAL
LOGICAL                       :: LCA_SMOOTH
LOGICAL                       :: LCA_TEST,LCA_ADVTEST
LOGICAL                       :: LCA_RANTROP
LOGICAL                       :: LCA_NBDEBUG
LOGICAL                       :: LCA_EXTRACT
CHARACTER(LEN=10)             :: CA_FORC
CHARACTER(LEN=10)             :: CA_WIND
CHARACTER(LEN=10)             :: CA_PROB
INTEGER(KIND=JPIM)            :: NFRCASEED
INTEGER(KIND=JPIM)            :: NTESTPROC, NTESTGP
INTEGER(KIND=JPIM)            :: NDXUNREAL, NDYUNREAL

REAL(KIND=JPRB)   ,ALLOCATABLE:: SLLONG(:), SLLAT(:)
REAL(KIND=JPRB)   ,ALLOCATABLE:: SLDLONG(:), SLDLAT(:), SLDDLAT(:)
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TECUCONVCA


!     ------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TECUCONVCA), INTENT(IN) :: SELF
INTEGER          , INTENT(IN) :: KDEPTH
INTEGER          , INTENT(IN) :: KOUTNO





















































END SUBROUTINE PRINT_CONFIGURATION

!     ------------------------------------------------------
!     CA-initialization
!       - read NAMELIST
!       - if CA active initialize constants and allocate arrays
!     ------------------------------------------------------
SUBROUTINE INI_CUCONVCA(YDGEOMETRY,YDDYNA,YDECUCONVCA,YDSL)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYNA      , ONLY : TDYNA
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0       , ONLY : LECMWF
USE YOMCST       , ONLY : RPI
USE YOMLUN       , ONLY : NULOUT, NULNAM
USE YOMGRIB      , ONLY : NENSFNB
USE EINT_MOD     , ONLY : SL_STRUCT
IMPLICIT NONE
TYPE(GEOMETRY)  , INTENT(IN)   :: YDGEOMETRY
TYPE(TDYNA)     , INTENT(IN)   :: YDDYNA
TYPE(TECUCONVCA), INTENT(INOUT),TARGET :: YDECUCONVCA
TYPE(SL_STRUCT) , INTENT(INOUT):: YDSL
































































END SUBROUTINE INI_CUCONVCA


!     ------------------------------------------------------
!     prepare arrays containing latitude, longitude and 
!     gridbox size for every gridbox in SL-buffer
!     ------------------------------------------------------
SUBROUTINE SETUP_LATLONGHELP(YDGEOMETRY,YDECUCONVCA,YDSL)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : RPI

USE YOMMP0   , ONLY : MY_REGION_NS, MY_REGION_EW
USE EINT_MOD , ONLY : SL_STRUCT
IMPLICIT NONE
TYPE(GEOMETRY)  , INTENT(IN)    :: YDGEOMETRY
TYPE(TECUCONVCA), INTENT(INOUT) :: YDECUCONVCA
TYPE(SL_STRUCT) , INTENT(INOUT) :: YDSL














END SUBROUTINE SETUP_LATLONGHELP


!     ------------------------------------------------------
!     calculate neighbours on reduced gaussian grid
!     for each gridcell on this processor
!       - neighbour given in terms of index in the SL-buffer
!     ------------------------------------------------------
SUBROUTINE CALCULATE_NEIGHBOURS(YDECUCONVCA,YDGEOMETRY,YDSL)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0  , ONLY : MY_REGION_NS, MY_REGION_EW
USE YOMCST  , ONLY : RPI
USE EINT_MOD, ONLY : SL_STRUCT
USE YOMCT0  , ONLY : LRPLANE
IMPLICIT NONE
TYPE(TECUCONVCA),INTENT(INOUT) :: YDECUCONVCA
TYPE(GEOMETRY)  , INTENT(IN)    :: YDGEOMETRY
TYPE(SL_STRUCT) , INTENT(INOUT) :: YDSL



















END SUBROUTINE CALCULATE_NEIGHBOURS


!     ------------------------------------------------------
!     initialize CA state
!       - set KCELL=KLIVES and KFERT=NFERTYRS where KINI == 1
!       - used for initialization at first timestep and for
!         seeding the CA with "new" cells later on
!     ------------------------------------------------------
SUBROUTINE INITIALIZE_CELLS(YDECUCONVCA,YDGEOMETRY,KINI,KLIVES,KCELL,KFERT)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
IMPLICIT NONE
TYPE(TECUCONVCA), INTENT(INOUT):: YDECUCONVCA
TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN)    :: KLIVES(YDGEOMETRY%YRGEM%NGPTOT)
INTEGER(KIND=JPIM), INTENT(IN)    :: KINI(YDGEOMETRY%YRGEM%NGPTOT)
INTEGER(KIND=JPIM), INTENT(INOUT) :: KCELL(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)
INTEGER(KIND=JPIM), INTENT(INOUT) :: KFERT(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)






END SUBROUTINE INITIALIZE_CELLS


!     ------------------------------------------------------
!     update CA
!       - evolve CA according to the probabilities set in the
!         initial setup
!       - advect CA if LD_ADVECT is set
!       - average to coarse CA grid
!     ------------------------------------------------------
SUBROUTINE UPDCELAUT_RGG(YDGEOMETRY,YDECUCONVCA,YDECUMF,YDSL, &
&                       KLIVE,KDX,KDY,KFERT,KCELL,PFERTIN,PCELLIN,PWGHT,PRAND1D,LD_ADVECT)
USE YOECUMF      , ONLY : TECUMF
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT3       , ONLY : NSTEP
USE EINT_MOD     , ONLY : SL_STRUCT
IMPLICIT NONE
TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
TYPE(TECUCONVCA)  , INTENT(INOUT) :: YDECUCONVCA
TYPE(TECUMF)      , INTENT(INOUT) :: YDECUMF
TYPE(SL_STRUCT)   , INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM), INTENT(INOUT) :: KLIVE(YDGEOMETRY%YRGEM%NGPTOT)                
INTEGER(KIND=JPIM), INTENT(IN)    :: KDX(YDGEOMETRY%YRGEM%NGPTOT),KDY(YDGEOMETRY%YRGEM%NGPTOT)   
REAL(KIND=JPRB)   , INTENT(IN)    :: PRAND1D(NIJH*NIJH*YDGEOMETRY%YRGEM%NGPTOT)    
REAL(KIND=JPRB)   , INTENT(INOUT) :: PFERTIN(YDSL%NASLB1,NIJH*NIJH)    
REAL(KIND=JPRB)   , INTENT(INOUT) :: PCELLIN(YDSL%NASLB1,NIJH*NIJH)    
INTEGER(KIND=JPIM),INTENT(OUT) :: KFERT(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)      
INTEGER(KIND=JPIM),INTENT(OUT) :: KCELL(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)      
REAL(KIND=JPRB)   ,INTENT(OUT) :: PWGHT(YDGEOMETRY%YRGEM%NGPTOT)                
LOGICAL,INTENT(IN)               :: LD_ADVECT














END SUBROUTINE UPDCELAUT_RGG


!     ------------------------------------------------------
!     - average over NIJHxNIJH cell blocks, smooth and weight
!     ------------------------------------------------------
SUBROUTINE WEIGHTING_FIELD(YDGEOMETRY,YDECUCONVCA,YDSL,KCELL,PWGHT)



USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE EINT_MOD     , ONLY : SL_STRUCT
USE ORDER_INDEPENDENT_SUMMATION_MOD, ONLY : ORDER_INDEP_GLOBAL_SUM
IMPLICIT NONE
TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
TYPE(TECUCONVCA)  , INTENT(INOUT) :: YDECUCONVCA
TYPE(SL_STRUCT)   , INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM), INTENT(IN)    :: KCELL(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PWGHT(YDGEOMETRY%YRGEM%NGPTOT)







END SUBROUTINE WEIGHTING_FIELD


!     ------------------------------------------------------
!       apply smoothing
!     ------------------------------------------------------
SUBROUTINE SMOOTH_CA(YDGEOMETRY,YDECUCONVCA,YDSL,PFLD)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK

USE EINT_MOD , ONLY : SL_STRUCT
IMPLICIT NONE
TYPE(GEOMETRY)   , INTENT(IN)    :: YDGEOMETRY
TYPE(TECUCONVCA) , INTENT(INOUT) :: YDECUCONVCA
TYPE(SL_STRUCT)  , INTENT(INOUT) :: YDSL
REAL(KIND=JPRB)  , INTENT(INOUT) :: PFLD(YDGEOMETRY%YRGEM%NGPTOT)












END SUBROUTINE SMOOTH_CA


!     ------------------------------------------------------
!     perform advection (using departure points from SL-advection)
!     ------------------------------------------------------
SUBROUTINE ADVECT_CA_SOPH(YDGEOMETRY,YDECUCONVCA,YDECUMF,YDSL,KCELL,KFERT,PFERTIN,PCELLIN)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULERR
USE YOECUMF  , ONLY : TECUMF
USE YOMCST   , ONLY : RPI
USE EINT_MOD , ONLY : SL_STRUCT
IMPLICIT NONE
TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
TYPE(TECUCONVCA)  , INTENT(INOUT) :: YDECUCONVCA
TYPE(TECUMF)      , INTENT(INOUT) :: YDECUMF
TYPE(SL_STRUCT)   , INTENT(INOUT) :: YDSL
REAL(KIND=JPRB)   , INTENT(INOUT) :: PFERTIN(YDSL%NASLB1,NIJH*NIJH)    
REAL(KIND=JPRB)   , INTENT(INOUT) :: PCELLIN(YDSL%NASLB1,NIJH*NIJH)    
INTEGER(KIND=JPIM), INTENT(INOUT) :: KFERT(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)      
INTEGER(KIND=JPIM), INTENT(INOUT) :: KCELL(YDGEOMETRY%YRGEM%NGPTOT,NIJH*NIJH)










 














END SUBROUTINE ADVECT_CA_SOPH

!     ------------------------------------------------------
!     binary 2D output for debuging
!     ------------------------------------------------------
SUBROUTINE WRITE_FIELD(CDFILENAME, PFIELD, KPOINTS, KLEVS)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN  , ONLY : NULOUT
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)      :: KPOINTS, KLEVS
CHARACTER(LEN=*)  ,INTENT(IN)      :: CDFILENAME
REAL(KIND=JPRB)   ,INTENT(IN)      :: PFIELD(KPOINTS,KLEVS)









END SUBROUTINE WRITE_FIELD


END MODULE YOE_CUCONVCA
