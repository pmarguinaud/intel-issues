MODULE FIELD_GFL_WRAPPER

! New model fields  AJG 5/11/2012
!
! For initial implementation in the IFS, we will point at the old 
! GFL/GMV/GMVS and surface fields, to allow a progressive 
! changeover. This subroutine does that mapping.
!

!**   MODIFICATIONS
!     -------------

! B. Ingleby       2019-03-12 Add y2sh

USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0             , ONLY : NUNDEFLD
USE SURFACE_FIELDS_MIX , ONLY : TSURF,TYPE_SURF_MTL
USE YOMGMV             , ONLY : TGMV
USE YOMGFL             , ONLY : TGFL
USE YOM_YGFL           , ONLY : TYPE_GFLD,TYPE_GFL_COMP
USE YOMLUN             , ONLY : NULOUT,NULERR

USE FIELD_DEFINITIONS_BASE, ONLY : JP_T0, JP_T9, JP_T1, JP_PC,JP_PT,TYPE_FVAR ,JP_T0_DL,JP_T0_DM,JP_T9_DL,JP_T9_DM,RESET_FVAR
USE FIELD_DEFINITIONS, ONLY : FID

IMPLICIT NONE

PRIVATE

INTEGER(KIND=JPIM), PARAMETER :: JP_NOT_SET = 0, JP_NOT_ACTIVE=-1

TYPE TGFL_MAP_INTERNAL  
INTEGER(KIND=JPIM), POINTER :: MT9_MAP(:) => NULL()
INTEGER(KIND=JPIM), POINTER :: MT1_MAP(:) => NULL()
INTEGER(KIND=JPIM), POINTER :: MT0_MAP(:) => NULL()
INTEGER(KIND=JPIM), POINTER :: MPC_MAP(:) => NULL()
INTEGER(KIND=JPIM), POINTER :: MPT_MAP(:) => NULL()
INTEGER(KIND=JPIM), POINTER :: MT0_DL_MAP(:) => NULL()
INTEGER(KIND=JPIM), POINTER :: MT0_DM_MAP(:) => NULL()
INTEGER(KIND=JPIM), POINTER :: MT9_DL_MAP(:) => NULL()
INTEGER(KIND=JPIM), POINTER :: MT9_DM_MAP(:) => NULL()
INTEGER(KIND=JPIM), POINTER :: MT0_SIZE(:) => NULL()
END TYPE TGFL_MAP_INTERNAL

PUBLIC GFL_AS_A_WHOLE, GFL_ATTACH, GFL_FID_SETUP, GFL_MAPPING_SETUP,GFL_MAP,&
 & TGFL_MAP_INTERNAL
#include "abor1.intfb.h"

contains

! ----------------------------------------------------------------------
! Attaches to the GFL blocks as a whole. I.e., this just gives out direct
! access to the GFL!
! ----------------------------------------------------------------------
SUBROUTINE GFL_AS_A_WHOLE(KBLOCK, KATTACH, YDGFL,PFPTR)
INTEGER(KIND=JPIM),       INTENT(IN)    :: KBLOCK     
INTEGER(KIND=JPIM),       INTENT(IN)    :: KATTACH    
TYPE(TGFL),      TARGET,  INTENT(INOUT) :: YDGFL
REAL(KIND=JPRB), POINTER, INTENT(OUT)   :: PFPTR(:,:,:)




END SUBROUTINE GFL_AS_A_WHOLE

! ----------------------------------------------------------------------
! Sets up the the mapping between new fields and old
! NFTBC - not deallocated ever ....
! ----------------------------------------------------------------------
SUBROUTINE GFL_MAPPING_SETUP(YDMAP,YDGMV,YDGFL,YDGFLDESC, YDSURF,FVAR)
TYPE(TGFL_MAP_INTERNAL), INTENT(INOUT) :: YDMAP
TYPE(TGMV)             , INTENT(INOUT) :: YDGMV
TYPE(TGFL)             , INTENT(INOUT) :: YDGFL
TYPE(TYPE_GFLD)        , INTENT(IN)    :: YDGFLDESC
TYPE(TSURF)            , INTENT(INOUT) :: YDSURF
TYPE(TYPE_FVAR)        , INTENT(INOUT) :: FVAR(:)






END SUBROUTINE GFL_MAPPING_SETUP

! ----------------------------------------------------------------------
! Maps from KID to the old GFL/GMV/surface pointer
!
! The GFL poiner can be +ve (== a "valid" ID)
!                        0  (== not_set, i.e. the mapping needs to be added above)
!                -99999999  (The value of an unused pointer in GFL world)       
! ----------------------------------------------------------------------
FUNCTION GFL_MAP(YDMAP,KID, KATTACH, LD_IGNORE_NOT_SET)
TYPE(TGFL_MAP_INTERNAL),INTENT(IN) :: YDMAP
INTEGER(KIND=JPIM), INTENT(IN) :: KID     
INTEGER(KIND=JPIM), INTENT(IN) :: KATTACH 
LOGICAL, OPTIONAL,  INTENT(IN) :: LD_IGNORE_NOT_SET 
INTEGER(KIND=JPIM)             :: GFL_MAP









END FUNCTION GFL_MAP

! ----------------------------------------------------------------------
! Sets up the list of valid FIDS when attaching to a GFL. 
! ----------------------------------------------------------------------
SUBROUTINE GFL_FID_SETUP(YDMAP,KFIELDS_ALL, KATTACH, KFIELDS_USED, LDON, LDACTIVE,KFIDS_USED,KFIDS)
TYPE(TGFL_MAP_INTERNAL),INTENT(IN) :: YDMAP
INTEGER(KIND=JPIM), INTENT(IN)    :: KFIELDS_ALL    
INTEGER(KIND=JPIM), INTENT(IN)    :: KATTACH        
INTEGER(KIND=JPIM), INTENT(OUT)   :: KFIELDS_USED   
LOGICAL,            INTENT(INOUT) :: LDON(:)        
LOGICAL,            INTENT(INOUT) :: LDACTIVE(:)        
INTEGER(KIND=JPIM), ALLOCATABLE, INTENT(OUT) :: KFIDS_USED(:)  
INTEGER(KIND=JPIM), OPTIONAL,INTENT(IN) :: KFIDS(:) 










END SUBROUTINE GFL_FID_SETUP

! ----------------------------------------------------------------------
! Attaches to individual members in the GFL, GMV, and surface structures
! ----------------------------------------------------------------------
SUBROUTINE GFL_ATTACH(FVAR,YDMAP,KLON, KLEV, KID, KBLOCK, KATTACH, YDGMV,YDGFL,YDSURF,FPTR_1D, FPTR_2D, &
 & FPTR_3D,FPTR_UNASSIGNED,KDIM3 )
TYPE(TGFL_MAP_INTERNAL),INTENT(IN) :: YDMAP
TYPE(TYPE_FVAR),    INTENT(IN) :: FVAR(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KLON, KLEV 
INTEGER(KIND=JPIM), INTENT(IN) :: KID        
INTEGER(KIND=JPIM), INTENT(IN) :: KBLOCK     
INTEGER(KIND=JPIM), INTENT(IN) :: KATTACH    
TYPE(TGMV),TARGET , INTENT(IN) :: YDGMV
TYPE(TGFL),TARGET , INTENT(IN) :: YDGFL
TYPE(TSURF),TARGET, INTENT(IN) :: YDSURF
REAL(KIND=JPRB), POINTER, OPTIONAL :: FPTR_1D(:), FPTR_2D(:,:), FPTR_3D(:,:,:)
INTEGER(KIND=JPIM), INTENT(IN),OPTIONAL :: KDIM3
REAL(KIND=JPRB), TARGET, INTENT(IN) :: FPTR_UNASSIGNED(:,:,:)




END SUBROUTINE GFL_ATTACH


SUBROUTINE SET_MAP(FVAR,KMAP,KFID,KTO,LDFIRST,YDGFLC,YDSURF,KFIELDS,KMAPFIELDS)
TYPE(TYPE_FVAR), INTENT(IN) :: FVAR(:)
INTEGER(KIND=JPIM), INTENT(OUT) :: KMAP(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KFID
INTEGER(KIND=JPIM), INTENT(IN) :: KTO
LOGICAL, OPTIONAL, INTENT(IN)   :: LDFIRST
TYPE(TYPE_GFL_COMP),OPTIONAL,INTENT(IN) :: YDGFLC
CLASS(TYPE_SURF_MTL),OPTIONAL,INTENT(IN):: YDSURF
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM),OPTIONAL :: KMAPFIELDS(:)



END SUBROUTINE SET_MAP

END MODULE FIELD_GFL_WRAPPER
