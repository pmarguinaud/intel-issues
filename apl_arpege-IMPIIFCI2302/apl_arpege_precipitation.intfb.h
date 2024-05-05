INTERFACE
SUBROUTINE APL_ARPEGE_PRECIPITATION (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS,&
 & YDMF_PHYS_SURF, YDVARS, YDMODEL, PFLU_NEIJ, PNEBS, PNEB_CVPP, PQC_DET_PCMT, PQI, PQL, PQLIS,&
 & PQLI_CVPP, PQR, PQS, PQV, PSEDIQI, PSEDIQL, PTENHA, PTENQVA, YDSTA) 
USE PARKIND1, ONLY : JPIM, JPRB
USE MF_PHYS_BASE_STATE_TYPE_MOD&
 & , ONLY : MF_PHYS_BASE_STATE_TYPE 
USE CPG_OPTS_TYPE_MOD , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_TYPE_MOD , ONLY : MF_PHYS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD&
 & , ONLY : MF_PHYS_SURF_TYPE 
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE TYPE_MODEL , ONLY : MODEL
USE YOMSTA , ONLY : TSTA
TYPE(MF_PHYS_BASE_STATE_TYPE), INTENT(IN) :: YDMF_PHYS_BASE_STATE
TYPE(CPG_BNDS_TYPE), INTENT(IN) :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE), INTENT(IN) :: YDCPG_OPTS
TYPE(MF_PHYS_TYPE), INTENT(INOUT) :: YDMF_PHYS
TYPE(MF_PHYS_SURF_TYPE), INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES), INTENT(INOUT) :: YDVARS
TYPE(MODEL), INTENT(IN) :: YDMODEL
REAL(KIND=JPRB), INTENT(IN) :: PFLU_NEIJ (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PNEBS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PNEB_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PQC_DET_PCMT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PQI(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PQL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PQLIS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PQLI_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PQR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PQS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PQV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PSEDIQI(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PSEDIQL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(OUT) :: PTENHA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(OUT) :: PTENQVA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
TYPE(TSTA), INTENT(IN) :: YDSTA
END SUBROUTINE APL_ARPEGE_PRECIPITATION
END INTERFACE
