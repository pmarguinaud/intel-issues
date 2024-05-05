INTERFACE
SUBROUTINE APL_ARPEGE_HYDRO_BUDGET (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDCPG_GPAR,&
 & YDMF_PHYS, YDMF_PHYS_SURF, YDMODEL, PC3, PCN, PDSA_C1, PDSA_C2, PFLU_FEVI, PFLU_NEIJ, PFLU_VEG,&
 & PFPLSL, PFPLSN, PWFC, PWLMX, PWPMX, PWSEQ, PWSMX) 
USE PARKIND1, ONLY : JPIM, JPRB
USE MF_PHYS_BASE_STATE_TYPE_MOD&
 & , ONLY : MF_PHYS_BASE_STATE_TYPE 
USE CPG_OPTS_TYPE_MOD , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE CPG_TYPE_MOD , ONLY : CPG_GPAR_TYPE
USE MF_PHYS_TYPE_MOD , ONLY : MF_PHYS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD&
 & , ONLY : MF_PHYS_SURF_TYPE 
USE TYPE_MODEL , ONLY : MODEL
TYPE (MF_PHYS_BASE_STATE_TYPE), INTENT(IN) :: YDMF_PHYS_BASE_STATE
TYPE(CPG_BNDS_TYPE), INTENT(IN) :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE), INTENT(IN) :: YDCPG_OPTS
TYPE(CPG_GPAR_TYPE), INTENT(INOUT) :: YDCPG_GPAR
TYPE(MF_PHYS_TYPE), INTENT(INOUT) :: YDMF_PHYS
TYPE(MF_PHYS_SURF_TYPE), INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(MODEL), INTENT(IN) :: YDMODEL
REAL(KIND=JPRB), INTENT(IN) :: PC3(YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PCN(YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PDSA_C1 (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PDSA_C2 (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PFLU_FEVI (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KTSSG+1)
REAL(KIND=JPRB), INTENT(IN) :: PFLU_NEIJ (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PFLU_VEG (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PFPLSL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PFPLSN(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PWFC(YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PWLMX(YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PWPMX(YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PWSEQ(YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PWSMX(YDCPG_OPTS%KLON)
END SUBROUTINE APL_ARPEGE_HYDRO_BUDGET
END INTERFACE
