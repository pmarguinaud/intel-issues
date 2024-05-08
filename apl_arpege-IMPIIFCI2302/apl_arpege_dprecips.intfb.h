INTERFACE
SUBROUTINE APL_ARPEGE_DPRECIPS (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS, &
&YDVARS, YDMODEL, PPRC_DPRECIPS, PPRC_DPRECIPS2)
USE PARKIND1, ONLY : JPIM, JPRB
USE MF_PHYS_BASE_STATE_TYPE_MOD&
 & , ONLY : MF_PHYS_BASE_STATE_TYPE 
USE CPG_OPTS_TYPE_MOD , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_TYPE_MOD , ONLY : MF_PHYS_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE TYPE_MODEL , ONLY : MODEL
TYPE(MF_PHYS_BASE_STATE_TYPE), INTENT(IN) :: YDMF_PHYS_BASE_STATE
TYPE(CPG_BNDS_TYPE), INTENT(IN) :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE), INTENT(IN) :: YDCPG_OPTS
TYPE(MF_PHYS_TYPE), INTENT(INOUT) :: YDMF_PHYS
TYPE(FIELD_VARIABLES), INTENT(INOUT) :: YDVARS
TYPE(MODEL), INTENT(IN) :: YDMODEL
REAL(KIND=JPRB), INTENT(OUT) :: PPRC_DPRECIPS (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%NDTPREC)
REAL(KIND=JPRB), INTENT(OUT) :: PPRC_DPRECIPS2 (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%NDTPREC2)
END SUBROUTINE APL_ARPEGE_DPRECIPS
END INTERFACE
