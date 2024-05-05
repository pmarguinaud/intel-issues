INTERFACE
SUBROUTINE APL_ARPEGE_ATMOSPHERE_UPDATE (YDMF_PHYS_BASE_STATE, YDMF_PHYS_NEXT_STATE, YDGEOMETRY,&
 & YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDMF_PHYS, YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS, YDMODEL,&
 & YDDDH, PDIFEXT, PDSA_CPS, PMSC_FRMQ, PPFL_FEFB1, PPFL_FEFB2, PPFL_FEFB3, PPFL_FTKE, PPFL_FTKEI, PFCHQ) 
USE PARKIND1, ONLY : JPIM, JPRB
USE MF_PHYS_BASE_STATE_TYPE_MOD&
 & , ONLY : MF_PHYS_BASE_STATE_TYPE 
USE MF_PHYS_NEXT_STATE_TYPE_MOD&
 & , ONLY : MF_PHYS_NEXT_STATE_TYPE 
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE CPG_OPTS_TYPE_MOD , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE CPG_TYPE_MOD , ONLY : CPG_MISC_TYPE, CPG_DYN_TYPE
USE MF_PHYS_TYPE_MOD , ONLY : MF_PHYS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD&
 & , ONLY : MF_PHYS_SURF_TYPE 
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE TYPE_MODEL , ONLY : MODEL
USE DDH_MIX , ONLY : TYP_DDH
TYPE (MF_PHYS_BASE_STATE_TYPE), INTENT(IN) :: YDMF_PHYS_BASE_STATE
TYPE (MF_PHYS_NEXT_STATE_TYPE), INTENT(INOUT) :: YDMF_PHYS_NEXT_STATE
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE), INTENT(IN) :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE), INTENT(IN) :: YDCPG_OPTS
TYPE(CPG_MISC_TYPE), INTENT(INOUT) :: YDCPG_MISC
TYPE(MF_PHYS_TYPE), INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE), INTENT(IN) :: YDCPG_DYN0
TYPE(MF_PHYS_SURF_TYPE), INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES), INTENT(INOUT) :: YDVARS
TYPE(MODEL), INTENT(IN) :: YDMODEL
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
REAL(KIND=JPRB), INTENT(IN) :: PDIFEXT (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB), INTENT(IN) :: PDSA_CPS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(IN) :: PMSC_FRMQ (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PPFL_FEFB1 (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PPFL_FEFB2 (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PPFL_FEFB3 (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PPFL_FTKE (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PPFL_FTKEI (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(IN) :: PFCHQ (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
END SUBROUTINE APL_ARPEGE_ATMOSPHERE_UPDATE
END INTERFACE
