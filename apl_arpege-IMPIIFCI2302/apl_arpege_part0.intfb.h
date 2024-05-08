INTERFACE
SUBROUTINE APL_ARPEGE_PART0(YDMF_PHYS_BASE_STATE, YDMF_PHYS_NEXT_STATE, YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS,&
 & YDCPG_MISC, YDCPG_GPAR, YDCPG_PHY0, YDMF_PHYS, YDCPG_DYN0, YDMF_PHYS_SURF, YDCPG_SL2, YDVARS, YDMODEL, YDDDH,&
 & YDPHYSMWAVE, KMOC_CLPH, PDSA_C1, PDSA_C2, PDSA_CPS, PDSA_LHS, PDSA_RS, PFLU_CDN, PFLU_CD, PFLU_CHN, PFLU_CH,&
 & PFLU_EMIS, PFLU_FEVI, PFLU_NEIJ, PFLU_QSATS, PFLU_QSAT, PFLU_VEG, PMSC_FRMQ, PFCHQ, PMSC_LH, PMSC_LSCPE,&
 & PMSC_QW, PMSC_TW, PPFL_FEFB1, PPFL_FEFB2, PPFL_FEFB3, PPFL_FPLCH, PPFL_FPLSH, PPFL_FTKEI, PPFL_FTKE,&
 & PRDG_CVGQ, PRDG_CVGT, PRDG_LCVQ, PRDG_MU0LU, PRDG_MU0M, PRDG_MU0N, PRDG_MU0, PRDG_SOLO, PSAV_DDAL, PSAV_DDOM,&
 & PSAV_ENTCH, PSAV_FHPS, PSAV_GZ0F, PSAV_GZ0HF, PSAV_HV, PSAV_PBLH, PSAV_QSH , PSAV_UDAL, PSAV_UDGRO, PSAV_UDOM,&
 & PSAV_UNEBH, PQCO2) 
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD , ONLY : MF_PHYS_TYPE
USE CPG_TYPE_MOD , ONLY : CPG_MISC_TYPE, CPG_DYN_TYPE,&
 & CPG_GPAR_TYPE, CPG_PHY_TYPE 
USE CPG_SL2_TYPE_MOD , ONLY : CPG_SL2_TYPE
USE CPG_OPTS_TYPE_MOD , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD&
 & , ONLY : MF_PHYS_SURF_TYPE 
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE MF_PHYS_BASE_STATE_TYPE_MOD&
 & , ONLY : MF_PHYS_BASE_STATE_TYPE 
USE MF_PHYS_NEXT_STATE_TYPE_MOD&
 & , ONLY : MF_PHYS_NEXT_STATE_TYPE 
USE TYPE_MODEL , ONLY : MODEL
USE PARKIND1 , ONLY : JPIM ,JPRB
USE DDH_MIX , ONLY : TYP_DDH
USE YOE_PHYS_MWAVE , ONLY : TEPHYSMWAVE
TYPE(MF_PHYS_BASE_STATE_TYPE), INTENT(IN) :: YDMF_PHYS_BASE_STATE
TYPE(MF_PHYS_NEXT_STATE_TYPE), INTENT(INOUT) :: YDMF_PHYS_NEXT_STATE
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE), INTENT(IN) :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE), INTENT(IN) :: YDCPG_OPTS
TYPE(CPG_MISC_TYPE), INTENT(INOUT) :: YDCPG_MISC
TYPE(CPG_GPAR_TYPE), INTENT(INOUT) :: YDCPG_GPAR
TYPE(CPG_PHY_TYPE), INTENT(IN) :: YDCPG_PHY0
TYPE(MF_PHYS_TYPE), INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE), INTENT(IN) :: YDCPG_DYN0
TYPE(MF_PHYS_SURF_TYPE), INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(CPG_SL2_TYPE), INTENT(INOUT) :: YDCPG_SL2
TYPE(FIELD_VARIABLES), INTENT(INOUT) :: YDVARS
TYPE(MODEL), INTENT(IN) :: YDMODEL
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TEPHYSMWAVE), INTENT(INOUT) :: YDPHYSMWAVE
INTEGER(KIND=JPIM), INTENT(INOUT) :: KMOC_CLPH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PDSA_C1 (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PDSA_C2 (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PDSA_CPS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PDSA_LHS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PDSA_RS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PFLU_CDN (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PFLU_CD (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PFLU_CHN (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PFLU_CH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PFLU_EMIS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PFLU_FEVI (YDCPG_OPTS%KLON, 1:YDMODEL%YRML_PHY_G%YRDPHY%NTSSG+1)
REAL(KIND=JPRB), INTENT(INOUT) :: PFLU_NEIJ (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PFLU_QSATS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PFLU_QSAT (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PFLU_VEG (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PMSC_FRMQ (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PFCHQ (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PMSC_LH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PMSC_LSCPE (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PMSC_QW (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PMSC_TW (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PPFL_FEFB1 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PPFL_FEFB2 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PPFL_FEFB3 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PPFL_FPLCH (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PPFL_FPLSH (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PPFL_FTKEI (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PPFL_FTKE (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PRDG_CVGQ (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PRDG_CVGT (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PRDG_LCVQ (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PRDG_MU0LU (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PRDG_MU0M (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PRDG_MU0N (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PRDG_MU0 (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PRDG_SOLO(YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_DDAL (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_DDOM (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_ENTCH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_FHPS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_GZ0F (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_GZ0HF (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_HV (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_PBLH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_QSH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_UDAL (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_UDGRO (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_UDOM (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PSAV_UNEBH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(INOUT) :: PQCO2(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
END SUBROUTINE APL_ARPEGE_PART0
END INTERFACE
