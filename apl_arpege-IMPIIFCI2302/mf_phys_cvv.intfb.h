INTERFACE
SUBROUTINE MF_PHYS_CVV (YDCPG_BNDS, YDCPG_OPTS, PCVV_T0, PCVV_T1)
USE PARKIND1, ONLY : JPIM, JPRB
USE CPG_OPTS_TYPE_MOD , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
TYPE(CPG_BNDS_TYPE),INTENT(IN) :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),INTENT(IN) :: YDCPG_OPTS
REAL (KIND=JPRB), INTENT(IN) :: PCVV_T0 (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(OUT) :: PCVV_T1 (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
END SUBROUTINE MF_PHYS_CVV
END INTERFACE