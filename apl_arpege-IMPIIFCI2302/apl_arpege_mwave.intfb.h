INTERFACE
SUBROUTINE APL_ARPEGE_MWAVE(YDCPG_MISC, YDMF_PHYS, YDCPG_BNDS, YDCPG_OPTS, YDPHYSMWAVE)
USE CPG_TYPE_MOD , ONLY : CPG_MISC_TYPE
USE MF_PHYS_TYPE_MOD , ONLY : MF_PHYS_TYPE
USE YOE_PHYS_MWAVE , ONLY : TEPHYSMWAVE
USE CPG_OPTS_TYPE_MOD, ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
TYPE(CPG_MISC_TYPE) ,INTENT(IN) :: YDCPG_MISC
TYPE(MF_PHYS_TYPE) ,INTENT(IN) :: YDMF_PHYS
TYPE(CPG_BNDS_TYPE) ,INTENT(IN) :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE) ,INTENT(IN) :: YDCPG_OPTS
TYPE(TEPHYSMWAVE) ,INTENT(INOUT) :: YDPHYSMWAVE
END SUBROUTINE APL_ARPEGE_MWAVE
END INTERFACE
