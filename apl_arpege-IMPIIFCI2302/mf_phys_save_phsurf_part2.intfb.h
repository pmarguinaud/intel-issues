INTERFACE
SUBROUTINE MF_PHYS_SAVE_PHSURF_PART2 (YDCPG_BNDS, YDCPG_OPTS, PSAV_DDAL, PSAV_DDOM, PSAV_ENTCH,&
 & PSAV_FHPS, PSAV_GZ0F, PSAV_GZ0HF, PSAV_HV, PSAV_PBLH, PSAV_QSH, PSAV_UDAL, PSAV_UDGRO, PSAV_UDOM,&
 & PSAV_UNEBH, PGSD_VF_PZ0F, PGSD_VH_PPBLH, PGSD_VH_PQSH, PGSD_VH_PSPSH, PGSD_VK_PUDGRO, PGSD_VV_PHV,&
 & PGSD_VV_PZ0H, PDAL_T0, PDOM_T0, PUAL_T0, PUEN_T0, PUNEBH_T0, PUOM_T0, YDMODEL) 
USE PARKIND1, ONLY : JPIM, JPRB
USE CPG_OPTS_TYPE_MOD , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE TYPE_MODEL , ONLY : MODEL
TYPE(CPG_BNDS_TYPE), INTENT(IN) :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE), INTENT(IN) :: YDCPG_OPTS
REAL (KIND=JPRB), INTENT(IN) :: PSAV_DDAL (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_DDOM (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_ENTCH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_FHPS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_GZ0F (YDCPG_OPTS%KLON)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_GZ0HF (YDCPG_OPTS%KLON)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_HV (YDCPG_OPTS%KLON)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_PBLH (YDCPG_OPTS%KLON)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_QSH (YDCPG_OPTS%KLON)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_UDAL (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_UDGRO (YDCPG_OPTS%KLON)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_UDOM (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(IN) :: PSAV_UNEBH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), INTENT(OUT) :: PGSD_VF_PZ0F (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(OUT) :: PGSD_VH_PPBLH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(OUT) :: PGSD_VH_PQSH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(OUT) :: PGSD_VH_PSPSH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(OUT) :: PGSD_VK_PUDGRO (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(OUT) :: PGSD_VV_PHV (YDCPG_OPTS%KLON)
REAL(KIND=JPRB), INTENT(OUT) :: PGSD_VV_PZ0H (YDCPG_OPTS%KLON)
REAL (KIND=JPRB), INTENT(OUT) :: PDAL_T0 (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(OUT) :: PDOM_T0 (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(OUT) :: PUAL_T0 (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(OUT) :: PUEN_T0 (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(OUT) :: PUNEBH_T0 (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB), INTENT(OUT) :: PUOM_T0 (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
TYPE(MODEL) ,INTENT(IN) :: YDMODEL
END SUBROUTINE MF_PHYS_SAVE_PHSURF_PART2
END INTERFACE
