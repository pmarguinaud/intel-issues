MODULE TYPE_FAOPH

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!========== HANDLING OF OUTPUT FIELDS IN FA FILES  ======

! CFPCA : names of output FA frames
! CFAMODEL : FA model name for GRIB2
! NHEADMAX : maximum size of the header of a record
! NGPSIZPK : size of a gridpoint field in FA file
! NSPSIZPK : size of a spectral field in FA file
! LFAEXTERN : .TRUE. to write fields in a separate GRIB file

TYPE TFAOPH

CHARACTER(LEN=16)  :: CFPCA
CHARACTER(LEN=64)  :: CFAMODEL
INTEGER(KIND=JPIM) :: NHEADMAX
INTEGER(KIND=JPIM) :: NGPSIZPK
INTEGER(KIND=JPIM) :: NSPSIZPK
LOGICAL            :: LFAEXTERN

END TYPE TFAOPH

!     ------------------------------------------------------------------
END MODULE TYPE_FAOPH
