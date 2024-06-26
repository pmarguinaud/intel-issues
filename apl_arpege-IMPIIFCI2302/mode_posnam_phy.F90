MODULE MODE_POSNAM_PHY
IMPLICIT NONE
CONTAINS
SUBROUTINE POSNAM_PHY(KUNITNML, CDNAML, LDNEEDNAM, LDFOUND, KLUOUT)
!Wrapper to call the AROME version of posnam

IMPLICIT NONE

INTEGER,          INTENT(IN)    :: KUNITNML  !< Logical unit to access the namelist
CHARACTER(LEN=*), INTENT(IN)    :: CDNAML    !< Namelist name
LOGICAL,          INTENT(IN)    :: LDNEEDNAM !< True to abort if namelist is absent
LOGICAL,          INTENT(OUT)   :: LDFOUND   !< True if namelist has been found
INTEGER,          INTENT(IN)    :: KLUOUT    !< Logical unit for output


END SUBROUTINE POSNAM_PHY
END MODULE MODE_POSNAM_PHY
