MODULE MM

TYPE TT
  REAL :: XX
  REAL :: YY
END TYPE

TYPE (TT), SAVE, TARGET :: SS

REAL, POINTER :: XX => SS%XX
REAL, POINTER :: YY => SS%YY

END MODULE

USE MM

PRINT *, LOC (XX)
PRINT *, LOC (YY)

IF (LOC (XX) == LOC (YY)) THEN
  PRINT *, " XX == YY "
ENDIF


END

