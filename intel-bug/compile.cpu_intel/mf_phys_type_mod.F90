MODULE MF_PHYS_TYPE_MOD

USE FIELD_MODULE

IMPLICIT NONE

TYPE MF_PHYS_OUT_TYPE
  REAL (KIND=8), POINTER, CONTIGUOUS :: DIFCQ (:, :) => NULL ()
  TYPE (FIELD_3D), POINTER :: F_DIFCQ => NULL ()
END TYPE MF_PHYS_OUT_TYPE

CONTAINS

SUBROUTINE MF_PHYS_OUT_TYPE_UPDATE_VIEW (SELF, BLOCK_INDEX)

TYPE (MF_PHYS_OUT_TYPE)       :: SELF
INTEGER(KIND=4), INTENT (IN)  :: BLOCK_INDEX

SELF%DIFCQ  => SELF%F_DIFCQ%PTR (:,:,BLOCK_INDEX)

WRITE (0, *) __FILE__, ':', __LINE__, " ASSOCIATED (SELF%DIFCQ) = ", ASSOCIATED (SELF%DIFCQ)
END SUBROUTINE MF_PHYS_OUT_TYPE_UPDATE_VIEW

END MODULE MF_PHYS_TYPE_MOD



