PROGRAM LATTEX

! ISO/IEC 1539-1 (2010-10-15)
! 12.5.2.3 Argument association
! 1 Except in references to intrinsic inquiry functions, a pointer actual argument that corresponds 
! to a nonoptional nonpointer dummy argument shall be pointer associated with a target.

IMPLICIT NONE

REAL(KIND=8)                 :: PMIXNL(1784,15)
REAL(KIND=8), TARGET         :: PGMV(1784,15,34)
REAL(KIND=8), TARGET         :: PGMVT1(1784,15,6)
REAL(KIND=8)                 :: PGMVTNDSI(1784,15,15)
REAL(KIND=8), TARGET         :: PB1(1784,460)
REAL(KIND=8), TARGET         :: PB2(1784,136)
REAL(KIND=8)                 :: ZMOY1VWV(1784,15)
REAL(KIND=8)                 :: ZT0NHY(1784,15)

LOGICAL :: LLVD5

REAL(KIND=8), POINTER :: ZT1SVD(:,:)
REAL(KIND=8), POINTER :: ZT9CVWVNL(:,:)
REAL(KIND=8), POINTER :: ZT9NHY(:,:)
REAL(KIND=8), POINTER :: ZT9VWVNL(:,:)
REAL(KIND=8), POINTER :: ZSLB1VD0(:), ZSLB1VD9(:)
REAL(KIND=8), POINTER :: ZSLB1VD9_NL(:), ZSLB1VD9_SI(:), ZSLB1VDF9(:)
REAL(KIND=8), POINTER :: ZSLB2VDSI(:)


INTERFACE
SUBROUTINE LATTEX_DNT(KLON,KLEV,KFLSA,KFLEN,&
 & PMOY1,PMIXNL,&
 & PXSI,PXNLT9,PXT1,PXL0,PXL9,PXLF9,PCXNLT9,&
 & PSIDDHXT1,PSIDDHXT9,PSIDDHXL0,PXLF0,LDNESC,&
 & LDVD5,PDYT0,PDYT9)

IMPLICIT NONE

INTEGER(KIND=4)       :: KLON,KLEV,KFLSA,KFLEN
REAL(KIND=8)          :: PMOY1(KLON,KLEV)
REAL(KIND=8)          :: PMIXNL(KLON,KLEV)
REAL(KIND=8)          :: PXSI(KLON,KLEV)
REAL(KIND=8)          :: PXNLT9(KLON,KLEV)
REAL(KIND=8)          :: PXT1(KLON,KLEV)
REAL(KIND=8)          :: PXL0(KLON,KFLSA:KFLEN)
REAL(KIND=8)          :: PXL9(KLON,KFLSA:KFLEN)
REAL(KIND=8)          :: PXLF9(KLON,KFLSA:KFLEN)
REAL(KIND=8)          :: PCXNLT9(KLON,KLEV)
REAL(KIND=8)          :: PSIDDHXT1(KLON,KLEV)
REAL(KIND=8)          :: PSIDDHXT9(KLON,KLEV)
REAL(KIND=8)          :: PSIDDHXL0(KLON,KFLSA:KFLEN)
REAL(KIND=8)          :: PXLF0(KLON,KFLSA:KFLEN)
LOGICAL,OPTIONAL      :: LDNESC
LOGICAL,OPTIONAL      :: LDVD5
REAL(KIND=8),OPTIONAL :: PDYT0(KLON,KLEV)
REAL(KIND=8),OPTIONAL :: PDYT9(KLON,KLEV)

END SUBROUTINE LATTEX_DNT

SUBROUTINE SC2PRG (ZSLB1VDF9, ZT9CVWVNL, ZSLB1VD9_NL, ZT9NHY, CDARG)

REAL(KIND=8), POINTER :: ZSLB1VDF9(:)
REAL(KIND=8), POINTER :: ZT9CVWVNL(:,:)
REAL(KIND=8), POINTER :: ZSLB1VD9_NL(:)
REAL(KIND=8), POINTER :: ZT9NHY(:,:)
CHARACTER*1 :: CDARG

END SUBROUTINE

END INTERFACE

CHARACTER*1 :: CLARG

CLARG = '0'

IF (IARGC () > 0) CALL GETARG (1, CLARG)

ZSLB1VD0    => PB1 (:, 427)
ZSLB1VD9    => PB1 (:, 155)
ZSLB1VD9_SI => PB1 (:, 257)
ZSLB2VDSI   => PB2 (:,  61)
ZT1SVD      => PGMVT1 (:,:,5)
ZT9VWVNL    => PGMV   (:,:,30)

CALL SC2PRG (ZSLB1VDF9, ZT9CVWVNL, ZSLB1VD9_NL, ZT9NHY, CLARG)

WRITE (0, *) "BEFORE LATTEX_DNT"
CALL LATTEX_DNT(&
 & 1784,15,0,16,&
 & ZMOY1VWV,PMIXNL,&
 & ZSLB2VDSI,ZT9VWVNL,ZT1SVD,&
 & ZSLB1VD0,ZSLB1VD9,ZSLB1VDF9,ZT9CVWVNL,&
 & PGMVTNDSI(1,1,1),PGMVTNDSI(1,1,1),ZSLB1VD9_SI,&
 & ZSLB1VD9_NL,LDVD5=LLVD5,PDYT0=ZT0NHY,PDYT9=ZT9NHY)
WRITE (0, *) "AFTER LATTEX_DNT"

STOP

END 

