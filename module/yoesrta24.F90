MODULE YOESRTA24


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA24* - SRTM COEFFICIENTS FOR INTERVAL 24
!     BAND 24: 12850-16000 cm-1 (low - H2O,O2; high - O2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG24 = 16

REAL_B :: KA(9,5,13,NG24) ,ABSA(585,NG24)
REAL_B :: KB(5,13:59,NG24),ABSB(235,NG24)
REAL_B :: SELFREF(10,NG24),FORREF(3,NG24)
REAL_B :: SFLUXREF(NG24,9)
REAL_B :: ABSO3A(NG24), ABSO3B(NG24), RAYLA(NG24,9), RAYLB(NG24)
REAL_B :: STRRAT
INTEGER_M :: LAYREFFR

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL 
! SFLUXREF: REAL
! ABSO3A  : REAL
! ABSO3B  : REAL
! RAYLA   : REAL
! RAYLB   : REAL   
! STRRAT  : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA24
