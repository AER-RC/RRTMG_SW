MODULE YOESRTA22


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA22* - SRTM COEFFICIENTS FOR INTERVAL 22
!     BAND 22:  7700-8050 cm-1 (low - H2O,O2; high - O2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG22 = 16

REAL_B :: KA(9,5,13,NG22) ,ABSA(585,NG22)
REAL_B :: KB(5,13:59,NG22),ABSB(235,NG22)
REAL_B :: SELFREF(10,NG22),FORREF(3,NG22)
REAL_B :: SFLUXREF(NG22,9)
REAL_B :: RAYL            ,STRRAT
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
! RAYL    : REAL
! STRRAT  : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA22
