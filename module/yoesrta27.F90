MODULE YOESRTA27


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA27* - SRTM COEFFICIENTS FOR INTERVAL 27
!     BAND 27: 29000-38000 cm-1 (low - O3; high - O3)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG27 = 16

REAL_B :: KA(5,13,NG27)   ,ABSA(65,NG27)
REAL_B :: KB(5,13:59,NG27),ABSB(235,NG27)
REAL_B :: SFLUXREF(NG27)  ,RAYL(NG27)
REAL_B :: SCALEKUR
INTEGER_M :: LAYREFFR

EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL 
! KB      : REAL    
! SFLUXREF: REAL 
! RAYL    : REAL    
! SCALEKUR: REAL
! LAYREFFR:INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA27
