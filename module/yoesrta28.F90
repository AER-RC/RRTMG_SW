MODULE YOESRTA28


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA28* - SRTM COEFFICIENTS FOR INTERVAL 28
!     BAND 28: 38000-50000 cm-1 (low - O3, O2; high - O3, O2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG28 = 16

REAL_B :: KA(9,5,13,NG28)   ,ABSA(585,NG28)
REAL_B :: KB(5,5,13:59,NG28),ABSB(1175,NG28)
REAL_B :: SFLUXREF(NG28,5)
REAL_B :: RAYL              ,STRRAT
INTEGER_M :: LAYREFFR

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SFLUXREF: REAL 
! RAYL    : REAL    
! STRRAT  : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA28
