MODULE YOESRTA29


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA29* - SRTM COEFFICIENTS FOR INTERVAL 29
!     BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG29 = 16

REAL_B :: KA(5,13,NG29)   ,ABSA(65,NG29)
REAL_B :: KB(5,13:59,NG29),ABSB(235,NG29)
REAL_B :: SELFREF(10,NG29),FORREF(4,NG29)
REAL_B :: SFLUXREF(NG29)  ,ABSH2O(NG29)  , ABSCO2(NG29)
REAL_B :: RAYL
INTEGER_M :: LAYREFFR

EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))

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
! ABSH2O  : REAL
! ABSCO2  : REAL   
! RAYL    : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA29
