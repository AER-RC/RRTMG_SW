MODULE YOESRTA19


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA19* - SRTM COEFFICIENTS FOR INTERVAL 19
!     BAND 19:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG19 = 16

REAL_B :: KA(9,5,13,NG19) ,ABSA(585,NG19)
REAL_B :: KB(5,13:59,NG19),ABSB(235,NG19)
REAL_B :: SELFREF(10,NG19),FORREF(3,NG19)
REAL_B :: SFLUXREF(NG19,9)
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
END MODULE YOESRTA19
