MODULE YOESRTA18


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA18* - SRTM COEFFICIENTS FOR INTERVAL 16
!     BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG18 = 16, NGS17=32

REAL_B :: KA(9,5,13,NG18) ,ABSA(585,NG18)
REAL_B :: KB(5,13:59,NG18),ABSB(235,NG18)
REAL_B :: SELFREF(10,NG18),FORREF(3,NG18)
REAL_B :: SFLUXREF(NG18,9)
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
END MODULE YOESRTA18
