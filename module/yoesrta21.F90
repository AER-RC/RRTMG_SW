MODULE YOESRTA21


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA21* - SRTM COEFFICIENTS FOR INTERVAL 21
!     BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG21 = 16

REAL_B :: KA(9,5,13,NG21)   ,ABSA(585,NG21)
REAL_B :: KB(5,5,13:59,NG21),ABSB(1175,NG21)
REAL_B :: SELFREF(10,NG21)  ,FORREF(4,NG21)
REAL_B :: SFLUXREF(NG21,9)
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
! SELFREF : REAL 
! FORREF  : REAL    
! SFLUXREF: REAL
! RAYL    : REAL
! STRRAT  : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA21
