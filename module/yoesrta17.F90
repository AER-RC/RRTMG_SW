MODULE YOESRTA17


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA17* - SRTM COEFFICIENTS FOR INTERVAL 17
!     BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG17 = 16, NGS16=16

REAL_B :: KA(9,5,13,NG17)   ,ABSA(585,NG17)
REAL_B :: KB(5,5,13:59,NG17),ABSB(1175,NG17)
REAL_B :: SELFREF(10,NG17)  ,FORREF(4,NG17)
REAL_B :: SFLUXREF(NG17,5)
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
END MODULE YOESRTA17
