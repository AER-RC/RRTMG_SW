MODULE YOESRTA20


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA20* - SRTM COEFFICIENTS FOR INTERVAL 20
!     BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG20 = 16

REAL_B :: KA(5,13,NG20)   ,ABSA(65,NG20)
REAL_B :: KB(5,13:59,NG20),ABSB(235,NG20)
REAL_B :: SELFREF(10,NG20),FORREF(4,NG20)
REAL_B :: SFLUXREF(NG20)  ,ABSCH4(NG20)
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
! ABSCH4  : REAL  
! RAYL    : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA20
