MODULE YOESRTA25


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA25* - SRTM COEFFICIENTS FOR INTERVAL 25
!     BAND 25: 16000-22650 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG25 = 16

REAL_B :: KA(5,13,NG25) ,ABSA(65,NG25)
REAL_B :: SFLUXREF(NG25)
REAL_B :: RAYL(NG25), ABSO3A(NG25), ABSO3B(NG25)
INTEGER_M :: LAYREFFR

EQUIVALENCE (KA(1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! SFLUXREF: REAL
! RAYL    : REAL
! ABSO3A  : REAL
! ABSO3B  : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA25
