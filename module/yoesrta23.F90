MODULE YOESRTA23


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA23* - SRTM COEFFICIENTS FOR INTERVAL 23
!     BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG23 = 16

REAL_B :: KA(5,13,NG23)   ,ABSA(65,NG23)
REAL_B :: SELFREF(10,NG23),FORREF(3,NG23)
REAL_B :: SFLUXREF(NG23)  ,RAYL(NG23)
REAL_B :: GIVFAC
INTEGER_M :: LAYREFFR

EQUIVALENCE (KA(1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! SELFREF : REAL 
! FORREF  : REAL    
! SFLUXREF: REAL
! RAYL    : REAL
! GIVFAC  : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA23
