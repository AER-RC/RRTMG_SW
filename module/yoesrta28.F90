MODULE YOESRTA28

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA28* - SRTM COEFFICIENTS FOR INTERVAL 28
!     BAND 28: 38000-50000 cm-1 (low - O3, O2; high - O3, O2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: JPG = 16, NG28 = 16

REAL_B :: KA(9,5,13,JPG)
REAL_B :: KB(5,5,13:59,JPG)
REAL_B :: SFLUXREF(JPG,5)
REAL_B :: RAYL,STRRAT
INTEGER_M :: LAYREFFR

REAL_B :: KAC(9,5,13,NG28)   ,ABSA(585,NG28)
REAL_B :: KBC(5,5,13:59,NG28),ABSB(1175,NG28)
REAL_B :: SFLUXREFC(NG28,5)

EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SFLUXREF: REAL 
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYL    : REAL    
! STRRAT  : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA28
