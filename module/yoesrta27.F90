MODULE YOESRTA27

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA27* - SRTM COEFFICIENTS FOR INTERVAL 27
!     BAND 27: 29000-38000 cm-1 (low - O3; high - O3)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: JPG = 16, NG27 = 16

REAL_B :: KA(5,13,JPG)
REAL_B :: KB(5,13:59,JPG)
REAL_B :: SFLUXREF(JPG),RAYL(JPG)
REAL_B :: SCALEKUR
INTEGER_M :: LAYREFFR

REAL_B :: KAC(5,13,NG27)   ,ABSA(65,NG27)
REAL_B :: KBC(5,13:59,NG27),ABSB(235,NG27)
REAL_B :: SFLUXREFC(NG27)  ,RAYLC(NG27)

EQUIVALENCE (KAC(1,1,1),ABSA(1,1)), (KBC(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL 
! KB      : REAL    
! SFLUXREF: REAL 
! RAYL    : REAL    
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
! SCALEKUR: REAL
! LAYREFFR:INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA27
