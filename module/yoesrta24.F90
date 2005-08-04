MODULE YOESRTA24

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA24* - SRTM COEFFICIENTS FOR INTERVAL 24
!     BAND 24: 12850-16000 cm-1 (low - H2O,O2; high - O2)
!     -----------------------------------------------------------------

INTEGER_B, PARAMETER :: JPG = 16, NG24 = 16

REAL_B :: KA(9,5,13,JPG)
REAL_B :: KB(5,13:59,JPG)
REAL_B :: SELFREF(10,JPG),FORREF(3,JPG)
REAL_B :: SFLUXREF(JPG,9)
REAL_B :: ABSO3A(JPG), ABSO3B(JPG), RAYLA(JPG,9), RAYLB(JPG)
REAL_B :: STRRAT
INTEGER_B :: LAYREFFR

REAL_B :: KAC(9,5,13,NG24) ,ABSA(585,NG24)
REAL_B :: KBC(5,13:59,NG24),ABSB(235,NG24)
REAL_B :: SELFREFC(10,NG24),FORREFC(3,NG24)
REAL_B :: SFLUXREFC(NG24,9)
REAL_B :: ABSO3AC(NG24), ABSO3BC(NG24), RAYLAC(NG24,9), RAYLBC(NG24)

EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL 
! SFLUXREF: REAL
! ABSO3A  : REAL
! ABSO3B  : REAL
! RAYLA   : REAL
! RAYLB   : REAL   
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! ABSO3AC : REAL     Reduced g-point array for ABSO3A
! ABSO3BC : REAL     Reduced g-point array for ABSO3B
! RAYLAC  : REAL     Reduced g-point array for RAYLA
! RAYLBC  : REAL     Reduced g-point array for RAYLB
! STRRAT  : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA24
