MODULE YOESRTA18

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA18* - SRTM COEFFICIENTS FOR INTERVAL 16
!     BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)
!     -----------------------------------------------------------------

INTEGER_B, PARAMETER :: JPG = 16, NG18 = 16

REAL_B :: KA(9,5,13,JPG)
REAL_B :: KB(5,13:59,JPG)
REAL_B :: SELFREF(10,JPG),FORREF(3,JPG)
REAL_B :: SFLUXREF(JPG,9)
REAL_B :: RAYL,STRRAT

INTEGER_B :: LAYREFFR

REAL_B :: KAC(9,5,13,NG18) ,ABSA(585,NG18)
REAL_B :: KBC(5,13:59,NG18),ABSB(235,NG18)
REAL_B :: SELFREFC(10,NG18),FORREFC(3,NG18)
REAL_B :: SFLUXREFC(NG18,9)

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
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYL    : REAL
! STRRAT  : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA18
