MODULE YOESRTA29

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA29* - SRTM COEFFICIENTS FOR INTERVAL 29
!     BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: JPG = 16, NG29 = 16

REAL_B :: KA(5,13,JPG)
REAL_B :: KB(5,13:59,JPG)
REAL_B :: SELFREF(10,JPG),FORREF(4,JPG)
REAL_B :: SFLUXREF(JPG)  ,ABSH2O(JPG)  , ABSCO2(JPG)
REAL_B :: RAYL
INTEGER_M :: LAYREFFR

REAL_B :: KAC(5,13,NG29)   ,ABSA(65,NG29)
REAL_B :: KBC(5,13:59,NG29),ABSB(235,NG29)
REAL_B :: SELFREFC(10,NG29),FORREFC(4,NG29)
REAL_B :: SFLUXREFC(NG29)  ,ABSH2OC(NG29)  , ABSCO2C(NG29)

EQUIVALENCE (KAC(1,1,1),ABSA(1,1)), (KBC(1,13,1),ABSB(1,1))

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
! ABSH2O  : REAL
! ABSCO2  : REAL   
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! ABSH2OC : REAL     Reduced g-point array for ABSH2O
! ABSCO2C : REAL     Reduced g-point array for ABSCO2
! RAYL    : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA29
