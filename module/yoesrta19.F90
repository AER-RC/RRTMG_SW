MODULE YOESRTA19

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA19* - SRTM COEFFICIENTS FOR INTERVAL 19
!     BAND 19:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)
!     -----------------------------------------------------------------

INTEGER_B, PARAMETER :: JPG = 16, NG19 = 16

REAL_B :: KA(9,5,13,JPG)
REAL_B :: KB(5,13:59,JPG)
REAL_B :: SELFREF(10,JPG),FORREF(3,JPG)
REAL_B :: SFLUXREF(JPG,9)
REAL_B :: RAYL,STRRAT
INTEGER_B :: LAYREFFR

REAL_B :: KAC(9,5,13,NG19) ,ABSA(585,NG19)
REAL_B :: KBC(5,13:59,NG19),ABSB(235,NG19)
REAL_B :: SELFREFC(10,NG19),FORREFC(3,NG19)
REAL_B :: SFLUXREFC(NG19,9)

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
END MODULE YOESRTA19
