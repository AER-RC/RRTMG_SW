MODULE YOESRTA17

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA17* - SRTM COEFFICIENTS FOR INTERVAL 17
!     BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: JPG = 16, NG17 = 16

REAL_B :: KA(9,5,13,JPG)
REAL_B :: KB(5,5,13:59,JPG)
REAL_B :: SELFREF(10,JPG),FORREF(4,JPG)
REAL_B :: SFLUXREF(JPG,5)
REAL_B :: RAYL,STRRAT
INTEGER_M :: LAYREFFR

REAL_B :: KAC(9,5,13,NG17),ABSA(585,NG17)
REAL_B :: KBC(5,5,13:59,NG17),ABSB(1175,NG17)
REAL_B :: SELFREFC(10,NG17),FORREFC(4,NG17)
REAL_B :: SFLUXREFC(NG17,5)

EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,1,13,1),ABSB(1,1))

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
END MODULE YOESRTA17
