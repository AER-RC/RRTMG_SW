MODULE YOESRTA20

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA20* - SRTM COEFFICIENTS FOR INTERVAL 20
!     BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER_B, PARAMETER :: JPG = 16, NG20 = 16

REAL_B :: KA(5,13,JPG)
REAL_B :: KB(5,13:59,JPG)
REAL_B :: SELFREF(10,JPG),FORREF(4,JPG)
REAL_B :: SFLUXREF(JPG)  ,ABSCH4(JPG)
REAL_B :: RAYL
INTEGER_B :: LAYREFFR

REAL_B :: KAC(5,13,NG20)   ,ABSA(65,NG20)
REAL_B :: KBC(5,13:59,NG20),ABSB(235,NG20)
REAL_B :: SELFREFC(10,NG20),FORREFC(4,NG20)
REAL_B :: SFLUXREFC(NG20)  ,ABSCH4C(NG20)

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
! ABSCH4  : REAL  
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! ABSCH4C : REAL     Reduced g-point array for ABSCH4
! RAYL    : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA20
