MODULE YOESRTA25

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA25* - SRTM COEFFICIENTS FOR INTERVAL 25
!     BAND 25: 16000-22650 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: JPG = 16, NG25 = 16

REAL_B :: KA(5,13,JPG)
REAL_B :: SFLUXREF(JPG)
REAL_B :: RAYL(JPG), ABSO3A(JPG), ABSO3B(JPG)
INTEGER_M :: LAYREFFR

REAL_B :: KAC(5,13,NG25) ,ABSA(65,NG25)
REAL_B :: SFLUXREFC(NG25)
REAL_B :: RAYLC(NG25), ABSO3AC(NG25), ABSO3BC(NG25)

EQUIVALENCE (KAC(1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! SFLUXREF: REAL
! RAYL    : REAL
! ABSO3A  : REAL
! ABSO3B  : REAL
! KAC     : REAL     Reduced g-point array for KA
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
! ABSO3AC : REAL     Reduced g-point array for ABSO3A
! ABSO3BC : REAL     Reduced g-point array for ABSO3B
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA25

