MODULE YOESRTA23

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA23* - SRTM COEFFICIENTS FOR INTERVAL 23
!     BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: JPG = 16, NG23 = 16

REAL_B :: KA(5,13,JPG)
REAL_B :: SELFREF(10,JPG),FORREF(3,JPG)
REAL_B :: SFLUXREF(JPG)  ,RAYL(JPG)
REAL_B :: GIVFAC
INTEGER_M :: LAYREFFR

REAL_B :: KAC(5,13,NG23)   ,ABSA(65,NG23)
REAL_B :: SELFREFC(10,NG23),FORREFC(3,NG23)
REAL_B :: SFLUXREFC(NG23)  ,RAYLC(NG23)

EQUIVALENCE (KAC(1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! SELFREF : REAL 
! FORREF  : REAL    
! SFLUXREF: REAL
! RAYL    : REAL
! KAC     : REAL     Reduced g-point array for KA
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
! GIVFAC  : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA23
