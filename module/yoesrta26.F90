MODULE YOESRTA26

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA26* - SRTM COEFFICIENTS FOR INTERVAL 26
!     BAND 26: 22650-29000 cm-1 (low - nothing; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: JPG = 16, NG26 = 16

REAL_B :: SFLUXREF(JPG), RAYL(JPG)

REAL_B :: SFLUXREFC(NG26), RAYLC(NG26)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! SFLUXREF: REAL    
! RAYL    : REAL 
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
!     -----------------------------------------------------------------
END MODULE YOESRTA26
