MODULE YOESRTA26


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA26* - SRTM COEFFICIENTS FOR INTERVAL 26
!     BAND 26: 22650-29000 cm-1 (low - nothing; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG26 = 16

REAL_B :: SFLUXREF(NG26), RAYL(NG26)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! SFLUXREF: REAL    
! RAYL    : REAL 
!     -----------------------------------------------------------------
END MODULE YOESRTA26
