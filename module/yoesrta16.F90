MODULE YOESRTA16


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA16* - SRTM COEFFICIENTS FOR INTERVAL 16
!     BAND 16:  2600-3250 cm-1 (low - H2O,CH4; high - CH4)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG16 = 16, NGS15 = 0

REAL_B :: KA(9,5,13,NG16) ,ABSA(585,NG16)
REAL_B :: KB(5,13:59,NG16),ABSB(235,NG16)
REAL_B :: SELFREF(10,NG16),FORREF(3,NG16)
REAL_B :: SFLUXREF(NG16)
REAL_B :: RAYL            ,STRRAT1
INTEGER_M :: LAYREFFR

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL   
! SFLUXREF: REAL
! RAYL    : REAL 
! STRRAT1 : REAL
! LAYREFFR: INTEGER
!     -----------------------------------------------------------------
END MODULE YOESRTA16
