MODULE YOESRTAER


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOESRTAER* - AEROSOL OPTICAL PROPERTIES FOR SRTM 
!     ------------------------------------------------------------------


!-- AEROSOL OPTICAL PROPERTIES

!-- coefficients and parameters related to SRTM_SW 16 sp.int.

REAL_B :: RSRTAUA(14,6)
REAL_B :: RSRPIZA(14,6)
REAL_B :: RSRASYA(14,6)

!     -----------------------------------------------------------------

!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      03/03/07

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : -------
!RSRTAUA :  REAL     S.W. NORMALIZED OPTICAL THICKNESS AT 0.55 MICRON
!RSRPIZA :  REAL     S.W. SINGLE SCATTERING ALBEDO
!RSRASYA :  REAL     S.W. ASSYMETRY FACTOR

!     -----------------------------------------------------------------
END MODULE YOESRTAER
