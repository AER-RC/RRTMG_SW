MODULE YOESRTWN


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

INTEGER_M , DIMENSION(16:29) :: NG
INTEGER_M , DIMENSION(16:29) :: NSPA
INTEGER_M , DIMENSION(16:29) :: NSPB

REAL_B , DIMENSION(16:29) :: WAVENUM1
REAL_B , DIMENSION(16:29) :: WAVENUM2
REAL_B , DIMENSION(16:29) :: DELWAVE

REAL_B, DIMENSION(59) :: PREF
REAL_B, DIMENSION(59) :: PREFLOG
REAL_B, DIMENSION(59) :: TREF

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----    : -------
!  NG     : INTEGER : Number of k-coefficients in spectral intervals
!  NSPA   : INTEGER :
!  NSPB   : INTEGER :
! WAVENUM1: REAL    : Lower wavenumber spectral limit
! WAVENUM2: REAL    : Higher wavenumber spectral limit
! DELWAVE : REAL    : Spectral interval width
! PREF    : REAL    : Reference pressure
! PREFLOG : REAL    : Log reference pressure
! TREF    : REAL    : Reference temperature
!     -----------------------------------------------------------------
END MODULE YOESRTWN
