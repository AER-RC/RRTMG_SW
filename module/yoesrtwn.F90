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

INTEGER_M , DIMENSION(224) :: NGM
INTEGER_M , DIMENSION(14) :: NGC, NGS
! Use for 112 g-points
INTEGER_M , DIMENSION(112) :: NGB, NGN
! Use for 224 g-points
!INTEGER_M , DIMENSION(224) :: NGB, NGN

REAL_B , DIMENSION(16) :: WT, WTSM
REAL_B , DIMENSION(224) :: RWGT

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15
!     M. J. IACONO          AER             12/09/03

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
!  NGC    : INTEGER : The number of new g-points in each band
!  NGS    : INTEGER : The cumulative sum of new g-points for each band
!  NGM    : INTEGER : The index of each new g-point relative to the
!                     original 16 g-points for each band.
!  NGN    : INTEGER : The number of original g-points that are combined 
!                     to make each new g-point in each band.
!  NGB    : INTEGER : The band index for each new g-point.
!  WT     : REAL    : RRTM weights for 16 g-points.
!  WTSUM  : REAL    : Sum of the weights
!  RWGT   : REAL    : 
!     -----------------------------------------------------------------
END MODULE YOESRTWN
