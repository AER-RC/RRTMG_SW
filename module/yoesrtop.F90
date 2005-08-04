MODULE YOESRTOP


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTOP* - SRTM CLOUD OPTICAL PROPERTIES
!     OVER BANDS 16 to 29:  2600-50000, then 800-2600 cm-1
!     -----------------------------------------------------------------

REAL_B :: EXTLIQ1(58,16:29), SSALIQ1(58,16:29), ASYLIQ1(58,16:29)
REAL_B :: EXTICE2(43,16:29), SSAICE2(43,16:29), ASYICE2(43,16:29)
REAL_B :: EXTICE3(46,16:29), SSAICE3(46,16:29), ASYICE3(46,16:29)
REAL_B :: FDLICE3(46,16:29), FDELTA(16:29)
REAL_B :: ABSCLD1

REAL_B :: ABSCOICE(16:29), EXTCOICE(16:29), GICE(16:29), SSACOICE(16:29), FORWICE(16:29)
REAL_B :: ABSCOLIQ(16:29), EXTCOLIQ(16:29), GLIQ(16:29), SSACOLIQ(16:29), FORWLIQ(16:29)


!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     J. DELAMERE/M. IACONO AER             08/01/2005

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! xxxLIQ1 : REAL   : OPTICAL PROPERTIES (EXTINCTION COEFFICIENT, SINGLE 
!                    SCATTERING ALBEDO, ASSYMETRY FACTOR) FROM
!                    HU & STAMNES, 1993, J. CLIM., 6, 728-742  
! xxxICE2 : REAL   : OPTICAL PROPERTIES (EXTINCTION COEFFICIENT, SINGLE 
!                    SCATTERING ALBEDO, ASSYMETRY FACTOR) FROM STREAMER V3.0,
!                    KEY, J., STREAMER USER'S GUIDE, COOPERATIVE INSTITUDE 
!                    FOR METEOROLOGICAL STUDIES, 2001, 95 PP.
! xxxICE3 : REAL   : OPTICAL PROPERTIES (EXTINCTION COEFFICIENT, SINGLE 
!                    SCATTERING ALBEDO, ASSYMETRY FACTOR) FROM
!                    FU, 1996, J. CLIM., 9, 
!     -----------------------------------------------------------------
END MODULE YOESRTOP
