!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2005, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!

SUBROUTINE RRTMG_SW_CLDPROP &
  &( KLEV, ICLDATM, INFLAG, ICEFLAG, LIQFLAG, NSTR &
  &, CLDFRAC, CLDDAT1, CLDDAT2, CLDDAT3, CLDDAT4, CLDDATMOM &
  &, TAUCLDORIG, TAUCLOUD, SSACLOUD, XMOM &
  &)

!  PURPOSE:  COMPUTE THE CLOUD OPTICAL DEPTH(S) FOR EACH CLOUDY
!  LAYER.  NOTE:  ONLY INFLAG = 0 AND INFLAG=2/LIQFLAG=1/ICEFLAG=2,3
!  (HU & STAMNES, KEY, Q. FU) ARE IMPLEMENTED.


#include "tsmbkind.h"

USE PARSRTM ,  ONLY : JPLAY, JPBAND, JPB1, JPB2
USE YOESRTOP,  ONLY : EXTLIQ1, SSALIQ1, ASYLIQ1 &
                   &, EXTICE2, SSAICE2, ASYICE2 &
                   &, EXTICE3, SSAICE3, ASYICE3, FDLICE3 &
                   &, FDELTA &
		   &, ABSCOICE, EXTCOICE, SSACOICE, GICE, FORWICE &
		   &, ABSCOLIQ, EXTCOLIQ, SSACOLIQ, GLIQ, FORWLIQ

!     ------------------------------------------------------------------

IMPLICIT NONE

!-- real arguments

REAL_B :: CLDFRAC(JPLAY), TAUCLOUD(JPLAY,JPBAND), SSACLOUD(JPLAY,JPBAND)
REAL_B :: TAUCLDORIG(JPLAY,JPBAND), XMOM(0:16,JPLAY,JPBAND)
REAL_B :: CLDDAT1(JPLAY), CLDDAT2(JPLAY), CLDDAT3(JPLAY), CLDDAT4(JPLAY)
REAL_B :: CLDDATMOM(0:16,JPLAY)

!-- integer arguments

INTEGER_B :: KLEV
INTEGER_B :: ICLDATM, INFLAG, ICEFLAG, LIQFLAG, NSTR

!-- real locals

REAL_B :: EPS
REAL_B :: TAUCLDORIG_A, FFP, FFP1, FFPSSA, SSACLOUD_A, TAUCLOUD_A
REAL_B :: CWP, FICE, RADICE, FACTOR, FINT, FLIQ, RADLIQ
REAL_B :: TAUICEORIG, SCATICE, SSAICE, TAUICE &
       &, TAULIQORIG, SCATLIQ, SSALIQ, TAULIQ

!-- integer locals
INTEGER_B :: NCBANDS, NLAYERS
INTEGER_B :: IB, IB1, IB2, LAY, ISTR , INDEX


! HVRCLD = '$Revision$'

! Initialize

EPS = 1.E-06_JPRB

ICLDATM = 0
NCBANDS = 29
NLAYERS = KLEV
IB1     = JPB1
IB2     = JPB2

DO LAY = 1, NLAYERS
   DO IB = IB1 , IB2
      TAUCLOUD(LAY,IB) = _ZERO_
      SSACLOUD(LAY,IB) = _ZERO_
      DO ISTR = 0,NSTR
        XMOM(ISTR,LAY,IB) = _ZERO_
      END DO
   END DO
END DO

! Main layer loop
DO LAY = 1, NLAYERS

  IF (CLDFRAC(LAY) .GE. EPS) THEN
    ICLDATM = 1

!   Ice clouds and water clouds combined.
    IF (INFLAG .EQ. 0) THEN
       TAUCLDORIG_A = CLDDAT1(LAY)
       FFP = CLDDATMOM(NSTR,LAY)
       FFP1 = 1.0 - FFP
       FFPSSA = 1.0 - FFP * CLDDAT2(LAY)
       SSACLOUD_A = FFP1*CLDDAT2(LAY)/FFPSSA
       TAUCLOUD_A = FFPSSA*TAUCLDORIG_A

       DO IB = IB1 , IB2
         TAUCLDORIG(LAY,IB) = TAUCLDORIG_A
         SSACLOUD(LAY,IB) = SSACLOUD_A
         TAUCLOUD(LAY,IB) = TAUCLOUD_A

         DO ISTR = 0,NSTR
           XMOM(ISTR,LAY,IB) = (CLDDATMOM(ISTR,LAY) - FFP)/ &
     &     (FFP1)
         END DO
       END DO

!   Separate treatement of ice clouds and water clouds.
     ELSE IF(INFLAG .EQ. 2) THEN       
       CWP = CLDDAT1(LAY)
       FICE = CLDDAT2(LAY)
       RADICE = CLDDAT3(LAY)

!   Calculation of absorption coefficients due to ice clouds.
       IF (FICE .EQ. 0.0) THEN
         DO IB = IB1 , IB2
           EXTCOICE(IB) = 0.0_JPRB
           SSACOICE(IB) = 0.0_JPRB
           GICE(IB)     = 0.0_JPRB
           FORWICE(IB)  = 0.0_JPRB
         END DO

       ELSEIF (ICEFLAG .EQ. 2) THEN
         IF (RADICE .LT. 5.0 .OR. RADICE .GT. 131.) STOP 'ICE RADIUS OUT OF BOUNDS'
         FACTOR = (RADICE - 2.)/3.
         INDEX = INT(FACTOR)
         IF (INDEX .EQ. 43) INDEX = 42
         FINT = FACTOR - FLOAT(INDEX)
         DO IB = IB1, IB2
            EXTCOICE(IB) = FICE * (EXTICE2(INDEX,IB) + FINT * &
     &                     (EXTICE2(INDEX+1,IB) -  EXTICE2(INDEX,IB)))
            SSACOICE(IB) = SSAICE2(INDEX,IB) + FINT * &
     &                     (SSAICE2(INDEX+1,IB) -  SSAICE2(INDEX,IB))
            GICE(IB) = ASYICE2(INDEX,IB) + FINT * &
     &                     (ASYICE2(INDEX+1,IB) -  ASYICE2(INDEX,IB))
            FORWICE(IB) = GICE(IB)**NSTR
         END DO

       ELSE IF (ICEFLAG .EQ. 3) THEN
         IF (RADICE .LT. 5.0 .OR. RADICE .GT. 140.0) STOP 'ICE EFFECTIVE SIZE OUT OF BOUNDS'
         FACTOR = (RADICE - 2._JPRB)/3._JPRB
         INDEX = INT(FACTOR)
         IF (INDEX .EQ. 46) INDEX = 45
         FINT = FACTOR - FLOAT(INDEX)

         DO IB = IB1 , IB2
           EXTCOICE(IB) = FICE * (EXTICE3(INDEX,IB) + FINT * &
     &                    (EXTICE3(INDEX+1,IB) - EXTICE3(INDEX,IB)))
           SSACOICE(IB) = SSAICE3(INDEX,IB) + FINT * &
     &                    (SSAICE3(INDEX+1,IB) - SSAICE3(INDEX,IB))
           GICE(IB) = ASYICE3(INDEX,IB) + FINT * &
     &                    (ASYICE3(INDEX+1,IB) - ASYICE3(INDEX,IB))
           FDELTA(IB) = FDLICE3(INDEX,IB) + FINT * &
     &                    (FDLICE3(INDEX+1,IB) - FDLICE3(INDEX,IB))
           if (fdelta(ib) .lt. 0.0) STOP 'FDELTA LESS THAN 0.0'
           if (fdelta(ib) .gt. 1.0) STOP 'FDELTA GT THAN 1.0'                     
           FORWICE(IB) = FDELTA(IB) + 0.5 / SSACOICE(IB)
! See Fu 1996 p. 2067 
           IF (FORWICE(IB) .GT. GICE(IB)) FORWICE(IB) = GICE(IB)
! Check to ensure all calculated quantities are within physical limits.
           if (extcoice(ib) .lt. 0.0_JPRB) STOP 'ICE EXTINCTION LESS THAN 0.0'
           if (ssacoice(ib) .gt. 1.0_JPRB) STOP 'ICE SSA GRTR THAN 1.0'
           if (ssacoice(ib) .lt. 0.0_JPRB) STOP 'ICE SSA LESS THAN 0.0'
           if (gice(ib) .gt. 1.0_JPRB) STOP 'ICE ASYM GRTR THAN 1.0'
           if (gice(ib) .lt. 0.0_JPRB) STOP 'ICE ASYM LESS THAN 0.0'
         END DO

       ENDIF
                  
!  Calculation of absorption coefficients due to water clouds.
       FLIQ = 1. - FICE
       IF (FLIQ .EQ. 0.0) THEN
         DO IB = IB1 , IB2
           EXTCOLIQ(IB) = 0.0
           SSACOLIQ(IB) = 0.0
           GLIQ(IB) = 0.0
           FORWLIQ(IB) = 0.0
         END DO

       ELSE IF (LIQFLAG .EQ. 1) THEN
         RADLIQ = CLDDAT4(LAY)
         IF (RADLIQ .LT. 1.5 .OR. RADLIQ .GT. 60.) STOP 'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
         INDEX = INT(RADLIQ - 1.5)
         IF (INDEX .EQ. 0) INDEX = 1
         IF (INDEX .EQ. 58) INDEX = 57
         FINT = RADLIQ - 1.5 - FLOAT(INDEX)

         DO IB = IB1 , IB2
           EXTCOLIQ(IB) = FLIQ * (EXTLIQ1(INDEX,IB) + FINT * &
     &                    (EXTLIQ1(INDEX+1,IB) - (EXTLIQ1(INDEX,IB))))
           SSACOLIQ(IB) = SSALIQ1(INDEX,IB) + FINT * &
     &                    (SSALIQ1(INDEX+1,IB) - SSALIQ1(INDEX,IB))
           GLIQ(IB) = ASYLIQ1(INDEX,IB) + FINT * &
     &                    (ASYLIQ1(INDEX+1,IB) - ASYLIQ1(INDEX,IB))
           FORWLIQ(IB) = GLIQ(IB)**NSTR
! Check to ensure all calculated quantities are within physical limits.
           if (extcoliq(ib) .lt. 0.0_JPRB) STOP 'LIQUID EXTINCTION LESS THAN 0.0'
           if (ssacoliq(ib) .gt. 1.0_JPRB) STOP 'LIQUID SSA GRTR THAN 1.0'
           if (ssacoliq(ib) .lt. 0.0_JPRB) STOP 'LIQUID SSA LESS THAN 0.0'
           if (gliq(ib) .gt. 1.0_JPRB) STOP 'LIQUID ASYM GRTR THAN 1.0'
           if (gliq(ib) .lt. 0.0_JPRB) STOP 'LIQUID ASYM LESS THAN 0.0'
         END DO

       END IF

       DO IB = IB1 , IB2
         TAULIQORIG = CWP * EXTCOLIQ(IB)
         TAUICEORIG = CWP * EXTCOICE(IB)
         TAUCLDORIG(LAY,IB) = TAULIQORIG + TAUICEORIG

         SSALIQ = SSACOLIQ(IB) * (1. - FORWLIQ(IB)) / &
     &                 (1. - FORWLIQ(IB) * SSACOLIQ(IB))
         TAULIQ =  (1. - FORWLIQ(IB) * SSACOLIQ(IB)) * &
     &                 TAULIQORIG
         SSAICE = SSACOICE(IB) * (1. - FORWICE(IB)) / &
     &                 (1. - FORWICE(IB) * SSACOICE(IB))
         TAUICE =  (1. - FORWICE(IB) * SSACOICE(IB)) * &
     &                 TAUICEORIG
         SCATLIQ = SSALIQ * TAULIQ
         SCATICE = SSAICE * TAUICE
         TAUCLOUD(LAY,IB) = TAULIQ + TAUICE
         SSACLOUD(LAY,IB) = (SCATLIQ + SCATICE) / &
     &                 TAUCLOUD(LAY,IB)
         XMOM(0,LAY,IB) = 1.0

         IF (ICEFLAG .EQ. 3) THEN
!     In accordance with the 1996 Fu paper, equation A.3, 
!     the moments for ice were calculated depending on whether using spheres
!     or hexagonal ice crystals.
             DO ISTR = 1, NSTR
                XMOM(ISTR,LAY,IB) = (1.0/(SCATLIQ+SCATICE))* &
     &              (SCATLIQ*(GLIQ(IB)**ISTR - FORWLIQ(IB)) / &
     &              (1. - FORWLIQ(IB)) + SCATICE * ((GICE(IB)-FORWICE(IB))/ &
     &              (1.0 - FORWICE(IB)))**ISTR)
             END DO
         ELSE 
!     This code is the standard method for delta-m scaling. 
             DO ISTR = 1, NSTR
                XMOM(ISTR,LAY,IB) = (SCATLIQ *  &
     &              (GLIQ(IB)**ISTR - FORWLIQ(IB)) / &
     &              (1. - FORWLIQ(IB)) + SCATICE * (GICE(IB)**ISTR - FORWICE(IB)) / &
     &              (1. - FORWICE(IB)))/(SCATLIQ + SCATICE)
             END DO
         ENDIF 
       END DO

     ENDIF

   ENDIF

 END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE RRTMG_SW_CLDPROP

