C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

SUBROUTINE RRTM_SW_CLDPROP &
  &( KLEV, ICLDATM, INFLAG, ICEFLAG, LIQFLAG, NSTR &
  &, CLDFRAC, CLDDAT1, CLDDAT2, CLDDAT3, CLDDAT4, CLDDATMOM &
  &, TAUCLDORIG, TAUCLOUD, SSACLOUD, XMOM &
  &)

!  path:      $Source$
!  author:    $Author$
!  revision:  $Revision$
!  created:   $Date$

!  PURPOSE:  COMPUTE THE CLOUD OPTICAL DEPTH(S) FOR EACH CLOUDY
!  LAYER.  NOTE:  ONLY INFLAG = 0 AND INFLAG=2/LIQFLAG=1/ICEFLAG=3
!  (HU & STAMNES, Q. FU) ARE IMPLEMENTED.


#include "tsmbkind.h"

USE PARSRTM ,  ONLY : JPLAY, JPBAND, JPB1, JPB2
USE YOESRTOP,  ONLY : EXTLIQ1, SSALIQ1, ASYLIQ1 &
                   &, EXTICE3, SSAICE3, ASYICE3, FDLICE3 &
                   &, FDELTA , ABSCLD1 &
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

INTEGER_M :: KLEV
INTEGER_M :: ICLDATM, INFLAG, ICEFLAG, LIQFLAG, NSTR

!-- real locals

REAL_B :: EPS
REAL_B :: TAUCLDORIG_A, FFP, FFP1, FFPSSA, SSACLOUD_A, TAUCLOUD_A
REAL_B :: CWP, FICE, RADICE, FACTOR, FINT, FLIQ, RADLIQ
REAL_B :: TAUICEORIG, SCATICE, SSAICE, TAUICE &
       &, TAULIQORIG, SCATLIQ, SSALIQ, TAULIQ

!-- integer locals
INTEGER_M :: NCBANDS, NLAYERS
INTEGER_M :: IB, IB1, IB2, LAY, ISTR , INDEX
INTEGER_M :: NDBUG
 

! HVRCLD = '$Revision$'

EPS = 1.E-06_JPRB
NDBUG = 4

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

!IF (NDBUG.LE.2) THEN
!  print *,'cldprop before loop INFLAG, ICEFLAG, LIQFLAG, IB1, IB2:',INFLAG,ICEFLAG,LIQFLAG,IB1,IB2
!END IF

DO LAY = 1, NLAYERS

  IF (CLDFRAC(LAY) .GE. EPS) THEN
    ICLDATM = 1

!    IF (NDBUG.LE.2) THEN
!      print 9101,LAY,ICLDATM,CLDFRAC(LAY),CLDDAT1(LAY),CLDDAT2(LAY),CLDDAT3(LAY)&
!       &,CLDDAT4(LAY),(CLDDATMOM(ISTR,LAY),ISTR=0,NSTR)
!9101  format(1x,'Cld :',2I3,f7.4,7E12.5)
!    END IF

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

!       IF (NDBUG.LE.1) THEN
!         print 9102,LAY,CWP,FICE,RADICE
!9102     format(1x,'A',I3,3E13.6)
!       END IF

!   Calculation of absorption coefficients due to ice clouds.
       IF (FICE .EQ. 0.0) THEN
         DO IB = IB1 , IB2
           EXTCOICE(IB) = 0.0_JPRB
           SSACOICE(IB) = 1.0_JPRB
           GICE(IB)     = 1.0_JPRB
           FORWICE(IB)  = 0.0_JPRB

!	   IF (NDBUG.LE.1) THEN
!             print 9103,LAY,FICE,CWP,RADICE,IB,EXTCOICE(IB),SSACOICE(IB),GICE(IB),FORWICE(IB)
!9103         format(1x,'B',I3,F6.3,2E13.6,I3,4E12.5)
!	   END IF

         END DO

       ELSE IF (ICEFLAG .EQ. 3) THEN
         IF (RADICE .LT. 10.0 .OR. RADICE .GT. 140.0) STOP 'ICE EFFECTIVE SIZE OUT OF BOUNDS'
         FACTOR = (RADICE - 5._JPRB)/5._JPRB
         INDEX = INT(FACTOR)
         IF (INDEX .EQ. 27) INDEX = 26
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

!   	   IF (NDBUG.LE.1) THEN
!             print 9104,LAY,FICE,CWP,RADICE,IB,EXTCOICE(IB),SSACOICE(IB),GICE(IB),FORWICE(IB),FDELTA(IB)
!9104         format(1x,'C',I3,F5.3,2E13.6,I3,5E12.5)
!           END IF

         END DO
       ENDIF
!       IF (NDBUG.LE.1) THEN
!         print *,'end of ice computations for LAY=',LAY
!       ENDIF
                  
!  Calculation of absorption coefficients due to water clouds.
       FLIQ = 1. - FICE
       IF (FLIQ .EQ. 0.0) THEN
         DO IB = IB1 , IB2
           EXTCOLIQ(IB) = 0.0
           SSACOLIQ(IB) = 1.0
           GLIQ(IB) = 1.0
           FORWLIQ(IB) = 0.0

!	   IF (NDBUG.LE.1) THEN
!             print 9105,LAY,FLIQ,CWP,IB,EXTCOLIQ(IB),SSACOLIQ(IB),GLIQ(IB),FORWLIQ(IB)
!9105         format(1x,'D',I3,F5.3,1E13.6,I3,4E12.5)
!           END IF

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

!	   IF (NDBUG.LE.1) THEN
!             print 9106,LAY,FLIQ,CWP,RADLIQ,IB,EXTCOLIQ(IB),SSACOLIQ(IB),GLIQ(IB),FORWLIQ(IB)
!9106         format(1x,'E',I3,F5.3,2E13.6,I3,5E12.5)
!           END IF

         END DO
       END IF

!       IF (NDBUG.LE.1) THEN
!         print *,'end of liquid water computations for LAY=',LAY
!       END IF

               
       DO IB = IB1 , IB2
         TAULIQORIG = CWP * EXTCOLIQ(IB)
         TAUICEORIG = CWP * EXTCOICE(IB)
         TAUCLDORIG(LAY,IB) = TAULIQORIG + TAUICEORIG

!	 IF (NDBUG.LE.1) THEN
!           print 9107,IB,TAULIQORIG,TAUICEORIG,TAUCLDORIG(LAY,IB),CWP &
!            &,EXTCOLIQ(IB),EXTCOICE(IB),SSACOLIQ(IB),SSACOICE(IB) &
!	    &,FORWLIQ(IB),FORWICE(IB)
!9107       format(1x,'F',I3,10E12.5)
!         END IF

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

!	 IF (NDBUG.LE.1) THEN
!           print 9108,IB,TAULIQORIG,TAUICEORIG,SSALIQ,TAULIQ,SCATLIQ,SSAICE,TAUICE,SCATICE
!9108       format(1x,'G',I3,8E13.6)
!         END IF

         DO ISTR = 1, NSTR
!This commented code is the standard method for delta-m scaling. In accordance
!  with the 1996 Fu paper, equation A.3, the moments for ice were calculated
!  as in the uncommented code.
!                     XMOM(ISTR,LAY,IB) = (SCATLIQ *  &
!     &                    (GLIQ(IB)**ISTR - FORWLIQ(IB)) / &
!     &                    (1. - FORWLIQ(IB)) &
!     &                    + SCATICE * &
!     &                    (GICE(IB)**ISTR - FORWICE(IB)) /  &
!     &                    (1. - FORWICE(IB)))/(SCATLIQ + SCATICE)

           XMOM(ISTR,LAY,IB) = (1.0/(scatliq+scatice))* &
     &                    (SCATLIQ*(GLIQ(IB)**ISTR - FORWLIQ(IB)) / &
     &                    (1. - FORWLIQ(IB)) &
     &                    + SCATICE * &
     &               ((gice(ib)-forwice(ib))/(1.0-forwice(ib)))**ISTR)
         END DO

!	 IF (NDBUG.LE.1) THEN
!           print 9109,IB,TAUCLOUD(LAY,IB),SSACLOUD(LAY,IB),XMOM(1,LAY,IB)
!9109       format(1x,'H',I3,3E13.6)
!         END IF


       END DO

     ENDIF

   ENDIF

 END DO

! IF (NDBUG.LE.1) THEN
!   print *,'about to leave RRTM_SW_CLDPROP'
! END IF

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_SW_CLDPROP

