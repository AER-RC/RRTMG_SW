C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

SUBROUTINE RRTM_SW_REFTRA &
  &( KLEV  , KMODTS &
  &, LRTCHK &
  &, PGG   , PRMUZ, PTAU , PW    &
  &, PREF  , PREFD, PTRA , PTRAD &
  &)
  
!**** *RRTM_SW_REFTRA* - REFLECTIVITY AND TRANSMISSIVITY

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLEAR OR 
!     CLOUDY LAYER USING A CHOICE OF VARIOUS APPROXIMATIONS.
!
!**   INTERFACE.
!     ----------
!          *RRTM_SW_REFTRA* IS CALLED BY *RRTM_SW_SPCVRT*
!
!
!        EXPLICIT ARGUMENTS :
!        --------------------
! INPUTS
! ------ 
!      KMODTS  = 1 EDDINGTON (JOSEPH ET AL., 1976)
!              = 2 PIFM (ZDUNKOWSKI ET AL., 1980)
!              = 3 DISCRETE ORDINATES (LIOU, 1973)
!      LRTCHK  = .T. IF CLOUDY
!              = .F. IF CLEAR-SKY
!      PGG     = ASSYMETRY FACTOR
!      PRMUZ   = COSINE SOLAR ZENITH ANGLE
!      PTAU    = OPTICAL THICKNESS
!      PW      = SINGLE SCATTERING ALBEDO
!
! OUTPUTS
! -------
!      PREF    : COLLIMATED BEAM REFLECTIVITY
!      PREFD   : DIFFUSE BEAM REFLECTIVITY 
!      PTRA    : COLLIMATED BEAM TRANSMISSIVITY
!      PTRAD   : DIFFUSE BEAM TRANSMISSIVITY
!
!
!     METHOD.
!     -------
!          STANDARD DELTA-EDDINGTON, P.I.F.M., OR D.O.M. LAYER CALCULATIONS.
!
!
!     EXTERNALS.
!     ----------
!          NONE
!
!
!     REFERENCE.
!     ----------
!
!
!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 03-02-27
!
!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

#include "tsmbkind.h"

USE PARSRTM , ONLY : JPLAY
!USE YOESW   , ONLY : NDBUG

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER_M :: KIDIA, KFDIA, KLON
INTEGER_M :: KLEV
INTEGER_M :: KMODTS

LOGICAL :: LRTCHK(JPLAY)

REAL_B :: PGG(JPLAY), PRMUZ, PTAU(JPLAY), PW(JPLAY)

REAL_B :: PREF(JPLAY),PREFD(JPLAY), PTRA(JPLAY), PTRAD(JPLAY)

!     ------------------------------------------------------------------

!     LOCAL INTEGER SCALARS
INTEGER_M :: JK, JL
!INTEGER_M :: NDBUG

REAL_B :: ZA, ZA1, ZA2
REAL_B :: ZBETA, ZDEND, ZDENR, ZDENT
REAL_B :: ZE1, ZE2, ZEM1, ZEM2, ZEMM, ZEP1, ZEP2
REAL_B :: ZG, ZG3, ZGAMMA1, ZGAMMA2, ZGAMMA3, ZGAMMA4, ZGT
REAL_B :: ZR1, ZR2, ZR3, ZR4, ZR5, ZRK, ZRK2, ZRKG, ZRM1, ZRP, ZRP1, ZRPP
REAL_B :: ZSR3, ZT1, ZT2, ZT3, ZT4, ZT5, ZTO1
REAL_B :: ZW, ZWCRIT, ZG, ZTO1, ZW, ZWO

!     ------------------------------------------------------------------

ZSR3=SQRT(3._JPRB)
ZWCRIT=0.9995_JPRB
KMODTS=2

!NDBUG=3

DO JK=1,KLEV
!  if (NDBUG < 2) then
!    print 9000,JK,LRTCHK(JK),PTAU(JK),PW(JK),PGG(JK),PRMUZ
9000 format(1x,'RRTM_SW_REFTRA:inputs:',I3,L8,4E13.6)
!  end if
  IF (.NOT.LRTCHK(JK)) THEN
    PREF(JK) =_ZERO_
    PTRA(JK) =_ONE_
    PREFD(JK)=_ZERO_
    PTRAD(JK)=_ONE_
!    if (NDBUG < 2) then
!      print 9001,JL,JK,PREF(JK),PTRA(JK),PREFD(JK),PTRAD(JK)
9001  format(1x,'RRTM_SW_REFTRA:not.LRTCKH:',2I3,4F10.6)
!    end if
  ELSE
    ZTO1=PTAU(JK)
    ZW  =PW(JK)
    ZG  =PGG(JK)  
!
!-- GENERAL TWO-STREAM EXPRESSIONS
!
    ZG3= 3._JPRB * ZG
    IF (KMODTS == 1) THEN
      ZGAMMA1= (7._JPRB - ZW * (4._JPRB + ZG3)) * 0.25_JPRB
      ZGAMMA2=-(1._JPRB - ZW * (4._JPRB - ZG3)) * 0.25_JPRB
      ZGAMMA3= (2._JPRB - ZG3 * PRMUZ ) * 0.25_JPRB
    ELSE IF (KMODTS == 2) THEN  
      ZGAMMA1= (8._JPRB - ZW * (5._JPRB + ZG3)) * 0.25_JPRB
      ZGAMMA2=  3._JPRB *(ZW * (1._JPRB - ZG )) * 0.25_JPRB
      ZGAMMA3= (2._JPRB - ZG3 * PRMUZ ) * 0.25_JPRB
    ELSE IF (KMODTS == 3) THEN  
      ZGAMMA1= ZSR3 * (2._JPRB - ZW * (1._JPRB + ZG)) * 0.5_JPRB
      ZGAMMA2= ZSR3 * ZW * (1._JPRB - ZG ) * 0.5_JPRB
      ZGAMMA3= (1._JPRB - ZSR3 * ZG * PRMUZ ) * 0.5_JPRB
    END IF
    ZGAMMA4= 1._JPRB - ZGAMMA3
    
!-- RECOMPUTE ORIGINAL S.S.A. TO TEST FOR CONSERVATIVE SOLUTION
    ZWO= ZW / (1._JPRB - (1._JPRB - ZW) * (ZG / (1._JPRB - ZG))**2)
    
    IF (ZWO >= ZWCRIT) THEN
!!-- conservative scattering
!      
      ZA  = ZGAMMA1 * PRMUZ 
      ZA1 = ZA - ZGAMMA3
      ZGT = ZGAMMA1 * ZTO1
        
!-- Homogeneous reflectance and transmittance
!
! collimated beam
!     
      ZE1 = MIN ( ZTO1 / PRMUZ , 500._JPRB)
      ZE2 = EXP ( - ZE1 )
      PREF(JK) = (ZGT - ZA1 * (1._JPRB - ZE2)) / (1._JPRB + ZGT)
      PTRA(JK) = 1._JPRB - PREF(JK)
!
! isotropic incidence
!
      PREFD(JK) = ZGT / (1._JPRB + ZGT)
      PTRAD(JK) = 1._JPRB - PREFD(JK)        
!    
!      if (NDBUG < 2) then
!        print 9002,JL,JK,PREF(JK),PTRA(JK),PREFD(JK),PTRAD(JK)
9002  format(1x,'RRTM_SW_REFTRA: consrv: LRTCHK:',2I3,4F10.6)
!      end if
!        
    ELSE
!
!-- non-conservative scattering
!
      ZA1 = ZGAMMA1 * ZGAMMA4 + ZGAMMA2 * ZGAMMA3
      ZA2 = ZGAMMA1 * ZGAMMA3 + ZGAMMA2 * ZGAMMA4
      ZRK = SQRT ( ZGAMMA1**2 - ZGAMMA2**2)
      ZRP = ZRK * PRMUZ               
      ZRP1 = 1._JPRB + ZRP
      ZRM1 = 1._JPRB - ZRP
      ZRK2 = 2._JPRB * ZRK
      ZRPP = 1._JPRB - ZRP*ZRP
      ZRKG = ZRK + ZGAMMA1
      ZR1  = ZRM1 * (ZA2 + ZRK * ZGAMMA3)
      ZR2  = ZRP1 * (ZA2 - ZRK * ZGAMMA3)
      ZR3  = ZRK2 * (ZGAMMA3 - ZA2 * PRMUZ )
      ZR4  = ZRPP * ZRKG
      ZR5  = ZRPP * (ZRK - ZGAMMA1)
      ZT1  = ZRP1 * (ZA1 + ZRK * ZGAMMA4)
      ZT2  = ZRM1 * (ZA1 - ZRK * ZGAMMA4)
      ZT3  = ZRK2 * (ZGAMMA4 + ZA1 * PRMUZ )
      ZT4  = ZR4
      ZT5  = ZR5
      ZBETA = - ZR5 / ZR4
        
!-- Homogeneous reflectance and transmittance
!
      ZE1 = MIN ( ZRK * ZTO1, 500._JPRB)
      ZE2 = MIN ( ZTO1 / PRMUZ , 500._JPRB)
      ZEP1 = EXP( ZE1 )
      ZEM1 = EXP(-ZE1 )
      ZEP2 = EXP( ZE2 )
      ZEM2 = EXP(-ZE2 )
!        
! collimated beam
!     
      ZDENR = ZR4*ZEP1 + ZR5*ZEM1
      PREF(JK) = ZWO * (ZR1*ZEP1 - ZR2*ZEM1 - ZR3*ZEM2) / ZDENR
      ZDENT = ZT4*ZEP1 + ZT5*ZEM1
      PTRA(JK) = ZEM2 * (1._JPRB - ZWO * (ZT1*ZEP1 - ZT2*ZEM1 - ZT3*ZEP2) / ZDENT)
!
! diffuse beam
!
      ZEMM = ZEM1*ZEM1
      ZDEND = 1._JPRB / ( (1._JPRB - ZBETA*ZEMM ) * ZRKG)
      PREFD(JK) =  ZGAMMA2 * (1._JPRB - ZEMM) * ZDEND
      PTRAD(JK) =  ZRK2*ZEM1*ZDEND
!
!      if (NDBUG < 2) then        
!        print 9003,JL,JK,PREF(JK),PTRA(JK),PREFD(JK),PTRAD(JK)
9003  format(1x,'RRTM_SW_REFTRA: OMG<1:  LRTCHK:',2I3,4F10.6)
!      end if
    END IF
!
  END IF         
!    
END DO    
!
!     ------------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_SW_REFTRA     
