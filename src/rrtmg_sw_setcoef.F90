C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

SUBROUTINE RRTM_SW_SETCOEF &
  &( KLEV   , NMOL    &
  &, PAVEL  , TAVEL   , PZ     , TZ     , TBOUND  &
  &, COLDRY , WKL     &
  &, LAYTROP, LAYSWTCH, LAYLOW &
  &, CO2MULT, COLCH4  , COLCO2 , COLH2O , COLMOL  , COLN2O  , COLO2 , COLO3 &
  &, FORFAC , FORFRAC , INDFOR , SELFFAC, SELFFRAC, INDSELF &
  &, FAC00  , FAC01   , FAC10  , FAC11  &
  &, JP     , JT      , JT1    &
  &)


!     J. Delamere, AER, Inc. (version 2.5, 02/04/01)
!
!     Modifications:
!     JJMorcrette 030224   rewritten / adapted to ECMWF F90 system


!     Purpose:  For a given atmosphere, calculate the indices and
!     fractions related to the pressure and temperature interpolations.


#include "tsmbkind.h"

USE PARSRTM , ONLY : JPLAY, JPG
USE YOESRTWN, ONLY : PREF, PREFLOG, TREF
!  USE YOESWN  , ONLY : NDBUG
				   
IMPLICIT NONE

!-- Input arguments

INTEGER_M :: KLEV, NMOL

REAL_B :: PAVEL(JPLAY) , TAVEL(JPLAY) , PZ(0:JPLAY)  , TZ(0:JPLAY)
REAL_B :: COLDRY(JPLAY), COLMOL(JPLAY), WKL(35,JPLAY)


!-- Output arguments

INTEGER_M :: LAYTROP       , LAYSWTCH     , LAYLOW
INTEGER_M :: JP(JPLAY)     , JT(JPLAY)    , JT1(JPLAY)
INTEGER_M :: INDSELF(JPLAY), INDFOR(JPLAY)

REAL_B :: COLH2O(JPLAY), COLCO2(JPLAY) , COLN2O(JPLAY)
REAL_B :: COLO2(JPLAY) , COLO3(JPLAY)  , COLCH4(JPLAY) , CO2MULT(JPLAY)
REAL_B :: FAC00(JPLAY) , FAC01(JPLAY)  , FAC10(JPLAY)  , FAC11(JPLAY)
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
REAL_B :: TBOUND


!-- local integers

INTEGER_M :: NLAYERS, INDBOUND, INDLEV0, LAY
INTEGER_M :: JP1

!-- local reals

REAL_B :: STPFAC, TBNDFRAC, T0FRAC, PLOG, FP, FT, FT1, WATER, SCALEFAC
REAL_B :: FACTOR, CO2REG, COMPFP

NLAYERS = KLEV

STPFAC = 296._JPRB/1013._JPRB

INDBOUND = TBOUND - 159._JPRB
TBNDFRAC = TBOUND - INT(TBOUND)
INDLEV0  = TZ(0) - 159._JPRB
T0FRAC   = TZ(0) - INT(TZ(0))

LAYTROP  = 0
LAYSWTCH = 0
LAYLOW   = 0

!IF (NDBUG.LE.3) THEN
!  print *,'-------- Computed in SETCOEF --------'
!  print 8990
!8990 format(18x,'  T     FAC00,    01,    10,    11  CO2MULT     MOL   &
!    &    CH4      CO2      H2O      N2O      O2      O3      SFAC  &
!    &    SFRAC    FFAC    FFRAC  ISLF IFOR')
!END IF

DO LAY = 1, NLAYERS
!        Find the two reference pressures on either side of the
!        layer pressure.  Store them in JP and JP1.  Store in FP the
!        fraction of the difference (in ln(pressure)) between these
!        two values that the layer pressure lies.

  PLOG = LOG(PAVEL(LAY))
  JP(LAY) = INT(36. - 5*(PLOG+0.04))
  IF (JP(LAY) .LT. 1) THEN
    JP(LAY) = 1
  ELSEIF (JP(LAY) .GT. 58) THEN
    JP(LAY) = 58
  ENDIF
  JP1 = JP(LAY) + 1
  FP = 5. * (PREFLOG(JP(LAY)) - PLOG)

!        Determine, for each reference pressure (JP and JP1), which
!        reference temperature (these are different for each  
!        reference pressure) is nearest the layer temperature but does
!        not exceed it.  Store these indices in JT and JT1, resp.
!        Store in FT (resp. FT1) the fraction of the way between JT
!        (JT1) and the next highest reference temperature that the 
!        layer temperature falls.

  JT(LAY) = INT(3. + (TAVEL(LAY)-TREF(JP(LAY)))/15.)
  IF (JT(LAY) .LT. 1) THEN
    JT(LAY) = 1
  ELSEIF (JT(LAY) .GT. 4) THEN
    JT(LAY) = 4
  ENDIF
  FT = ((TAVEL(LAY)-TREF(JP(LAY)))/15.) - FLOAT(JT(LAY)-3)
  JT1(LAY) = INT(3. + (TAVEL(LAY)-TREF(JP1))/15.)
  IF (JT1(LAY) .LT. 1) THEN
    JT1(LAY) = 1
  ELSEIF (JT1(LAY) .GT. 4) THEN
    JT1(LAY) = 4
  ENDIF
  FT1 = ((TAVEL(LAY)-TREF(JP1))/15.) - FLOAT(JT1(LAY)-3)

  WATER = WKL(1,LAY)/COLDRY(LAY)
  SCALEFAC = PAVEL(LAY) * STPFAC / TAVEL(LAY)

!        If the pressure is less than ~100mb, perform a different
!        set of species interpolations.

  IF (PLOG .LE. 4.56) GO TO 5300
  LAYTROP =  LAYTROP + 1
  IF (PLOG .GE. 6.62) LAYLOW = LAYLOW + 1

!        Set up factors needed to separately include the water vapor
!        foreign-continuum in the calculation of absorption coefficient.

  FORFAC(LAY) = SCALEFAC / (1.+WATER)
  FACTOR = (332.0-TAVEL(LAY))/36.0
  INDFOR(LAY) = MIN(2, MAX(1, INT(FACTOR)))
  FORFRAC(LAY) = FACTOR - FLOAT(INDFOR(LAY))

!
!        Set up factors needed to separately include the water vapor
!        self-continuum in the calculation of absorption coefficient.

  SELFFAC(LAY) = WATER * FORFAC(LAY)
  FACTOR = (TAVEL(LAY)-188.0)/7.2
  INDSELF(LAY) = MIN(9, MAX(1, INT(FACTOR)-7))
  SELFFRAC(LAY) = FACTOR - FLOAT(INDSELF(LAY) + 7)

!        Calculate needed column amounts.

  COLH2O(LAY) = 1.E-20 * WKL(1,LAY)
  COLCO2(LAY) = 1.E-20 * WKL(2,LAY)
  COLO3(LAY) = 1.E-20 * WKL(3,LAY)
!         COLO3(LAY) = 0.
!         COLO3(LAY) = colo3(lay)/1.16
  COLN2O(LAY) = 1.E-20 * WKL(4,LAY)
  COLCH4(LAY) = 1.E-20 * WKL(6,LAY)
  COLO2(LAY) = 1.E-20 * WKL(7,LAY)
  COLMOL(LAY) = 1.E-20 * COLDRY(LAY) + COLH2O(LAY)
!         colco2(lay) = 0.
!         colo3(lay) = 0.
!         coln2o(lay) = 0.
!         colch4(lay) = 0.
!         colo2(lay) = 0.
!         colmol(lay) = 0.
  IF (COLCO2(LAY) .EQ. 0.) COLCO2(LAY) = 1.E-32 * COLDRY(LAY)
  IF (COLN2O(LAY) .EQ. 0.) COLN2O(LAY) = 1.E-32 * COLDRY(LAY)
  IF (COLCH4(LAY) .EQ. 0.) COLCH4(LAY) = 1.E-32 * COLDRY(LAY)
  IF (COLO2(LAY) .EQ. 0.) COLO2(LAY) = 1.E-32 * COLDRY(LAY)
!        Using E = 1334.2 cm-1.
  CO2REG = 3.55E-24 * COLDRY(LAY)
  CO2MULT(LAY)= (COLCO2(LAY) - CO2REG) * &
    &        272.63*EXP(-1919.4/TAVEL(LAY))/(8.7604E-4*TAVEL(LAY))
  GO TO 5400

!        Above LAYTROP.
 5300    CONTINUE

!        Set up factors needed to separately include the water vapor
!        foreign-continuum in the calculation of absorption coefficient.

  FORFAC(LAY) = SCALEFAC / (1.+WATER)
  FACTOR = (TAVEL(LAY)-188.0)/36.0
  INDFOR(LAY) = 3
  FORFRAC(LAY) = FACTOR - 1.0
!
!        Calculate needed column amounts.

  COLH2O(LAY) = 1.E-20 * WKL(1,LAY)
  COLCO2(LAY) = 1.E-20 * WKL(2,LAY)
  COLO3(LAY)  = 1.E-20 * WKL(3,LAY)
  COLN2O(LAY) = 1.E-20 * WKL(4,LAY)
  COLCH4(LAY) = 1.E-20 * WKL(6,LAY)
  COLO2(LAY)  = 1.E-20 * WKL(7,LAY)
  COLMOL(LAY) = 1.E-20 * COLDRY(LAY) + COLH2O(LAY)
  IF (COLCO2(LAY) .EQ. 0.) COLCO2(LAY) = 1.E-32 * COLDRY(LAY)
  IF (COLN2O(LAY) .EQ. 0.) COLN2O(LAY) = 1.E-32 * COLDRY(LAY)
  IF (COLCH4(LAY) .EQ. 0.) COLCH4(LAY) = 1.E-32 * COLDRY(LAY)
  IF (COLO2(LAY)  .EQ. 0.) COLO2(LAY)  = 1.E-32 * COLDRY(LAY)
  CO2REG = 3.55E-24 * COLDRY(LAY)
  CO2MULT(LAY)= (COLCO2(LAY) - CO2REG) * &
    &        272.63*EXP(-1919.4/TAVEL(LAY))/(8.7604E-4*TAVEL(LAY))

  SELFFAC(LAY) =_ZERO_
  SELFFRAC(LAY)=_ZERO_
  INDSELF(LAY) = 0

 5400    CONTINUE

!        We have now isolated the layer ln pressure and temperature,
!        between two reference pressures and two reference temperatures 
!        (for each reference pressure).  We multiply the pressure 
!        fraction FP with the appropriate temperature fractions to get 
!        the factors that will be needed for the interpolation that yields
!        the optical depths (performed in routines TAUGBn for band n).

  COMPFP = 1. - FP
  FAC10(LAY) = COMPFP * FT
  FAC00(LAY) = COMPFP * (1. - FT)
  FAC11(LAY) = FP * FT1
  FAC01(LAY) = FP * (1. - FT1)

!  IF (NDBUG.LE.3) THEN
!    print 9000,LAY,LAYTROP,JP(LAY),JT(LAY),JT1(LAY),TAVEL(LAY) &
!      &,FAC00(LAY),FAC01(LAY),FAC10(LAY),FAC11(LAY) &
!      &,CO2MULT(LAY),COLMOL(LAY),COLCH4(LAY),COLCO2(LAY),COLH2O(LAY) &
!      &,COLN2O(LAY),COLO2(LAY),COLO3(LAY),SELFFAC(LAY),SELFFRAC(LAY) &
!      &,FORFAC(LAY),FORFRAC(LAY),INDSELF(LAY),INDFOR(LAY)
!9000 format(1x,2I3,3I4,F6.1,4F7.2,12E9.2,2I5)
!  END IF

END DO

!----------------------------------------------------------------------- 
RETURN
END SUBROUTINE RRTM_SW_SETCOEF

