C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

SUBROUTINE TAUMOL16 &
  &( KLEV    &
  &, FAC00   , FAC01   , FAC10   , FAC11   &
  &, JP      , JT      , JT1     , ONEMINUS&
  &, COLH2O  , COLCH4  , COLMOL  &
  &, LAYTROP , SELFFAC , SELFFRAC, INDSELF , FORFAC  , FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG    , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 16:  2600-3250 cm-1 (low - H2O,CH4; high - CH4)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

#include "tsmbkind.h"


USE PARSRTM  , ONLY : JPLAY, JPBAND, NG16
USE YOESRTA16, ONLY : ABSA, ABSB, FORREF, SELFREF &
  &                 , SFLUXREF, RAYL &
  &                 , LAYREFFR, STRRAT1
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE

!-- Output
REAL_B :: TAUG(JPLAY,NG16), TAUR(JPLAY,NG16), SSA(JPLAY,NG16), SFLUXZEN(NG16)

!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLCH4(JPLAY), COLCO2(JPLAY), COLMOL(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

DO LAY = 1, LAYTROP
  SPECCOMB = COLH2O(LAY) + STRRAT1*COLCH4(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 8._JPRB*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (_ONE_ - FS) * FAC00(LAY)
  FAC010 = (_ONE_ - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (_ONE_ - FS) * FAC01(LAY)
  FAC011 = (_ONE_ - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(16) + JS
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(16) + JS
  INDS = INDSELF(LAY)
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  print 9001,LAY,IND0,IND1,INDS,INDF,FAC000,FAC010,FAC100,FAC110,FAC001,FAC011,FAC101,FAC111 &
!   &,TAURAY,SELFFAC(LAY),SELFFRAC(LAY),FORFAC(LAY),FORFRAC(LAY)
9001 format(1x,'T16 ',5I4,13E12.3)


!  DO IG = 1, NG(16)
  DO IG = 1, NG16
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSA(IND0   ,IG) + &
     &           FAC100 * ABSA(IND0 +1,IG) + &
     &           FAC010 * ABSA(IND0 +9,IG) + &
     &           FAC110 * ABSA(IND0+10,IG) + &
     &           FAC001 * ABSA(IND1   ,IG) + &
     &           FAC101 * ABSA(IND1 +1,IG) + &
     &           FAC011 * ABSA(IND1 +9,IG) + &
     &           FAC111 * ABSA(IND1+10,IG)) + &
     &           COLH2O(LAY) * &
     &           (SELFFAC(LAY) * (SELFREF(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG)))) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
!    print 9002,LAY,IG,ABSA(IND0,IG),ABSA(IND0+1,IG),ABSA(IND0+9,IG),ABSA(IND0+10,IG) & 
!      &, ABSA(IND1,IG),ABSA(IND1+1,IG),ABSA(IND1+9,IG),ABSA(IND1+10,IG) & 
!      &, SELFREF(INDS+1,IG),SELFREF(INDS,IG),FORREF(INDF+1,IG),FORREF(INDF,IG)
9002 format(1x,'U16 ',2I3,12E12.3)
  END DO
END DO

LAYSOLFR = NLAYERS

DO LAY = LAYTROP+1, NLAYERS
  IF (JP(LAY-1) .LT. LAYREFFR .AND. JP(LAY) .GE. LAYREFFR) &
   &        LAYSOLFR = LAY
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(16) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(16) + 1
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(16)
  DO IG = 1, NG16
    TAUG(LAY,IG) = COLCH4(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0  ,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1  ,IG) + &
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG) 
    TAUR(LAY,IG) = TAURAY  
  END DO
END DO

!DO LAY=1,NLAYERS
!  print 9003,LAY,(TAUG(LAY,IG),IG=1,NG16)
9003 format(1x,'O16 ',I3,16E13.5)
!END DO

!----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL16

SUBROUTINE TAUMOL17 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCO2 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

#include "tsmbkind.h"


USE PARSRTM  , ONLY : JPLAY, JPBAND, NG17
USE YOESRTA17, ONLY : ABSA, ABSB, FORREF, SELFREF &
  &                 , SFLUXREF, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG17), TAUR(JPLAY,NG17), SSA(JPLAY,NG17), SFLUXZEN(NG17)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLCO2(JPLAY), COLMOL(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  
DO LAY = 1, LAYTROP
  SPECCOMB = COLH2O(LAY) + STRRAT*COLCO2(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 8.*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (1. - FS) * FAC00(LAY)
  FAC010 = (1. - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (1. - FS) * FAC01(LAY)
  FAC011 = (1. - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(17) + JS
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(17) + JS
  INDS = INDSELF(LAY)
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(17)
  DO IG = 1, NG17
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSA(IND0,IG) + &
     &           FAC100 * ABSA(IND0+1,IG) + &
     &           FAC010 * ABSA(IND0+9,IG) + &
     &           FAC110 * ABSA(IND0+10,IG) + &
     &           FAC001 * ABSA(IND1,IG) + &
     &           FAC101 * ABSA(IND1+1,IG) + &
     &           FAC011 * ABSA(IND1+9,IG) + &
     &           FAC111 * ABSA(IND1+10,IG)) + &
     &           COLH2O(LAY) * &
     &           (SELFFAC(LAY) * (SELFREF(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG)))) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

LAYSOLFR = NLAYERS

DO LAY = LAYTROP+1, NLAYERS
  IF (JP(LAY-1) .LT. LAYREFFR .AND. JP(LAY) .GE. LAYREFFR) &
    &        LAYSOLFR = LAY
  SPECCOMB = COLH2O(LAY) + STRRAT*COLCO2(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 4.*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (1. - FS) * FAC00(LAY)
  FAC010 = (1. - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (1. - FS) * FAC01(LAY)
  FAC011 = (1. - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(17) + JS
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(17) + JS
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(17)
  DO IG = 1, NG17
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSB(IND0,IG) + &
     &           FAC100 * ABSB(IND0+1,IG) + &
     &           FAC010 * ABSB(IND0+5,IG) + &
     &           FAC110 * ABSB(IND0+6,IG) + &
     &           FAC001 * ABSB(IND1,IG) + &
     &           FAC101 * ABSB(IND1+1,IG) + &
     &           FAC011 * ABSB(IND1+5,IG) + &
     &           FAC111 * ABSB(IND1+6,IG)) + &
     &           COLH2O(LAY) * &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG,JS) &
     &           + FS * (SFLUXREF(IG,JS+1) - SFLUXREF(IG,JS))
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL17

SUBROUTINE TAUMOL18 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCH4 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG18
USE YOESRTA18, ONLY : ABSA, ABSB, FORREF, SELFREF &
  &                 , SFLUXREF, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG18), TAUR(JPLAY,NG18), SSA(JPLAY,NG18), SFLUXZEN(NG18)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLCH4(JPLAY), COLMOL(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

      LAYSOLFR = LAYTROP
      
DO LAY = 1, LAYTROP
  IF (JP(LAY) .LT. LAYREFFR .AND. JP(LAY+1) .GE. LAYREFFR) &
   &        LAYSOLFR = MIN(LAY+1,LAYTROP)
  SPECCOMB = COLH2O(LAY) + STRRAT*COLCH4(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 8.*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (1. - FS) * FAC00(LAY)
  FAC010 = (1. - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (1. - FS) * FAC01(LAY)
  FAC011 = (1. - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(18) + JS
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(18) + JS
  INDS = INDSELF(LAY)
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(18)
  DO IG = 1, NG18
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSA(IND0,IG) + &
     &           FAC100 * ABSA(IND0+1,IG) + &
     &           FAC010 * ABSA(IND0+9,IG) + &
     &           FAC110 * ABSA(IND0+10,IG) + &
     &           FAC001 * ABSA(IND1,IG) + &
     &           FAC101 * ABSA(IND1+1,IG) + &
     &           FAC011 * ABSA(IND1+9,IG) + &
     &           FAC111 * ABSA(IND1+10,IG)) + &
     &           COLH2O(LAY) * &
     &           (SELFFAC(LAY) * (SELFREF(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG)))) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG,JS)  &
     &           + FS * (SFLUXREF(IG,JS+1) - SFLUXREF(IG,JS))
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(18) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(18) + 1
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(18)
  DO IG = 1, NG18
    TAUG(LAY,IG) = COLCH4(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + &	  
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL18

SUBROUTINE TAUMOL19 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCO2 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.


!     BAND 19:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG19
USE YOESRTA19, ONLY : ABSA, ABSB, FORREF, SELFREF &
  &                 , SFLUXREF, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG19), TAUR(JPLAY,NG19), SSA(JPLAY,NG19), SFLUXZEN(NG19)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLCO2(JPLAY), COLMOL(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

LAYSOLFR = LAYTROP
      
DO LAY = 1, LAYTROP
  IF (JP(LAY) .LT. LAYREFFR .AND. JP(LAY+1) .GE. LAYREFFR) &
   &        LAYSOLFR = MIN(LAY+1,LAYTROP)
  SPECCOMB = COLH2O(LAY) + STRRAT*COLCO2(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 8.*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (1. - FS) * FAC00(LAY)
  FAC010 = (1. - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (1. - FS) * FAC01(LAY)
  FAC011 = (1. - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(19) + JS
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(19) + JS
  INDS = INDSELF(LAY)
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(19)
  DO IG = 1 , NG19
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSA(IND0,IG) + &
     &           FAC100 * ABSA(IND0+1,IG) + &
     &           FAC010 * ABSA(IND0+9,IG) + &
     &           FAC110 * ABSA(IND0+10,IG) + &
     &           FAC001 * ABSA(IND1,IG) + &
     &           FAC101 * ABSA(IND1+1,IG) + &
     &           FAC011 * ABSA(IND1+9,IG) + &
     &           FAC111 * ABSA(IND1+10,IG)) + &
     &           COLH2O(LAY) * &
     &           (SELFFAC(LAY) * (SELFREF(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) + & 
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG)))) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG,JS) &
     &           + FS * (SFLUXREF(IG,JS+1) - SFLUXREF(IG,JS))
    TAUR(LAY,IG) = TAURAY   
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(19) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(19) + 1
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(19)
  DO IG = 1 , NG19
    TAUG(LAY,IG) = COLCO2(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + &
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG) 
    TAUR(LAY,IG) = TAURAY   
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL19

SUBROUTINE TAUMOL20 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCH4 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.


!     BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)


! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG20
USE YOESRTA20, ONLY : ABSA, ABSB, FORREF, SELFREF &
  &                 , SFLUXREF, ABSCH4, RAYL &
  &                 , LAYREFFR
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG20), TAUR(JPLAY,NG20), SSA(JPLAY,NG20), SFLUXZEN(NG20)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLCH4(JPLAY), COLMOL(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV



!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

LAYSOLFR = LAYTROP

DO LAY = 1, LAYTROP
  IF (JP(LAY) .LT. LAYREFFR .AND. JP(LAY+1) .GE. LAYREFFR) &
     &        LAYSOLFR = MIN(LAY+1,LAYTROP)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(20) + 1
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(20) + 1
  INDS = INDSELF(LAY)
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(20)
  DO IG = 1 , NG20
    TAUG(LAY,IG) = COLH2O(LAY) * &
     &          ((FAC00(LAY) * ABSA(IND0,IG) + &
     &           FAC10(LAY) * ABSA(IND0+1,IG) + &
     &           FAC01(LAY) * ABSA(IND1,IG) + &
     &           FAC11(LAY) * ABSA(IND1+1,IG)) + &
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + & 
     &           SELFFRAC(LAY) * &
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG)))) &
     &           + COLCH4(LAY) * ABSCH4(IG)
!     &           + TAURAY &
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY 
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG) 
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(20) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(20) + 1
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(20)
  DO IG = 1 , NG20
    TAUG(LAY,IG) = COLH2O(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + &
     &           FAC11(LAY) * ABSB(IND1+1,IG) + &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG)))) + &
     &           COLCH4(LAY) * ABSCH4(IG)
!     &           TAURAY + &
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY 
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL20

SUBROUTINE TAUMOL21 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCO2 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG21
USE YOESRTA21, ONLY : ABSA, ABSB, FORREF, SELFREF &
  &                 , SFLUXREF, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG21), TAUR(JPLAY,NG21), SSA(JPLAY,NG21), SFLUXZEN(NG21)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLCO2(JPLAY), COLMOL(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV



!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  
      LAYSOLFR = LAYTROP
      
DO LAY = 1, LAYTROP
  IF (JP(LAY) .LT. LAYREFFR .AND. JP(LAY+1) .GE. LAYREFFR) &
   &        LAYSOLFR = MIN(LAY+1,LAYTROP)
  SPECCOMB = COLH2O(LAY) + STRRAT*COLCO2(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 8.*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (1. - FS) * FAC00(LAY)
  FAC010 = (1. - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (1. - FS) * FAC01(LAY)
  FAC011 = (1. - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(21) + JS
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(21) + JS
  INDS = INDSELF(LAY)
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(21)
  DO IG = 1 , NG21
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSA(IND0,IG) + &
     &           FAC100 * ABSA(IND0+1,IG) + &
     &           FAC010 * ABSA(IND0+9,IG) + &
     &           FAC110 * ABSA(IND0+10,IG) + &
     &           FAC001 * ABSA(IND1,IG) + &
     &           FAC101 * ABSA(IND1+1,IG) + &
     &           FAC011 * ABSA(IND1+9,IG) + &
     &           FAC111 * ABSA(IND1+10,IG)) + &
     &           COLH2O(LAY) * &
     &           (SELFFAC(LAY) * (SELFREF(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))))
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG,JS) &
     &           + FS * (SFLUXREF(IG,JS+1) - SFLUXREF(IG,JS))
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  SPECCOMB = COLH2O(LAY) + STRRAT*COLCO2(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 4.*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (1. - FS) * FAC00(LAY)
  FAC010 = (1. - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (1. - FS) * FAC01(LAY)
  FAC011 = (1. - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(21) + JS
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(21) + JS
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(21)
  DO IG = 1 , NG21
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSB(IND0,IG) + &
     &           FAC100 * ABSB(IND0+1,IG) + &
     &           FAC010 * ABSB(IND0+5,IG) + &
     &           FAC110 * ABSB(IND0+6,IG) + &
     &           FAC001 * ABSB(IND1,IG) + &
     &           FAC101 * ABSB(IND1+1,IG) + &
     &           FAC011 * ABSB(IND1+5,IG) + &
     &           FAC111 * ABSB(IND1+6,IG)) + &
     &           COLH2O(LAY) * &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG)))
!     &           + TAURAY 
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL21

SUBROUTINE TAUMOL22 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLMOL , COLO2   &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 22:  7700-8050 cm-1 (low - H2O,O2; high - O2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG22
USE YOESRTA22, ONLY : ABSA, ABSB, FORREF, SELFREF &
  &                 , SFLUXREF, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG22), TAUR(JPLAY,NG22), SSA(JPLAY,NG22), SFLUXZEN(NG22)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLMOL(JPLAY), COLO2(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY, O2ADJ , O2CONT

NLAYERS = KLEV



!     The following factor is the ratio of total O2 band intensity (lines 
!     and Mate continuum) to O2 band intensity (line only).  It is needed
!     to adjust the optical depths since the k's include only lines.
      O2ADJ = 1.6_JPRB
      
!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

      LAYSOLFR = LAYTROP

DO LAY = 1, LAYTROP
  IF (JP(LAY) .LT. LAYREFFR .AND. JP(LAY+1) .GE. LAYREFFR) &
   &        LAYSOLFR = MIN(LAY+1,LAYTROP)
  O2CONT = 4.35e-4*colo2(lay)/(350.0*2.0)
  SPECCOMB = COLH2O(LAY) + O2ADJ*STRRAT*COLO2(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 8.*(SPECPARM)
!         ODADJ = SPECPARM + O2ADJ * (1. - SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (1. - FS) * FAC00(LAY)
  FAC010 = (1. - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (1. - FS) * FAC01(LAY)
  FAC011 = (1. - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(22) + JS
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(22) + JS
  INDS = INDSELF(LAY)
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(22)
  DO IG = 1 , NG22
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSA(IND0,IG) + &
     &           FAC100 * ABSA(IND0+1,IG) + &
     &           FAC010 * ABSA(IND0+9,IG) + &
     &           FAC110 * ABSA(IND0+10,IG) + &
     &           FAC001 * ABSA(IND1,IG) + &
     &           FAC101 * ABSA(IND1+1,IG) + &
     &           FAC011 * ABSA(IND1+9,IG) + &
     &           FAC111 * ABSA(IND1+10,IG)) + &
     &           COLH2O(LAY) * &
     &           (SELFFAC(LAY) * (SELFREF(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG)))) &
     &           + O2CONT
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG,JS) &
     &           + FS * (SFLUXREF(IG,JS+1) - SFLUXREF(IG,JS))
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  O2CONT = 4.35e-4*colo2(lay)/(350.0*2.0)
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(22) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(22) + 1
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(22)
  DO IG = 1 , NG22
    TAUG(LAY,IG) = COLO2(LAY) * O2ADJ * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + &
     &           FAC11(LAY) * ABSB(IND1+1,IG)) + &
     &           O2CONT
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL22

SUBROUTINE TAUMOL23 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLMOL &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG23
USE YOESRTA23, ONLY : ABSA, FORREF, SELFREF &
  &                 , SFLUXREF, RAYL &
  &                 , LAYREFFR, GIVFAC 
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG23), TAUR(JPLAY,NG23), SSA(JPLAY,NG23), SFLUXZEN(NG23)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLCO2(JPLAY), COLMOL(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV



!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

LAYSOLFR = LAYTROP

DO LAY = 1, LAYTROP
  IF (JP(LAY) .LT. LAYREFFR .AND. JP(LAY+1) .GE. LAYREFFR) &
   &        LAYSOLFR = MIN(LAY+1,LAYTROP)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(23) + 1
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(23) + 1
  INDS = INDSELF(LAY)
  INDF = INDFOR(LAY)

!  DO IG = 1, NG(23)
  DO IG = 1 , NG23
    TAURAY = COLMOL(LAY) * RAYL(IG)
    TAUG(LAY,IG) = COLH2O(LAY) * &
     &          (GIVFAC * (FAC00(LAY) * ABSA(IND0,IG) + &
     &           FAC10(LAY) * ABSA(IND0+1,IG) + &
     &           FAC01(LAY) * ABSA(IND1,IG) + &
     &           FAC11(LAY) * ABSA(IND1+1,IG)) + &
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG)))) 
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG) 
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
!  DO IG = 1, NG(23)
  DO IG = 1 , NG23
!    TAUG(LAY,IG) = COLMOL(LAY) * RAYL(IG)
!    SSA(LAY,IG) = 1.0
    TAUG(LAY,IG) = _ZERO_
    TAUR(LAY,IG) = COLMOL(LAY) * RAYL(IG) 
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL23

SUBROUTINE TAUMOL24 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLMOL , COLO2   , COLO3    &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 24:  12850-16000 cm-1 (low - H2O,O2; high - O2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG24
USE YOESRTA24, ONLY : ABSA, ABSB, FORREF, SELFREF &
  &                 , SFLUXREF, ABSO3A, ABSO3B, RAYLA, RAYLB &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG24), TAUR(JPLAY,NG24), SSA(JPLAY,NG24), SFLUXZEN(NG24)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLMOL(JPLAY), COLO2(JPLAY), COLO3(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

LAYSOLFR = LAYTROP

DO LAY = 1, LAYTROP
  IF (JP(LAY) .LT. LAYREFFR .AND. JP(LAY+1) .GE. LAYREFFR) &
     &        LAYSOLFR = MIN(LAY+1,LAYTROP)
  SPECCOMB = COLH2O(LAY) + STRRAT*COLO2(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 8.*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (1. - FS) * FAC00(LAY)
  FAC010 = (1. - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (1. - FS) * FAC01(LAY)
  FAC011 = (1. - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(24) + JS
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(24) + JS
  INDS = INDSELF(LAY)
  INDF = INDFOR(LAY)

!  DO IG = 1, NG(24)
  DO IG = 1 , NG24
    TAURAY = COLMOL(LAY) * (RAYLA(IG,JS) + &
     &           FS * (RAYLA(IG,JS+1) - RAYLA(IG,JS)))
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSA(IND0,IG) + &
     &           FAC100 * ABSA(IND0+1,IG) + &
     &           FAC010 * ABSA(IND0+9,IG) + &
     &           FAC110 * ABSA(IND0+10,IG) + &
     &           FAC001 * ABSA(IND1,IG) + &
     &           FAC101 * ABSA(IND1+1,IG) + &
     &           FAC011 * ABSA(IND1+9,IG) + &
     &           FAC111 * ABSA(IND1+10,IG)) + &
     &           COLO3(LAY) * ABSO3A(IG) + &
     &           COLH2O(LAY) * & 
     &           (SELFFAC(LAY) * (SELFREF(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + & 
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))))
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG,JS) &
     &           + FS * (SFLUXREF(IG,JS+1) - SFLUXREF(IG,JS))
    TAUR(LAY,IG) =  TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(24) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(24) + 1

!  DO IG = 1, NG(24)
  DO IG = 1 , NG24
    TAURAY = COLMOL(LAY) * RAYLB(IG)
    TAUG(LAY,IG) = COLO2(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + &
     &           FAC11(LAY) * ABSB(IND1+1,IG)) + &
     &           COLO3(LAY) * ABSO3B(IG)
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL24

SUBROUTINE TAUMOL25 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLMOL , COLO3   &
  &, LAYTROP &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 25:  16000-22650 cm-1 (low - H2O; high - nothing)

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG25
USE YOESRTA25, ONLY : ABSA &
  &                 , SFLUXREF, ABSO3A, ABSO3B, RAYL &
  &                 , LAYREFFR
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG25), TAUR(JPLAY,NG25), SSA(JPLAY,NG25), SFLUXZEN(NG25)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLMOL(JPLAY), COLO3(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV



!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

LAYSOLFR = LAYTROP

DO LAY = 1, LAYTROP
  IF (JP(LAY) .LT. LAYREFFR .AND. JP(LAY+1) .GE. LAYREFFR) &
    &        LAYSOLFR = MIN(LAY+1,LAYTROP)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(25) + 1
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(25) + 1

!  DO IG = 1, NG(25)
  DO IG = 1 , NG25
    TAURAY = COLMOL(LAY) * RAYL(IG)
    TAUG(LAY,IG) = COLH2O(LAY) * &
     &          (FAC00(LAY) * ABSA(IND0,IG) + &
     &           FAC10(LAY) * ABSA(IND0+1,IG) + &
     &           FAC01(LAY) * ABSA(IND1,IG) + &
     &           FAC11(LAY) * ABSA(IND1+1,IG)) + &
     &           COLO3(LAY) * ABSO3A(IG) 
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG) 
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
!  DO IG = 1, NG(25)
  DO IG = 1 , NG25
    TAURAY = COLMOL(LAY) * RAYL(IG)
    TAUG(LAY,IG) = COLO3(LAY) * ABSO3B(IG) 
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL25

SUBROUTINE TAUMOL26 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCO2 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG26
USE YOESRTA26, ONLY : SFLUXREF, RAYL
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG26), TAUR(JPLAY,NG26), SSA(JPLAY,NG26), SFLUXZEN(NG26)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from AER
REAL_B :: TAUAERL(JPLAY,JPBAND)

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLCO2(JPLAY), COLMOL(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV



!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  
LAYSOLFR = LAYTROP

DO LAY = 1, LAYTROP
!  DO IG = 1, NG(26)
  DO IG = 1 , NG26 
!    TAUG(LAY,IG) = COLMOL(LAY) * RAYL(IG)
!    SSA(LAY,IG) = 1.0
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG) 
    TAUG(LAY,IG) = _ZERO_
    TAUR(LAY,IG) = COLMOL(LAY) * RAYL(IG) 
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
!  DO IG = 1, NG(26)
  DO IG = 1 , NG26
!    TAUG(LAY,IG) = COLMOL(LAY) * RAYL(IG)
!    SSA(LAY,IG) = 1.0
    TAUG(LAY,IG) = _ZERO_
    TAUR(LAY,IG) = COLMOL(LAY) * RAYL(IG) 
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL26

SUBROUTINE TAUMOL27 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLMOL  , COLO3  &
  &, LAYTROP &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 27:  29000-38000 cm-1 (low - O3; high - O3)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG27
USE YOESRTA27, ONLY : ABSA, ABSB &
  &                 , SFLUXREF, RAYL &
  &                 , LAYREFFR, SCALEKUR
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG27), TAUR(JPLAY,NG27), SSA(JPLAY,NG27), SFLUXZEN(NG27)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLMOL(JPLAY), COLO3(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV



!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  
DO LAY = 1, LAYTROP
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(27) + 1
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(27) + 1

!  DO IG = 1, NG(27)
  DO IG = 1 , NG27
    TAURAY = COLMOL(LAY) * RAYL(IG)
    TAUG(LAY,IG) = COLO3(LAY) * &
     &          (FAC00(LAY) * ABSA(IND0,IG) + &
     &           FAC10(LAY) * ABSA(IND0+1,IG) + &
     &           FAC01(LAY) * ABSA(IND1,IG) + &
     &           FAC11(LAY) * ABSA(IND1+1,IG))
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

LAYSOLFR = NLAYERS

DO LAY = LAYTROP+1, NLAYERS
  IF (JP(LAY-1) .LT. LAYREFFR .AND. JP(LAY) .GE. LAYREFFR) &
     &        LAYSOLFR = LAY
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(27) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(27) + 1

!  DO IG = 1, NG(27)
  DO IG = 1 , NG27
    TAURAY = COLMOL(LAY) * RAYL(IG)
    TAUG(LAY,IG) = COLO3(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + & 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY.EQ.LAYSOLFR) SFLUXZEN(IG) = SCALEKUR * SFLUXREF(IG) 
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL27

SUBROUTINE TAUMOL28 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLMOL  , COLO2  , COLO3   &
  &, LAYTROP &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 28:  38000-50000 cm-1 (low - O3,O2; high - O3,O2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG28
USE YOESRTA28, ONLY : ABSA, ABSB &
  &                 , SFLUXREF, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,NG28), TAUR(JPLAY,NG28), SSA(JPLAY,NG28), SFLUXZEN(NG28)


!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLMOL(JPLAY), COLO2(JPLAY), COLO3(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY

NLAYERS = KLEV

      
!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

DO LAY = 1, LAYTROP
  SPECCOMB = COLO3(LAY) + STRRAT*COLO2(LAY)
  SPECPARM = COLO3(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 8.*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (1. - FS) * FAC00(LAY)
  FAC010 = (1. - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (1. - FS) * FAC01(LAY)
  FAC011 = (1. - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(28) + JS
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(28) + JS
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(28)
  DO IG = 1 , NG28
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSA(IND0,IG) + &
     &           FAC100 * ABSA(IND0+1,IG) + &
     &           FAC010 * ABSA(IND0+9,IG) + &
     &           FAC110 * ABSA(IND0+10,IG) + &
     &           FAC001 * ABSA(IND1,IG) + &
     &           FAC101 * ABSA(IND1+1,IG) + &
     &           FAC011 * ABSA(IND1+9,IG) + &
     &           FAC111 * ABSA(IND1+10,IG)) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

LAYSOLFR = NLAYERS

DO LAY = LAYTROP+1, NLAYERS
  IF (JP(LAY-1) .LT. LAYREFFR .AND. JP(LAY) .GE. LAYREFFR) &
   &        LAYSOLFR = LAY
  SPECCOMB = COLO3(LAY) + STRRAT*COLO2(LAY)
  SPECPARM = COLO3(LAY)/SPECCOMB 
  IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
  SPECMULT = 4.*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT, _ONE_ )
  FAC000 = (1. - FS) * FAC00(LAY)
  FAC010 = (1. - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (1. - FS) * FAC01(LAY)
  FAC011 = (1. - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(28) + JS
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(28) + JS
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(28)
  DO IG = 1 , NG28
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSB(IND0,IG) + &
     &           FAC100 * ABSB(IND0+1,IG) + &
     &           FAC010 * ABSB(IND0+5,IG) + &
     &           FAC110 * ABSB(IND0+6,IG) + &
     &           FAC001 * ABSB(IND1,IG) + &
     &           FAC101 * ABSB(IND1+1,IG) + &
     &           FAC011 * ABSB(IND1+5,IG) + &
     &           FAC111 * ABSB(IND1+6,IG)) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG,JS) &
     &           + FS * (SFLUXREF(IG,JS+1) - SFLUXREF(IG,JS))
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL28

SUBROUTINE TAUMOL29 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    & 
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCO2 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)

! Modifications
!
!     JJMorcrette 2002-10-03 adapted to ECMWF environment

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, NG29
USE YOESRTA29, ONLY : ABSA, ABSB, FORREF, SELFREF &
  &                 , SFLUXREF, ABSH2O, ABSCO2, RAYL &
  &                 , LAYREFFR
USE YOESRTWN , ONLY : NG, NSPA, NSPB

  
IMPLICIT NONE

!-- Output
REAL_B :: TAUG(JPLAY,NG29), TAUR(JPLAY,NG29), SSA(JPLAY,NG29), SFLUXZEN(NG29)

!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY), FAC01(JPLAY), FAC10(JPLAY), FAC11(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY), COLCO2(JPLAY), COLMOL(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: SELFFAC(JPLAY), SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)

!-- from FOREIGN
REAL_B :: FORFAC(JPLAY), FORFRAC(JPLAY)
INTEGER_M :: INDFOR(JPLAY)



!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, INDF, JS, LAY, LAYSOLFR, NLAYERS

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
        & FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM, &
        & TAURAY



NLAYERS = KLEV

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  

DO LAY = 1, LAYTROP
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(29) + 1
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(29) + 1
  INDS = INDSELF(LAY)
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(29)
  DO IG = 1, NG29
    TAUG(LAY,IG) = COLH2O(LAY) * &
     &          ((FAC00(LAY) * ABSA(IND0,IG) + &
     &           FAC10(LAY) * ABSA(IND0+1,IG) + &
     &           FAC01(LAY) * ABSA(IND1,IG) + &
     &           FAC11(LAY) * ABSA(IND1+1,IG)) + &
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREF(INDF,IG) + & 
     &           FORFRAC(LAY) * &
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG)))) &
     &           + COLCO2(LAY) * ABSCO2(IG) 
!     &           + TAURAY &
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

LAYSOLFR = NLAYERS
DO LAY = LAYTROP+1, NLAYERS
  IF (JP(LAY-1) .LT. LAYREFFR .AND. JP(LAY) .GE. LAYREFFR) &
     &        LAYSOLFR = LAY
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(29) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(29) + 1
  TAURAY = COLMOL(LAY) * RAYL

!  DO IG = 1, NG(29)
  DO IG = 1 , NG29
    TAUG(LAY,IG) = COLCO2(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + &
     &           FAC11(LAY) * ABSB(IND1+1,IG)) &  
     &           + COLH2O(LAY) * ABSH2O(IG) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREF(IG) 
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE TAUMOL29

