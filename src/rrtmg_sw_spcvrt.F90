!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

SUBROUTINE RRTMG_SW_SPCVRT &
 &( KLEV   , KMOL    , KSW    , ONEMINUS,ISTART  , IEND &
 &, PAVEL  , TAVEL   , PZ     , TZ     , TBOUND  , PALBD   , PALBP &
 &, PCLFR  , PTAUC   , PASYC  , POMGC  , PTAUCORIG &
 &, PTAUA  , PASYA   , POMGA  , PRMU0   &
 &, COLDRY , WKL     , ADJFLUX  &
 &, LAYTROP, LAYSWTCH, LAYLOW &
 &, CO2MULT, COLCH4  , COLCO2 , COLH2O , COLMOL  , COLN2O  , COLO2 , COLO3 &
 &, FORFAC , FORFRAC , INDFOR , SELFFAC, SELFFRAC, INDSELF &
 &, FAC00  , FAC01   , FAC10  , FAC11  &
 &, JP     , JT      , JT1	 &
!-- output arrays 
! &, PBBFD, PBBFU, PUVFD, PUVFU, PVSFD, PVSFU , PNIFD , PNIFU &
! &, PBBCD, PBBCU, PUVCD, PUVCU, PVSCD, PVSCU , PNICD , PNICU &
 &, PBBFD, PBBFU &
 &, PBBCD, PBBCU &
 &, PBBFDdir, PBBCDdir &
 &)


!**** *RRTMG_SW_SPCVRT* - SPECTRAL LOOP TO COMPUTE THE SHORTWAVE RADIATION FLUXES.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE TWO-STREAM METHOD OF BARKER

!**   INTERFACE.
!     ----------

!          *RRTMG_SW_SPCVRT* IS CALLED FROM *SRTM_SRTM_224GP*


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------

!          *SWVRTQDR*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION
!     AUTHOR.
!     -------
!        from Howard Barker
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 03-02-27
!        Add adjustment for Earth/Sun distance : MJIacono, AER, October 2003
!        Bug fix for use of PALBP and PALBD: MJIacono, AER, November 2003
!        Bug fix to apply delta scaling to clear sky: AER, December 2004
!        Code modified so that delta scaling not done in cloudy profiles
!        if routine cldprop is used; delta scaling can be applied by
!        switching code below if cldprop is not used to get cloud 
!        properties:                                  AER, January 2005
!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPB1, JPB2, JPGPT
USE YOESRTWN , ONLY : NGC

IMPLICIT NONE

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER_M :: KAER   , KLEV    , KMOL  , KSW, ISTART, IEND
INTEGER_M :: LAYTROP, LAYSWTCH, LAYLOW

REAL_B :: ONEMINUS
REAL_B :: PALBD(KSW)     , PALBP(KSW)     , PRMU0 
REAL_B :: PCLFR(JPLAY)	 , PTAUC(JPLAY,KSW),PASYC(JPLAY,KSW), POMGC(JPLAY,KSW)
REAL_B :: PTAUCORIG(JPLAY,KSW)
REAL_B :: PTAUA(JPLAY,KSW),PASYA(JPLAY,KSW),POMGA(JPLAY,KSW)
REAL_B :: PAVEL(JPLAY)   , TAVEL(JPLAY)   , PZ(0:JPLAY)     , TZ(0:JPLAY)  , TBOUND
REAL_B :: COLDRY(JPLAY)  , COLMOL(JPLAY)  , WKL(35,JPLAY)
REAL_B :: CO2MULT(JPLAY) , COLCH4(JPLAY)  , COLCO2(JPLAY)   , COLH2O(JPLAY)
REAL_B :: COLN2O(JPLAY)  , COLO2(JPLAY)   , COLO3(JPLAY)
REAL_B :: FORFAC(JPLAY)  , FORFRAC(JPLAY) , SELFFAC(JPLAY)  , SELFFRAC(JPLAY)
REAL_B :: FAC00(JPLAY)   ,   FAC01(JPLAY) , FAC10(JPLAY)    , FAC11(JPLAY)
REAL_B :: ADJFLUX(JPBAND)

INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)
INTEGER_M :: JP(JPLAY)    , JT(JPLAY)     , JT1(JPLAY)

REAL_B :: &
  &   PBBCD(JPLAY+1)          , PBBCU(JPLAY+1) &
  &,  PUVCD(JPLAY+1)          , PUVCU(JPLAY+1) &
  &,  PVSCD(JPLAY+1)          , PVSCU(JPLAY+1) &
  &,  PNICD(JPLAY+1)          , PNICU(JPLAY+1) &
  &,  PBBFD(JPLAY+1)          , PBBFU(JPLAY+1) &
  &,  PBBFDdir(JPLAY+1)       , PBBCDdir(JPLAY+1) &
  &,  PUVFD(JPLAY+1)          , PUVFU(JPLAY+1) &
  &,  PVSFD(JPLAY+1)          , PVSFU(JPLAY+1) &
  &,  PNIFD(JPLAY+1)          , PNIFU(JPLAY+1)

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

LOGICAL :: LRTCHKCLR(JPLAY),LRTCHKCLD(JPLAY)

REAL_B :: &
  &   ZCLEAR      , ZCLOUD       &
  &,  ZDBT(JPLAY+1), ZDBT_nodel(JPLAY+1) &
  &,  ZGC(JPLAY)    , ZGCC(JPLAY)    , ZGCO(JPLAY)     &
  &,  ZOMC(JPLAY)   , ZOMCC(JPLAY)   , ZOMCO(JPLAY)    &
  &,  ZRDND(JPLAY+1), ZRDNDC(JPLAY+1)&
  &,  ZREF(JPLAY+1) , ZREFC(JPLAY+1) , ZREFO(JPLAY+1)  &
  &,  ZREFD(JPLAY+1), ZREFDC(JPLAY+1), ZREFDO(JPLAY+1) &
  &,  ZRUP(JPLAY+1) , ZRUPD(JPLAY+1) &
  &,  ZRUPC(JPLAY+1), ZRUPDC(JPLAY+1)&
  &,  ZS1(JPLAY+1)  &
  &,  ZTAUC(JPLAY)  , ZTAUO(JPLAY)    &
  &,  ZTDN(JPLAY+1) , ZTDND(JPLAY+1) , ZTDBT(JPLAY+1)  &
  &,  ZTOC(JPLAY)   , ZTOR(JPLAY)    & 
  &,  ZTRA(JPLAY+1) , ZTRAC(JPLAY+1) , ZTRAO(JPLAY+1)  &
  &,  ZTRAD(JPLAY+1), ZTRADC(JPLAY+1), ZTRADO(JPLAY+1) 
REAL_B :: &
  &   ZDBTC(JPLAY+1), ZTDBTC(JPLAY+1), ZINCFLX(JPGPT)  &
  &,  ZDBTC_nodel(JPLAY+1), ZTDBT_nodel(JPLAY+1), ZTDBTC_nodel(JPLAY+1)
  

!     LOCAL INTEGER SCALARS
INTEGER_M :: IB1, IB2, IBM, IGT, IKL, IKP, IKX, IW, JB, JG, JL, JK, KGS, KMODTS

!     LOCAL REAL SCALARS
REAL_B :: ZDBTMC, ZDBTMO, ZF, ZGW, ZINCFLUX, ZREFLECT, ZWF, TAUORIG
REAL_B :: REPCLC

!-- Output of RRTMG_SW_TAUMOLn routines

REAL_B :: ZTAUG(JPLAY,16), ZTAUR(JPLAY,16), ZSSA(JPLAY,16), ZSFLXZEN(16)

!-- Output of RRTMG_SW_VRTQDR routine
REAL_B :: &
  &   ZCD(JPLAY+1,JPGPT), ZCU(JPLAY+1,JPGPT) &
  &,  ZFD(JPLAY+1,JPGPT), ZFU(JPLAY+1,JPGPT)
REAL_B :: &
  &   ZBBCD(JPLAY+1)          , ZBBCU(JPLAY+1) &
  &,  ZBBFD(JPLAY+1)          , ZBBFU(JPLAY+1) &
  &,  ZBBFDdir(JPLAY+1)       , ZBBCDdir(JPLAY+1)
!     ------------------------------------------------------------------

!-- Two-stream model 1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
! KMODTS is set in SWREFTRA

IB1=ISTART
IB2=IEND

REPCLC=1.E-12_JPRB

IW=0
ZINCFLUX=_ZERO_

DO JK=1,KLEV+1
  PBBCD(JK)=_ZERO_
  PBBCU(JK)=_ZERO_
  PBBFD(JK)=_ZERO_
  PBBFU(JK)=_ZERO_
END DO

JB=IB1-1
DO JB = IB1, IB2
  IBM = JB-15
  IGT = NGC(IBM)

  IF (JB .EQ. 16) THEN
    CALL TAUMOL16 &
      &( KLEV    &
      &, FAC00   , FAC01   , FAC10   , FAC11   &
      &, JP      , JT      , JT1     , ONEMINUS &
      &, COLH2O  , COLCH4  , COLMOL  &
      &, LAYTROP , SELFFAC , SELFFRAC, INDSELF, FORFAC  , FORFRAC, INDFOR &
      &, ZSFLXZEN, ZTAUG    , ZTAUR    &
      &)

  ELSE IF (JB .EQ. 17) THEN
    CALL TAUMOL17 &
      &( KLEV    &
      &, FAC00   , FAC01  , FAC10   , FAC11    &
      &, JP      , JT     , JT1     , ONEMINUS &
      &, COLH2O  , COLCO2 , COLMOL  &
      &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
      &, ZSFLXZEN, ZTAUG   , ZTAUR    &
      &)

  ELSE IF (JB .EQ. 18) THEN
    CALL TAUMOL18 &
      &( KLEV    &
      &, FAC00   , FAC01  , FAC10   , FAC11    &
      &, JP      , JT     , JT1     , ONEMINUS &
      &, COLH2O  , COLCH4 , COLMOL  &
      &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
      &, ZSFLXZEN, ZTAUG   , ZTAUR    &
      &)

  ELSE IF (JB .EQ. 19) THEN
    CALL TAUMOL19 &
      &( KLEV    &
      &, FAC00   , FAC01  , FAC10   , FAC11    &
      &, JP      , JT     , JT1     , ONEMINUS &
      &, COLH2O  , COLCO2 , COLMOL  &
      &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
      &, ZSFLXZEN, ZTAUG   , ZTAUR    &
      &)

  ELSE IF (JB .EQ. 20) THEN
    CALL TAUMOL20 &
     &( KLEV    &
     &, FAC00   , FAC01  , FAC10   , FAC11    &
     &, JP      , JT     , JT1     , ONEMINUS &
     &, COLH2O  , COLCH4 , COLMOL  &
     &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
     &, ZSFLXZEN, ZTAUG   , ZTAUR    &
     &)

  ELSE IF (JB .EQ. 21) THEN
    CALL TAUMOL21 &
     &( KLEV    &
     &, FAC00   , FAC01  , FAC10   , FAC11    &
     &, JP      , JT     , JT1     , ONEMINUS &
     &, COLH2O  , COLCO2 , COLMOL  &
     &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
     &, ZSFLXZEN, ZTAUG   , ZTAUR    &
     &)

  ELSE IF (JB .EQ. 22) THEN
    CALL TAUMOL22 &
     &( KLEV    &
     &, FAC00   , FAC01  , FAC10   , FAC11    &
     &, JP      , JT     , JT1     , ONEMINUS &
     &, COLH2O  , COLMOL , COLO2   &
     &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
     &, ZSFLXZEN, ZTAUG   , ZTAUR    &
     &)

  ELSE IF (JB .EQ. 23) THEN
    CALL TAUMOL23 &
     &( KLEV    &
     &, FAC00   , FAC01  , FAC10   , FAC11    &
     &, JP      , JT     , JT1     , ONEMINUS &
     &, COLH2O  , COLMOL &
     &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
     &, ZSFLXZEN, ZTAUG   , ZTAUR    &
     &)

  ELSE IF (JB .EQ. 24) THEN
    CALL TAUMOL24 &
     &( KLEV    &
     &, FAC00   , FAC01  , FAC10   , FAC11    &
     &, JP      , JT     , JT1     , ONEMINUS &
     &, COLH2O  , COLMOL , COLO2   , COLO3    &
     &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
     &, ZSFLXZEN, ZTAUG   , ZTAUR    &
     &)

  ELSE IF (JB .EQ. 25) THEN
!--- visible 16000-22650 cm-1   0.4415 - 0.6250 um
    CALL TAUMOL25 &
     &( KLEV    &
     &, FAC00   , FAC01  , FAC10   , FAC11    &
     &, JP      , JT     , JT1     , ONEMINUS &
     &, COLH2O  , COLMOL , COLO3   &
     &, LAYTROP &
     &, ZSFLXZEN, ZTAUG   , ZTAUR    &
     &)

  ELSE IF (JB .EQ. 26) THEN
!--- UV-A 22650-29000 cm-1   0.3448 - 0.4415 um
    CALL TAUMOL26 &
     &( KLEV    &
     &, FAC00   , FAC01  , FAC10   , FAC11    &
     &, JP      , JT     , JT1     , ONEMINUS &
     &, COLH2O  , COLCO2 , COLMOL  &
     &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
     &, ZSFLXZEN, ZTAUG   , ZTAUR    &
     &)

  ELSE IF (JB .EQ. 27) THEN
!--- UV-B 29000-38000 cm-1   0.2632 - 0.3448 um
    CALL TAUMOL27 &
     &( KLEV    &
     &, FAC00   , FAC01  , FAC10   , FAC11    &
     &, JP      , JT     , JT1     , ONEMINUS &
     &, COLMOL  , COLO3  &
     &, LAYTROP &
     &, ZSFLXZEN, ZTAUG   , ZTAUR    &
     &)

  ELSE IF (JB .EQ. 28) THEN
!--- UV-C 38000-50000 cm-1   0.2000 - 0.2632 um
    CALL TAUMOL28 &
     &( KLEV    &
     &, FAC00   , FAC01  , FAC10   , FAC11    &
     &, JP      , JT     , JT1     , ONEMINUS &
     &, COLMOL  , COLO2  , COLO3   &
     &, LAYTROP &
     &, ZSFLXZEN, ZTAUG   , ZTAUR    &
     &)

  ELSE IF (JB .EQ. 29) THEN
    CALL TAUMOL29 &
     &( KLEV    &
     &, FAC00   , FAC01  , FAC10   , FAC11    & 
     &, JP      , JT     , JT1     , ONEMINUS &
     &, COLH2O  , COLCO2 , COLMOL  &
     &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
     &, ZSFLXZEN, ZTAUG   , ZTAUR    &
     &)

  ENDIF

  DO JK=1,KLEV+1
    ZBBCD(JK)=_ZERO_
    ZBBCU(JK)=_ZERO_
    ZBBFD(JK)=_ZERO_
    ZBBFU(JK)=_ZERO_
  END DO

  DO JG=1,IGT
    IW=IW+1

! mji - add adjustment for correct Earth/Sun distance
!    ZINCFLX(IW)=ZSFLXZEN(JG)*PRMU0
!    ZINCFLUX   =ZINCFLUX+ZSFLXZEN(JG)*PRMU0           
    ZINCFLX(IW)=ADJFLUX(JB)*ZSFLXZEN(JG)*PRMU0
    ZINCFLUX   =ZINCFLUX+ADJFLUX(JB)*ZSFLXZEN(JG)*PRMU0           

!-- CALL to compute layer reflectances and transmittances for direct 
!  and diffuse sources, first clear then cloudy

! ZREFC(JK)  direct albedo for clear
! ZREFO(JK)  direct albedo for cloud
! ZREFDC(JK) diffuse albedo for clear
! ZREFDO(JK) diffuse albedo for cloud
! ZTRAC(JK)  direct transmittance for clear
! ZTRAO(JK)  direct transmittance for cloudy
! ZTRADC(JK) diffuse transmittance for clear
! ZTRADO(JK) diffuse transmittance for cloudy
!  
! ZREF(JK)   direct reflectance
! ZREFD(JK)  diffuse reflectance
! ZTRA(JK)   direct transmittance
! ZTRAD(JK)  diffuse transmittance
!
! ZDBTC(JK)  clear direct beam transmittance
! ZDBTO(JK)  cloudy direct beam transmittance
! ZDBT(JK)   layer mean direct beam transmittance
! ZTDBT(JK)  total direct beam transmittance at levels

!-- clear-sky    
!----- TOA direct beam    
    ZTDBTC(1)=1._JPRB
    ZTDBTC_nodel(1)=1._JPRB
!----- surface values
    ZDBTC(KLEV+1) =_ZERO_
    ZTRAC(KLEV+1) =_ZERO_
    ZTRADC(KLEV+1)=_ZERO_
    ZREFC(KLEV+1) =PALBP(IBM)
    ZREFDC(KLEV+1)=PALBD(IBM)
    ZRUPC(KLEV+1) =PALBP(IBM)
    ZRUPDC(KLEV+1)=PALBD(IBM)
           
!-- total sky    
!----- TOA direct beam    
    ZTDBT(1)=1._JPRB
    ZTDBT_nodel(1)=1._JPRB
!----- surface values
    ZDBT(KLEV+1) =_ZERO_
    ZTRA(KLEV+1) =_ZERO_
    ZTRAD(KLEV+1)=_ZERO_
    ZREF(KLEV+1) =PALBP(IBM)
    ZREFD(KLEV+1)=PALBD(IBM)
    ZRUP(KLEV+1) =PALBP(IBM)
    ZRUPD(KLEV+1)=PALBD(IBM)
    
    
    DO JK=1,KLEV

!-- NB: a two-stream calculations from top to bottom, but RRTMG_SW quantities 
!       are given bottom to top and are reverse here:

      IKL=KLEV+1-JK

!-- Set logical flag to do REFTRA calculation
!   Do REFTRA for all clear layers
       LRTCHKCLR(JK)=.TRUE.
!   Do REFTRA only for cloudy layers in profile, since already done for clear layers
       LRTCHKCLD(JK)=.FALSE.
       LRTCHKCLD(JK)=(PCLFR(IKL) > REPCLC)

!-- clear-sky optical parameters      
!-- original
!      ZTAUC(JK)=ZTAUR(IKL,JG)+ZTAUG(IKL,JG)
!      ZOMCC(JK)=ZTAUR(IKL,JG)/ZTAUC(JK)
!      ZGCC (JK)=0.0001_JPRB
!-- total sky optical parameters        
!      ZTAUO(JK)=ZTAUR(IKL,JG)+ZTAUG(IKL,JG)+PTAUC(IKL,IBM)
!      ZOMCO(JK)=PTAUC(IKL,IBM)*POMGC(IKL,IBM)+ZTAUR(IKL,JG)
!      ZGCO (JK)=(PTAUC(IKL,IBM)*POMGC(IKL,IBM)*PASYC(IKL,IBM) &
!        & +ZTAUR(IKL,JG)*0.0001_JPRB)/ZOMCO(JK)
!      ZOMCO(JK)=ZOMCO(JK)/ZTAUO(JK)

!-- clear-sky optical parameters including aerosols
      ZTAUC(JK) = ZTAUR(IKL,JG) + ZTAUG(IKL,JG) + PTAUA(IKL,IBM)
      ZOMCC(JK) = ZTAUR(IKL,JG)*_ONE_ + PTAUA(IKL,IBM)*POMGA(IKL,IBM)
      ZGCC (JK) = PASYA(IKL,IBM)*POMGA(IKL,IBM)*PTAUA(IKL,IBM) / ZOMCC(JK)
      ZOMCC(JK) = ZOMCC(JK) / ZTAUC(JK)

!-- Pre-delta-scaling clear and cloudy direct beam transmittance (must use 'orig', unscaled cloud OD)       
!   This block of code is only needed for direct beam calculation
!     
      ZCLEAR     = _ONE_ - PCLFR(IKL)
      ZCLOUD     = PCLFR(IKL)
      ZDBTMC     = EXP(-ZTAUC(JK)/PRMU0)
      TAUORIG    = ZTAUC(JK) + PTAUCORIG(IKL,IBM)
      ZDBTMO     = EXP(-TAUORIG/PRMU0)
      ZDBT_nodel(JK)   = ZCLEAR*ZDBTMC+ZCLOUD*ZDBTMO
      ZTDBT_nodel(JK+1)= ZDBT_nodel(JK)*ZTDBT_nodel(JK)
!-- clear-sky        
      ZDBTC_nodel(JK)=ZDBTMC
      ZTDBTC_nodel(JK+1)=ZDBTC_nodel(JK)*ZTDBTC_nodel(JK)
!  Only needed for direct beam calculation ^^^

!-- Delta scaling - clear   
      ZF=ZGCC(JK)*ZGCC(JK)
      ZWF=ZOMCC(JK)*ZF
      ZTAUC(JK)=(1._JPRB-ZWF)*ZTAUC(JK)
      ZOMCC(JK)=(ZOMCC(JK)-ZWF)/(1._JPRB-ZWF)
      ZGCC (JK)=(ZGCC(JK)-ZF)/(1._JPRB-ZF)

!-- total sky optical parameters (cloud properties already delta-scaled)
!   Use this code if cloud properties are derived in rrtmg_sw_cldprop       
      ZTAUO(JK) = ZTAUC(JK) + PTAUC(IKL,IBM)
      ZOMCO(JK) = ZTAUC(JK)*ZOMCC(JK) + PTAUC(IKL,IBM)*POMGC(IKL,IBM) 
      ZGCO (JK) = (PTAUC(IKL,IBM)*POMGC(IKL,IBM)*PASYC(IKL,IBM)  &
        &       +  ZTAUC(JK)*ZOMCC(JK)*PASYA(IKL,IBM)) /  ZOMCO(JK)
      ZOMCO(JK) = ZOMCO(JK) / ZTAUO(JK)

!-- total sky optical parameters (if cloud properties not delta scaled)
!   Use this code if cloud properties are not derived in rrtmg_sw_cldprop       
!      ZTAUO(JK) = ZTAUR(IKL,JG) + ZTAUG(IKL,JG) + PTAUA(IKL,IBM) + PTAUC(IKL,IBM)
!      ZOMCO(JK) = PTAUA(IKL,IBM)*POMGA(IKL,IBM) + PTAUC(IKL,IBM)*POMGC(IKL,IBM) &
!        &       + ZTAUR(IKL,JG)*_ONE_
!      ZGCO (JK) = (PTAUC(IKL,IBM)*POMGC(IKL,IBM)*PASYC(IKL,IBM)  &
!        &       +  PTAUA(IKL,IBM)*POMGA(IKL,IBM)*PASYA(IKL,IBM)) &
!        &       /  ZOMCO(JK)
!      ZOMCO(JK) = ZOMCO(JK) / ZTAUO(JK)
!
!-- Delta scaling - clouds; Use only if subroutine rrtmg_sw_cldprop 
!   is not used to get cloud properties and to apply delta scaling
!      ZF=ZGCO(JK)*ZGCO(JK)
!      ZWF=ZOMCO(JK)*ZF
!      ZTAUO(JK)=(1._JPRB-ZWF)*ZTAUO(JK)
!      ZOMCO(JK)=(ZOMCO(JK)-ZWF)/(1._JPRB-ZWF)
!      ZGCO (JK)=(ZGCO(JK)-ZF)/(1._JPRB-ZF)

    END DO    

! Clear sky reflectivities
    CALL RRTMG_SW_REFTRA ( KLEV, KMODTS &
      &, LRTCHKCLR, ZGCC  , PRMU0, ZTAUC , ZOMCC &
      &, ZREFC , ZREFDC, ZTRAC, ZTRADC )

! Total sky reflectivities      
    CALL RRTMG_SW_REFTRA ( KLEV, KMODTS &
      &, LRTCHKCLD, ZGCO  , PRMU0, ZTAUO , ZOMCO &
      &, ZREFO , ZREFDO, ZTRAO, ZTRADO )

!
    DO JK=1,KLEV
!      
!-- combine clear and cloudy contributions for total sky
!
      IKL=KLEV+1-JK 
      ZCLEAR   = _ONE_ - PCLFR(IKL)
      ZCLOUD   = PCLFR(IKL)

      ZREF(JK) = ZCLEAR*ZREFC(JK) +ZCLOUD*ZREFO(JK)
      ZREFD(JK)= ZCLEAR*ZREFDC(JK)+ZCLOUD*ZREFDO(JK)
      ZTRA(JK) = ZCLEAR*ZTRAC(JK) +ZCLOUD*ZTRAO(JK)
      ZTRAD(JK)= ZCLEAR*ZTRADC(JK)+ZCLOUD*ZTRADO(JK)
!
!-- direct beam transmittance        
!     
      ZDBTMC     = EXP(-ZTAUC(JK)/PRMU0)
      ZDBTMO     = EXP(-ZTAUO(JK)/PRMU0)
      ZDBT(JK)   = ZCLEAR*ZDBTMC+ZCLOUD*ZDBTMO
      ZTDBT(JK+1)= ZDBT(JK)*ZTDBT(JK)
        
!-- clear-sky        
      ZDBTC(JK)=ZDBTMC
      ZTDBTC(JK+1)=ZDBTC(JK)*ZTDBTC(JK)

    END DO           
                 
 
!-- vertical quadrature producing clear-sky fluxes

    CALL RRTMG_SW_VRTQDR ( KLEV, IW &
      &, ZREFC, ZREFDC, ZTRAC , ZTRADC &
      &, ZDBTC, ZRDNDC, ZRUPC , ZRUPDC, ZTDBTC &
      &, ZCD  , ZCU   )
      

!-- vertical quadrature producing cloudy fluxes

    CALL RRTMG_SW_VRTQDR ( KLEV, IW &
      &, ZREF , ZREFD , ZTRA , ZTRAD &
      &, ZDBT , ZRDND , ZRUP , ZRUPD , ZTDBT &
      &, ZFD  , ZFU   )
 
                                              
!-- up and down-welling fluxes at levels
!-- mji: two-stream calculations go from top to bottom; reverse layer
!    indexing to go bottom to top for output arrays

    DO JK=1,KLEV+1
      IKL=KLEV+2-JK

!-- accumulation of spectral fluxes over band         
      ZBBFU(IKL) = ZBBFU(IKL) + ZINCFLX(IW)*ZFU(JK,IW)  
      ZBBFD(IKL) = ZBBFD(IKL) + ZINCFLX(IW)*ZFD(JK,IW)
      ZBBCU(IKL) = ZBBCU(IKL) + ZINCFLX(IW)*ZCU(JK,IW)
      ZBBCD(IKL) = ZBBCD(IKL) + ZINCFLX(IW)*ZCD(JK,IW)
      ZBBFDdir(IKL) = ZBBFDdir(IKL) + ZINCFLX(IW)*ztdbt_nodel(jk)
      ZBBCDdir(IKL) = ZBBCDdir(IKL) + ZINCFLX(IW)*ztdbtc_nodel(jk)

!-- accumulation of spectral fluxes over whole spectrum  
      PBBFU(IKL) = PBBFU(IKL) + ZINCFLX(IW)*ZFU(JK,IW)
      PBBFD(IKL) = PBBFD(IKL) + ZINCFLX(IW)*ZFD(JK,IW)
      PBBCU(IKL) = PBBCU(IKL) + ZINCFLX(IW)*ZCU(JK,IW)
      PBBCD(IKL) = PBBCD(IKL) + ZINCFLX(IW)*ZCD(JK,IW)
      PBBFDdir(IKL) = PBBFDdir(IKL) + ZINCFLX(IW)*ztdbt_nodel(jk)
      PBBCDdir(IKL) = PBBCDdir(IKL) + ZINCFLX(IW)*ztdbtc_nodel(jk)

!      PBBFU(JK)=PBBFU(JK)+RWGT(IW)*ZFU(JK,IW)
!      PBBFD(JK)=PBBFD(JK)+RWGT(IW)*ZFD(JK,IW)
!      PBBCU(JK)=PBBCU(JK)+RWGT(IW)*ZCU(JK,IW)
!      PBBCD(JK)=PBBCD(JK)+RWGT(IW)*ZCD(JK,IW)
!      IF (IW <= NUV) THEN
!        PUVFD(JK)=PUVFD(JK)+RWGT(IW)*ZFD(JK,IW)
!        PUVFU(JK)=PUVFU(JK)+RWGT(IW)*ZFU(JK,IW)
!        PUVCD(JK)=PUVCD(JK)+RWGT(IW)*ZCD(JK,IW)
!        PUVCU(JK)=PUVCU(JK)+RWGT(IW)*ZCU(JK,IW)
!      ELSE IF (IW == NUV+1 .AND. IW <= NVS) THEN  
!        PVSFD(JK)=PVSFD(JK)+RWGT(IW)*ZFD(JK,IW)
!        PVSFU(JK)=PVSFU(JK)+RWGT(IW)*ZFU(JK,IW)
!        PVSCD(JK)=PVSCD(JK)+RWGT(IW)*ZCD(JK,IW)
!        PVSCU(JK)=PVSCU(JK)+RWGT(IW)*ZCU(JK,IW)
!      ELSE IF (IW > NVS) THEN  
!        PNIFD(JK)=PNIFD(JK)+RWGT(IW)*ZFD(JK,IW)
!        PNIFU(JK)=PNIFU(JK)+RWGT(IW)*ZFU(JK,IW)
!        PNICD(JK)=PNICD(JK)+RWGT(IW)*ZCD(JK,IW)
!        PNICU(JK)=PNICU(JK)+RWGT(IW)*ZCU(JK,IW)
!      END IF  
    END DO

  END DO             
!-- end loop on JG

END DO                    
!-- end loop on JB

!     ------------------------------------------------------------------
RETURN
END SUBROUTINE RRTMG_SW_SPCVRT

