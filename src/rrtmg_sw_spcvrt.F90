C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

SUBROUTINE RRTM_SW_SPCVRT &
 &( KLEV   , KMOL    , KSW    , ONEMINUS,ISTART  , IEND &
 &, PAVEL  , TAVEL   , PZ     , TZ     , TBOUND  , PALBD   , PALBP &
 &, PCLFR  , PTAUC   , PASYC  , POMGC  , PTAUA   , PASYA   , POMGA , PRMU0   &
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
 &)


!**** *RRTM_SW_SPCVRT* - SPECTRAL LOOP TO COMPUTE THE SHORTWAVE RADIATION FLUXES.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE TWO-STREAM METHOD OF BARKER

!**   INTERFACE.
!     ----------

!          *RRTM_SW_SPCVRT* IS CALLED FROM *SRTM_SRTM_224GP*


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
!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPB1, JPB2, JPGPT
USE YOESRTWN , ONLY : NG

IMPLICIT NONE

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER_M :: KAER   , KLEV    , KMOL  , KSW, ISTART, IEND
INTEGER_M :: LAYTROP, LAYSWTCH, LAYLOW

REAL_B :: ONEMINUS
REAL_B :: PALBD(KSW)     , PALBP(KSW)     , PRMU0 
REAL_B :: PCLFR(JPLAY)	 , PTAUC(JPLAY,KSW),PASYC(JPLAY,KSW), POMGC(JPLAY,KSW)
REAL_B :: PTAUA(JPLAY,KSW),PASYA(JPLAY,KSW),POMGA(JPLAY,KSW)
REAL_B :: PAVEL(JPLAY)   , TAVEL(JPLAY)   , PZ(0:JPLAY)     , TZ(0:JPLAY)  , TBOUND
REAL_B :: COLDRY(JPLAY)  , COLMOL(JPLAY)  , WKL(35,JPLAY)
REAL_B :: CO2MULT(JPLAY) , COLCH4(JPLAY)  , COLCO2(JPLAY)   , COLH2O(JPLAY)
REAL_B :: COLN2O(JPLAY)  , COLO2(JPLAY)   , COLO3(JPLAY)
REAL_B :: FORFAC(JPLAY)  , FORFRAC(JPLAY) , SELFFAC(JPLAY)  , SELFFRAC(JPLAY)
REAL_B :: FAC00(JPLAY)   ,   FAC01(JPLAY) , FAC10(JPLAY)    , FAC11(JPLAY)
REAL_B :: ADJFLUX

INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)
INTEGER_M :: JP(JPLAY)    , JT(JPLAY)     , JT1(JPLAY)

REAL_B :: &
  &   PBBCD(JPLAY+1)          , PBBCU(JPLAY+1) &
  &,  PUVCD(JPLAY+1)          , PUVCU(JPLAY+1) &
  &,  PVSCD(JPLAY+1)          , PVSCU(JPLAY+1) &
  &,  PNICD(JPLAY+1)          , PNICU(JPLAY+1) &
  &,  PBBFD(JPLAY+1)          , PBBFU(JPLAY+1) &
  &,  PUVFD(JPLAY+1)          , PUVFU(JPLAY+1) &
  &,  PVSFD(JPLAY+1)          , PVSFU(JPLAY+1) &
  &,  PNIFD(JPLAY+1)          , PNIFU(JPLAY+1)

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

LOGICAL :: LRTCHK(JPLAY)

REAL_B :: &
  &   ZCLEAR      , ZCLOUD       &
  &,  ZDBT(JPLAY+1) &
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
  &   ZDBTC(JPLAY+1), ZTDBTC(JPLAY+1), ZINCFLX(JPGPT)
  

!     LOCAL INTEGER SCALARS
INTEGER_M :: IB1, IB2, IBM, IGT, IKL, IKP, IKX, IW, JB, JG, JL, JK, KGS, KMODTS
!INTEGER_M :: NDBUG

!     LOCAL REAL SCALARS
REAL_B :: ZDBTMC, ZDBTMO, ZF, ZGW, ZINCFLUX, ZREFLECT, ZWF
REAL_B :: REPCLC

!-- Output of RRTM_SW_TAUMOLn routines

REAL_B :: ZTAUG(JPLAY,16), ZTAUR(JPLAY,16), ZSSA(JPLAY,16), ZSFLXZEN(16)

!-- Output of RRTM_SW_VRTQDR routine
REAL_B :: &
  &   ZCD(JPLAY+1,JPGPT), ZCU(JPLAY+1,JPGPT) &
  &,  ZFD(JPLAY+1,JPGPT), ZFU(JPLAY+1,JPGPT)
REAL_B :: &
  &   ZBBCD(JPLAY+1)          , ZBBCU(JPLAY+1) &
  &,  ZBBFD(JPLAY+1)          , ZBBFU(JPLAY+1)
!     ------------------------------------------------------------------

!-- Two-stream model 1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
! KMODTS is set in SWREFTRA
!NDBUG=1

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
  IGT = NG(JB)

! mji out
!  print *,'=== spectral band === JB= ',JB,' ====== i.e. IBM= ',IBM,' with IGT= ',IGT
!  print *,'    ClearUpw   ClearDnw  TotalUpw  TotalDnw for band= ',JB            
!-- for each band, computes the gaseous and Rayleigh optical thickness 
!  for each g-point

!  print*, 'JB = ', JB
!  print*, 'KLEV, FAC00, FAC01, FAC10, FAC11 = ',KLEV, FAC00, FAC01, FAC10, FAC11
!  print*, 'JP, JT, JT1, ONEMINUS = ', JP, JT, JT1, ONEMINUS
!  print*, 'COLH2O, COLCH4, COLMOL = ', COLH2O, COLCH4, COLMOL
!  print*, 'LAYTROP, SELFFAC, SELFFRAC, INDSELF = ', LAYTROP, SELFFAC, SELFFRAC, INDSELF
!  print*, 'FORFAC, FORFRAC, INDFOR = ', FORFAC, FORFRAC, INDFOR

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

!  IF (NDBUG.LE.3) THEN
!    print *,'Incident Solar Flux'
!    PRINT 9010,(ZSFLXZEN(JG),JG=1,16)
!9010 format(1x,'SolFlx ',16F8.4)
!    print *,'Optical thickness for molecular absorption for JB= ',JB 
!    DO JK=1,KLEV
!      PRINT 9011,JK,(ZTAUG(JK,JG),JG=1,16)
!9011  format(1x,'TauGas ',I3,16E9.2)
!    END DO
!    print *,'Optical thickness for Rayleigh scattering for JB= ',JB 
!    DO JK=1,KLEV
!      PRINT 9012,JK,(ZTAUR(JK,JG),JG=1,16)
!9012  format(1x,'TauRay ',I3,16E9.2)
!    END DO
!    print *,'Cloud optical properties for JB= ',JB
!    DO JK=1,KLEV
!      PRINT 9013,JK,PCLFR(JK),PTAUC(JK,IBM),POMGC(JK,IBM),PASYC(JK,IBM)
!9013  format(1x,'Cloud optprop ',I3,f8.4,f8.3,2f8.5)
!    END DO
!  END IF

  DO JK=1,KLEV+1
    ZBBCD(JK)=_ZERO_
    ZBBCU(JK)=_ZERO_
    ZBBFD(JK)=_ZERO_
    ZBBFU(JK)=_ZERO_
  END DO

  DO JG=1,IGT
    IW=IW+1

!    IF (NDBUG.LE.1) THEN
!      print *,' === JG= ',JG,' === for JB= ',JB,' with IW, IBM, JPLAY, KLEV=',IW,IBM,JPLAY,KLEV
!    END IF

! mji - add adjustment for correct Earth/Sun distance
!    ZINCFLX(IW)=ZSFLXZEN(JG)*PRMU0
!    ZINCFLUX   =ZINCFLUX+ZSFLXZEN(JG)*PRMU0           
    ZINCFLX(IW)=ADJFLUX*ZSFLXZEN(JG)*PRMU0
    ZINCFLUX   =ZINCFLUX+ADJFLUX*ZSFLXZEN(JG)*PRMU0           

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
!----- surface values
    ZDBTC(KLEV+1) =_ZERO_
    ZTRAC(KLEV+1) =_ZERO_
    ZTRADC(KLEV+1)=_ZERO_
    ZREFC(KLEV+1) =PALBP(IBM)
    ZREFDC(KLEV+1)=PALBP(IBM)
    ZRUPC(KLEV+1) =PALBP(IBM)
    ZRUPDC(KLEV+1)=PALBP(IBM)
           
!-- total sky    
!----- TOA direct beam    
    ZTDBT(1)=1._JPRB
!----- surface values
    ZDBT(KLEV+1) =_ZERO_
    ZTRA(KLEV+1) =_ZERO_
    ZTRAD(KLEV+1)=_ZERO_
    ZREF(KLEV+1) =PALBD(IBM)
    ZREFD(KLEV+1)=PALBD(IBM)
    ZRUP(KLEV+1) =PALBD(IBM)
    ZRUPD(KLEV+1)=PALBD(IBM)
!    if (NDBUG < 2) print *,'SWSPCTRL after 1 with JB,JG,IBM and IW= ',JB,JG,IBM,IW
    
    
    DO JK=1,KLEV

!-- NB: a two-stream calculations from top to bottom, but RRTM_SW quantities 
!       are given bottom to top and are reverse here:

      IKL=KLEV+1-JK

!-- clear-sky optical parameters      
      LRTCHK(JK)=.TRUE.

!      print 9000,JK,JG,IKL,ZTAUR(IKL,JG),ZTAUG(IKL,JG),PTAUC(IKL,IBM)
!9000  format(1x,'Cloud quantities ',3I4,3E12.5)

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

!-- total sky optical parameters        
      ZTAUO(JK) = ZTAUR(IKL,JG) + ZTAUG(IKL,JG) + PTAUA(IKL,IBM) + PTAUC(IKL,IBM)
      ZOMCO(JK) = PTAUA(IKL,IBM)*POMGA(IKL,IBM) + PTAUC(IKL,IBM)*POMGC(IKL,IBM) &
        &       + ZTAUR(IKL,JG)*_ONE_
      ZGCO (JK) = (PTAUC(IKL,IBM)*POMGC(IKL,IBM)*PASYC(IKL,IBM)  &
        &       +  PTAUA(IKL,IBM)*POMGA(IKL,IBM)*PASYA(IKL,IBM)) &
        &       /  ZOMCO(JK)
      ZOMCO(JK) = ZOMCO(JK) / ZTAUO(JK)

!      if (NDBUG <2) THEN
!        print 9001,JK,JG,LRTCHK(JK),0.00,ZTAUC(JK),ZOMCC(JK),ZGCC(JK),ZTAUR(IKL,JG),ZTAUG(IKL,JG)
!9001    format(1x,'SPCVRT; clear :',2I3,L4,7(1x,E13.6))
!        print 9002,JK,JG,LRTCHK(JK),PCLFR(JK),ZTAUO(JK),ZOMCO(JK),ZGCO(JK) &
!          &,PTAUC(IKL,IBM),POMGC(IKL,IBM),PASYC(IKL,IBM)
!9002    format(1x,'SPCVRT; total0:',2I3,L4,7(1x,E13.6))
!      end if

    END DO    
!    if (NDBUG < 2) print *,'SWSPCTRL after 2'
   
    CALL RRTM_SW_REFTRA ( KLEV, KMODTS &
      &, LRTCHK, ZGCC  , PRMU0, ZTAUC , ZOMCC &
      &, ZREFC , ZREFDC, ZTRAC, ZTRADC )

!    if (NDBUG < 2) print *,'SWSPCTR after SWREFTRA for clear-sky'
    
      
!-- Delta scaling    
    DO JK=1,KLEV
      IKL=KLEV+1-JK
      LRTCHK(JK)=.FALSE.
      ZF=ZGCO(JK)*ZGCO(JK)
      ZWF=ZOMCO(JK)*ZF
      ZTAUO(JK)=(1._JPRB-ZWF)*ZTAUO(JK)
      ZOMCO(JK)=(ZOMCO(JK)-ZWF)/(1._JPRB-ZWF)
      ZGCO (JK)=(ZGCO(JK)-ZF)/(1._JPRB-ZF)
      LRTCHK(JK)=(PCLFR(IKL) > REPCLC)

!      if (NDBUG < 2) THEN 
!        print 9003,JK,LRTCHK(JK),PCLFR(IKL),ZTAUO(JK),ZOMCO(JK),ZGCO(JK) &
!          &,PTAUC(IKL,IBM),POMGC(IKL,IBM),PASYC(IKL,IBM)
9003    format(1x,'totalD:',I3,L4,7(1x,E13.6))
!      end if

    END DO
!    if (NDBUG < 2) print *,'SWSPCTR after Delta scaling'
    
    CALL RRTM_SW_REFTRA ( KLEV, KMODTS &
      &, LRTCHK, ZGCO  , PRMU0, ZTAUO , ZOMCO &
      &, ZREFO , ZREFDO, ZTRAO, ZTRADO )
!    if (NDBUG < 2) print *,'SWSPCTR after SWREFTRA for cloudy'

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

!      if (NDBUG < 2) print 9200,JK,ZREFC(JK),ZREFDC(JK),ZTRAC(JK),ZTRADC(JK),ZDBTC(JK),ZTDBTC(JK+1)
!      if (NDBUG < 2) print 9199,JK,ZREF(JK),ZREFD(JK),ZTRA(JK),ZTRAD(JK),ZDBT(JK),ZTDBT(JK+1)
9199  format(1x,'Comb total:',I3,6E13.6)
9200  format(1x,'Comb clear:',I3,6E13.6)

    END DO           
!    if (NDBUG < 2) print *,'RRTM_SW_SPCVRT after combining clear and cloudy'
                 
 
!-- vertical quadrature producing clear-sky fluxes

!    print *,'RRTM_SW_SPCVRT after 3 before RRTM_SW_VRTQDR clear'
    
    CALL RRTM_SW_VRTQDR ( KLEV, IW &
      &, ZREFC, ZREFDC, ZTRAC , ZTRADC &
      &, ZDBTC, ZRDNDC, ZRUPC , ZRUPDC, ZTDBTC &
      &, ZCD  , ZCU   )
      
!    IF (NDBUG < 2) THEN
!      print *,'RRTM_SW_SPCVRT out of RRTM_SW_VRTQDR for clear IW=',IW  
!      DO JK=1,KLEV+1
!        print 9201,JK,ZCD(JK,IW),ZCU(JK,IW)
9201    format(1x,'clear-sky contrib to fluxes',I3,2F12.4)
!      END DO      
!    END IF

!-- vertical quadrature producing cloudy fluxes

!    print *,'RRTM_SW_SPCVRT after 4 before RRTM_SW_VRTQDR cloudy'
    
    CALL RRTM_SW_VRTQDR ( KLEV, IW &
      &, ZREF , ZREFD , ZTRA , ZTRAD &
      &, ZDBT , ZRDND , ZRUP , ZRUPD , ZTDBT &
      &, ZFD  , ZFU   )
 
!    IF (NDBUG < 2) THEN     
!      print *,'RRTM_SW_SPCVRT out of RRTM_SW_VRTQDR for cloudy IW=',IW
!      DO JK=1,KLEV+1
!        print 9202,JK,ZFD(JK,IW),ZFU(JK,IW)
9202    format(1x,'cloudy sky contrib to fluxes',I3,2F12.4)
!      END DO      
!    END IF

                                              
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

!-- accumulation of spectral fluxes over whole spectrum  
      PBBFU(IKL) = PBBFU(IKL) + ZINCFLX(IW)*ZFU(JK,IW)
      PBBFD(IKL) = PBBFD(IKL) + ZINCFLX(IW)*ZFD(JK,IW)
      PBBCU(IKL) = PBBCU(IKL) + ZINCFLX(IW)*ZCU(JK,IW)
      PBBCD(IKL) = PBBCD(IKL) + ZINCFLX(IW)*ZCD(JK,IW)

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
!      if (NDBUG < 2) then
!      if (JG.EQ.IGT) THEN 
!          print 9206,JB,JG,JK,IW,PBBCU(JK),PBBCD(JK),PBBFU(JK),PBBFD(JK)
9206      format(1x,'fluxes up to:',3I3,I4,6E13.6)       
!      end if
    END DO

!    if (NDBUG < 2) print *,'RRTM_SW_SPCVRT end of JG=',JG,' for JB=',JB,' i.e. IW=',IW

  END DO             
!-- end loop on JG

! mji out
!  DO JK=1,KLEV+1  
!    print 9207,JK,ZBBCU(JK),ZBBCD(JK),ZBBFU(JK),ZBBFD(JK)
!9207 format(1x,I3,4F10.3)
!  END DO
!!  print *,' --- JB= ',JB,' with IB1, IB2= ',IB1,IB2

END DO                    
!-- end loop on JB

!if (NDBUG < 2) print *,'RRTM_SW_SPCVRT about to come out'
!     ------------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_SW_SPCVRT

