C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

!***************************************************************************
!                                                                          *
!                               RRTM_SW                                    *
!                                                                          *
!                                                                          *
!                                                                          *
!                   A RAPID RADIATIVE TRANSFER MODEL                       *
!                    FOR THE SOLAR SPECTRAL REGION                         *
!                                                                          *
!            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                  *
!                        840 MEMORIAL DRIVE                                *
!                        CAMBRIDGE, MA 02139                               *
!                                                                          *
!                                                                          *
!                           ELI J. MLAWER                                  *   
!                         JENNIFER DELAMERE                                *
!                         STEVEN J. TAUBMAN~                               *
!                         SHEPARD A. CLOUGH                                *
!                                                                          *
!                                                                          *
!                         ~currently at GFDL                               *
!                                                                          *
!                                                                          *
!                                                                          *
!                       email:  mlawer@aer.com                             *
!                                                                          *
!        The authors wish to acknowledge the contributions of the          *
!        following people:  Patrick D. Brown, Michael J. Iacono,           *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                    *
!                                                                          *
!***************************************************************************

PROGRAM RRTM_SW_224

!-- Interface to RRTM_SW; conversion to F90 formatting; addition of 
!   2-stream radiative transfer
!     J.-J. Morcrette, ECMWF, 030225
!-- Additional modifications for GCM application
!     This stand-alone version uses RRTATM for input
!     M. J. Iacono, AER, Inc., August 2003

!- Input from calling routine:
!- NOTE: In RRTM_SW, the layer dimension goes from bottom to top
!  KLON, JPLON = number of longitudes (1 in single-column mode)
!  JPLAY = maximum number of layers
!  KLEV, NLAYERS = number of layers
!  KSW, JPSW = number of sw bands (14)
!  NSTR = Number of streams for radiative transfer (= 2)
!  KOVLP = cloud overlap method (1=Maximum-random, 2=Maximum, or 3=Random)
!  ZENITH = Cosine of solar zenith angle
!  ADJFLUX = Incoming solar flux adjustment for current Earth/Sun distance
!  PAER = Aerosol optical thickness (6 types; see susrtaer.F90)
!  SEMISS = surface emissivity by band
!  PAPH = pressure levels (mb)
!  PAP = pressure layers (mb)
!  PTS = surface temp (K)
!  PTH = temperature levels (K)
!  PT = temperature layers (K)
!-NOTE: Water vapor, CO2, ozone, CH4, and N2O are needed in units of vmr
!- If passed in as mmr, use molecular weights below to convert to vmr
!  PQ = water vapor 
!  PCCO2 = CO2 
!  POZN = ozone
!  PCH4 = methane
!  PN2O = nitrous ozide
!  PCLFR = cloud fraction

!-Variables related to the cloud optical properties (original RRTM_SW)
!- (See susrtop.F90 for further description):
!  INFLAG = Flag for cloud optical properties
!  ICEFLAG = Flag for ice particle specification
!  LIQFLAG = Flag for liquid droplet specification
!-For INFLAG=0/LIQFLAG=0/ICEFLAG=0:
!  CLDDAT1 = Total (ice and water) cloud optical depth for layer
!  CLDDAT2 = Singe-scattering albedo for layer
!  CLDDATMOM= Moments of the phase function (from 0 to NSTR)
!-For INFLAG=2/LIQFLAG=1/ICEFLAG=3:
!  CLDDAT1 = Cloud water path for the layer (g/m2)
!  CLDDAT2 = Fraction of the cloud layer's water path in the form 
!           of ice particles
!  CLDDAT3 = Generalized effective size of the ice crystals, dge (microns)
!           (see Q. Fu, 1996 for a definition of this parameter).
!           Valid sizes are 5.0 to 140.0 microns.
!  CLDDAT4 = Liquid droplet effective radius, re, (microns)
!           of ice particles
!           Valid sizes are 2.5 to 60.0 microns.

!-Return to calling routine:
!  PFUP = total sky flux (W/m2)
!  PFDOWN = total sky flux (W/m2)
!  PCUP = clear sky up flux (W/m2)
!  PCDOWN = clear sky down flux (W/m2)
!  PHEAT = total sky heating rate (K/day)
!  PHEAC = clear sky heating rate (K/day)


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLON, JPLAY, JPAER, JPSW, JPBAND, JPB1, JPB2
USE YOESRTAER, ONLY : RSRTAUA, RSRPIZA, RSRASYA 
USE YOESRTWN,  ONLY : WAVENUM1, WAVENUM2

IMPLICIT NONE

!-- Input arguments

INTEGER_M :: KLON, KLEV, KSW, KOVLP

REAL_B :: PAER(JPLON,6,JPLAY)
REAL_B :: PAP(JPLON,JPLAY)     , PAPH(JPLON,JPLAY+1) , PDP(JPLON,JPLAY)
REAL_B :: PTS(JPLON)           , PT(JPLON,JPLAY)     , PTH(JPLON,JPLAY+1)
REAL_B :: PQ(JPLON,JPLAY)      , POZN(JPLON,JPLAY)   , PCLFR(JPLON,JPLAY)
REAL_B :: PCCO2, PCH4, PN2O

!-- Output arguments

REAL_B :: PFCS(JPLON,JPLAY+1), PFLS(JPLON,JPLAY+1)
REAL_B :: PHEAC(JPLON,JPLAY), PHEAT(JPLON,JPLAY)
REAL_B :: PFDOWN(JPLON,JPLAY+1), PCDOWN(JPLON,JPLAY+1)
REAL_B :: PFUP(JPLON,JPLAY+1), PCUP(JPLON,JPLAY+1)

!-----------------------------------------------------------------------

!-- dummy integers
INTEGER_M :: JB

INTEGER_M :: ICLDATM, INFLAG, ICEFLAG, LIQFLAG, NMOL, NSTR, IWR
INTEGER_M :: I, IK, IMOL, J1, J2, JA, JAE, JL, JK, JMOM, JSW
INTEGER_M :: NLAYERS, ISTART, IEND, IFLAG, IOUT, ICLD, IAER, ISCCOS, IDELM, INDFORM

INTEGER_M :: LAYTROP, LAYSWTCH, LAYLOW
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!-- dummy reals

REAL_B :: CLDFRAC(JPLAY), CLDDAT1(JPLAY), CLDDAT2(JPLAY), CLDDAT3(JPLAY), CLDDAT4(JPLAY)
REAL_B :: CLDDATMOM(0:16,JPLAY)
REAL_B :: TAUCLDORIG(JPLAY), TAUCLOUD(JPLAY,JPBAND), SSACLOUD(JPLAY,JPBAND)
REAL_B :: XMOM(0:16,JPLAY,JPBAND)
REAL_B :: TAUAER(JPLAY,JPBAND), SSAAER(JPLAY,JPBAND), PHASE(32,JPLAY,JPBAND)
REAL_B :: ZENITH, ADJFLUX, SEMISS(JPBAND)

REAL_B :: PZ(0:JPLAY)   , TZ(0:JPLAY)   , PAVEL(JPLAY)  , TAVEL(JPLAY)   , HTR(0:JPLAY)
REAL_B :: COLDRY(JPLAY) , COLMOL(JPLAY) , WKL(35,JPLAY)
REAL_B :: CO2MULT(JPLAY), COLCH4(JPLAY) , COLCO2(JPLAY) , COLH2O(JPLAY)
REAL_B :: COLN2O(JPLAY) , COLO2(JPLAY)  , COLO3(JPLAY)
REAL_B :: FORFAC(JPLAY) , FORFRAC(JPLAY), SELFFAC(JPLAY), SELFFRAC(JPLAY)
REAL_B :: FAC00(JPLAY)  , FAC01(JPLAY)  , FAC10(JPLAY)  , FAC11(JPLAY)
REAL_B :: TBOUND        , ONEMINUS	, ZRMU0
REAL_B :: ZALBD(JPSW)    , ZALBP(JPSW)    , ZCLFR(JPLAY)
REAL_B :: ZTAUC(JPLAY,JPSW), ZASYC(JPLAY,JPSW), ZOMGC(JPLAY,JPSW)
REAL_B :: ZTAUA(JPLAY,JPSW), ZASYA(JPLAY,JPSW), ZOMGA(JPLAY,JPSW)

REAL_B :: ZBBCD(JPLAY+1), ZBBCU(JPLAY+1), ZBBFD(JPLAY+1), ZBBFU(JPLAY+1)
!REAL_B :: ZUVCD(JPLAY+1), ZUVCU(JPLAY+1), ZUVFD(JPLAY+1), ZUVFU(JPLAY+1)
!REAL_B :: ZVSCD(JPLAY+1), ZVSCU(JPLAY+1), ZVSFD(JPLAY+1), ZVSFU(JPLAY+1)
!REAL_B :: ZNICD(JPLAY+1), ZNICU(JPLAY+1), ZNIFD(JPLAY+1), ZNIFU(JPLAY+1)

REAL_B :: ZCLEAR, ZCLOUD, ZEPSEC, ZTOTCC, ZDPGCP
REAL_B :: DIRDOWNFLUX, DIFDOWNFLUX

REAL_B :: RDAY, RG, RMD, RKBOL, RNAVO, R, RD, RCPD, RCDAY

!-- character strings for RRTM output
CHARACTER PAGE

CHARACTER*50 OUTFORM(7)

!     Setup format statements for output

DATA OUTFORM &
   &/'(1X,I3,3X,F7.6,4X,4(F10.4,4X),F10.4,4X,F10.5)',&
   & '(1X,I3,4X,F6.5,4X,4(F10.4,4X),F10.4,4X,F10.5)',&
   & '(1X,I3,5X,F5.4,4X,4(F10.4,4X),F10.4,4X,F10.5)',&
   & '(1X,I3,5X,F5.3,4X,4(F10.4,4X),F10.4,4X,F10.5)',&
   & '(1X,I3,5X,F5.2,4X,4(F10.4,4X),F10.4,4X,F10.5)',&
   & '(1X,I3,5X,F5.1,4X,4(F10.4,4X),F10.4,4X,F10.5)',&
   & '(1X,I3,4X,F6.1,4X,4(F10.4,4X),F10.4,4X,F10.5)'/

!      HVRRTM = '$Revision$'
!      CHARACTER*15 HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
!     *            HVDUM1,HVRUTL,HVREXT,
!     *            HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
!     *            HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR
!
!      CHARACTER*15 HVRKG16,HVRKG17,HVRKG18,HVRKG19,
!     *             HVRKG20,HVRKG21,HVRKG22,HVRKG23,HVRKG24,
!     *             HVRKG25,HVRKG27,HVRKG28,HVRKG29

!-----------------------------------------------------------------------
!-- calculate information needed by the radiative transfer routine 

ZEPSEC  = 1.E-06_JPRB
ONEMINUS=_ONE_ -  ZEPSEC

NSTR	= 2
NMOL	= 6
KLON    = JPLON
KLEV    = JPLAY
KSW     = JPSW
KOVLP   = 1

IWR = 10
PAGE = CHAR(12)

!     RCDAY is the factor by which one must multiply delta-flux/
!     delta-pressure, with flux in w/m-2 and pressure in mbar, to get
!     the heating rate in units of degrees/day.  It is the equivalent
!     to the constant HEATFAC used in RRTM_LW, and it is equal to
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(86400)(1e-5)/(1.004)
!     Constants are set for consistency between RRTM_LW and RRTM_SW.
RDAY = 86400.
RG   = 9.8066
RMD  = 28.9644
! RKBOL= 1.380658E-23
! RNAVO= 6.0221367E+23
RKBOL= 1.3806503E-23
RNAVO= 6.02214199E+23
R    = RNAVO*RKBOL
RD   = 1000.*R/RMD
RCPD = 3.5*RD
! RCDAY= RDAY * RG / RCPD
RCDAY = 8.4391

! mji - Initialization routine for RRTM_SW (called once only).
! This calls SUSRTM, KGBn, SUSRTOP and should be placed in the
! GCM's initialization section, if RRTM_SW is used in a GCM.

CALL RRTM_SW_INIT


! Main longitude loop
DO JL = 1, KLON

!- Input atmospheric profile from INPUT_RRTM.
   CALL READPROF(NLAYERS,IOUT,ICLD,IAER,ISCCOS,IDELM, &
  &     PAVEL,TAVEL,PZ,TZ,ZENITH,ADJFLUX,SEMISS, &
  &     WKL,COLDRY,INFLAG,ICEFLAG,LIQFLAG, &
  &     CLDFRAC,CLDDAT1,CLDDAT2,CLDDAT3,CLDDAT4,CLDDATMOM, &
  &     TAUAER,SSAAER,PHASE)

   KLEV = NLAYERS
   ISTART = JPB1
   IEND = JPB2
   IFLAG = IOUT

1000 CONTINUE
   IF (IFLAG .GT. 0 .AND. IFLAG .LE. JPB2) THEN
      ISTART = IFLAG
      IEND = IFLAG
   ENDIF

!  DO JK=1,KLEV 
!    print 9101,JK,CLDFRAC(JK),CLDDAT1(JK),CLDDAT2(JK),CLDDAT3(JK)&
!    &,CLDDAT4(JK),(CLDDATMOM(JMOM,JK),JMOM=0,NSTR)
!9101 format(1x,'RRTM_SW_EC224GP Cld :',I3,f7.4,7E12.5)
!  END DO

   IF (ICLD.EQ.1) THEN
      CALL RRTM_SW_CLDPROP &
     &( KLEV, ICLDATM, INFLAG, ICEFLAG, LIQFLAG, NSTR &
     &, CLDFRAC, CLDDAT1, CLDDAT2, CLDDAT3, CLDDAT4, CLDDATMOM &
     &, TAUCLDORIG, TAUCLOUD, SSACLOUD, XMOM &
     &)
   ENDIF

!- coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients
! Passed in from READPROF in stand-alone version

   TBOUND=TZ(0)

   DO JK=1,KLEV 
      PCLFR(JL,JK) = CLDFRAC(JK)
!     print 9102,JK,PCLFR(JL,JK)
!9102 format(1x,'RRTM_SW_EC224GP PCLFR :',I3,f7.4)
   END DO

   ZCLEAR=_ONE_
   ZCLOUD=_ZERO_
   ZTOTCC=_ZERO_
   DO JK = 1, KLEV
      PDP(JL,JK)= (PZ(JK-1)-PZ(JK))

! Set up for cloud overlap
      IF (KOVLP == 1) THEN
         ZCLEAR=ZCLEAR*(_ONE_-MAX(PCLFR(JL,JK),ZCLOUD)) &
        &   /(_ONE_-MIN(ZCLOUD,_ONE_-ZEPSEC))
         ZCLOUD=PCLFR(JL,JK)
         ZTOTCC=_ONE_-ZCLEAR
      ELSE IF (KOVLP == 2) THEN
         ZCLOUD=MAX(ZCLOUD,PCLFR(JL,JK))
         ZCLEAR=_ONE_-ZCLOUD
         ZTOTCC=ZCLOUD
      ELSE IF (KOVLP == 3) THEN
         ZCLEAR=ZCLEAR*(_ONE_-PCLFR(JL,JK))
         ZCLOUD=_ONE_-ZCLEAR
         ZTOTCC=ZCLOUD
      END IF

!    print 9202,JK,KOVLP,ZCLEAR,ZCLOUD,ZTOTCC
!9202 format(1x,'RRTM_SW, OVLP: ',2I3,3F8.1)

   END DO

! DO JK = 1, KLEV
!    print 9200,JK,PAVEL(JK),TAVEL(JK),(WKL(JA,JK),JA=1,4),WKL(6,JK),COLDRY(JK)
!9200 format(1x,'RRTM_SW ',I3,2F8.1,6E13.5)
!    print 9201,JK,PDP(JL,JK),PZ(JK)
!9201 format(1x,'RRTM_SW ',I3,2F8.1)
! END DO

   IF (ZTOTCC == _ZERO_) THEN
      DO JK=1,KLEV
         ZCLFR(JK)=_ZERO_   
      END DO
   ELSE
      DO JK=1,KLEV
         ZCLFR(JK)=PCLFR(JL,JK)/ZTOTCC
      END DO
   END IF

!  print *,'just before RRTM_SW_SETCOEF'

   CALL RRTM_SW_SETCOEF &
     &( KLEV   , NMOL &
     &, PAVEL  , TAVEL   , PZ     , TZ     , TBOUND &
     &, COLDRY , WKL     &
     &, LAYTROP, LAYSWTCH, LAYLOW &
     &, CO2MULT, COLCH4  , COLCO2 , COLH2O , COLMOL  , COLN2O  , COLO2 , COLO3 &
     &, FORFAC , FORFRAC , INDFOR , SELFFAC, SELFFRAC, INDSELF &
     &, FAC00  , FAC01   , FAC10  , FAC11  &
     &, JP     , JT      , JT1	 &
     &)
  

!- Call the radiation transfer routine

!- Put cloud optical properties in arrays for 2-stream 
  
  DO JSW=1,KSW
     ZALBD(JSW)=1.-SEMISS(JPB1-1+JSW)
     ZALBP(JSW)=1.-SEMISS(JPB1-1+JSW)
     DO JK=1,KLEV
        ZTAUC(JK,JSW)=_ZERO_
        ZASYC(JK,JSW)=_ZERO_
        ZOMGC(JK,JSW)=_ZERO_
        IF (CLDFRAC(JK) .GE. ZEPSEC) THEN
           ZTAUC(JK,JSW) = TAUCLOUD(JK,JPB1-1+JSW)
           ZASYC(JK,JSW) = XMOM(1,JK,JPB1-1+JSW)
           ZOMGC(JK,JSW) = SSACLOUD(JK,JPB1-1+JSW)
        ENDIF
!      print 9002,JSW,JK,ZCLFR(JK),ZTAUC(JK,JSW),ZASYC(JK,JSW),ZOMGC(JK,JSW)
!!9002  format(1x,'srtm_srtm_224gp ClOPropECmodel ',2I3,f8.4,3E12.5)
!9002  format(1x,'RRTM_SW_EC224gp CldpropRRTM_SW ',2I3,f8.4,3E12.5)
     END DO
  END DO

!- mixing of aerosols
!- Put aerosol optical properties in arrays for 2-stream 

! IAER = 0: no aerosols
   IF (IAER.EQ.0) THEN

     DO JSW=1,KSW
        DO JK=1,KLEV
           ZTAUA(JK,JSW)=_ZERO_
           ZASYA(JK,JSW)=_ZERO_
           ZOMGA(JK,JSW)=_ZERO_
       ENDDO
     ENDDO

! IAER = 6: Use ECMWF six aerosol types. See susrtaer.F90 for details.
! Set aerosol optical thickness (PAER) here manually for each aerosol 
! and layer.
   ELSEIF (IAER.EQ.6) THEN

     DO JK = 1, KLEV
       DO JAE = 1, JPAER
         PAER(JL,JAE,JK) = 1.0E-15
!         PAER(JL,JAE,JK) = 1.0E-04
       ENDDO
     ENDDO
     DO JSW=1,KSW
        DO JK=1,KLEV
           ZTAUA(JK,JSW)=_ZERO_
           ZASYA(JK,JSW)=_ZERO_
           ZOMGA(JK,JSW)=_ZERO_
           DO JAE=1,JPAER
              ZTAUA(JK,JSW)=ZTAUA(JK,JSW)+RSRTAUA(JSW,JAE)*PAER(JL,JAE,JK)
              ZOMGA(JK,JSW)=ZOMGA(JK,JSW)+RSRTAUA(JSW,JAE)*PAER(JL,JAE,JK) &
             & *RSRPIZA(JSW,JAE)
              ZASYA(JK,JSW)=ZASYA(JK,JSW)+RSRTAUA(JSW,JAE)*PAER(JL,JAE,JK) &
             &*RSRPIZA(JSW,JAE)*RSRASYA(JSW,JAE)
           END DO
           IF (ZOMGA(JK,JSW) /= _ZERO_) THEN
              ZASYA(JK,JSW)=ZASYA(JK,JSW)/ZOMGA(JK,JSW)
           END IF
           IF (ZTAUA(JK,JSW) /= _ZERO_) THEN
              ZOMGA(JK,JSW)=ZOMGA(JK,JSW)/ZTAUA(JK,JSW)
           END IF
        END DO
     END DO

! IAER=10: Direct specification of aerosol properties from IN_AER_RRTM
   ELSEIF (IAER.EQ.10) THEN

     DO JSW=1,KSW
        DO JK=1,KLEV
           ZTAUA(JK,JSW)=TAUAER(JK,JPB1-1+JSW)
           ZASYA(JK,JSW)=PHASE(1,JK,JPB1-1+JSW)
           ZOMGA(JK,JSW)=SSAAER(JK,JPB1-1+JSW)
        END DO
     END DO

   ENDIF

!  DO JK=1,KLEV
!    print 9013,JK,(PAER(JL,JAE,JK),JAE=1,6)
!  END DO
!9013 format(1x,'Aerosol opt.thickness: ',I3,6F10.5)
!  DO JSW=1,KSW
!    print 9003,JSW,(RSRTAUA(JSW,JAE),JAE=1,6)
!    print 9003,JSW,(RSRPIZA(JSW,JAE),JAE=1,6)
!    print 9003,JSW,(RSRASYA(JSW,JAE),JAE=1,6)
!9003 format(1x,'Aerosol opt.prop: ',I3,6F10.5)
!  END DO
 
!  print *,'just after setting aerosol properties'

   ZRMU0=ZENITH

!  print*, 'ZRMU0, ZENITH = ', ZRMU0, ZENITH

   DO JK=1,KLEV+1
      ZBBCU(JK)=_ZERO_
      ZBBCD(JK)=_ZERO_
      ZBBFU(JK)=_ZERO_
      ZBBFD(JK)=_ZERO_
!      ZUVCU(JK)=_ZERO_
!      ZUVCD(JK)=_ZERO_
!      ZUVFU(JK)=_ZERO_
!      ZUVFD(JK)=_ZERO_
!      ZVSCU(JK)=_ZERO_
!      ZVSCD(JK)=_ZERO_
!      ZVSFU(JK)=_ZERO_
!      ZVSFD(JK)=_ZERO_
!      ZNICU(JK)=_ZERO_
!      ZNICD(JK)=_ZERO_
!      ZNIFU(JK)=_ZERO_
!      ZNIFD(JK)=_ZERO_
   END DO

!  print *,'just before calling STRM_SPCVRT for JL=',JL

   CALL RRTM_SW_SPCVRT &
     &( KLEV   , NMOL    , KSW    ,ONEMINUS, ISTART  , IEND &
     &, PAVEL  , TAVEL   , PZ     , TZ     , TBOUND  , ZALBD   , ZALBP &
     &, ZCLFR  , ZTAUC   , ZASYC  , ZOMGC  , ZTAUA   , ZASYA   , ZOMGA , ZRMU0   &
     &, COLDRY , WKL     , ADJFLUX  &	 
     &, LAYTROP, LAYSWTCH, LAYLOW &
     &, CO2MULT, COLCH4  , COLCO2 , COLH2O , COLMOL  , COLN2O  , COLO2 , COLO3 &
     &, FORFAC , FORFRAC , INDFOR , SELFFAC, SELFFRAC, INDSELF &
     &, FAC00  , FAC01   , FAC10  , FAC11  &
     &, JP     , JT      , JT1	 &
!     &, ZBBFD  , ZBBFU   , ZUVFD  , ZUVFU  , ZVSFD   , ZVSFU   , ZNIFD , ZNIFU &
!     &, ZBBCD  , ZBBCU   , ZUVCD  , ZUVCU  , ZVSCD   , ZVSCU   , ZNICD , ZNICU &
     &, ZBBFD  , ZBBFU &
     &, ZBBCD  , ZBBCU &
     &)


!  print *,'SRTM_SRTM_224GP before potential scaling'
!  DO JK=1,KLEV+1
!    print 9004,JK,ZBBCU(JK),ZBBCD(JK),ZBBFU(JK),ZBBFD(JK)
9004 format(1x,'Clear-sky and total fluxes U & D ',I3,4F10.3)
!  END DO


!   print *,'Clear-sky and total fluxes over whole SW spectrum'
!   print *,'    ClearUpw   ClearDnw  TotalUpw  TotalDnw'
   DO JK=1,KLEV+1
      PCUP(JL,JK)=ZBBCU(JK)
      PCDOWN(JL,JK)=ZBBCD(JK)
      PFUP(JL,JK)=(_ONE_-ZCLEAR)*ZBBFU(JK)+ZCLEAR*ZBBCU(JK)
      PFDOWN(JL,JK)=(_ONE_-ZCLEAR)*ZBBFD(JK)+ZCLEAR*ZBBCD(JK)
   END DO
!   DO JK=1,KLEV+1
!     print 9005,JK,PCUP(JL,JK),PCDOWN(JL,JK),PFUP(JL,JK),PFDOWN(JL,JK)
!9005 format(1x,I3,4F10.3)
!   END DO

   DO JK = 1 , KLEV+1
      PFLS(JL,JK) = PFDOWN(JL,JK) - PFUP(JL,JK)
      PFCS(JL,JK) = PCDOWN(JL,JK) - PCUP(JL,JK)
!  print 9504,JK,ZFLUC(JL,2,JK),ZFLUC(JL,1,JK),PFCS(JL,JK) &
!    &,          ZFLUX(JL,2,JK),ZFLUX(JL,1,JK),PFLS(JL,JK)
!9504 format(1x,'SRTM-FLX',1I3,2(2x,f10.3,1x,f10.3,1x,f10.3))
   ENDDO

   HTR(KLEV) = 0.
   DO JK=1,KLEV
      ZDPGCP=RCDAY/PDP(JL,JK)
      PHEAC(JL,JK)=(PFCS(JL,JK+1)-PFCS(JL,JK))*ZDPGCP
      PHEAT(JL,JK)=(PFLS(JL,JK+1)-PFLS(JL,JK))*ZDPGCP
      HTR(JK-1) = PHEAT(JL,JK)
!  print 9505,JK,ZCEAT(JL,JK),ZHEAT(JL,JK),PCLFR(JL,JK)
!9505 format(1x,'SRTM-HR',I3,3(10x,f10.3,14x))
   END DO

! Stand-alone version output
   IF (IOUT .LT. 0) GO TO 4000
   DIRDOWNFLUX = 0.
   DIFDOWNFLUX = 0.

! ***    Process output for this atmosphere.
   OPEN (IWR,FILE='OUTPUT_RRTM',FORM='FORMATTED')
   IF ( ISTART .NE. IEND) THEN
      IEND = IEND-1
      ISTART = JPB2
   ENDIF

   if (isccos .eq. 1) then 
      write(iwr,9880) 
   elseif (isccos .eq. 2) then
      write(iwr,9881)
   else
      write(iwr,9879)
   endif

   if (idelm .eq. 0) then
      write(iwr,9883)
   else
      write(iwr,9882)
   endif

   WRITE(IWR,9899)WAVENUM1(ISTART),WAVENUM2(IEND)
   WRITE(IWR,9900)
   WRITE(IWR,9901)

   DO 3000 I = KLEV, 0, -1
      IF (PZ(I) .LT. 1.E-2) THEN
         INDFORM = 1
      ELSEIF (PZ(I) .LT. 1.E-1) THEN
         INDFORM = 2
      ELSEIF (PZ(I) .LT. 1.) THEN
         INDFORM = 3
      ELSEIF (PZ(I) .LT. 10.) THEN
         INDFORM = 4
      ELSEIF (PZ(I) .LT. 100.) THEN
         INDFORM = 5
      ELSEIF (PZ(I) .LT. 1000.) THEN
         INDFORM = 6
      ELSE
         INDFORM = 7
      ENDIF
      WRITE(IWR,OUTFORM(INDFORM)) I, PZ(I), PFUP(1,I+1), &
     &     DIFDOWNFLUX, DIRDOWNFLUX, &
     &     PFDOWN(1,I+1), PFLS(1,I+1), HTR(I)
 3000 CONTINUE
      WRITE(IWR,9903)PAGE
 3001 CONTINUE

   IF (IOUT .GE. 0 .AND. IOUT .LE. JPB2) GO TO 3500
   IF (IFLAG .EQ. 98) THEN
      IFLAG = JPB1
   ELSEIF (IFLAG .LT. JPB2) THEN
      IFLAG = IFLAG + 1
   ELSE
      GO TO 3500
   ENDIF
   GO TO 1000
3500 CONTINUE

! ***    Output module version numbers

!         WRITE(IWR,9910) HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,&
!     &                   HVRUTL,&
!     &                   HVRD1M,HVRR1M,HVREPK,HVRLPK,&
!     &                   HVRCLD,HVRDIS,HVREXT,&
!     &                   HVRKG16,HVRKG17,HVRKG18,HVRKG19,&
!     &                   HVRKG20,HVRKG21,HVRKG22,HVRKG23,HVRKG24,&
!     &                   HVRKG25,HVRKG27,HVRKG28,HVRKG29

   CLOSE(IWR)

4000 CONTINUE

 9879 format(1x)
 9880 format(1x,'All output fluxes have been adjusted to account for ins&
     &trumental cosine response.') 
 9881 format(1x,'The output diffuse fluxes have been adjusted to account&
     & for instrumental cosine response.') 

 9882 format(1x,'The downwelling direct and diffuse fluxes have been com&
     &puted using the delta-M scaling approximation.') 
 9883 format(1x)

 9899 FORMAT(1X,'Wavenumbers: ',F6.0,' - ',F6.0,' cm-1')
 9900 FORMAT(1X,'LEVEL PRESSURE   UPWARD FLUX   DIFDOWN FLUX  DIRDOWN FL&
     &UX  DOWNWARD FLUX   NET FLUX    HEATING RATE')
 9901 FORMAT(1X,'         mb          W/m2          W/m2          W/m2&
     &        W/m2          W/m2       degree/day')
 9902 FORMAT(1X,I3,3X,F11.6,4X,1P,2(G12.6,2X),G13.6,3X,G16.9,0P)
 9903 FORMAT(A)
 9910 FORMAT('  Modules and versions used in this calculation:',/,/,5X,&
     &       '       rrtm.f: ',4X,A15,10X,'    rtrdis.f: ',6X,A15,/,5X,&
     &       '     rrtatm.f: ',4X,A15,10X,'   setcoef.f: ',6X,A15,/,5X,&
     &       '  taumoldis.f: ',4X,A15,10X,'   util_xx.f: ',6X,A15,/,5X,&
     &       '     D1MACH.f: ',4X,A15,10X,'    R1MACH.f: ',6X,A15,/,5X,&
     &       '    ErrPack.f: ',4X,A15,10X,'    LINPAK.f: ',6X,A15,/,5X,&
     &       ' cldprop_sw.f: ',4X,A15,10X,'    disort.f: ',6X,A15,/,5X,&
     &       '      extra.f: ',4X,A15,10X,'    k_gB16.f: ',6X,A15,/,5X,&
     &       '     k_gB17.f: ',4X,A15,10X,'    k_gB18.f: ',6X,A15,/,5X,&
     &       '     k_gB19.f: ',4X,A15,10X,'    k_gB20.f: ',6X,A15,/,5X,&
     &       '     k_gB21.f: ',4X,A15,10X,'    k_gB22.f: ',6X,A15,/,5X,&
     &       '     k_gB23.f: ',4X,A15,10X,'    k_gB24.f: ',6X,A15,/,5X,&
     &       '     k_gB25.f: ',4X,A15,10X,'    k_gB27.f: ',6X,A15,/,5X,&
     &       '     k_gB28.f: ',4X,A15,10X,'    k_gB29.f: ',6X,A15,/)

END DO

END PROGRAM RRTM_SW_224

!*****************************************************************
FUNCTION EARTH_SUN(IDN)

!  Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  IDN        : Day of the year
!  EARTH_SUN  : square of the ratio of mean to actual Earth-Sun distance

PI   =  3.141592654
GAMMA = 2.*PI*(IDN-1)/365.

!  Use Iqbal's equation 1.2.1

EARTH_SUN = 1.000110 + .034221 * COS(GAMMA) + .001289 * SIN(GAMMA) + &
  &          .000719 * COS(2.*GAMMA) + .000077 * SIN(2.*GAMMA)

RETURN
END FUNCTION EARTH_SUN
!*****************************************************************
