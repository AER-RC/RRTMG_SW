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

SUBROUTINE RRTM_SW_224 &
  &( PTS   , PTH    , PT  &
  &, PQ    , PCCO2  , POZON , PCH4  , PN2O  &
  &, PAPH  , PAP    , PAER  , PALBD , PALBP &
  &, PRMU0 , PCLFR  , KOVLP , DYOFYR &
  &, INFLAG, ICEFLAG, LIQFLAG &
  &, CLDD1SW,CLDD2SW, CLDD3SW,CLDD4SW,CLDSWMOM &
  &, PCUP  , PCDOWN , PHEAC &
  &, PFUP  , PFDOWN , PHEAT &
  &)

!-- Interface to RRTM_SW; conversion to F90 formatting; addition of 
!   2-stream radiative transfer
!     J.-J. Morcrette, ECMWF, 030225
!-- Additional modifications for GCM application
!     M. J. Iacono, AER Inc., August 2003

!- Input from calling routine:
!- Note: In RRTM_SW, the level dimension goes from bottom to top
!  KLON = number of longitudes (1 in single-column mode)
!  JPLAY = maximum number of layers
!  KLEV, NLAYERS = number of layers
!  NSTR = Number of streams for radiative transfer (= 2)
!  KOVLP = cloud overlap method (1=Maximum-random, 2=Maximum, or 3=Random)
!  DYOFYR = day of year (1 to 365; set to 0 for calculation at 1 AU)
!  PRMU0 = Cosine of solar zenith angle
!  ADJFLUX = Incoming solar flux adjustment for current Earth/Sun distance
!  PAER = Aerosol optical thickness (six aerosol types from ECMWF model)
!  PALBD = sw surface albedo for diffuse radiation
!  PALBP = sw surface albdeo for parallel radiation (usually = PALBD)
!  PAPH = pressure levels (Pa)
!  PAP = pressure layers (Pa)
!  PTS = surface temp (K)
!  PTH = temperature levels (K)
!  PT = temperature layers (K)
!- Water vapor, CO2, ozone, CH4, and N2O are needed in units of vmr
!- If passed in as mmr, use molecular weights below to convert to vmr
!  PQ = water vapor 
!  PCCO2 = CO2 
!  POZON = ozone
!  PCH4 = methane
!  PN2O = nitrous ozide
!  PCLFR = cloud fraction

!-Variables related to the cloud optical properties (original RRTM_SW)
!- (See susrtop.F90 for further description):
!  INFLAG = Flag for cloud optical properties
!  ICEFLAG = Flag for ice particle specification
!  LIQFLAG = Flag for liquid droplet specification
!-For INFLAG=0/LIQFLAG=0/ICEFLAG=0:
!  CLDD1SW = Total (ice and water) cloud optical depth for layer
!  CLDD2SW = Singe-scattering albedo for layer
!  CLDSWMOM= Moments of the phase function (from 0 to NSTR)
!-For INFLAG=2/LIQFLAG=1/ICEFLAG=3:
!  CLDD1SW = Cloud water path for the layer (g/m2)
!  CLDD2SW = Fraction of the cloud layer's water path in the form 
!            of ice particles
!  CLDD3SW = Generalized effective size of the ice crystals, dge (microns)
!            (see Q. Fu, 1996 for a definition of this parameter).
!            Valid sizes are 5.0 to 140.0 microns.
!  CLDD4SW = Liquid droplet effective radius, re, (microns)
!            of ice particles
!            Valid sizes are 2.5 to 60.0 microns.

!- Return to calling routine:
!  PFUP = total sky flux (W/m2)
!  PFDOWN = total sky flux (W/m2)
!  PCUP = clear sky up flux (W/m2)
!  PCDOWN = clear sky down flux (W/m2)
!  PHEAT = total sky heating rate (K/day)
!  PHEAC = clear sky heating rate (K/day)


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLON, JPLAY, JPAER, JPSW, JPBAND, JPB1, JPB2
USE YOESRTAER, ONLY : RSRTAUA, RSRPIZA, RSRASYA 

IMPLICIT NONE

!-- Input arguments

INTEGER_M :: KSW, KOVLP, DYOFYR

REAL_B :: PAER(JPLON,JPAER,JPLAY)
REAL_B :: PALBD(JPLON,JPSW)   , PALBP(JPLON,JPSW)
REAL_B :: PAP(JPLON,JPLAY)     , PAPH(JPLON,JPLAY+1) , PDP(JPLON,JPLAY)
REAL_B :: CLDD1SW(JPLON,JPLAY) , CLDD2SW(JPLON,JPLAY), CLDD3SW(JPLON,JPLAY)
REAL_B :: CLDD4SW(JPLON,JPLAY) , CLDSWMOM(JPLON,0:2,JPLAY)
REAL_B :: PTS(JPLON)          , PT(JPLON,JPLAY)     , PTH(JPLON,JPLAY+1)
REAL_B :: PQ(JPLON,JPLAY)      , POZON(JPLON,JPLAY)  , PCLFR(JPLON,JPLAY)
REAL_B :: PCCO2, PCH4, PN2O

!-- Output arguments

REAL_B :: PFCS(JPLON,JPLAY+1), PFLS(JPLON,JPLAY+1)
REAL_B :: PHEAC(JPLON,JPLAY), PHEAT(JPLON,JPLAY)
REAL_B :: PFDOWN(JPLON,JPLAY+1), PCDOWN(JPLON,JPLAY+1)
REAL_B :: PFUP(JPLON,JPLAY+1), PCUP(JPLON,JPLAY+1)
 
!-----------------------------------------------------------------------

!-- dummy integers

INTEGER_M :: ICLDATM, INFLAG, ICEFLAG, LIQFLAG, NMOL, NSTR
INTEGER_M :: IK, IMOL, J1, J2, JA, JAE, JL, JK, JMOM, JSW

INTEGER_M :: LAYTROP, LAYSWTCH, LAYLOW
INTEGER_M :: INDFOR(JPLAY), INDSELF(JPLAY)
INTEGER_M :: JP(JPLAY), JT(JPLAY), JT1(JPLAY)

!-- dummy reals

REAL_B :: CLDFRAC(JPLAY), CLDDAT1(JPLAY), CLDDAT2(JPLAY), CLDDAT3(JPLAY), CLDDAT4(JPLAY)
REAL_B :: CLDDATMOM(0:16,JPLAY)
REAL_B :: TAUCLDORIG(JPLAY), TAUCLOUD(JPLAY,JPBAND), SSACLOUD(JPLAY,JPBAND)
REAL_B :: XMOM(0:16,JPLAY,JPBAND)
REAL_B :: PRMU0(JPLON)
REAL_B :: ADJFLUX, EARTH_SUN

REAL_B :: PZ(0:JPLAY)   , TZ(0:JPLAY)   , PAVEL(JPLAY)  , TAVEL(JPLAY)
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
REAL_B :: ZUVCD(JPLAY+1), ZUVCU(JPLAY+1), ZUVFD(JPLAY+1), ZUVFU(JPLAY+1)
REAL_B :: ZVSCD(JPLAY+1), ZVSCU(JPLAY+1), ZVSFD(JPLAY+1), ZVSFU(JPLAY+1)
REAL_B :: ZNICD(JPLAY+1), ZNICU(JPLAY+1), ZNIFD(JPLAY+1), ZNIFU(JPLAY+1)

REAL_B :: amd                  ! Effective molecular weight of dry air (g/mol)
REAL_B :: amw                  ! Molecular weight of water vapor (g/mol)
REAL_B :: amco2                ! Molecular weight of carbon dioxide (g/mol)
REAL_B :: amo                  ! Molecular weight of ozone (g/mol)
REAL_B :: amo2                 ! Molecular weight of oxygen (g/mol)
REAL_B :: amch4                ! Molecular weight of methane (g/mol)
REAL_B :: amn2o                ! Molecular weight of nitrous oxide (g/mol)
REAL_B :: amc11                ! Molecular weight of CFC11 (g/mol) - CFCL3
REAL_B :: amc12                ! Molecular weight of CFC12 (g/mol) - CF2CL2
REAL_B :: avgdro               ! Avogadro's number (molecules/mole)
REAL_B :: gravit               ! Gravitational acceleration (cm/sec2)
REAL_B :: amm                  ! Molecular weight of moist air
REAL_B :: o2mmr                ! Oxygen mass mixing ratio

! Atomic weights for conversion from mass to volume mixing ratios; these
!  are the same values used in ECRT to assure accurate conversion to vmr
data amd   /  28.970_JPRB    /
data amw   /  18.0154_JPRB   /
data amco2 /  44.011_JPRB    /
data amo   /  47.9982_JPRB   /
data amo2  /  31.9999_JPRB   /
data amch4 /  16.043_JPRB    /
data amn2o /  44.013_JPRB    /
data amc11 / 137.3686_JPRB   /
data amc12 / 120.9140_JPRB   /
data avgdro/ 6.02214E23_JPRB /
data gravit/ 9.8066E02_JPRB  /

! Oxygen mass mixing ratio 
data o2mmr /  0.23143_JPRB   /

REAL_B :: ZCLEAR, ZCLOUD, ZEPSEC, ZTOTCC, ZDPGCP, ZEPZEN

REAL_B :: RDAY, RG, RMD, RKBOL, RNAVO, R, RD, RCPD, RCDAY

!-----------------------------------------------------------------------
!-- calculate information needed by the radiative transfer routine 

ZEPSEC  = 1.E-06_JPRB
ZEPZEN  = 1.E-10_JPRB
ONEMINUS=_ONE_ -  ZEPSEC

NSTR	= 2
NMOL	= 6
KSW     = JPSW
KOVLP   = 1
KLON    = JPLON
KLEV    = JPLAY

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
RD   = 1.E5*R/RMD
RCPD = 3.5*RD
!RCDAY= RDAY * RG / RCPD
RCDAY= 8.4391

! mji - Earth/Sun distance
! Input DYOFYR, the cumulative day of the year, or
! set to 0 for calculation at 1 AU
! DYOFYR = 0
IF (DYOFYR .EQ. 0) THEN
   ADJFLUX = 1.
ELSE
   ADJFLUX = EARTH_SUN(DYOFYR)
ENDIF
!print*, 'RRTM_SW; DYOFYR, ADJFLUX: ', DYOFYR, ADJFLUX

! Main longitude loop
DO JL = 1, KLON

   DO JK=1,KLEV 
      CLDFRAC(JK) = PCLFR (JL,JK)
      CLDDAT1(JK) = CLDD1SW(JL,JK)
      CLDDAT2(JK) = CLDD2SW(JL,JK)
      CLDDAT3(JK) = CLDD3SW(JL,JK)
      CLDDAT4(JK) = CLDD4SW(JL,JK)
      DO JMOM=0,NSTR
         CLDDATMOM(JMOM,JK)=CLDSWMOM(JL,JMOM,JK)
      END DO

!    print 9100,INFLAG,ICEFLAG,LIQFLAG
!9100 format(1x,'RRTM_SW_224 Flags :',I3)
!    print 9101,JK,CLDFRAC(JK),CLDDAT1(JK),CLDDAT2(JK),CLDDAT3(JK)&
!    &,CLDDAT4(JK),(CLDDATMOM(JMOM,JK),JMOM=0,NSTR)
!9101 format(1x,'RRTM_SW_224 Cld :',I3,F7.4,1P7E12.5)
   END DO


   CALL RRTM_SW_CLDPROP &
     &( KLEV, ICLDATM, INFLAG, ICEFLAG, LIQFLAG, NSTR &
     &, CLDFRAC, CLDDAT1, CLDDAT2, CLDDAT3, CLDDAT4, CLDDATMOM &
     &, TAUCLDORIG, TAUCLOUD, SSACLOUD, XMOM &
     &)

! mji out
! DO JSW = 16, 29
!  DO JK = 1, KLEV
!    print 9111,JSW,JK,TAUCLOUD(JK,JSW),SSACLOUD(JK,JSW),(XMOM(JMOM,JK,JSW),JMOM=0,NSTR)
!9111 format(1x,'RRTM_SW_224 Cldprop out :',2I3,1P7E12.5)
!  END DO
! END DO

!  DO JK = 1, KLEV+1
!    print 9199,JK,PAPH(JL,JK),PTH(JL,JK)
!  END DO
!  DO JK = 1, KLEV
!    print 9199,JK,PAP(JL,JK),PT(JL,JK),PQ(JL,JK),PCCO2,POZON(JL,JK),PN2O,PCH4
!  END DO
!9199 format(1x,'RRTM_SW_IN ',I3,2F9.1,1P6E13.5)

!- coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients

   DO J1=1,35
      DO J2=1,KLEV
         WKL(J1,J2)=_ZERO_ 
      ENDDO
   ENDDO

! Pass input arrays into RRTM arrays.
! Convert pressure from Pa to mb, and convert molecular amounts
! from mmr to vmr.
   TBOUND=PTS(JL)
   PZ(0) = paph(JL,1)*1.e-2
   TZ(0) = pth (JL,1)

   ZCLEAR=_ONE_
   ZCLOUD=_ZERO_
   ZTOTCC=_ZERO_
   DO JK = 1, KLEV
      PAVEL(JK) = pap(JL,JK)*1.e-2
      TAVEL(JK) = pt(JL,JK)
      PZ(JK)    = paph(JL,JK+1)*1.e-2
      TZ(JK)    = pth(JL,JK+1)
      PDP(JL,JK)= (paph(JL,JK)-paph(JL,JK+1))*1.e-2
      WKL(1,JK) = pq(JL,JK)*amd/amw
      WKL(2,JK) = pcco2*amd/amco2
      WKL(3,JK) = pozon(JL,JK)*amd/amo
      WKL(4,JK) = pn2o*amd/amn2o
      WKL(6,JK) = pch4*amd/amch4
      WKL(7,JK) = o2mmr*amd/amo2
      amm = (1-WKL(1,JK))*amd + WKL(1,JK)*amw
      COLDRY(JK) = (PZ(JK-1)-PZ(JK))*1.E3_JPRB*avgdro/(gravit*amm*(1+WKL(1,JK)))

!    print 9200,JK,PAVEL(JK),TAVEL(JK),(WKL(JA,JK),JA=1,4),WKL(6,JK),COLDRY(JK)
!    print 9200,JK,PZ(JK),TZ(JK),PQ(JL,JK),POZON(JL,JK)
!9200 format(1x,'RRTM_SW ',I3,2F9.1,1P6E13.5)
!    print*, 'PRMU0 = ', PRMU0(JL)

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

   END DO

! Convert molecular amount from vmr to molecules/cm2
   DO IMOL=1,NMOL
      DO JK=1,KLEV
         WKL(IMOL,JK)=COLDRY(JK)* WKL(IMOL,JK)
      END DO
   END DO

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
  
!  print *,'just after RRTM_SW_SETCOEF'

!- Call the radiation transfer routine
  
!- Put cloud optical properties in arrays for 2-stream

   DO JSW=1,KSW
      ZALBD(JSW)=PALBD(JL,JSW)
      ZALBP(JSW)=PALBP(JL,JSW)
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
!9002  format(1x,'RRTM_SW_EC224gp CldpropRRTM_SW ',2I3,f8.4,1P3E12.5)
      END DO
   END DO

!- mixing of aerosols
!- Put aerosol optical properties in arrays for 2-stream

!  DO JK=1,KLEV
!    print 9013,JK,(PAER(JL,JAE,JK),JAE=1,JPAER)
!9013 format(1x,'Aerosol opt.thickness: ',I3,1P6E12.5)
!  END DO
!  DO JSW=1,KSW
!    print 9003,JSW,(RSRTAUA(JSW,JAE),JAE=1,JPAER)
!    print 9003,JSW,(RSRPIZA(JSW,JAE),JAE=1,JPAER)
!    print 9003,JSW,(RSRASYA(JSW,JAE),JAE=1,JPAER)
!9003 format(1x,'Aerosol opt.prop: ',I3,6F10.5)
!  END DO
 
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

   ZRMU0=PRMU0(JL)
   IF (ZRMU0.EQ.0.) ZRMU0=ZEPZEN

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
     &( KLEV   , NMOL    , KSW    ,ONEMINUS, JPB1    , JPB2 &
     &, PAVEL  , TAVEL   , PZ     , TZ     , TBOUND  , ZALBD   , ZALBP &
     &, ZCLFR  , ZTAUC   , ZASYC  , ZOMGC  , ZTAUA   , ZASYA   , ZOMGA , ZRMU0   &
     &, COLDRY , WKL     , ADJFLUX &	 
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

!  print *,'Clear-sky and total fluxes over whole SW spectrum'
!  print *,'    ClearUpw   ClearDnw  TotalUpw  TotalDnw'

   DO JK=1,KLEV+1
      PCUP(JL,JK)=ZBBCU(JK)
      PCDOWN(JL,JK)=ZBBCD(JK)
      PFUP(JL,JK)=(_ONE_-ZCLEAR)*ZBBFU(JK)+ZCLEAR*ZBBCU(JK)
      PFDOWN(JL,JK)=(_ONE_-ZCLEAR)*ZBBFD(JK)+ZCLEAR*ZBBCD(JK)
   END DO
!  DO JK=1,KLEV+1
!    print 9005,JK,PCUP(JL,JK),PCDOWN(JL,JK),PFUP(JL,JK),PFDOWN(JL,JK)
!9005 format(1x,I3,4F10.3)
!  END DO

   DO JK = 1 , KLEV+1
      PFLS(JL,JK) = PFDOWN(JL,JK) - PFUP(JL,JK)
      PFCS(JL,JK) = PCDOWN(JL,JK) - PCUP(JL,JK)
!  print 9504,JK,ZFLUC(JL,2,JK),ZFLUC(JL,1,JK),PFCS(JL,JK) &
!    &,          ZFLUX(JL,2,JK),ZFLUX(JL,1,JK),PFLS(JL,JK)
!9504 format(1x,'SRTM-FLX',1I3,2(2x,f10.3,1x,f10.3,1x,f10.3))
   ENDDO

   DO JK=1,KLEV
      ZDPGCP=RCDAY/PDP(JL,JK)
      PHEAC(JL,JK)=(PFCS(JL,JK+1)-PFCS(JL,JK))*ZDPGCP
      PHEAT(JL,JK)=(PFLS(JL,JK+1)-PFLS(JL,JK))*ZDPGCP
!  print 9505,JK,ZCEAT(JL,JK),ZHEAT(JL,JK),PCLFR(JL,JK)
!9505 format(1x,'SRTM-HR',I3,3(10x,f10.3,14x))
   END DO

END DO

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_SW_224

!*****************************************************************
REAL*8 FUNCTION EARTH_SUN(IDN)

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

