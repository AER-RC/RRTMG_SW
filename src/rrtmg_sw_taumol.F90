!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

!******************************************************************************
!                                                                             *
!                  Optical depths developed for the                           *
!                                                                             *
!                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
!                                                                             *
!                                                                             *
!            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
!                        840 MEMORIAL DRIVE                                   *
!                        CAMBRIDGE, MA 02139                                  *
!                                                                             *
!                                                                             *
!                           ELI J. MLAWER                                     *
!                         STEVEN J. TAUBMAN                                   *
!                         SHEPARD A. CLOUGH                                   *
!                                                                             *
!                                                                             *
!                                                                             *
!                                                                             *
!                       email:  mlawer@aer.com                                *
!                                                                             *
!        The authors wish to acknowledge the contributions of the             *
!        following people:  Patrick D. Brown, Michael J. Iacono,              *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                       *
!                                                                             *
!******************************************************************************
!     TAUMOL                                                                  *
!                                                                             *
!     This file contains the subroutines TAUGBn (where n goes from            *
!     1 to 28).  TAUGBn calculates the optical depths and Planck fractions    *
!     per g-value and layer for band n.                                       *
!                                                                             *
!  Output:  optical depths (unitless)                                         *
!           fractions needed to compute Planck functions at every layer       *
!               and g-value                                                   *
!                                                                             *
!     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
!     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
!                                                                             *
!  Input                                                                      *
!                                                                             *
!     PARAMETER (MG=16, MXLAY=203, NBANDS=14)                                 *
!                                                                             *
!     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
!     COMMON /PRECISE/  ONEMINUS                                              *
!     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
!    &                  PZ(0:MXLAY),TZ(0:MXLAY),TBOUND                        *
!     COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,                              *
!    &                  COLH2O(MXLAY),COLCO2(MXLAY),                          *
!    &                  COLO3(MXLAY),COLN2O(MXLAY),COLCH4(MXLAY),             *
!    &                  COLO2(MXLAY),CO2MULT(MXLAY)                           *
!     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
!    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
!     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
!     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
!                                                                             *
!     Description:                                                            *
!     NG(IBAND) - number of g-values in band IBAND                            *
!     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
!                   atmospheres that are stored for band IBAND per            *
!                   pressure level and temperature.  Each of these            *
!                   atmospheres has different relative amounts of the         *
!                   key species for the band (i.e. different binary           *
!                   species parameters).                                      *
!     NSPB(IBAND) - same for upper atmosphere                                 *
!     ONEMINUS - since problems are caused in some cases by interpolation     *
!                parameters equal to or greater than 1, for these cases       *
!                these parameters are set to this value, slightly < 1.        *
!     PAVEL - layer pressures (mb)                                            *
!     TAVEL - layer temperatures (degrees K)                                  *
!     PZ - level pressures (mb)                                               *
!     TZ - level temperatures (degrees K)                                     *
!     LAYTROP - layer at which switch is made from one combination of         *
!               key species to another                                        *
!     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
!               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
!               respectively (molecules/cm**2)                                *
!     CO2MULT - for bands in which carbon dioxide is implemented as a         *
!               trace species, this is the factor used to multiply the        *
!               band's average CO2 absorption coefficient to get the added    *
!               contribution to the optical depth relative to 355 ppm.        *
!     FACij(LAY) - for layer LAY, these are factors that are needed to        *
!                  compute the interpolation factors that multiply the        *
!                  appropriate reference k-values.  A value of 0 (1) for      *
!                  i,j indicates that the corresponding factor multiplies     *
!                  reference k-value for the lower (higher) of the two        *
!                  appropriate temperatures, and altitudes, respectively.     *
!     JP - the index of the lower (in altitude) of the two appropriate        *
!          reference pressure levels needed for interpolation                 *
!     JT, JT1 - the indices of the lower of the two appropriate reference     *
!               temperatures needed for interpolation (for pressure           *
!               levels JP and JP+1, respectively)                             *
!     SELFFAC - scale factor needed to water vapor self-continuum, equals     *
!               (water vapor density)/(atmospheric density at 296K and        *
!               1013 mb)                                                      *
!     SELFFRAC - factor needed for temperature interpolation of reference     *
!                water vapor self-continuum data                              *
!     INDSELF - index of the lower of the two appropriate reference           *
!               temperatures needed for the self-continuum interpolation      *
!                                                                             *
!  Data input                                                                 *
!     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG) *
!        (note:  n is the band number)                                        *
!                                                                             *
!     Description:                                                            *
!     KA - k-values for low reference atmospheres (no water vapor             *
!          self-continuum) (units: cm**2/molecule)                            *
!     KB - k-values for high reference atmospheres (all sources)              *
!          (units: cm**2/molecule)                                            *
!     SELFREF - k-values for water vapor self-continuum for reference         *
!               atmospheres (used below LAYTROP)                              *
!               (units: cm**2/molecule)                                       *
!                                                                             *
!     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
!     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
!                                                                             *
!******************************************************************************

SUBROUTINE TAUMOL16 &
  &( KLEV    &
  &, FAC00   , FAC01   , FAC10   , FAC11   &
  &, JP      , JT      , JT1     , ONEMINUS&
  &, COLH2O  , COLCH4  , COLMOL  &
  &, LAYTROP , SELFFAC , SELFFRAC, INDSELF , FORFAC  , FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG    , TAUR    &
  &)

!     BAND 16:  2600-3250 cm-1 (low - H2O,CH4; high - CH4)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG16
USE YOESRTA16, ONLY : ABSA, ABSB, FORREFC, SELFREFC &
  &                 , SFLUXREFC, RAYL &
  &                 , LAYREFFR, STRRAT1
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE

!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)

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
     &           (SELFFAC(LAY) * (SELFREFC(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

LAYSOLFR = NLAYERS

DO LAY = LAYTROP+1, NLAYERS
  IF (JP(LAY-1) .LT. LAYREFFR .AND. JP(LAY) .GE. LAYREFFR) &
   &        LAYSOLFR = LAY
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(16) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(16) + 1
  TAURAY = COLMOL(LAY) * RAYL

  DO IG = 1, NG16
    TAUG(LAY,IG) = COLCH4(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0  ,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1  ,IG) + &
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG) 
    TAUR(LAY,IG) = TAURAY  
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL16
!----------------------------------------------------------------------

SUBROUTINE TAUMOL17 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCO2 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment

#include "tsmbkind.h"


USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG17
USE YOESRTA17, ONLY : ABSA, ABSB, FORREFC, SELFREFC &
  &                 , SFLUXREFC, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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
     &           (SELFFAC(LAY) * (SELFREFC(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))) 
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
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG))) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG,JS) &
     &           + FS * (SFLUXREFC(IG,JS+1) - SFLUXREFC(IG,JS))
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL17
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL18 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCH4 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG18
USE YOESRTA18, ONLY : ABSA, ABSB, FORREFC, SELFREFC &
  &                 , SFLUXREFC, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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
     &           (SELFFAC(LAY) * (SELFREFC(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG,JS)  &
     &           + FS * (SFLUXREFC(IG,JS+1) - SFLUXREFC(IG,JS))
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(18) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(18) + 1
  TAURAY = COLMOL(LAY) * RAYL

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

RETURN
END SUBROUTINE TAUMOL18
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL19 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCO2 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 19:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG19
USE YOESRTA19, ONLY : ABSA, ABSB, FORREFC, SELFREFC &
  &                 , SFLUXREFC, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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
     &           (SELFFAC(LAY) * (SELFREFC(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + & 
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG,JS) &
     &           + FS * (SFLUXREFC(IG,JS+1) - SFLUXREFC(IG,JS))
    TAUR(LAY,IG) = TAURAY   
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(19) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(19) + 1
  TAURAY = COLMOL(LAY) * RAYL

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

RETURN
END SUBROUTINE TAUMOL19
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL20 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCH4 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG20
USE YOESRTA20, ONLY : ABSA, ABSB, FORREFC, SELFREFC &
  &                 , SFLUXREFC, ABSCH4C, RAYL &
  &                 , LAYREFFR
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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

  DO IG = 1 , NG20
    TAUG(LAY,IG) = COLH2O(LAY) * &
     &          ((FAC00(LAY) * ABSA(IND0,IG) + &
     &           FAC10(LAY) * ABSA(IND0+1,IG) + &
     &           FAC01(LAY) * ABSA(IND1,IG) + &
     &           FAC11(LAY) * ABSA(IND1+1,IG)) + &
     &           SELFFAC(LAY) * (SELFREFC(INDS,IG) + & 
     &           SELFFRAC(LAY) * &
     &           (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))) &
     &           + COLCH4(LAY) * ABSCH4C(IG)
!     &           + TAURAY &
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY 
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG) 
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(20) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(20) + 1
  INDF = INDFOR(LAY)
  TAURAY = COLMOL(LAY) * RAYL

  DO IG = 1 , NG20
    TAUG(LAY,IG) = COLH2O(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + &
     &           FAC11(LAY) * ABSB(IND1+1,IG) + &
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))) + &
     &           COLCH4(LAY) * ABSCH4C(IG)
!     &           TAURAY + &
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY 
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL20
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL21 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCO2 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG21
USE YOESRTA21, ONLY : ABSA, ABSB, FORREFC, SELFREFC &
  &                 , SFLUXREFC, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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
     &           (SELFFAC(LAY) * (SELFREFC(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG))))
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG,JS) &
     &           + FS * (SFLUXREFC(IG,JS+1) - SFLUXREFC(IG,JS))
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
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))
!     &           + TAURAY 
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL21
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL22 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLMOL , COLO2   &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 22:  7700-8050 cm-1 (low - H2O,O2; high - O2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG22
USE YOESRTA22, ONLY : ABSA, ABSB, FORREFC, SELFREFC &
  &                 , SFLUXREFC, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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
     &           (SELFFAC(LAY) * (SELFREFC(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))) &
     &           + O2CONT
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG,JS) &
     &           + FS * (SFLUXREFC(IG,JS+1) - SFLUXREFC(IG,JS))
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  O2CONT = 4.35e-4*colo2(lay)/(350.0*2.0)
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(22) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(22) + 1
  TAURAY = COLMOL(LAY) * RAYL

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

RETURN
END SUBROUTINE TAUMOL22
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL23 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLMOL &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG23
USE YOESRTA23, ONLY : ABSA, FORREFC, SELFREFC &
  &                 , SFLUXREFC, RAYLC &
  &                 , LAYREFFR, GIVFAC 
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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

  DO IG = 1 , NG23
    TAURAY = COLMOL(LAY) * RAYLC(IG)
    TAUG(LAY,IG) = COLH2O(LAY) * &
     &          (GIVFAC * (FAC00(LAY) * ABSA(IND0,IG) + &
     &           FAC10(LAY) * ABSA(IND0+1,IG) + &
     &           FAC01(LAY) * ABSA(IND1,IG) + &
     &           FAC11(LAY) * ABSA(IND1+1,IG)) + &
     &           SELFFAC(LAY) * (SELFREFC(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + &
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))) 
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG) 
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  DO IG = 1 , NG23
!    TAUG(LAY,IG) = COLMOL(LAY) * RAYLC(IG)
!    SSA(LAY,IG) = 1.0
    TAUG(LAY,IG) = _ZERO_
    TAUR(LAY,IG) = COLMOL(LAY) * RAYLC(IG) 
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL23
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL24 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLMOL , COLO2   , COLO3    &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 24:  12850-16000 cm-1 (low - H2O,O2; high - O2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG24
USE YOESRTA24, ONLY : ABSA, ABSB, FORREFC, SELFREFC &
  &                 , SFLUXREFC, ABSO3AC, ABSO3BC, RAYLAC, RAYLBC &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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

  DO IG = 1 , NG24
    TAURAY = COLMOL(LAY) * (RAYLAC(IG,JS) + &
     &           FS * (RAYLAC(IG,JS+1) - RAYLAC(IG,JS)))
    TAUG(LAY,IG) = SPECCOMB * &
     &          (FAC000 * ABSA(IND0,IG) + &
     &           FAC100 * ABSA(IND0+1,IG) + &
     &           FAC010 * ABSA(IND0+9,IG) + &
     &           FAC110 * ABSA(IND0+10,IG) + &
     &           FAC001 * ABSA(IND1,IG) + &
     &           FAC101 * ABSA(IND1+1,IG) + &
     &           FAC011 * ABSA(IND1+9,IG) + &
     &           FAC111 * ABSA(IND1+10,IG)) + &
     &           COLO3(LAY) * ABSO3AC(IG) + &
     &           COLH2O(LAY) * & 
     &           (SELFFAC(LAY) * (SELFREFC(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + & 
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG))))
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG,JS) &
     &           + FS * (SFLUXREFC(IG,JS+1) - SFLUXREFC(IG,JS))
    TAUR(LAY,IG) =  TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(24) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(24) + 1

  DO IG = 1 , NG24
    TAURAY = COLMOL(LAY) * RAYLBC(IG)
    TAUG(LAY,IG) = COLO2(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + &
     &           FAC11(LAY) * ABSB(IND1+1,IG)) + &
     &           COLO3(LAY) * ABSO3BC(IG)
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL24
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL25 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLMOL , COLO3   &
  &, LAYTROP &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 25:  16000-22650 cm-1 (low - H2O; high - nothing)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG25
USE YOESRTA25, ONLY : ABSA &
  &                 , SFLUXREFC, ABSO3AC, ABSO3BC, RAYLC &
  &                 , LAYREFFR
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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

  DO IG = 1 , NG25
    TAURAY = COLMOL(LAY) * RAYLC(IG)
    TAUG(LAY,IG) = COLH2O(LAY) * &
     &          (FAC00(LAY) * ABSA(IND0,IG) + &
     &           FAC10(LAY) * ABSA(IND0+1,IG) + &
     &           FAC01(LAY) * ABSA(IND1,IG) + &
     &           FAC11(LAY) * ABSA(IND1+1,IG)) + &
     &           COLO3(LAY) * ABSO3AC(IG) 
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG) 
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  DO IG = 1 , NG25
    TAURAY = COLMOL(LAY) * RAYLC(IG)
    TAUG(LAY,IG) = COLO3(LAY) * ABSO3BC(IG) 
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL25
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL26 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCO2 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG26
USE YOESRTA26, ONLY : SFLUXREFC, RAYLC
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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
  DO IG = 1 , NG26 
!    TAUG(LAY,IG) = COLMOL(LAY) * RAYLC(IG)
!    SSA(LAY,IG) = 1.0
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG) 
    TAUG(LAY,IG) = _ZERO_
    TAUR(LAY,IG) = COLMOL(LAY) * RAYLC(IG) 
  END DO
END DO

DO LAY = LAYTROP+1, NLAYERS
  DO IG = 1 , NG26
!    TAUG(LAY,IG) = COLMOL(LAY) * RAYLC(IG)
!    SSA(LAY,IG) = 1.0
    TAUG(LAY,IG) = _ZERO_
    TAUR(LAY,IG) = COLMOL(LAY) * RAYLC(IG) 
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL26
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL27 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLMOL  , COLO3  &
  &, LAYTROP &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 27:  29000-38000 cm-1 (low - O3; high - O3)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG27
USE YOESRTA27, ONLY : ABSA, ABSB &
  &                 , SFLUXREFC, RAYLC &
  &                 , LAYREFFR, SCALEKUR
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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

  DO IG = 1 , NG27
    TAURAY = COLMOL(LAY) * RAYLC(IG)
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

  DO IG = 1 , NG27
    TAURAY = COLMOL(LAY) * RAYLC(IG)
    TAUG(LAY,IG) = COLO3(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + & 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
!     &          + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY.EQ.LAYSOLFR) SFLUXZEN(IG) = SCALEKUR * SFLUXREFC(IG) 
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL27
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL28 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    &
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLMOL  , COLO2  , COLO3   &
  &, LAYTROP &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 28:  38000-50000 cm-1 (low - O3,O2; high - O3,O2)

! Modifications
!
!     JJMorcrette 2003-02-24 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG28
USE YOESRTA28, ONLY : ABSA, ABSB &
  &                 , SFLUXREFC, RAYL &
  &                 , LAYREFFR, STRRAT
USE YOESRTWN , ONLY : NG, NSPA, NSPB

IMPLICIT NONE


!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)


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
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG,JS) &
     &           + FS * (SFLUXREFC(IG,JS+1) - SFLUXREFC(IG,JS))
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL28
!-----------------------------------------------------------------------

SUBROUTINE TAUMOL29 &
  &( KLEV    &
  &, FAC00   , FAC01  , FAC10   , FAC11    & 
  &, JP      , JT     , JT1     , ONEMINUS &
  &, COLH2O  , COLCO2 , COLMOL  &
  &, LAYTROP , SELFFAC, SELFFRAC, INDSELF  , FORFAC, FORFRAC, INDFOR &
  &, SFLUXZEN, TAUG   , TAUR    &
  &)

!     BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)

! Modifications
!
!     JJMorcrette 2002-10-03 adapted to ECMWF environment


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPBAND, JPG, NG29
USE YOESRTA29, ONLY : ABSA, ABSB, FORREFC, SELFREFC &
  &                 , SFLUXREFC, ABSH2OC, ABSCO2C, RAYL &
  &                 , LAYREFFR
USE YOESRTWN , ONLY : NG, NSPA, NSPB

  
IMPLICIT NONE

!-- Output
REAL_B :: TAUG(JPLAY,JPG), TAUR(JPLAY,JPG), SSA(JPLAY,JPG), SFLUXZEN(JPG)

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

  DO IG = 1, NG29
    TAUG(LAY,IG) = COLH2O(LAY) * &
     &          ((FAC00(LAY) * ABSA(IND0,IG) + &
     &           FAC10(LAY) * ABSA(IND0+1,IG) + &
     &           FAC01(LAY) * ABSA(IND1,IG) + &
     &           FAC11(LAY) * ABSA(IND1+1,IG)) + &
     &           SELFFAC(LAY) * (SELFREFC(INDS,IG) + &
     &           SELFFRAC(LAY) * &
     &           (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
     &           FORFAC(LAY) * (FORREFC(INDF,IG) + & 
     &           FORFRAC(LAY) * &
     &           (FORREFC(INDF+1,IG) - FORREFC(INDF,IG)))) &
     &           + COLCO2(LAY) * ABSCO2C(IG) 
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

  DO IG = 1 , NG29
    TAUG(LAY,IG) = COLCO2(LAY) * &
     &          (FAC00(LAY) * ABSB(IND0,IG) + &
     &           FAC10(LAY) * ABSB(IND0+1,IG) + &
     &           FAC01(LAY) * ABSB(IND1,IG) + &
     &           FAC11(LAY) * ABSB(IND1+1,IG)) &  
     &           + COLH2O(LAY) * ABSH2OC(IG) 
!     &           + TAURAY
!    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
    IF (LAY .EQ. LAYSOLFR) SFLUXZEN(IG) = SFLUXREFC(IG) 
    TAUR(LAY,IG) = TAURAY
  END DO
END DO

RETURN
END SUBROUTINE TAUMOL29



