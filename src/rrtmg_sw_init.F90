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

SUBROUTINE RRTMG_SW_INIT

#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPG, JPSW, JPGPT
USE YOESRTWN , ONLY : NG, NSPA, NSPB, &
                    & NGM, WT, NGC, NGS, NGN, NGB, RWGT, WTSM

IMPLICIT NONE

! Local variables
INTEGER_B :: IGC, IGCSM, IBND, IG, IND, IPR, IPRSM
REAL_B :: WTSUM


!-- Read in the basic coefficients to configure RRTMG_SW.
!-- Creates module YOESRTWN with BG, NSPA, NSPB, WAVENUM1, WAVENUM2,
!   DELWAVE, PREF, PREFLOG, TREF

CALL SUSRTM


!-- Read in the molecular absorption coefficients
!-- Band 16 module renamed to avoid conflict with RRTM_LW band 16.

CALL KGB16S
CALL KGB17
CALL KGB18
CALL KGB19
CALL KGB20
CALL KGB21
CALL KGB22
CALL KGB23
CALL KGB24
CALL KGB25
CALL KGB26
CALL KGB27
CALL KGB28
CALL KGB29


!-- Perform g-point reduction from 16 per band (224 total points) to
!-- a band dependent number (112 total points) for all absorption
!-- coefficient input data and Planck fraction input data.
!-- Compute relative weighting for new g-point combinations.

      IGCSM = 0
      DO 500 IBND = 1,JPSW
         IPRSM = 0
         IF (NGC(IBND).LT.JPG) THEN
            DO 450 IGC = 1,NGC(IBND)
               IGCSM = IGCSM + 1
               WTSUM = 0.
               DO 420 IPR = 1, NGN(IGCSM)
                  IPRSM = IPRSM + 1
                  WTSUM = WTSUM + WT(IPRSM)
 420           CONTINUE
               WTSM(IGC) = WTSUM
 450        CONTINUE
            DO 400 IG = 1,NG(IBND+15)
               IND = (IBND-1)*JPG + IG
               RWGT(IND) = WT(IG)/WTSM(NGM(IND))
 400        CONTINUE
         ELSE
            DO 300 IG = 1,NG(IBND+15)
               IGCSM = IGCSM + 1
               IND = (IBND-1)*JPG + IG
               RWGT(IND) = 1.0
 300        CONTINUE
         ENDIF
 500  CONTINUE

CALL CMBGB16S
CALL CMBGB17
CALL CMBGB18
CALL CMBGB19
CALL CMBGB20
CALL CMBGB21
CALL CMBGB22
CALL CMBGB23
CALL CMBGB24
CALL CMBGB25
CALL CMBGB26
CALL CMBGB27
CALL CMBGB28
CALL CMBGB29


!-- Read in the cloud optical properties for RRTMG_SW_CLDPROP.
!-- Creates module YOESRTOP with EXTLIQ1, SSALIQ1, ASYLIQ1, 
!-- EXTICE3, SSAICE3, ASYICE3, FDLICE3  

CALL SUSRTOP

!-- Read in the aerosol optical properties for six aerosol
!-- types derived from ECMWF SW model.
!-- Creates module YOESRTAER with RSRTAUA, RSRPIZA, RSRASYA.
!-- The six aerosol types are respectively:
!--  1/ continental average                 2/ maritime
!--  3/ desert                              4/ urban
!--  5/ volcanic active                     6/ stratospheric background

CALL SUSRTAER

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE RRTMG_SW_INIT


!-----------------------------------------------------------------------
SUBROUTINE CMBGB16S

!
!  Original version:       Michael J. Iacono; July, 1998
!  Revision for RRTM_SW:   Michael J. Iacono; November, 2002
!  Revision for RRTMG_SW:  Michael J. Iacono; December, 2003
!
!  The subroutines CMBGB16->CMBGB29 input the absorption coefficient
!  data for each band, which are defined for 16 g-points and 14 spectral
!  bands. The data are combined with appropriate weighting following the
!  g-point mapping arrays specified in RRTMG_SW_INIT.  Solar source 
!  function data in array SFLUXREF are combined without weighting.  All
!  g-point reduced data are put into new arrays for use in RRTMG_SW.
!
!  BAND 16:  2600-3250 cm-1 (low key- H2O,CH4; high key - CH4)
!
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"

USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA16, ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, &
                    & KAC, KBC, SELFREFC, FORREFC, SFLUXREFC

IMPLICIT NONE

! Local variables
INTEGER_B :: JN, JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF

      DO 2000 JN = 1,9
         DO 2000 JT = 1,5
            DO 2200 JP = 1,13
               IPRSM = 0
               DO 2400 IGC = 1,NGC(1)
                  SUMK = 0.
                  DO 2600 IPR = 1, NGN(IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KA(JN,JT,JP,IPRSM)*RWGT(IPRSM)
 2600             CONTINUE
                  KAC(JN,JT,JP,IGC) = SUMK
 2400          CONTINUE
 2200       CONTINUE
 2000 CONTINUE

      DO 3000 JT = 1,5
         DO 3200 JP = 13,59
            IPRSM = 0
            DO 3400 IGC = 1,NGC(1)
               SUMK = 0.
               DO 3600 IPR = 1, NGN(IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KB(JT,JP,IPRSM)*RWGT(IPRSM)
 3600          CONTINUE
               KBC(JT,JP,IGC) = SUMK
 3400       CONTINUE
 3200    CONTINUE
 3000 CONTINUE

      DO 4000 JT = 1,10
         IPRSM = 0
         DO 4400 IGC = 1,NGC(1)
            SUMK = 0.
            DO 4600 IPR = 1, NGN(IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM)
 4600       CONTINUE
            SELFREFC(JT,IGC) = SUMK
 4400    CONTINUE
 4000 CONTINUE

      DO 5000 JT = 1,3
         IPRSM = 0
         DO 5400 IGC = 1,NGC(1)
            SUMK = 0.
            DO 5600 IPR = 1, NGN(IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + FORREF(JT,IPRSM)*RWGT(IPRSM)
 5600       CONTINUE
            FORREFC(JT,IGC) = SUMK
 5400    CONTINUE
 5000 CONTINUE

      IPRSM = 0
      DO 6000 IGC = 1,NGC(1)
         SUMF = 0.
         DO 6200 IPR = 1, NGN(IGC)
            IPRSM = IPRSM + 1
            SUMF = SUMF + SFLUXREF(IPRSM)
 6200    CONTINUE
         SFLUXREFC(IGC) = SUMF
 6000 CONTINUE

RETURN
END SUBROUTINE CMBGB16S

!-----------------------------------------------------------------------
SUBROUTINE CMBGB17

!     BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA17, ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, &
                    & KAC, KBC, SELFREFC, FORREFC, SFLUXREFC

IMPLICIT NONE

! Local variables
INTEGER_B :: JN, JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF

      DO 2000 JN = 1,9
         DO 2000 JT = 1,5
            DO 2200 JP = 1,13
               IPRSM = 0
               DO 2400 IGC = 1,NGC(2)
                  SUMK = 0.
                  DO 2600 IPR = 1, NGN(NGS(1)+IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KA(JN,JT,JP,IPRSM)*RWGT(IPRSM+16)
 2600             CONTINUE
                  KAC(JN,JT,JP,IGC) = SUMK
 2400          CONTINUE
 2200       CONTINUE
 2000 CONTINUE

      DO 3000 JN = 1,5
         DO 3000 JT = 1,5
            DO 3200 JP = 13,59
               IPRSM = 0
               DO 3400 IGC = 1,NGC(2)
                  SUMK = 0.
                  DO 3600 IPR = 1, NGN(NGS(1)+IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KB(JN,JT,JP,IPRSM)*RWGT(IPRSM+16)
 3600             CONTINUE
                  KBC(JN,JT,JP,IGC) = SUMK
 3400          CONTINUE
 3200       CONTINUE
 3000 CONTINUE

      DO 4000 JT = 1,10
         IPRSM = 0
         DO 4400 IGC = 1,NGC(2)
            SUMK = 0.
            DO 4600 IPR = 1, NGN(NGS(1)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+16)
 4600       CONTINUE
            SELFREFC(JT,IGC) = SUMK
 4400    CONTINUE
 4000 CONTINUE

      DO 5000 JT = 1,4
         IPRSM = 0
         DO 5400 IGC = 1,NGC(2)
            SUMK = 0.
            DO 5600 IPR = 1, NGN(NGS(1)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+16)
 5600       CONTINUE
            FORREFC(JT,IGC) = SUMK
 5400    CONTINUE
 5000 CONTINUE

      DO 7000 JP = 1,5
         IPRSM = 0
         DO 7400 IGC = 1,NGC(2)
            SUMF = 0.
            DO 7600 IPR = 1, NGN(NGS(1)+IGC)
               IPRSM = IPRSM + 1
               SUMF = SUMF + SFLUXREF(IPRSM,JP)
 7600       CONTINUE
            SFLUXREFC(IGC,JP) = SUMF
 7400    CONTINUE
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB17

!-----------------------------------------------------------------------
SUBROUTINE CMBGB18

!     BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA18, ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, &
                    & KAC, KBC, SELFREFC, FORREFC, SFLUXREFC

IMPLICIT NONE

! Local variables
INTEGER_B :: JN, JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF

      DO 2000 JN = 1,9
         DO 2000 JT = 1,5
            DO 2200 JP = 1,13
               IPRSM = 0
               DO 2400 IGC = 1,NGC(3)
                  SUMK = 0.
                  DO 2600 IPR = 1, NGN(NGS(2)+IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KA(JN,JT,JP,IPRSM)*RWGT(IPRSM+32)
 2600             CONTINUE
                  KAC(JN,JT,JP,IGC) = SUMK
 2400          CONTINUE
 2200       CONTINUE
 2000 CONTINUE

      DO 3000 JT = 1,5
         DO 3200 JP = 13,59
            IPRSM = 0
            DO 3400 IGC = 1,NGC(3)
               SUMK = 0.
               DO 3600 IPR = 1, NGN(NGS(2)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KB(JT,JP,IPRSM)*RWGT(IPRSM+32)
 3600          CONTINUE
               KBC(JT,JP,IGC) = SUMK
 3400       CONTINUE
 3200    CONTINUE
 3000 CONTINUE

      DO 4000 JT = 1,10
         IPRSM = 0
         DO 4400 IGC = 1,NGC(3)
            SUMK = 0.
            DO 4600 IPR = 1, NGN(NGS(2)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+32)
 4600       CONTINUE
            SELFREFC(JT,IGC) = SUMK
 4400    CONTINUE
 4000 CONTINUE

      DO 5000 JT = 1,3
         IPRSM = 0
         DO 5400 IGC = 1,NGC(3)
            SUMK = 0.
            DO 5600 IPR = 1, NGN(NGS(2)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+32)
 5600       CONTINUE
            FORREFC(JT,IGC) = SUMK
 5400    CONTINUE
 5000 CONTINUE

      DO 7000 JP = 1,9
         IPRSM = 0
         DO 7400 IGC = 1,NGC(3)
            SUMF = 0.
            DO 7600 IPR = 1, NGN(NGS(2)+IGC)
               IPRSM = IPRSM + 1
               SUMF = SUMF + SFLUXREF(IPRSM,JP)
 7600       CONTINUE
            SFLUXREFC(IGC,JP) = SUMF
 7400    CONTINUE
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB18

!-----------------------------------------------------------------------
SUBROUTINE CMBGB19

!     BAND 19:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA19, ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, &
                    & KAC, KBC, SELFREFC, FORREFC, SFLUXREFC

IMPLICIT NONE

! Local variables
INTEGER_B :: JN, JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF

      DO 2000 JN = 1,9
         DO 2000 JT = 1,5
            DO 2200 JP = 1,13
               IPRSM = 0
               DO 2400 IGC = 1,NGC(4)
                  SUMK = 0.
                  DO 2600 IPR = 1, NGN(NGS(3)+IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KA(JN,JT,JP,IPRSM)*RWGT(IPRSM+48)
 2600             CONTINUE
                  KAC(JN,JT,JP,IGC) = SUMK
 2400          CONTINUE
 2200       CONTINUE
 2000 CONTINUE

      DO 3000 JT = 1,5
         DO 3200 JP = 13,59
            IPRSM = 0
            DO 3400 IGC = 1,NGC(4)
               SUMK = 0.
               DO 3600 IPR = 1, NGN(NGS(3)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KB(JT,JP,IPRSM)*RWGT(IPRSM+48)
 3600          CONTINUE
               KBC(JT,JP,IGC) = SUMK
 3400       CONTINUE
 3200    CONTINUE
 3000 CONTINUE

      DO 4000 JT = 1,10
         IPRSM = 0
         DO 4400 IGC = 1,NGC(4)
            SUMK = 0.
            DO 4600 IPR = 1, NGN(NGS(3)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+48)
 4600       CONTINUE
            SELFREFC(JT,IGC) = SUMK
 4400    CONTINUE
 4000 CONTINUE

      DO 5000 JT = 1,3
         IPRSM = 0
         DO 5400 IGC = 1,NGC(4)
            SUMK = 0.
            DO 5600 IPR = 1, NGN(NGS(3)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+48)
 5600       CONTINUE
            FORREFC(JT,IGC) = SUMK
 5400    CONTINUE
 5000 CONTINUE

      DO 7000 JP = 1,9
         IPRSM = 0
         DO 7400 IGC = 1,NGC(4)
            SUMF = 0.
            DO 7600 IPR = 1, NGN(NGS(3)+IGC)
               IPRSM = IPRSM + 1
               SUMF = SUMF + SFLUXREF(IPRSM,JP)
 7600       CONTINUE
            SFLUXREFC(IGC,JP) = SUMF
 7400    CONTINUE
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB19

!-----------------------------------------------------------------------
SUBROUTINE CMBGB20

!     BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA20, ONLY : KA, KB, SELFREF, FORREF, ABSCH4, SFLUXREF, &
                    & KAC, KBC, SELFREFC, FORREFC, ABSCH4C, SFLUXREFC

IMPLICIT NONE

! Local variables
INTEGER_B :: JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF1, SUMF2

      DO 2000 JT = 1,5
         DO 2200 JP = 1,13
            IPRSM = 0
            DO 2400 IGC = 1,NGC(5)
               SUMK = 0.
               DO 2600 IPR = 1, NGN(NGS(4)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KA(JT,JP,IPRSM)*RWGT(IPRSM+64)
 2600          CONTINUE
               KAC(JT,JP,IGC) = SUMK
 2400       CONTINUE
 2200    CONTINUE
         DO 3200 JP = 13,59
            IPRSM = 0
            DO 3400 IGC = 1,NGC(5)
               SUMK = 0.
               DO 3600 IPR = 1, NGN(NGS(4)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KB(JT,JP,IPRSM)*RWGT(IPRSM+64)
 3600          CONTINUE
               KBC(JT,JP,IGC) = SUMK
 3400       CONTINUE
 3200    CONTINUE
 2000 CONTINUE

      DO 4000 JT = 1,10
         IPRSM = 0
         DO 4400 IGC = 1,NGC(5)
            SUMK = 0.
            DO 4600 IPR = 1, NGN(NGS(4)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+64)
 4600       CONTINUE
            SELFREFC(JT,IGC) = SUMK
 4400    CONTINUE
 4000 CONTINUE

      DO 5000 JT = 1,4
         IPRSM = 0
         DO 5400 IGC = 1,NGC(5)
            SUMK = 0.
            DO 5600 IPR = 1, NGN(NGS(4)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+64)
 5600       CONTINUE
            FORREFC(JT,IGC) = SUMK
 5400    CONTINUE
 5000 CONTINUE

      IPRSM = 0
      DO 7000 IGC = 1,NGC(5)
         SUMF1 = 0.
         SUMF2 = 0.
         DO 7200 IPR = 1, NGN(NGS(4)+IGC)
            IPRSM = IPRSM + 1
            SUMF1 = SUMF1 + SFLUXREF(IPRSM)
            SUMF2 = SUMF2 + ABSCH4(IPRSM)*RWGT(IPRSM+64)
 7200    CONTINUE
         SFLUXREFC(IGC) = SUMF1
         ABSCH4C(IGC) = SUMF2
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB20

!-----------------------------------------------------------------------
SUBROUTINE CMBGB21

!     BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA21, ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, &
                    & KAC, KBC, SELFREFC, FORREFC, SFLUXREFC

IMPLICIT NONE

! Local variables
INTEGER_B :: JN, JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF

      DO 2000 JN = 1,9
         DO 2000 JT = 1,5
            DO 2200 JP = 1,13
               IPRSM = 0
               DO 2400 IGC = 1,NGC(6)
                  SUMK = 0.
                  DO 2600 IPR = 1, NGN(NGS(5)+IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KA(JN,JT,JP,IPRSM)*RWGT(IPRSM+80)
 2600             CONTINUE
                  KAC(JN,JT,JP,IGC) = SUMK
 2400          CONTINUE
 2200       CONTINUE
 2000 CONTINUE

      DO 3000 JN = 1,5
         DO 3000 JT = 1,5
            DO 3200 JP = 13,59
               IPRSM = 0
               DO 3400 IGC = 1,NGC(6)
                  SUMK = 0.
                  DO 3600 IPR = 1, NGN(NGS(5)+IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KB(JN,JT,JP,IPRSM)*RWGT(IPRSM+80)
 3600             CONTINUE
                  KBC(JN,JT,JP,IGC) = SUMK
 3400          CONTINUE
 3200       CONTINUE
 3000 CONTINUE

      DO 4000 JT = 1,10
         IPRSM = 0
         DO 4400 IGC = 1,NGC(6)
            SUMK = 0.
            DO 4600 IPR = 1, NGN(NGS(5)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+80)
 4600       CONTINUE
            SELFREFC(JT,IGC) = SUMK
 4400    CONTINUE
 4000 CONTINUE

      DO 5000 JT = 1,4
         IPRSM = 0
         DO 5400 IGC = 1,NGC(6)
            SUMK = 0.
            DO 5600 IPR = 1, NGN(NGS(5)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+80)
 5600       CONTINUE
            FORREFC(JT,IGC) = SUMK
 5400    CONTINUE
 5000 CONTINUE

      DO 7000 JP = 1,9
         IPRSM = 0
         DO 7400 IGC = 1,NGC(6)
            SUMF = 0.
            DO 7600 IPR = 1, NGN(NGS(5)+IGC)
               IPRSM = IPRSM + 1
               SUMF = SUMF + SFLUXREF(IPRSM,JP)
 7600       CONTINUE
            SFLUXREFC(IGC,JP) = SUMF
 7400    CONTINUE
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB21

!-----------------------------------------------------------------------
SUBROUTINE CMBGB22

!     BAND 22:  7700-8050 cm-1 (low - H2O,O2; high - O2)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA22, ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, &
                    & KAC, KBC, SELFREFC, FORREFC, SFLUXREFC

IMPLICIT NONE

! Local variables
INTEGER_B :: JN, JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF

      DO 2000 JN = 1,9
         DO 2000 JT = 1,5
            DO 2200 JP = 1,13
               IPRSM = 0
               DO 2400 IGC = 1,NGC(7)
                  SUMK = 0.
                  DO 2600 IPR = 1, NGN(NGS(6)+IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KA(JN,JT,JP,IPRSM)*RWGT(IPRSM+96)
 2600             CONTINUE
                  KAC(JN,JT,JP,IGC) = SUMK
 2400          CONTINUE
 2200       CONTINUE
 2000 CONTINUE

      DO 3000 JT = 1,5
         DO 3200 JP = 13,59
            IPRSM = 0
            DO 3400 IGC = 1,NGC(7)
               SUMK = 0.
               DO 3600 IPR = 1, NGN(NGS(6)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KB(JT,JP,IPRSM)*RWGT(IPRSM+96)
 3600          CONTINUE
               KBC(JT,JP,IGC) = SUMK
 3400       CONTINUE
 3200    CONTINUE
 3000 CONTINUE

      DO 4000 JT = 1,10
         IPRSM = 0
         DO 4400 IGC = 1,NGC(7)
            SUMK = 0.
            DO 4600 IPR = 1, NGN(NGS(6)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+96)
 4600       CONTINUE
            SELFREFC(JT,IGC) = SUMK
 4400    CONTINUE
 4000 CONTINUE

      DO 5000 JT = 1,3
         IPRSM = 0
         DO 5400 IGC = 1,NGC(7)
            SUMK = 0.
            DO 5600 IPR = 1, NGN(NGS(6)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+96)
 5600       CONTINUE
            FORREFC(JT,IGC) = SUMK
 5400    CONTINUE
 5000 CONTINUE

      DO 7000 JP = 1,9
         IPRSM = 0
         DO 7400 IGC = 1,NGC(7)
            SUMF = 0.
            DO 7600 IPR = 1, NGN(NGS(6)+IGC)
               IPRSM = IPRSM + 1
               SUMF = SUMF + SFLUXREF(IPRSM,JP)
 7600       CONTINUE
            SFLUXREFC(IGC,JP) = SUMF
 7400    CONTINUE
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB22

!-----------------------------------------------------------------------
SUBROUTINE CMBGB23

!     BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA23, ONLY : KA, SELFREF, FORREF, SFLUXREF, RAYL, &
                    & KAC, SELFREFC, FORREFC, SFLUXREFC, RAYLC

IMPLICIT NONE

! Local variables
INTEGER_B :: JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF1, SUMF2

      DO 2000 JT = 1,5
         DO 2200 JP = 1,13
            IPRSM = 0
            DO 2400 IGC = 1,NGC(8)
               SUMK = 0.
               DO 2600 IPR = 1, NGN(NGS(7)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KA(JT,JP,IPRSM)*RWGT(IPRSM+112)
 2600          CONTINUE
               KAC(JT,JP,IGC) = SUMK
 2400       CONTINUE
 2200    CONTINUE
 2000 CONTINUE

      DO 4000 JT = 1,10
         IPRSM = 0
         DO 4400 IGC = 1,NGC(8)
            SUMK = 0.
            DO 4600 IPR = 1, NGN(NGS(7)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+112)
 4600       CONTINUE
            SELFREFC(JT,IGC) = SUMK
 4400    CONTINUE
 4000 CONTINUE

      DO 5000 JT = 1,3
         IPRSM = 0
         DO 5400 IGC = 1,NGC(8)
            SUMK = 0.
            DO 5600 IPR = 1, NGN(NGS(7)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+112)
 5600       CONTINUE
            FORREFC(JT,IGC) = SUMK
 5400    CONTINUE
 5000 CONTINUE

      IPRSM = 0
      DO 7000 IGC = 1,NGC(8)
         SUMF1 = 0.
         SUMF2 = 0.
         DO 7200 IPR = 1, NGN(NGS(7)+IGC)
            IPRSM = IPRSM + 1
            SUMF1 = SUMF1 + SFLUXREF(IPRSM)
            SUMF2 = SUMF2 + RAYL(IPRSM)*RWGT(IPRSM+112)
 7200    CONTINUE
         SFLUXREFC(IGC) = SUMF1
         RAYLC(IGC) = SUMF2
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB23

!-----------------------------------------------------------------------
SUBROUTINE CMBGB24

!     BAND 24:  12850-16000 cm-1 (low - H2O,O2; high - O2)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA24, ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, &
                    & ABSO3A, ABSO3B, RAYLA, RAYLB, &
                    & KAC, KBC, SELFREFC, FORREFC, SFLUXREFC, &
                    & ABSO3AC, ABSO3BC, RAYLAC, RAYLBC

IMPLICIT NONE

! Local variables
INTEGER_B :: JN, JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF1, SUMF2, SUMF3

      DO 2000 JN = 1,9
         DO 2000 JT = 1,5
            DO 2200 JP = 1,13
               IPRSM = 0
               DO 2400 IGC = 1,NGC(9)
                  SUMK = 0.
                  DO 2600 IPR = 1, NGN(NGS(8)+IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KA(JN,JT,JP,IPRSM)*RWGT(IPRSM+128)
 2600             CONTINUE
                  KAC(JN,JT,JP,IGC) = SUMK
 2400          CONTINUE
 2200       CONTINUE
 2000 CONTINUE

      DO 3000 JT = 1,5
         DO 3200 JP = 13,59
            IPRSM = 0
            DO 3400 IGC = 1,NGC(9)
               SUMK = 0.
               DO 3600 IPR = 1, NGN(NGS(8)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KB(JT,JP,IPRSM)*RWGT(IPRSM+128)
 3600          CONTINUE
               KBC(JT,JP,IGC) = SUMK
 3400       CONTINUE
 3200    CONTINUE
 3000 CONTINUE

      DO 4000 JT = 1,10
         IPRSM = 0
         DO 4400 IGC = 1,NGC(9)
            SUMK = 0.
            DO 4600 IPR = 1, NGN(NGS(8)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+128)
 4600       CONTINUE
            SELFREFC(JT,IGC) = SUMK
 4400    CONTINUE
 4000 CONTINUE

      DO 5000 JT = 1,3
         IPRSM = 0
         DO 5400 IGC = 1,NGC(9)
            SUMK = 0.
            DO 5600 IPR = 1, NGN(NGS(8)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+128)
 5600       CONTINUE
            FORREFC(JT,IGC) = SUMK
 5400    CONTINUE
 5000 CONTINUE

      IPRSM = 0
      DO 6000 IGC = 1,NGC(9)
         SUMF1 = 0.
         SUMF2 = 0.
         SUMF3 = 0.
         DO 6200 IPR = 1, NGN(NGS(8)+IGC)
            IPRSM = IPRSM + 1
            SUMF1 = SUMF1 + RAYLB(IPRSM)*RWGT(IPRSM+128)
            SUMF2 = SUMF2 + ABSO3A(IPRSM)*RWGT(IPRSM+128)
            SUMF3 = SUMF3 + ABSO3B(IPRSM)*RWGT(IPRSM+128)
 6200    CONTINUE
         RAYLBC(IGC) = SUMF1
         ABSO3AC(IGC) = SUMF2
         ABSO3BC(IGC) = SUMF3
 6000 CONTINUE

      DO 7000 JP = 1,9
         IPRSM = 0
         DO 7400 IGC = 1,NGC(9)
            SUMF1 = 0.
            SUMF2 = 0.
            DO 7600 IPR = 1, NGN(NGS(8)+IGC)
               IPRSM = IPRSM + 1
               SUMF1 = SUMF1 + SFLUXREF(IPRSM,JP)
               SUMF2 = SUMF2 + RAYLA(IPRSM,JP)*RWGT(IPRSM+128)
 7600       CONTINUE
            SFLUXREFC(IGC,JP) = SUMF1
            RAYLAC(IGC,JP) = SUMF2
 7400    CONTINUE
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB24

!-----------------------------------------------------------------------
SUBROUTINE CMBGB25

!     BAND 25:  16000-22650 cm-1 (low - H2O; high - nothing)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA25, ONLY : KA, SFLUXREF, ABSO3A, ABSO3B, RAYL, &
                    & KAC, SFLUXREFC, ABSO3AC, ABSO3BC, RAYLC

IMPLICIT NONE

! Local variables
INTEGER_B :: JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF1, SUMF2, SUMF3, SUMF4

      DO 2000 JT = 1,5
         DO 2200 JP = 1,13
            IPRSM = 0
            DO 2400 IGC = 1,NGC(10)
               SUMK = 0.
               DO 2600 IPR = 1, NGN(NGS(9)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KA(JT,JP,IPRSM)*RWGT(IPRSM+144)
 2600          CONTINUE
               KAC(JT,JP,IGC) = SUMK
 2400       CONTINUE
 2200    CONTINUE
 2000 CONTINUE

      IPRSM = 0
      DO 6000 IGC = 1,NGC(10)
         SUMF1 = 0.
         SUMF2 = 0.
         SUMF3 = 0.
         SUMF4 = 0.
         DO 6200 IPR = 1, NGN(NGS(9)+IGC)
            IPRSM = IPRSM + 1
            SUMF1 = SUMF1 + SFLUXREF(IPRSM)
            SUMF2 = SUMF2 + ABSO3A(IPRSM)*RWGT(IPRSM+144)
            SUMF3 = SUMF3 + ABSO3B(IPRSM)*RWGT(IPRSM+144)
            SUMF4 = SUMF4 + RAYL(IPRSM)*RWGT(IPRSM+144)
 6200    CONTINUE
         SFLUXREFC(IGC) = SUMF1
         ABSO3AC(IGC) = SUMF2
         ABSO3BC(IGC) = SUMF3
         RAYLC(IGC) = SUMF4
 6000 CONTINUE

RETURN
END SUBROUTINE CMBGB25

!-----------------------------------------------------------------------
SUBROUTINE CMBGB26

!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA26, ONLY : SFLUXREF, RAYL, &
                    & SFLUXREFC, RAYLC

IMPLICIT NONE

! Local variables
INTEGER_B :: IGC, IPR, IPRSM
REAL_B :: SUMF1, SUMF2

      IPRSM = 0
      DO 6000 IGC = 1,NGC(11)
         SUMF1 = 0.
         SUMF2 = 0.
         DO 6200 IPR = 1, NGN(NGS(10)+IGC)
            IPRSM = IPRSM + 1
            SUMF1 = SUMF1 + RAYL(IPRSM)*RWGT(IPRSM+160)
            SUMF2 = SUMF2 + SFLUXREF(IPRSM)
 6200    CONTINUE
         RAYLC(IGC) = SUMF1
         SFLUXREFC(IGC) = SUMF2
 6000 CONTINUE

RETURN
END SUBROUTINE CMBGB26

!-----------------------------------------------------------------------
SUBROUTINE CMBGB27

!     BAND 27:  29000-38000 cm-1 (low - O3; high - O3)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA27, ONLY : KA, KB, SFLUXREF, RAYL, &
                    & KAC, KBC, SFLUXREFC, RAYLC

IMPLICIT NONE

! Local variables
INTEGER_B :: JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF1, SUMF2

      DO 2000 JT = 1,5
         DO 2200 JP = 1,13
            IPRSM = 0
            DO 2400 IGC = 1,NGC(12)
               SUMK = 0.
               DO 2600 IPR = 1, NGN(NGS(11)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KA(JT,JP,IPRSM)*RWGT(IPRSM+176)
 2600          CONTINUE
               KAC(JT,JP,IGC) = SUMK
 2400       CONTINUE
 2200    CONTINUE
         DO 3200 JP = 13,59
            IPRSM = 0
            DO 3400 IGC = 1,NGC(12)
               SUMK = 0.
               DO 3600 IPR = 1, NGN(NGS(11)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KB(JT,JP,IPRSM)*RWGT(IPRSM+176)
 3600          CONTINUE
               KBC(JT,JP,IGC) = SUMK
 3400       CONTINUE
 3200    CONTINUE
 2000 CONTINUE

      IPRSM = 0
      DO 7000 IGC = 1,NGC(12)
         SUMF1 = 0.
         SUMF2 = 0.
         DO 7200 IPR = 1, NGN(NGS(11)+IGC)
            IPRSM = IPRSM + 1
            SUMF1 = SUMF1 + SFLUXREF(IPRSM)
            SUMF2 = SUMF2 + RAYL(IPRSM)*RWGT(IPRSM+176)
 7200    CONTINUE
         SFLUXREFC(IGC) = SUMF1
         RAYLC(IGC) = SUMF2
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB27

!-----------------------------------------------------------------------
SUBROUTINE CMBGB28

!     BAND 28:  38000-50000 cm-1 (low - O3,O2; high - O3,O2)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA28, ONLY : KA, KB, SFLUXREF, &
                    & KAC, KBC, SFLUXREFC

IMPLICIT NONE

! Local variables
INTEGER_B :: JN, JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF

      DO 2000 JN = 1,9
         DO 2000 JT = 1,5
            DO 2200 JP = 1,13
               IPRSM = 0
               DO 2400 IGC = 1,NGC(13)
                  SUMK = 0.
                  DO 2600 IPR = 1, NGN(NGS(12)+IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KA(JN,JT,JP,IPRSM)*RWGT(IPRSM+192)
 2600             CONTINUE
                  KAC(JN,JT,JP,IGC) = SUMK
 2400          CONTINUE
 2200       CONTINUE
 2000 CONTINUE

      DO 3000 JN = 1,5
         DO 3000 JT = 1,5
            DO 3200 JP = 13,59
               IPRSM = 0
               DO 3400 IGC = 1,NGC(13)
                  SUMK = 0.
                  DO 3600 IPR = 1, NGN(NGS(12)+IGC)
                     IPRSM = IPRSM + 1
                     SUMK = SUMK + KB(JN,JT,JP,IPRSM)*RWGT(IPRSM+192)
 3600             CONTINUE
                  KBC(JN,JT,JP,IGC) = SUMK
 3400          CONTINUE
 3200       CONTINUE
 3000 CONTINUE

      DO 7000 JP = 1,5
         IPRSM = 0
         DO 7400 IGC = 1,NGC(13)
            SUMF = 0.
            DO 7600 IPR = 1, NGN(NGS(12)+IGC)
               IPRSM = IPRSM + 1
               SUMF = SUMF + SFLUXREF(IPRSM,JP)
 7600       CONTINUE
            SFLUXREFC(IGC,JP) = SUMF
 7400    CONTINUE
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB28

!-----------------------------------------------------------------------
SUBROUTINE CMBGB29

!     BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)
!-----------------------------------------------------------------------

! Modules
#include "tsmbkind.h"
USE YOESRTWN , ONLY : NGC, NGS, NGN, WT, WTSM, RWGT
USE YOESRTA29, ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, &
                    & ABSH2O, ABSCO2, &
                    & KAC, KBC, SELFREFC, FORREFC, SFLUXREFC, &
                    & ABSH2OC, ABSCO2C

IMPLICIT NONE

! Local variables
INTEGER_B :: JT, JP, IGC, IPR, IPRSM
REAL_B :: SUMK, SUMF1, SUMF2, SUMF3

      DO 2000 JT = 1,5
         DO 2200 JP = 1,13
            IPRSM = 0
            DO 2400 IGC = 1,NGC(14)
               SUMK = 0.
               DO 2600 IPR = 1, NGN(NGS(13)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KA(JT,JP,IPRSM)*RWGT(IPRSM+208)
 2600          CONTINUE
               KAC(JT,JP,IGC) = SUMK
 2400       CONTINUE
 2200    CONTINUE
         DO 3200 JP = 13,59
            IPRSM = 0
            DO 3400 IGC = 1,NGC(14)
               SUMK = 0.
               DO 3600 IPR = 1, NGN(NGS(13)+IGC)
                  IPRSM = IPRSM + 1
                  SUMK = SUMK + KB(JT,JP,IPRSM)*RWGT(IPRSM+208)
 3600          CONTINUE
               KBC(JT,JP,IGC) = SUMK
 3400       CONTINUE
 3200    CONTINUE
 2000 CONTINUE

      DO 4000 JT = 1,10
         IPRSM = 0
         DO 4400 IGC = 1,NGC(14)
            SUMK = 0.
            DO 4600 IPR = 1, NGN(NGS(13)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+208)
 4600       CONTINUE
            SELFREFC(JT,IGC) = SUMK
 4400    CONTINUE
 4000 CONTINUE

      DO 5000 JT = 1,4
         IPRSM = 0
         DO 5400 IGC = 1,NGC(14)
            SUMK = 0.
            DO 5600 IPR = 1, NGN(NGS(13)+IGC)
               IPRSM = IPRSM + 1
               SUMK = SUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+208)
 5600       CONTINUE
            FORREFC(JT,IGC) = SUMK
 5400    CONTINUE
 5000 CONTINUE

      IPRSM = 0
      DO 7000 IGC = 1,NGC(14)
         SUMF1 = 0.
         SUMF2 = 0.
         SUMF3 = 0.
         DO 7200 IPR = 1, NGN(NGS(13)+IGC)
            IPRSM = IPRSM + 1
            SUMF1 = SUMF1 + SFLUXREF(IPRSM)
            SUMF2 = SUMF2 + ABSCO2(IPRSM)*RWGT(IPRSM+208)
            SUMF3 = SUMF3 + ABSH2O(IPRSM)*RWGT(IPRSM+208)
 7200    CONTINUE
         SFLUXREFC(IGC) = SUMF1
         ABSCO2C(IGC) = SUMF2
         ABSH2OC(IGC) = SUMF3
 7000 CONTINUE

RETURN
END SUBROUTINE CMBGB29
