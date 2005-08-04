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

SUBROUTINE RRTMG_SW_VRTQDR &
 &( KLEV , KW &
 &, PREF , PREFD, PTRA , PTRAD &
 &, PDBT , PRDND, PRUP , PRUPD , PTDBT &
 &, PFD  , PFU  &
 &)
 
!**** *RRTMG_SW_VRTQDR* - VERTICAL QUADRATURE

!     PURPOSE.
!     --------

!          THIS ROUTINE PERFORMS THE VERTICAL INTEGRATION

!**   INTERFACE.
!     ----------

!          *RRTMG_SW_VRTQDR* IS CALLED FROM *RRTMG_SW_SPCVRT*


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------
!          NONE


!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        from Howard Barker
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 02-10-04
!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE PARSRTM  , ONLY : JPLAY, JPGPT 

IMPLICIT NONE

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER_B :: KLEV, KW

REAL_B :: PREF(JPLAY+1), PREFD(JPLAY+1), PTRA(JPLAY+1), PTRAD(JPLAY+1)
REAL_B :: PDBT(JPLAY+1), PRDND(JPLAY+1), PRUP(JPLAY+1), PRUPD(JPLAY+1), PTDBT(JPLAY+1)

REAL_B :: PFD(JPLAY+1,JPGPT), PFU(JPLAY+1,JPGPT) 

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL_B :: ZTDN(JPLAY+1)  

!     LOCAL INTEGER SCALARS
INTEGER_B :: IKP, IKX, JK

!     LOCAL REAL SCALARS
REAL_B :: ZREFLECT

!     ------------------------------------------------------------------

! PREF(JK)   direct reflectance
! PREFD(JK)  diffuse reflectance
! PTRA(JK)   direct transmittance
! PTRAD(JK)  diffuse transmittance
!
! PDBT(JK)   layer mean direct beam transmittance
! PTDBT(JK)  total direct beam transmittance at levels
                   
!-- link lowest layer with surface
             
ZREFLECT=_ONE_ / (_ONE_ -PREFD(KLEV+1)*PREFD(KLEV))
PRUP(KLEV)=PREF(KLEV)+(PTRAD(KLEV)* &
  &          ((PTRA(KLEV)-PDBT(KLEV))*PREFD(KLEV+1)+ &
  &           PDBT(KLEV)*PREF(KLEV+1)))*ZREFLECT
PRUPD(KLEV)=PREFD(KLEV)+PTRAD(KLEV)* &
  &          PTRAD(KLEV)*PREFD(KLEV+1)*ZREFLECT

    
!-- pass from bottom to top 

DO JK=1,KLEV-1
  IKP=KLEV+1-JK                       
  IKX=IKP-1
  ZREFLECT=_ONE_ / (_ONE_ -PRUPD(IKP)*PREFD(IKX))
  PRUP(IKX)=PREF(IKX)+(PTRAD(IKX)* &
    &          ((PTRA(IKX)-PDBT(IKX))*PRUPD(IKP)+ &
    &           PDBT(IKX)*PRUP(IKP)))*ZREFLECT
  PRUPD(IKX)=PREFD(IKX)+PTRAD(IKX)* &
    &          PTRAD(IKX)*PRUPD(IKP)*ZREFLECT
END DO
    
!-- upper boundary conditions

ZTDN(1)=_ONE_
PRDND(1)=_ZERO_
ZTDN(2)=PTRA(1)
PRDND(2)=PREFD(1)

!-- pass from top to bottom

DO JK=2,KLEV
  IKP=JK+1
  ZREFLECT=_ONE_ / (_ONE_ -PREFD(JK)*PRDND(JK))
  ZTDN(IKP)=PTDBT(JK)*PTRA(JK)+ &
    & (PTRAD(JK)*((ZTDN(JK)-PTDBT(JK))+ &
    &  PTDBT(JK)*PREF(JK)*PRDND(JK))) * ZREFLECT
  PRDND(IKP)=PREFD(JK)+PTRAD(JK)*PTRAD(JK) &
    & *PRDND(JK)*ZREFLECT
END DO
    
!-- up and down-welling fluxes at levels

DO JK=1,KLEV+1
  ZREFLECT=_ONE_ / (_ONE_ - PRDND(JK)*PRUPD(JK))
  PFU(JK,KW)=(PTDBT(JK)*PRUP(JK) + &
    &        (ZTDN(JK)-PTDBT(JK))*PRUPD(JK))*ZREFLECT
  PFD(JK,KW)=PTDBT(JK) + (ZTDN(JK)-PTDBT(JK)+ &
    &        PTDBT(JK)*PRUP(JK)*PRDND(JK))*ZREFLECT
END DO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE RRTMG_SW_VRTQDR

