C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

SUBROUTINE RRTM_SW_VRTQDR &
 &( KLEV , KW &
 &, PREF , PREFD, PTRA , PTRAD &
 &, PDBT , PRDND, PRUP , PRUPD , PTDBT &
 &, PFD  , PFU  &
 &)
 
!**** *RRTM_SW_VRTQDR* - VERTICAL QUADRATURE

!     PURPOSE.
!     --------

!          THIS ROUTINE PERFORMS THE VERTICAL INTEGRATION

!**   INTERFACE.
!     ----------

!          *RRTM_SW_VRTQDR* IS CALLED FROM *RRTM_SW_SPCVRT*


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

!USE YOESWN   , ONLY : NDBUG


IMPLICIT NONE

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER_M :: KLEV, KW

REAL_B :: PREF(JPLAY+1), PREFD(JPLAY+1), PTRA(JPLAY+1), PTRAD(JPLAY+1)
REAL_B :: PDBT(JPLAY+1), PRDND(JPLAY+1), PRUP(JPLAY+1), PRUPD(JPLAY+1), PTDBT(JPLAY+1)

REAL_B :: PFD(JPLAY+1,JPGPT), PFU(JPLAY+1,JPGPT) 

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL_B :: ZTDN(JPLAY+1)  

!     LOCAL INTEGER SCALARS
INTEGER_M :: IKP, IKX, JK
!INTEGER_M :: NDBUG

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
                   
!NDBUG=3
    
!-- link lowest layer with surface
             
ZREFLECT=_ONE_ / (_ONE_ -PREFD(KLEV+1)*PREFD(KLEV))
PRUP(KLEV)=PREF(KLEV)+(PTRAD(KLEV)* &
  &          ((PTRA(KLEV)-PDBT(KLEV))*PREFD(KLEV+1)+ &
  &           PDBT(KLEV)*PREF(KLEV+1)))*ZREFLECT
PRUPD(KLEV)=PREFD(KLEV)+PTRAD(KLEV)* &
  &          PTRAD(KLEV)*PREFD(KLEV+1)*ZREFLECT

!IF (NDBUG.LE.1) THEN
!  print 9201,PRUP(KLEV),PRUPD(KLEV)
9201 format(1x,'link surf:',6E13.6)
!  print *,'RRTM_SW_VRTQDR after linking with surface layer'
!END IF
    
    
!-- pass from bottom to top 

DO JK=1,KLEV-1
  IKP=KLEV+1-JK                       
  IKX=IKP-1
!  print 9202,JK,IKP,IKX
9202  format(1x,'Pass from bottom to top:',3I3)      
  ZREFLECT=_ONE_ / (_ONE_ -PRUPD(IKP)*PREFD(IKX))
  PRUP(IKX)=PREF(IKX)+(PTRAD(IKX)* &
    &          ((PTRA(IKX)-PDBT(IKX))*PRUPD(IKP)+ &
    &           PDBT(IKX)*PRUP(IKP)))*ZREFLECT
  PRUPD(IKX)=PREFD(IKX)+PTRAD(IKX)* &
    &          PTRAD(IKX)*PRUPD(IKP)*ZREFLECT

!  print 9203,PRUP(IKX),PRUPD(IKX)
9203 format(1x,'bot2top:',6E13.6)
END DO
!print *,'RRTM_SW_VRTQDR after passing from bottom to top'
    
    
!-- upper boundary conditions

ZTDN(1)=_ONE_
PRDND(1)=_ZERO_
ZTDN(2)=PTRA(1)
PRDND(2)=PREFD(1)

!IF (NDBUG.LE.1) THEN
!  print 9204,ZTDN(1),PRDND(1),ZTDN(2),PRDND(2)
9204 format(1x,'link upper bound:',6E13.6)
!  print *,'RRTM_SW_VRTQDR after upper boundary conditions'
!END IF
    
!-- pass from top to bottom

DO JK=2,KLEV
  IKP=JK+1
  ZREFLECT=_ONE_ / (_ONE_ -PREFD(JK)*PRDND(JK))
  ZTDN(IKP)=PTDBT(JK)*PTRA(JK)+ &
    & (PTRAD(JK)*((ZTDN(JK)-PTDBT(JK))+ &
    &  PTDBT(JK)*PREF(JK)*PRDND(JK))) * ZREFLECT
  PRDND(IKP)=PREFD(JK)+PTRAD(JK)*PTRAD(JK) &
    & *PRDND(JK)*ZREFLECT

!  IF (NDBUG.LE.1) THEN
!    print 9205,ZTDN(IKP),PRDND(IKP)
9205 format(1x,'top2bot2:',6E13.6)
!  END IF

END DO
!print *,'RRTM_SW_VRTQDR after passing from top to bottom'
                                              
    
!-- up and down-welling fluxes at levels

DO JK=1,KLEV+1
!  IF (NDBUG.LE.1) THEN
!    print 9207,JK,PRDND(JK),PRUPD(JK)
!    print 9208,JK,PTDBT(JK),PRUP(JK),ZTDN(JK)
9207 format(1x,'A',I3,4E13.6)      
9208 format(1x,'B',I3,4E13.6)      
!  END IF 

  ZREFLECT=_ONE_ / (_ONE_ - PRDND(JK)*PRUPD(JK))
  PFU(JK,KW)=(PTDBT(JK)*PRUP(JK) + &
    &        (ZTDN(JK)-PTDBT(JK))*PRUPD(JK))*ZREFLECT
  PFD(JK,KW)=PTDBT(JK) + (ZTDN(JK)-PTDBT(JK)+ &
    &        PTDBT(JK)*PRUP(JK)*PRDND(JK))*ZREFLECT

!  IF (NDBUG.LE.2) THEN
!    print 9206,JK,PFU(JK,KW),PFD(JK,KW)
9206 format(1x,'fluxes:',I3,6E13.6)
!  END IF

END DO
!print *,'RRTM_SW_VRTQDR after up and down flux'
   
!print *,'RRTM_SW_VRTQDR about to come out'
!     ------------------------------------------------------------------

RETURN
END SUBROUTINE RRTM_SW_VRTQDR

