C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

C     path:      Source: /storm/rc1/cvsroot/rc/rrtm_sw/src/param.f,v
C     author:    Author: jdelamer
C     revision:  Revision: 2.3
C     created:   Date: 2001/10/11 17:53:22

	parameter (mxlay = 203, nbands = 29)
	parameter (ib1 = 16, ib2 = 29)
        parameter (mg = 16)
	parameter (mxstr = 16)

      COMMON /BANDS/     WAVENUM1(IB1:IB2),
     &                   WAVENUM2(IB1:IB2),
     &                   DELWAVE(IB1:IB2)
