C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

SUBROUTINE RRTM_SW_INIT

!-- Read in the basic coefficients to configure RRTM_SW.
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

!-- Read in the cloud optical properties for RRTM_SW_CLDPROP.
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
END SUBROUTINE RRTM_SW_INIT

