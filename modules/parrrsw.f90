
      module parrrsw

      use parkind ,only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw main parameters
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! nlon   :  integer: number of columns or longitudes
! mxlay  :  integer: maximum number of layers
! mg     :  integer: number of original g-intervals per spectral band
! nbndsw :  integer: number of spectral bands
! naer   :  integer: number of aerosols
! ngpt   :  integer: total number of reduced g-intervals for rrtmg_lw
! ngNN   :  integer: number of reduced g-intervals per spectral band
! ngsNN  :  integer: cumulative number of g-intervals per band
!------------------------------------------------------------------

! Settings for single column mode.
! For GCM use, set nlon to number of longitudes, and
! mxlay to number of model layers
      integer(kind=jpim), parameter :: nlon   = 1      !jplon, klon
      integer(kind=jpim), parameter :: mxlay  = 203    !jplay, klev
      integer(kind=jpim), parameter :: mg     = 16     !jpg
      integer(kind=jpim), parameter :: nbndsw = 14     !jpsw, ksw
      integer(kind=jpim), parameter :: naer   = 6      !jpaer
      integer(kind=jpim), parameter :: mxmol  = 38
      integer(kind=jpim), parameter :: nstr   = 2
      integer(kind=jpim), parameter :: nmol   = 7
! Use for 112 g-point model   
      integer(kind=jpim), parameter :: ngpt   = 112    !jpgpt
! Use for 224 g-point model   
!      integer(kind=jpim), parameter :: ngpt   = 224   !jpgpt

! may need to rename these - from v2.6
      integer(kind=jpim), parameter :: jpband   = 29
      integer(kind=jpim), parameter :: jpb1     = 16   !istart
      integer(kind=jpim), parameter :: jpb2     = 29   !iend

      integer(kind=jpim), parameter :: jmcmu    = 32
      integer(kind=jpim), parameter :: jmumu    = 32
      integer(kind=jpim), parameter :: jmphi    = 3
      integer(kind=jpim), parameter :: jmxang   = 4
      integer(kind=jpim), parameter :: jmxstr   = 16
! ^

! Use for 112 g-point model   
      integer(kind=jpim), parameter :: ng16 = 6
      integer(kind=jpim), parameter :: ng17 = 12
      integer(kind=jpim), parameter :: ng18 = 8
      integer(kind=jpim), parameter :: ng19 = 8
      integer(kind=jpim), parameter :: ng20 = 10
      integer(kind=jpim), parameter :: ng21 = 10
      integer(kind=jpim), parameter :: ng22 = 2
      integer(kind=jpim), parameter :: ng23 = 10
      integer(kind=jpim), parameter :: ng24 = 8
      integer(kind=jpim), parameter :: ng25 = 6
      integer(kind=jpim), parameter :: ng26 = 6
      integer(kind=jpim), parameter :: ng27 = 8
      integer(kind=jpim), parameter :: ng28 = 6
      integer(kind=jpim), parameter :: ng29 = 12

      integer(kind=jpim), parameter :: ngs16 = 6
      integer(kind=jpim), parameter :: ngs17 = 18
      integer(kind=jpim), parameter :: ngs18 = 26
      integer(kind=jpim), parameter :: ngs19 = 34
      integer(kind=jpim), parameter :: ngs20 = 44
      integer(kind=jpim), parameter :: ngs21 = 54
      integer(kind=jpim), parameter :: ngs22 = 56
      integer(kind=jpim), parameter :: ngs23 = 66
      integer(kind=jpim), parameter :: ngs24 = 74
      integer(kind=jpim), parameter :: ngs25 = 80
      integer(kind=jpim), parameter :: ngs26 = 86
      integer(kind=jpim), parameter :: ngs27 = 94
      integer(kind=jpim), parameter :: ngs28 = 100
      integer(kind=jpim), parameter :: ngs29 = 112

! Use for 224 g-point model   
!      integer(kind=jpim), parameter :: ng16 = 16
!      integer(kind=jpim), parameter :: ng17 = 16
!      integer(kind=jpim), parameter :: ng18 = 16
!      integer(kind=jpim), parameter :: ng19 = 16
!      integer(kind=jpim), parameter :: ng20 = 16
!      integer(kind=jpim), parameter :: ng21 = 16
!      integer(kind=jpim), parameter :: ng22 = 16
!      integer(kind=jpim), parameter :: ng23 = 16
!      integer(kind=jpim), parameter :: ng24 = 16
!      integer(kind=jpim), parameter :: ng25 = 16
!      integer(kind=jpim), parameter :: ng26 = 16
!      integer(kind=jpim), parameter :: ng27 = 16
!      integer(kind=jpim), parameter :: ng28 = 16
!      integer(kind=jpim), parameter :: ng29 = 16

!      integer(kind=jpim), parameter :: ngs16 = 16
!      integer(kind=jpim), parameter :: ngs17 = 32
!      integer(kind=jpim), parameter :: ngs18 = 48
!      integer(kind=jpim), parameter :: ngs19 = 64
!      integer(kind=jpim), parameter :: ngs20 = 80
!      integer(kind=jpim), parameter :: ngs21 = 96
!      integer(kind=jpim), parameter :: ngs22 = 112
!      integer(kind=jpim), parameter :: ngs23 = 128
!      integer(kind=jpim), parameter :: ngs24 = 144
!      integer(kind=jpim), parameter :: ngs25 = 160
!      integer(kind=jpim), parameter :: ngs26 = 176
!      integer(kind=jpim), parameter :: ngs27 = 192
!      integer(kind=jpim), parameter :: ngs28 = 208
!      integer(kind=jpim), parameter :: ngs29 = 224

      end module parrrsw


