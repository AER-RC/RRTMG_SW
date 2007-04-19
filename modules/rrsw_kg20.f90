      module rrsw_kg20

      use parkind ,only : jpim, jprb
      use parrrsw, only : ng20

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 20
! band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! kao     : real     
! kbo     : real     
! selfrefo: real     
! forrefo : real
!sfluxrefo: real     
! absch4o : real     
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: no20 = 16

      real(kind=jprb) :: kao(5,13,no20)
      real(kind=jprb) :: kbo(5,13:59,no20)
      real(kind=jprb) :: selfrefo(10,no20), forrefo(4,no20)
      real(kind=jprb) :: sfluxrefo(no20)
      real(kind=jprb) :: absch4o(no20)

      integer(kind=jpim) :: layreffr
      real(kind=jprb) :: rayl 

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 20
! band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! ka      : real     
! kb      : real     
! absa    : real
! absb    : real
! selfref : real     
! forref  : real
! sfluxref: real     
! absch4  : real     
!-----------------------------------------------------------------

      real(kind=jprb) :: ka(5,13,ng20), absa(65,ng20)
      real(kind=jprb) :: kb(5,13:59,ng20), absb(235,ng20)
      real(kind=jprb) :: selfref(10,ng20), forref(4,ng20)
      real(kind=jprb) :: sfluxref(ng20)
      real(kind=jprb) :: absch4(ng20)

      equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg20

