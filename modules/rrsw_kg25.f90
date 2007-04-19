      module rrsw_kg25

      use parkind ,only : jpim, jprb
      use parrrsw, only : ng25

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 25
! band 25: 16000-22650 cm-1 (low - h2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! kao     : real     
!sfluxrefo: real     
! abso3ao : real     
! abso3bo : real     
! raylo   : real     
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: no25 = 16

      real(kind=jprb) :: kao(5,13,no25)
      real(kind=jprb) :: sfluxrefo(no25)
      real(kind=jprb) :: abso3ao(no25), abso3bo(no25)
      real(kind=jprb) :: raylo(no25)

      integer(kind=jpim) :: layreffr

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 25
! band 25: 16000-22650 cm-1 (low - h2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! ka      : real     
! absa    : real
! sfluxref: real     
! abso3a  : real     
! abso3b  : real     
! rayl    : real     
!-----------------------------------------------------------------

      real(kind=jprb) :: ka(5,13,ng25), absa(65,ng25)
      real(kind=jprb) :: sfluxref(ng25)
      real(kind=jprb) :: abso3a(ng25), abso3b(ng25)
      real(kind=jprb) :: rayl(ng25)

      equivalence (ka(1,1,1),absa(1,1))

      end module rrsw_kg25

