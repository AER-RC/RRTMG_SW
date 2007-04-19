      module parkind

      implicit none
      save

!------------------------------------------------------------------
! rrtmg kinds
! Define usual integer and real kinds for strong typing.
! Only jpim, jprb are currently active.
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!------------------------------------------------------------------

!
!     integer kinds
!     -------------
!     jpit = integer range +/-10**2
!     jpis = integer range +/-10**4
!     jpim = integer range +/-10**9
!     jpib = integer range +/-10**12
!
!      integer, parameter :: jpit = selected_int_kind(2)
!      integer, parameter :: jpis = selected_int_kind(4)
      integer, parameter :: jpim = selected_int_kind(9)
!      integer, parameter :: jpib = selected_int_kind(12)

! Special integer type to be used for sensitive address calculations;
! should be *8 for a machine with 8byte addressing for optimum performance
!#ifdef address64
!      integer, parameter :: jpia = jpib
!#else
!      integer, parameter :: jpia = jpim
!#endif

!
!     real kinds
!     ----------
!     jprt = real range +/-10**1, 2 decimal place precision
!     jprs = real range +/-10**2, 4 decimal place precision
!     jprm = real range +/-10**37, 6 decimal place precision
!     jprb = real range +/-10**300, 13 decimal place precision
!
!      integer, parameter :: jprt = selected_real_kind(2,1)
!      integer, parameter :: jprs = selected_real_kind(4,2)
!      integer, parameter :: jprm = selected_real_kind(6,37)
      integer, parameter :: jprb = selected_real_kind(13,300)
!
      end module parkind
