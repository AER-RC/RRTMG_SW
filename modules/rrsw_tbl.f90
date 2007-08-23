      module rrsw_tbl

      use parkind, only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw lookup table arrays

! Initial version: MJIacono, AER, may2007
! Revised: MJIacono, AER, aug2007
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ntbl   :  integer: Lookup table dimension
! tblint :  real   : Lookup table conversion factor
! tau_tbl:  real   : Clear-sky optical depth 
! exp_tbl:  real   : Exponential lookup table for transmittance
! od_lo  :  real   : Value of tau below which expansion is used
!                  : in place of lookup table
! pade   :  real   : Pade approximation constant
! bpade  :  real   : Inverse of Pade constant
!------------------------------------------------------------------

      integer(kind=jpim), parameter :: ntbl = 10000

      real(kind=jprb), parameter :: tblint = 10000.0

      real(kind=jprb), parameter :: od_lo = 0.06

      real(kind=jprb) :: tau_tbl
      real(kind=jprb) , dimension(0:ntbl) :: exp_tbl

      real(kind=jprb), parameter :: pade = 0.278_jprb
      real(kind=jprb) :: bpade

      end module rrsw_tbl

