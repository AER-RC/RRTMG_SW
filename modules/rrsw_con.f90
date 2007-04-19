      module rrsw_con

      use parkind, only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw constants

! Initial version: MJIacono, AER, jun2006
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! fluxfac:  real   : radiance to flux conversion factor 
! heatfac:  real   : flux to heating rate conversion factor
!oneminus:  real   : 1.-1.e-6
! pi     :  real   : pi
! grav   :  real   : acceleration of gravity (m/s2)
! planck :  real   : planck constant
! boltz  :  real   : boltzman constant
! clight :  real   : speed of light
! avogad :  real   : avogadro's constant 
! alosmt :  real   : 
! gascon :  real   : gas constant
! radcn1 :  real   : 
! radcn2 :  real   : 
!------------------------------------------------------------------

      real(kind=jprb) :: fluxfac, heatfac
      real(kind=jprb) :: oneminus, pi, grav
      real(kind=jprb) :: planck, boltz, clight
      real(kind=jprb) :: avogad, alosmt, gascon
      real(kind=jprb) :: radcn1, radcn2

      end module rrsw_con

