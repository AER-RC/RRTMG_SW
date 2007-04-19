!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

      subroutine vrtqdr_sw(klev, kw, &
                           pref, prefd, ptra, ptrad, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
                           pfd, pfu)
 
! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
! 
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
!
!-----------------------------------------------------------------------

! ------- Modules -------

      use parkind, only: jpim, jprb
      use parrrsw, only: mxlay, ngpt

      implicit none

! ------- Declarations -------

! Input

      integer(kind=jpim), intent (in) :: klev, kw

      real(kind=jprb), intent(in) :: pref(mxlay+1)
      real(kind=jprb), intent(in) :: prefd(mxlay+1)
      real(kind=jprb), intent(in) :: ptra(mxlay+1)
      real(kind=jprb), intent(in) :: ptrad(mxlay+1)

      real(kind=jprb), intent(in) :: pdbt(mxlay+1)
      real(kind=jprb), intent(in) :: ptdbt(mxlay+1)

      real(kind=jprb), intent(inout) :: prdnd(mxlay+1)
      real(kind=jprb), intent(inout) :: prup(mxlay+1)
      real(kind=jprb), intent(inout) :: prupd(mxlay+1)

! Output
      real(kind=jprb), intent(out) :: pfd(mxlay+1,ngpt)
      real(kind=jprb), intent(out) :: pfu(mxlay+1,ngpt) 

! Local

      integer(kind=jpim) :: ikp, ikx, jk

      real(kind=jprb) :: zreflect
      real(kind=jprb) :: ztdn(mxlay+1)  

! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------
                   
! Link lowest layer with surface
             
      zreflect = 1._jprb / (1._jprb - prefd(klev+1) * prefd(klev))
      prup(klev) = pref(klev) + (ptrad(klev) * &
                 ((ptra(klev) - pdbt(klev)) * prefd(klev+1) + &
                   pdbt(klev) * pref(klev+1))) * zreflect
      prupd(klev) = prefd(klev) + ptrad(klev) * ptrad(klev) * &
                    prefd(klev+1) * zreflect

! Pass from bottom to top 

      do jk = 1,klev-1
         ikp = klev+1-jk                       
         ikx = ikp-1
         zreflect = 1._jprb / (1._jprb -prupd(ikp) * prefd(ikx))
         prup(ikx) = pref(ikx) + (ptrad(ikx) * &
                   ((ptra(ikx) - pdbt(ikx)) * prupd(ikp) + &
                     pdbt(ikx) * prup(ikp))) * zreflect
         prupd(ikx) = prefd(ikx) + ptrad(ikx) * ptrad(ikx) * &
                      prupd(ikp) * zreflect
      enddo
    
! Upper boundary conditions

      ztdn(1) = 1._jprb
      prdnd(1) = 0._jprb
      ztdn(2) = ptra(1)
      prdnd(2) = prefd(1)

! Pass from top to bottom

      do jk = 2,klev
         ikp = jk+1
         zreflect = 1._jprb / (1._jprb - prefd(jk) * prdnd(jk))
         ztdn(ikp) = ptdbt(jk) * ptra(jk) + &
                    (ptrad(jk) * ((ztdn(jk) - ptdbt(jk)) + &
                     ptdbt(jk) * pref(jk) * prdnd(jk))) * zreflect
         prdnd(ikp) = prefd(jk) + ptrad(jk) * ptrad(jk) * &
                      prdnd(jk) * zreflect
      enddo
    
! Up and down-welling fluxes at levels

      do jk = 1,klev+1
         zreflect = 1._jprb / (1._jprb - prdnd(jk) * prupd(jk))
         pfu(jk,kw) = (ptdbt(jk) * prup(jk) + &
                      (ztdn(jk) - ptdbt(jk)) * prupd(jk)) * zreflect
         pfd(jk,kw) = ptdbt(jk) + (ztdn(jk) - ptdbt(jk)+ &
                      ptdbt(jk) * prup(jk) * prdnd(jk)) * zreflect
      enddo

      return
      end



