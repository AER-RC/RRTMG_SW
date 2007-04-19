!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

      subroutine cldprmc_sw(nlayers, inflag, iceflag, liqflag, cldfmc, &
                            ciwpmc, clwpmc, reicmc, relqmc, &
                            taormc, taucmc, ssacmc, asmcmc)

! Purpose: Compute the cloud optical properties for each cloudy layer
! and g-point interval for use by the McICA method.  
! Note: Only inflag = 0 and inflag=2/liqflag=1/iceflag=2,3 are available;
! (Hu & Stamnes, Key, and Fu) are implemented. 

! ------- Modules -------

      use parkind, only : jpim, jprb
      use parrrsw, only : mxlay, ngpt, nstr, jpband, jpb1, jpb2
      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                           extice2, ssaice2, asyice2, &
                           extice3, ssaice3, asyice3, fdlice3, &
                           abari, bbari, cbari, dbari, ebari, fbari
      use rrsw_wvn, only : wavenum1, wavenum2, ngb
      use rrsw_vsn, only : hvrclc, hnamclc

      implicit none


! ------- Declarations -------

! Input

      integer(kind=jpim), intent(in) :: nlayers
      integer(kind=jpim), intent(in) :: inflag
      integer(kind=jpim), intent(in) :: iceflag
      integer(kind=jpim), intent(in) :: liqflag

      real(kind=jprb), intent(in) :: cldfmc(ngpt,mxlay)
      real(kind=jprb), intent(in) :: ciwpmc(ngpt,mxlay)
      real(kind=jprb), intent(in) :: clwpmc(ngpt,mxlay)
      real(kind=jprb), intent(in) :: relqmc(mxlay)
      real(kind=jprb), intent(in) :: reicmc(mxlay)

! Output

      real(kind=jprb), intent(inout) :: taucmc(ngpt,mxlay)
      real(kind=jprb), intent(inout) :: ssacmc(ngpt,mxlay)
      real(kind=jprb), intent(inout) :: asmcmc(ngpt,mxlay)
      real(kind=jprb), intent(out) :: taormc(ngpt,mxlay)


! Local

      integer(kind=jpim) :: ncbands
      integer(kind=jpim) :: ib, lay, istr, index, icx, ig

      real(kind=jprb) :: eps
      real(kind=jprb) :: taucldorig_a, taucloud_a, ssacloud_a, ffp, ffp1, ffpssa
      real(kind=jprb) :: cwp, radice, factor, fint, radliq, dgeice
      real(kind=jprb) :: tauiceorig, scatice, ssaice, tauice, tauliqorig, scatliq, ssaliq, tauliq

      real(kind=jprb) :: fdelta(ngpt)
      real(kind=jprb) :: extcoice(ngpt), gice(ngpt), ssacoice(ngpt), forwice(ngpt)
      real(kind=jprb) :: extcoliq(ngpt), gliq(ngpt), ssacoliq(ngpt), forwliq(ngpt)

! Initialize

      hvrclc = '$Revision$'

! Initialize

      eps = 1.e-06_jprb

! Some of these initializations are done in rrtmg_sw_subcol.F90.
      do lay = 1, nlayers
         do ig = 1, ngpt
            taormc(ig,lay) = taucmc(ig,lay)
!            taucmc(ig,lay) = 0.0_jprb
!            ssacmc(ig,lay) = 0.0_jprb
!            asmcmc(ig,lay) = 0.0_jprb
         enddo
      enddo

! Main layer loop
      do lay = 1, nlayers

! Main g-point interval loop
         do ig = 1, ngpt 
            cwp = ciwpmc(ig,lay) + clwpmc(ig,lay)
            if (cldfmc(ig,lay) .ge. eps .and. &
               (cwp .ge. eps .or. taucmc(ig,lay) .ge. eps)) then

! (inflag=0): Cloud optical properties input directly
               if (inflag .eq. 0) then
! Cloud optical properties already defined in taucmc, ssacmc, asmcmc are unscaled;
! Apply delta-M scaling here (using Henyey-Greenstein approximation)
                  taucldorig_a = taucmc(ig,lay)
                  ffp = asmcmc(ig,lay)**2
                  ffp1 = 1.0_jprb - ffp
                  ffpssa = 1.0_jprb - ffp * ssacmc(ig,lay)
                  ssacloud_a = ffp1 * ssacmc(ig,lay) / ffpssa
                  taucloud_a = ffpssa * taucldorig_a

                  taormc(ig,lay) = taucldorig_a
                  ssacmc(ig,lay) = ssacloud_a
                  taucmc(ig,lay) = taucloud_a
                  asmcmc(ig,lay) = (asmcmc(ig,lay) - ffp) / (ffp1)

               elseif (inflag .eq. 1) then 
                  stop 'INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA'

! (inflag=2): Separate treatement of ice clouds and water clouds.
               elseif (inflag .eq. 2) then       
                  radice = reicmc(lay)

! Calculation of absorption coefficients due to ice clouds.
                  if (ciwpmc(ig,lay) .eq. 0.0) then
                     extcoice(ig) = 0.0_jprb
                     ssacoice(ig) = 0.0_jprb
                     gice(ig)     = 0.0_jprb
                     forwice(ig)  = 0.0_jprb

! (iceflag = 1): 
! Note: This option uses Ebert and Curry approach for all particle sizes similar to
! CAM3 implementation, though this is somewhat unjustified for large ice particles
                  elseif (iceflag .eq. 1) then
                     ib = ngb(ig)
                     if (wavenum2(ib) .gt. 1.43e04_jprb) then
                        icx = 1
                     elseif (wavenum2(ib) .gt. 7.7e03_jprb) then
                        icx = 2
                     elseif (wavenum2(ib) .gt. 5.3e03_jprb) then
                        icx = 3
                     elseif (wavenum2(ib) .gt. 4.0e03_jprb) then
                        icx = 4
                     elseif (wavenum2(ib) .ge. 2.5e03_jprb) then
                        icx = 5
                     endif
                     extcoice(ig) = (abari(icx) + bbari(icx)/radice)
                     ssacoice(ig) = 1._jprb - cbari(icx) - dbari(icx) * radice
                     gice(ig) = ebari(icx) + fbari(icx) * radice
! Check to ensure upper limit of gice is within physical limits for large particles
                     if (gice(ig).ge.1._jprb) gice(ig) = 1._jprb - eps
                     forwice(ig) = gice(ig)*gice(ig)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoice(ig) .lt. 0.0_jprb) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ig) .gt. 1.0_jprb) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ig) .lt. 0.0_jprb) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ig) .gt. 1.0_jprb) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ig) .lt. 0.0_jprb) stop 'ICE ASYM LESS THAN 0.0'

! For iceflag=2 option, combine with iceflag=0 option to handle large particle sizes.
! Use iceflag=2 option for ice particle effective radii from 5.0 to 131.0 microns
! and use iceflag=0 option for ice particles greater than 131.0 microns.
! *** NOTE: Transition between two methods has not been smoothed.

                  elseif (iceflag .eq. 2) then
                     if (radice .lt. 5.0_jprb) stop 'ICE RADIUS OUT OF BOUNDS'
                     if (radice .ge. 5.0_jprb .and. radice .le. 131._jprb) then
                        factor = (radice - 2._jprb)/3._jprb
                        index = int(factor)
                        if (index .eq. 43) index = 42
                        fint = factor - float(index)
                        ib = ngb(ig)
                        extcoice(ig) = extice2(index,ib) + fint * &
                                      (extice2(index+1,ib) -  extice2(index,ib))
                        ssacoice(ig) = ssaice2(index,ib) + fint * &
                                      (ssaice2(index+1,ib) -  ssaice2(index,ib))
                        gice(ig) = asyice2(index,ib) + fint * &
                                      (asyice2(index+1,ib) -  asyice2(index,ib))
                        forwice(ig) = gice(ig)*gice(ig)
! Check to ensure all calculated quantities are within physical limits.
                        if (extcoice(ig) .lt. 0.0_jprb) stop 'ICE EXTINCTION LESS THAN 0.0'
                        if (ssacoice(ig) .gt. 1.0_jprb) stop 'ICE SSA GRTR THAN 1.0'
                        if (ssacoice(ig) .lt. 0.0_jprb) stop 'ICE SSA LESS THAN 0.0'
                        if (gice(ig) .gt. 1.0_jprb) stop 'ICE ASYM GRTR THAN 1.0'
                        if (gice(ig) .lt. 0.0_jprb) stop 'ICE ASYM LESS THAN 0.0'
                     elseif (radice .gt. 131._jprb) then
                        ib = ngb(ig)
                        if (wavenum2(ib) .gt. 1.43e04_jprb) then
                           icx = 1
                        elseif (wavenum2(ib) .gt. 7.7e03_jprb) then
                           icx = 2
                        elseif (wavenum2(ib) .gt. 5.3e03_jprb) then
                           icx = 3
                        elseif (wavenum2(ib) .gt. 4.0e03_jprb) then
                           icx = 4
                        elseif (wavenum2(ib) .ge. 2.5e03_jprb) then
                           icx = 5
                        endif
                        extcoice(ig) = (abari(icx) + bbari(icx)/radice)
                        ssacoice(ig) = 1._jprb - cbari(icx) - dbari(icx) * radice
                        gice(ig) = ebari(icx) + fbari(icx) * radice
! Check to ensure upper limit of gice is within physical limits for large particles
                        if (gice(ig).ge.1.0_jprb) gice(ig) = 1.0_jprb-eps
                        forwice(ig) = gice(ig)*gice(ig)
! Check to ensure all calculated quantities are within physical limits.
                        if (extcoice(ig) .lt. 0.0_jprb) stop 'ICE EXTINCTION LESS THAN 0.0'
                        if (ssacoice(ig) .gt. 1.0_jprb) stop 'ICE SSA GRTR THAN 1.0'
                        if (ssacoice(ig) .lt. 0.0_jprb) stop 'ICE SSA LESS THAN 0.0'
                        if (gice(ig) .gt. 1.0_jprb) stop 'ICE ASYM GRTR THAN 1.0'
                        if (gice(ig) .lt. 0.0_jprb) stop 'ICE ASYM LESS THAN 0.0'
                     endif

! For iceflag=3 option, combine with iceflag=0 option to handle large particle sizes
! Use iceflag=3 option for ice particle effective radii from 5.0 to 91.0 microns
! (generalized effective size, dge, from 8 to 140 microns), and use iceflag=0 option
! for ice particle effective radii greater than 91.0 microns (dge = 140 microns).
! *** NOTE: Fu parameterization requires particle size in generalized effective size.
! *** NOTE: Transition between two methods has not been smoothed. 

                  elseif (iceflag .eq. 3) then
                     dgeice = radice
                     if (dgeice .lt. 5.0_jprb) stop 'ICE EFFECTIVE SIZE OUT OF BOUNDS'
                     if (dgeice .ge. 5.0_jprb .and. dgeice .le. 140._jprb) then
                        factor = (dgeice - 2._jprb)/3._jprb
                        index = int(factor)
                        if (index .eq. 46) index = 45
                        fint = factor - float(index)
                        ib = ngb(ig)
                        extcoice(ig) = extice3(index,ib) + fint * &
                                      (extice3(index+1,ib) - extice3(index,ib))
                        ssacoice(ig) = ssaice3(index,ib) + fint * &
                                      (ssaice3(index+1,ib) - ssaice3(index,ib))
                        gice(ig) = asyice3(index,ib) + fint * &
                                  (asyice3(index+1,ib) - asyice3(index,ib))
                        fdelta(ig) = fdlice3(index,ib) + fint * &
                                    (fdlice3(index+1,ib) - fdlice3(index,ib))
                        if (fdelta(ig) .lt. 0.0_jprb) stop 'FDELTA LESS THAN 0.0'
                        if (fdelta(ig) .gt. 1.0_jprb) stop 'FDELTA GT THAN 1.0'
                        forwice(ig) = fdelta(ig) + 0.5_jprb / ssacoice(ig)
! See Fu 1996 p. 2067 
                        if (forwice(ig) .gt. gice(ig)) forwice(ig) = gice(ig)
! Check to ensure all calculated quantities are within physical limits.  
                        if (extcoice(ig) .lt. 0.0_jprb) stop 'ICE EXTINCTION LESS THAN 0.0'
                        if (ssacoice(ig) .gt. 1.0_jprb) stop 'ICE SSA GRTR THAN 1.0'
                        if (ssacoice(ig) .lt. 0.0_jprb) stop 'ICE SSA LESS THAN 0.0'
                        if (gice(ig) .gt. 1.0_jprb) stop 'ICE ASYM GRTR THAN 1.0'
                        if (gice(ig) .lt. 0.0_jprb) stop 'ICE ASYM LESS THAN 0.0'
                     elseif (dgeice .gt. 140._jprb) then
                        ib = ngb(ig)
                        if (wavenum2(ib) .gt. 1.43e04_jprb) then
                           icx = 1
                        elseif (wavenum2(ib) .gt. 7.7e03_jprb) then
                           icx = 2
                        elseif (wavenum2(ib) .gt. 5.3e03_jprb) then
                           icx = 3
                        elseif (wavenum2(ib) .gt. 4.0e03_jprb) then
                           icx = 4
                        elseif (wavenum2(ib) .ge. 2.5e03_jprb) then
                           icx = 5
                        endif
                        extcoice(ig) = (abari(icx) + bbari(icx)/radice)
                        ssacoice(ig) = 1._jprb - cbari(icx) - dbari(icx) * radice
                        gice(ig) = ebari(icx) + fbari(icx) * radice
! Check to ensure upper limit of gice is within physical limits for large particles
                        if (gice(ig).ge.1._jprb) gice(ig) = 1._jprb - eps
                        forwice(ig) = gice(ig)*gice(ig)
! Check to ensure all calculated quantities are within physical limits.
                        if (extcoice(ig) .lt. 0.0_jprb) stop 'ICE EXTINCTION LESS THAN 0.0'
                        if (ssacoice(ig) .gt. 1.0_jprb) stop 'ICE SSA GRTR THAN 1.0'
                        if (ssacoice(ig) .lt. 0.0_jprb) stop 'ICE SSA LESS THAN 0.0'
                        if (gice(ig) .gt. 1.0_jprb) stop 'ICE ASYM GRTR THAN 1.0'
                        if (gice(ig) .lt. 0.0_jprb) stop 'ICE ASYM LESS THAN 0.0'
                     endif
                  endif

! Calculation of absorption coefficients due to water clouds.
                  if (clwpmc(ig,lay) .eq. 0.0_jprb) then
                     extcoliq(ig) = 0.0_jprb
                     ssacoliq(ig) = 0.0_jprb
                     gliq(ig) = 0.0_jprb
                     forwliq(ig) = 0.0_jprb

                  elseif (liqflag .eq. 1) then
                     radliq = relqmc(lay)
                     if (radliq .lt. 1.5_jprb .or. radliq .gt. 60._jprb) stop &
                        'liquid effective radius out of bounds'
                     index = int(radliq - 1.5_jprb)
                     if (index .eq. 0) index = 1
                     if (index .eq. 58) index = 57
                     fint = radliq - 1.5_jprb - float(index)
                     ib = ngb(ig)
                     extcoliq(ig) = extliq1(index,ib) + fint * &
                                   (extliq1(index+1,ib) - extliq1(index,ib))
                     ssacoliq(ig) = ssaliq1(index,ib) + fint * &
                                   (ssaliq1(index+1,ib) - ssaliq1(index,ib))
                     if (fint .lt. 0._jprb .and. ssacoliq(ig) .gt. 1._jprb) &
                                    ssacoliq(ig) = ssaliq1(index,ib)
                     gliq(ig) = asyliq1(index,ib) + fint * &
                               (asyliq1(index+1,ib) - asyliq1(index,ib))
                     forwliq(ig) = gliq(ig)*gliq(ig)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoliq(ig) .lt. 0.0_jprb) stop 'LIQUID EXTINCTION LESS THAN 0.0'
                     if (ssacoliq(ig) .gt. 1.0_jprb) stop 'LIQUID SSA GRTR THAN 1.0'
                     if (ssacoliq(ig) .lt. 0.0_jprb) stop 'LIQUID SSA LESS THAN 0.0'
                     if (gliq(ig) .gt. 1.0_jprb) stop 'LIQUID ASYM GRTR THAN 1.0'
                     if (gliq(ig) .lt. 0.0_jprb) stop 'LIQUID ASYM LESS THAN 0.0'
                  endif
   
                  tauliqorig = clwpmc(ig,lay) * extcoliq(ig)
                  tauiceorig = ciwpmc(ig,lay) * extcoice(ig)
                  taormc(ig,lay) = tauliqorig + tauiceorig

                  ssaliq = ssacoliq(ig) * (1._jprb - forwliq(ig)) / &
                          (1._jprb - forwliq(ig) * ssacoliq(ig))
                  tauliq = (1._jprb - forwliq(ig) * ssacoliq(ig)) * tauliqorig
                  ssaice = ssacoice(ig) * (1._jprb - forwice(ig)) / &
                          (1._jprb - forwice(ig) * ssacoice(ig))
                  tauice = (1._jprb - forwice(ig) * ssacoice(ig)) * tauiceorig

                  scatliq = ssaliq * tauliq
                  scatice = ssaice * tauice
                  taucmc(ig,lay) = tauliq + tauice

! Ensure non-zero taucmc and scatice
                  if(taucmc(ig,lay).eq.0.) taucmc(ig,lay) = eps
                  if(scatice.eq.0.) scatice = eps

                  ssacmc(ig,lay) = (scatliq + scatice) / taucmc(ig,lay)

                  if (iceflag .eq. 3) then
! In accordance with the 1996 Fu paper, equation A.3, 
! the moments for ice were calculated depending on whether using spheres
! or hexagonal ice crystals.
! Set asymetry parameter to first moment (istr=1)
                     istr = 1
                     asmcmc(ig,lay) = (1.0_jprb/(scatliq+scatice))* &
                        (scatliq*(gliq(ig)**istr - forwliq(ig)) / &
                        (1.0_jprb - forwliq(ig)) + scatice * ((gice(ig)-forwice(ig))/ &
                        (1.0_jprb - forwice(ig)))**istr)

                  else 
! This code is the standard method for delta-m scaling. 
! Set asymetry parameter to first moment (istr=1)
                     istr = 1
                     asmcmc(ig,lay) = (scatliq *  &
                        (gliq(ig)**istr - forwliq(ig)) / &
                        (1.0_jprb - forwliq(ig)) + scatice * (gice(ig)**istr - forwice(ig)) / &
                        (1.0_jprb - forwice(ig)))/(scatliq + scatice)
                  endif 

               endif

            endif

! End g-point interval loop
         enddo

! End layer loop
      enddo

      return
      end














