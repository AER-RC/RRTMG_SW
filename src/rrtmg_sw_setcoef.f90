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

      subroutine setcoef_sw(nlayers, pavel, tavel, pz, tz, tbound, coldry, wkl, &
                            laytrop, layswtch, laylow, jp, jt, jt1, &
                            co2mult, colch4, colco2, colh2o, colmol, coln2o, &
                            colo2, colo3, fac00, fac01, fac10, fac11, &
                            selffac, selffrac, indself, forfac, forfrac, indfor)

! Purpose:  For a given atmosphere, calculate the indices and
! fractions related to the pressure and temperature interpolations.

! Modifications:
! Original: J. Delamere, AER, Inc. (version 2.5, 02/04/01)
! Revised: Rewritten and adapted to ECMWF F90, JJMorcrette 030224
! Revised: For uniform rrtmg formatting, MJIacono, Jul 2006

! ------- Modules -------

      use parkind, only : jpim, jprb
      use parrrsw, only : mxlay, mxmol
      use rrsw_ref, only : pref, preflog, tref
      use rrsw_vsn, only : hvrset, hnamset

      implicit none

! ------ Declarations -------

! Input 

      integer(kind=jpim), intent(in) :: nlayers

      real(kind=jprb), intent(in) :: pavel(mxlay)
      real(kind=jprb), intent(in) :: tavel(mxlay)
      real(kind=jprb), intent(in) :: pz(0:mxlay)
      real(kind=jprb), intent(in) :: tz(0:mxlay)
      real(kind=jprb), intent(in) :: coldry(mxlay)
      real(kind=jprb), intent(in) :: wkl(mxmol,mxlay)

! Output

      integer(kind=jpim), intent(out) :: laytrop
      integer(kind=jpim), intent(out) :: layswtch
      integer(kind=jpim), intent(out) :: laylow
      integer(kind=jpim), intent(out) :: jp(mxlay)
      integer(kind=jpim), intent(out) :: jt(mxlay)
      integer(kind=jpim), intent(out) :: jt1(mxlay)
      integer(kind=jpim), intent(out) :: indself(mxlay)
      integer(kind=jpim), intent(out) :: indfor(mxlay)

      real(kind=jprb), intent(out) :: colh2o(mxlay)
      real(kind=jprb), intent(out) :: colco2(mxlay)
      real(kind=jprb), intent(out) :: coln2o(mxlay)
      real(kind=jprb), intent(out) :: colo2(mxlay)
      real(kind=jprb), intent(out) :: colo3(mxlay)
      real(kind=jprb), intent(out) :: colch4(mxlay)
      real(kind=jprb), intent(out) :: colmol(mxlay)
      real(kind=jprb), intent(out) :: co2mult(mxlay)

      real(kind=jprb), intent(out) :: fac00(mxlay)
      real(kind=jprb), intent(out) :: fac01(mxlay)
      real(kind=jprb), intent(out) :: fac10(mxlay)
      real(kind=jprb), intent(out) :: fac11(mxlay)
      real(kind=jprb), intent(out) :: forfac(mxlay)
      real(kind=jprb), intent(out) :: forfrac(mxlay)
      real(kind=jprb), intent(out) :: selffac(mxlay)
      real(kind=jprb), intent(out) :: selffrac(mxlay)
      real(kind=jprb), intent(out) :: tbound

! Local

      integer(kind=jpim) :: indbound
      integer(kind=jpim) :: indlev0
      integer(kind=jpim) :: lay
      integer(kind=jpim) :: jp1

      real(kind=jprb) :: stpfac
      real(kind=jprb) :: tbndfrac
      real(kind=jprb) :: t0frac
      real(kind=jprb) :: plog
      real(kind=jprb) :: fp
      real(kind=jprb) :: ft
      real(kind=jprb) :: ft1
      real(kind=jprb) :: water
      real(kind=jprb) :: scalefac
      real(kind=jprb) :: factor
      real(kind=jprb) :: co2reg
      real(kind=jprb) :: compfp


! Initializations
      stpfac = 296._jprb/1013._jprb

      indbound = tbound - 159._jprb
      tbndfrac = tbound - int(tbound)
      indlev0  = tz(0) - 159._jprb
      t0frac   = tz(0) - int(tz(0))

      laytrop  = 0
      layswtch = 0
      laylow   = 0

! Begin layer loop
      do lay = 1, nlayers
! Find the two reference pressures on either side of the
! layer pressure.  Store them in JP and JP1.  Store in FP the
! fraction of the difference (in ln(pressure)) between these
! two values that the layer pressure lies.

         plog = log(pavel(lay))
         jp(lay) = int(36._jprb - 5*(plog+0.04_jprb))
         if (jp(lay) .lt. 1) then
            jp(lay) = 1
         elseif (jp(lay) .gt. 58) then
            jp(lay) = 58
         endif
         jp1 = jp(lay) + 1
         fp = 5._jprb * (preflog(jp(lay)) - plog)

! Determine, for each reference pressure (JP and JP1), which
! reference temperature (these are different for each  
! reference pressure) is nearest the layer temperature but does
! not exceed it.  Store these indices in JT and JT1, resp.
! Store in FT (resp. FT1) the fraction of the way between JT
! (JT1) and the next highest reference temperature that the 
! layer temperature falls.

         jt(lay) = int(3._jprb + (tavel(lay)-tref(jp(lay)))/15._jprb)
         if (jt(lay) .lt. 1) then
            jt(lay) = 1
         elseif (jt(lay) .gt. 4) then
            jt(lay) = 4
         endif
         ft = ((tavel(lay)-tref(jp(lay)))/15._jprb) - float(jt(lay)-3)
         jt1(lay) = int(3._jprb + (tavel(lay)-tref(jp1))/15._jprb)
         if (jt1(lay) .lt. 1) then
            jt1(lay) = 1
         elseif (jt1(lay) .gt. 4) then
            jt1(lay) = 4
         endif
         ft1 = ((tavel(lay)-tref(jp1))/15._jprb) - float(jt1(lay)-3)

         water = wkl(1,lay)/coldry(lay)
         scalefac = pavel(lay) * stpfac / tavel(lay)

! If the pressure is less than ~100mb, perform a different
! set of species interpolations.

         if (plog .le. 4.56_jprb) go to 5300
         laytrop =  laytrop + 1
         if (plog .ge. 6.62_jprb) laylow = laylow + 1

! Set up factors needed to separately include the water vapor
! foreign-continuum in the calculation of absorption coefficient.

         forfac(lay) = scalefac / (1.+water)
         factor = (332.0_jprb-tavel(lay))/36.0_jprb
         indfor(lay) = min(2, max(1, int(factor)))
         forfrac(lay) = factor - float(indfor(lay))

! Set up factors needed to separately include the water vapor
! self-continuum in the calculation of absorption coefficient.

         selffac(lay) = water * forfac(lay)
         factor = (tavel(lay)-188.0_jprb)/7.2_jprb
         indself(lay) = min(9, max(1, int(factor)-7))
         selffrac(lay) = factor - float(indself(lay) + 7)

! Calculate needed column amounts.

         colh2o(lay) = 1.e-20_jprb * wkl(1,lay)
         colco2(lay) = 1.e-20_jprb * wkl(2,lay)
         colo3(lay) = 1.e-20_jprb * wkl(3,lay)
!           colo3(lay) = 0._jprb
!           colo3(lay) = colo3(lay)/1.16_jprb
         coln2o(lay) = 1.e-20_jprb * wkl(4,lay)
         colch4(lay) = 1.e-20_jprb * wkl(6,lay)
         colo2(lay) = 1.e-20_jprb * wkl(7,lay)
         colmol(lay) = 1.e-20_jprb * coldry(lay) + colh2o(lay)
!           colco2(lay) = 0._jprb
!           colo3(lay) = 0._jprb
!           coln2o(lay) = 0._jprb
!           colch4(lay) = 0._jprb
!           colo2(lay) = 0._jprb
!           colmol(lay) = 0._jprb
         if (colco2(lay) .eq. 0._jprb) colco2(lay) = 1.e-32_jprb * coldry(lay)
         if (coln2o(lay) .eq. 0._jprb) coln2o(lay) = 1.e-32_jprb * coldry(lay)
         if (colch4(lay) .eq. 0._jprb) colch4(lay) = 1.e-32_jprb * coldry(lay)
         if (colo2(lay) .eq. 0._jprb) colo2(lay) = 1.e-32_jprb * coldry(lay)
! Using E = 1334.2 cm-1.
         co2reg = 3.55e-24_jprb * coldry(lay)
         co2mult(lay)= (colco2(lay) - co2reg) * &
               272.63_jprb*exp(-1919.4_jprb/tavel(lay))/(8.7604e-4_jprb*tavel(lay))
         goto 5400

! Above laytrop.
 5300    continue

! Set up factors needed to separately include the water vapor
! foreign-continuum in the calculation of absorption coefficient.

         forfac(lay) = scalefac / (1.+water)
         factor = (tavel(lay)-188.0_jprb)/36.0_jprb
         indfor(lay) = 3
         forfrac(lay) = factor - 1.0_jprb

! Calculate needed column amounts.

         colh2o(lay) = 1.e-20_jprb * wkl(1,lay)
         colco2(lay) = 1.e-20_jprb * wkl(2,lay)
         colo3(lay)  = 1.e-20_jprb * wkl(3,lay)
         coln2o(lay) = 1.e-20_jprb * wkl(4,lay)
         colch4(lay) = 1.e-20_jprb * wkl(6,lay)
         colo2(lay)  = 1.e-20_jprb * wkl(7,lay)
         colmol(lay) = 1.e-20_jprb * coldry(lay) + colh2o(lay)
         if (colco2(lay) .eq. 0._jprb) colco2(lay) = 1.e-32_jprb * coldry(lay)
         if (coln2o(lay) .eq. 0._jprb) coln2o(lay) = 1.e-32_jprb * coldry(lay)
         if (colch4(lay) .eq. 0._jprb) colch4(lay) = 1.e-32_jprb * coldry(lay)
         if (colo2(lay)  .eq. 0._jprb) colo2(lay)  = 1.e-32_jprb * coldry(lay)
         co2reg = 3.55e-24_jprb * coldry(lay)
         co2mult(lay)= (colco2(lay) - co2reg) * &
               272.63_jprb*exp(-1919.4_jprb/tavel(lay))/(8.7604e-4_jprb*tavel(lay))

         selffac(lay) = 0._jprb
         selffrac(lay)= 0._jprb
         indself(lay) = 0

 5400    continue

! We have now isolated the layer ln pressure and temperature,
! between two reference pressures and two reference temperatures 
! (for each reference pressure).  We multiply the pressure 
! fraction FP with the appropriate temperature fractions to get 
! the factors that will be needed for the interpolation that yields
! the optical depths (performed in routines TAUGBn for band n).

         compfp = 1._jprb - fp
         fac10(lay) = compfp * ft
         fac00(lay) = compfp * (1._jprb - ft)
         fac11(lay) = fp * ft1
         fac01(lay) = fp * (1._jprb - ft1)

! End layer loop
      enddo

      return
      end

!***************************************************************************
      subroutine swatmref
!***************************************************************************

      use parkind, only : jpim, jprb
      use rrsw_ref, only : pref, preflog, tref

      implicit none
      save
 
! These pressures are chosen such that the ln of the first pressure
! has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
! each subsequent ln(pressure) differs from the previous one by 0.2.

      pref(:) = (/ &
          1.05363e+03_jprb,8.62642e+02_jprb,7.06272e+02_jprb,5.78246e+02_jprb,4.73428e+02_jprb, &
          3.87610e+02_jprb,3.17348e+02_jprb,2.59823e+02_jprb,2.12725e+02_jprb,1.74164e+02_jprb, &
          1.42594e+02_jprb,1.16746e+02_jprb,9.55835e+01_jprb,7.82571e+01_jprb,6.40715e+01_jprb, &
          5.24573e+01_jprb,4.29484e+01_jprb,3.51632e+01_jprb,2.87892e+01_jprb,2.35706e+01_jprb, &
          1.92980e+01_jprb,1.57998e+01_jprb,1.29358e+01_jprb,1.05910e+01_jprb,8.67114e+00_jprb, &
          7.09933e+00_jprb,5.81244e+00_jprb,4.75882e+00_jprb,3.89619e+00_jprb,3.18993e+00_jprb, &
          2.61170e+00_jprb,2.13828e+00_jprb,1.75067e+00_jprb,1.43333e+00_jprb,1.17351e+00_jprb, &
          9.60789e-01_jprb,7.86628e-01_jprb,6.44036e-01_jprb,5.27292e-01_jprb,4.31710e-01_jprb, &
          3.53455e-01_jprb,2.89384e-01_jprb,2.36928e-01_jprb,1.93980e-01_jprb,1.58817e-01_jprb, &
          1.30029e-01_jprb,1.06458e-01_jprb,8.71608e-02_jprb,7.13612e-02_jprb,5.84256e-02_jprb, &
          4.78349e-02_jprb,3.91639e-02_jprb,3.20647e-02_jprb,2.62523e-02_jprb,2.14936e-02_jprb, &
          1.75975e-02_jprb,1.44076e-02_jprb,1.17959e-02_jprb,9.65769e-03_jprb /)

      preflog(:) = (/ &
           6.9600e+00_jprb, 6.7600e+00_jprb, 6.5600e+00_jprb, 6.3600e+00_jprb, 6.1600e+00_jprb, &
           5.9600e+00_jprb, 5.7600e+00_jprb, 5.5600e+00_jprb, 5.3600e+00_jprb, 5.1600e+00_jprb, &
           4.9600e+00_jprb, 4.7600e+00_jprb, 4.5600e+00_jprb, 4.3600e+00_jprb, 4.1600e+00_jprb, &
           3.9600e+00_jprb, 3.7600e+00_jprb, 3.5600e+00_jprb, 3.3600e+00_jprb, 3.1600e+00_jprb, &
           2.9600e+00_jprb, 2.7600e+00_jprb, 2.5600e+00_jprb, 2.3600e+00_jprb, 2.1600e+00_jprb, &
           1.9600e+00_jprb, 1.7600e+00_jprb, 1.5600e+00_jprb, 1.3600e+00_jprb, 1.1600e+00_jprb, &
           9.6000e-01_jprb, 7.6000e-01_jprb, 5.6000e-01_jprb, 3.6000e-01_jprb, 1.6000e-01_jprb, &
          -4.0000e-02_jprb,-2.4000e-01_jprb,-4.4000e-01_jprb,-6.4000e-01_jprb,-8.4000e-01_jprb, &
          -1.0400e+00_jprb,-1.2400e+00_jprb,-1.4400e+00_jprb,-1.6400e+00_jprb,-1.8400e+00_jprb, &
          -2.0400e+00_jprb,-2.2400e+00_jprb,-2.4400e+00_jprb,-2.6400e+00_jprb,-2.8400e+00_jprb, &
          -3.0400e+00_jprb,-3.2400e+00_jprb,-3.4400e+00_jprb,-3.6400e+00_jprb,-3.8400e+00_jprb, &
          -4.0400e+00_jprb,-4.2400e+00_jprb,-4.4400e+00_jprb,-4.6400e+00_jprb /)

! These are the temperatures associated with the respective 
! pressures for the MLS standard atmosphere. 

      tref(:) = (/ &
           2.9420e+02_jprb, 2.8799e+02_jprb, 2.7894e+02_jprb, 2.6925e+02_jprb, 2.5983e+02_jprb, &
           2.5017e+02_jprb, 2.4077e+02_jprb, 2.3179e+02_jprb, 2.2306e+02_jprb, 2.1578e+02_jprb, &
           2.1570e+02_jprb, 2.1570e+02_jprb, 2.1570e+02_jprb, 2.1706e+02_jprb, 2.1858e+02_jprb, &
           2.2018e+02_jprb, 2.2174e+02_jprb, 2.2328e+02_jprb, 2.2479e+02_jprb, 2.2655e+02_jprb, &
           2.2834e+02_jprb, 2.3113e+02_jprb, 2.3401e+02_jprb, 2.3703e+02_jprb, 2.4022e+02_jprb, &
           2.4371e+02_jprb, 2.4726e+02_jprb, 2.5085e+02_jprb, 2.5457e+02_jprb, 2.5832e+02_jprb, &
           2.6216e+02_jprb, 2.6606e+02_jprb, 2.6999e+02_jprb, 2.7340e+02_jprb, 2.7536e+02_jprb, &
           2.7568e+02_jprb, 2.7372e+02_jprb, 2.7163e+02_jprb, 2.6955e+02_jprb, 2.6593e+02_jprb, &
           2.6211e+02_jprb, 2.5828e+02_jprb, 2.5360e+02_jprb, 2.4854e+02_jprb, 2.4348e+02_jprb, & 
           2.3809e+02_jprb, 2.3206e+02_jprb, 2.2603e+02_jprb, 2.2000e+02_jprb, 2.1435e+02_jprb, &
           2.0887e+02_jprb, 2.0340e+02_jprb, 1.9792e+02_jprb, 1.9290e+02_jprb, 1.8809e+02_jprb, &
           1.8329e+02_jprb, 1.7849e+02_jprb, 1.7394e+02_jprb, 1.7212e+02_jprb /)

      return
      end

