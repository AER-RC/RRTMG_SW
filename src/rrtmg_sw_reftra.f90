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

      subroutine reftra_sw(nlayers, lrtchk, pgg, prmuz, ptau, pw, &
                           pref, prefd, ptra, ptrad)
  
! Purpose: computes the reflectivity and transmissivity of a clear or 
!   cloudy layer using a choice of various approximations.
!
! Interface:  *rrtmg_sw_reftra* is called by *rrtmg_sw_spcvrt*
!
! Description:
! explicit arguments :
! --------------------
! inputs
! ------ 
!      lrtchk  = .t. if cloudy
!              = .f. if clear-sky
!      pgg     = assymetry factor
!      prmuz   = cosine solar zenith angle
!      ptau    = optical thickness
!      pw      = single scattering albedo
!
! outputs
! -------
!      pref    : collimated beam reflectivity
!      prefd   : diffuse beam reflectivity 
!      ptra    : collimated beam transmissivity
!      ptrad   : diffuse beam transmissivity
!
!
! Method:
! -------
!      standard delta-eddington, p.i.f.m., or d.o.m. layer calculations.
!      kmodts  = 1 eddington (joseph et al., 1976)
!              = 2 pifm (zdunkowski et al., 1980)
!              = 3 discrete ordinates (liou, 1973)
!
!
! Modifications:
! --------------
! Original: J-JMorcrette, ECMWF, Feb 2003
! Revised for F90 reformatting: MJIacono, AER, Jul 2006
!
! ------------------------------------------------------------------

! ------- Modules -------

      use parkind, only : jpim, jprb
      use parrrsw, only : mxlay
      use rrsw_vsn, only : hvrrft, hnamrft

      implicit none

! ------- Declarations ------

! Input

      integer(kind=jpim), intent(in) :: nlayers

      logical, intent(in) :: lrtchk(mxlay)

      real(kind=jprb), intent(in) :: pgg(mxlay)
      real(kind=jprb), intent(in) :: prmuz
      real(kind=jprb), intent(in) :: ptau(mxlay)
      real(kind=jprb), intent(in) :: pw(mxlay)

! Output

      real(kind=jprb), intent(out) :: pref(mxlay)
      real(kind=jprb), intent(out) :: prefd(mxlay)
      real(kind=jprb), intent(out) :: ptra(mxlay)
      real(kind=jprb), intent(out) :: ptrad(mxlay)

! Local

      integer(kind=jpim) :: jk, jl, kmodts

      real(kind=jprb) :: za, za1, za2
      real(kind=jprb) :: zbeta, zdend, zdenr, zdent
      real(kind=jprb) :: ze1, ze2, zem1, zem2, zemm, zep1, zep2
      real(kind=jprb) :: zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4, zgt
      real(kind=jprb) :: zr1, zr2, zr3, zr4, zr5
      real(kind=jprb) :: zrk, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
      real(kind=jprb) :: zsr3, zt1, zt2, zt3, zt4, zt5, zto1
      real(kind=jprb) :: zw, zwcrit, zwo

!     ------------------------------------------------------------------

! Initialize

      hvrrft = '$Revision$'

      zsr3=sqrt(3._jprb)
      zwcrit=0.9995_jprb
      kmodts=2

      do jk=1, nlayers
         if (.not.lrtchk(jk)) then
            pref(jk) =0._jprb
            ptra(jk) =1._jprb
            prefd(jk)=0._jprb
            ptrad(jk)=1._jprb
         else
            zto1=ptau(jk)
            zw  =pw(jk)
            zg  =pgg(jk)  

! General two-stream expressions

            zg3= 3._jprb * zg
            if (kmodts == 1) then
               zgamma1= (7._jprb - zw * (4._jprb + zg3)) * 0.25_jprb
               zgamma2=-(1._jprb - zw * (4._jprb - zg3)) * 0.25_jprb
               zgamma3= (2._jprb - zg3 * prmuz ) * 0.25_jprb
            else if (kmodts == 2) then  
               zgamma1= (8._jprb - zw * (5._jprb + zg3)) * 0.25_jprb
               zgamma2=  3._jprb *(zw * (1._jprb - zg )) * 0.25_jprb
               zgamma3= (2._jprb - zg3 * prmuz ) * 0.25_jprb
            else if (kmodts == 3) then  
               zgamma1= zsr3 * (2._jprb - zw * (1._jprb + zg)) * 0.5_jprb
               zgamma2= zsr3 * zw * (1._jprb - zg ) * 0.5_jprb
               zgamma3= (1._jprb - zsr3 * zg * prmuz ) * 0.5_jprb
            end if
            zgamma4= 1._jprb - zgamma3
    
! Recompute original s.s.a. to test for conservative solution

            zwo= zw / (1._jprb - (1._jprb - zw) * (zg / (1._jprb - zg))**2)
    
            if (zwo >= zwcrit) then
! Conservative scattering

               za  = zgamma1 * prmuz 
               za1 = za - zgamma3
               zgt = zgamma1 * zto1
        
! Homogeneous reflectance and transmittance,
! collimated beam

               ze1 = min ( zto1 / prmuz , 500._jprb)
               ze2 = exp ( - ze1 )

               pref(jk) = (zgt - za1 * (1._jprb - ze2)) / (1._jprb + zgt)
               ptra(jk) = 1._jprb - pref(jk)

! isotropic incidence

               prefd(jk) = zgt / (1._jprb + zgt)
               ptrad(jk) = 1._jprb - prefd(jk)        

            else
! Non-conservative scattering

               za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
               za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
               zrk = sqrt ( zgamma1**2 - zgamma2**2)
               zrp = zrk * prmuz               
               zrp1 = 1._jprb + zrp
               zrm1 = 1._jprb - zrp
               zrk2 = 2._jprb * zrk
               zrpp = 1._jprb - zrp*zrp
               zrkg = zrk + zgamma1
               zr1  = zrm1 * (za2 + zrk * zgamma3)
               zr2  = zrp1 * (za2 - zrk * zgamma3)
               zr3  = zrk2 * (zgamma3 - za2 * prmuz )
               zr4  = zrpp * zrkg
               zr5  = zrpp * (zrk - zgamma1)
               zt1  = zrp1 * (za1 + zrk * zgamma4)
               zt2  = zrm1 * (za1 - zrk * zgamma4)
               zt3  = zrk2 * (zgamma4 + za1 * prmuz )
               zt4  = zr4
               zt5  = zr5
               zbeta = - zr5 / zr4
        
! Homogeneous reflectance and transmittance

               ze1 = min ( zrk * zto1, 500._jprb)
               ze2 = min ( zto1 / prmuz , 500._jprb)

! jjm - original
!              zep1 = exp( ze1 )
!              zem1 = exp(-ze1 )
!              zep2 = exp( ze2 )
!              zem2 = exp(-ze2 )
!
! mji - optimization
               zep1 = exp( ze1 )
               zem1 = 1._jprb / zep1
               zep2 = exp( ze2 )
               zem2 = 1._jprb / zep2
!
        
! collimated beam

               zdenr = zr4*zep1 + zr5*zem1
               pref(jk) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
               zdent = zt4*zep1 + zt5*zem1
               ptra(jk) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent

! diffuse beam

               zemm = zem1*zem1
               zdend = 1._jprb / ( (1._jprb - zbeta*zemm ) * zrkg)
               prefd(jk) =  zgamma2 * (1._jprb - zemm) * zdend
               ptrad(jk) =  zrk2*zem1*zdend

            endif

         endif         

      enddo    

      return
      end


