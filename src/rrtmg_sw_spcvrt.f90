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

      subroutine spcvrt_sw &
            (nlayers, istart, iend, icpr, iout, &
             pavel, tavel, pz, tz, tbound, palbd, palbp, &
             pclfr, ptauc, pasyc, pomgc, ptaucorig, &
             ptaua, pasya, pomga, prmu0, coldry, wkl, adjflux, &
             laytrop, layswtch, laylow, jp, jt, jt1, &
             co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
             fac00, fac01, fac10, fac11, &
             selffac, selffrac, indself, forfac, forfrac, indfor, &
             pbbfd, pbbfu, pbbcd, pbbcu, puvfd, puvcd, pnifd, pnicd, &
             pbbfddir, pbbcddir, puvfddir, puvcddir, pnifddir, pnicddir)

! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
!          using the two-stream method of H. Barker. 
!
! Interface:  *spcvrt_sw* is called from *rrtmg_sw.F90* or rrtmg_sw.1col.F90*
!
! Method:
!    Adapted from two-stream model of H. Barker;
!    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
!        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
!
! Modifications:

! Original: H. Barker
! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
! Revision: Code modified so that delta scaling is not done in cloudy profiles
!           if routine cldprop is used; delta scaling can be applied by swithcing
!           code below if cldprop is not used to get cloud properties. 
!           AER, Jan 2005
! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
!
! ------------------------------------------------------------------

! ------- Modules -------

      use parkind, only : jpim, jprb
      use parrrsw, only : mxlay, nbndsw, ngpt, mxmol, jpband
      use rrsw_vsn, only : hvrspv, hnamspv
      use rrsw_wvn, only : ngc, ngs

      implicit none

! ------- Declarations ------

! Input

      integer(kind=jpim), intent(in) :: nlayers
      integer(kind=jpim), intent(in) :: istart
      integer(kind=jpim), intent(in) :: iend
      integer(kind=jpim), intent(in) :: icpr
      integer(kind=jpim), intent(in) :: iout
      integer(kind=jpim), intent(in) :: laytrop
      integer(kind=jpim), intent(in) :: layswtch
      integer(kind=jpim), intent(in) :: laylow

      integer(kind=jpim), intent(in) :: indfor(mxlay)
      integer(kind=jpim), intent(in) :: indself(mxlay)
      integer(kind=jpim), intent(in) :: jp(mxlay)
      integer(kind=jpim), intent(in) :: jt(mxlay)
      integer(kind=jpim), intent(in) :: jt1(mxlay)

      real(kind=jprb), intent(in) :: pavel(mxlay)
      real(kind=jprb), intent(in) :: tavel(mxlay)
      real(kind=jprb), intent(in) :: pz(0:mxlay)
      real(kind=jprb), intent(in) :: tz(0:mxlay)
      real(kind=jprb), intent(in) :: tbound
      real(kind=jprb), intent(in) :: wkl(mxmol,mxlay)
      real(kind=jprb), intent(in) :: coldry(mxlay)
      real(kind=jprb), intent(in) :: colmol(mxlay)
      real(kind=jprb), intent(in) :: adjflux(jpband)

      real(kind=jprb), intent(in) :: palbd(nbndsw)
      real(kind=jprb), intent(in) :: palbp(nbndsw)
      real(kind=jprb), intent(in) :: prmu0 
      real(kind=jprb), intent(in) :: pclfr(mxlay)
      real(kind=jprb), intent(in) :: ptauc(mxlay,nbndsw)
      real(kind=jprb), intent(in) :: pasyc(mxlay,nbndsw)
      real(kind=jprb), intent(in) :: pomgc(mxlay,nbndsw)
      real(kind=jprb), intent(in) :: ptaucorig(mxlay,nbndsw)
      real(kind=jprb), intent(in) :: ptaua(mxlay,nbndsw)
      real(kind=jprb), intent(in) :: pasya(mxlay,nbndsw)
      real(kind=jprb), intent(in) :: pomga(mxlay,nbndsw)

      real(kind=jprb), intent(in) :: colh2o(mxlay)
      real(kind=jprb), intent(in) :: colco2(mxlay)
      real(kind=jprb), intent(in) :: colch4(mxlay)
      real(kind=jprb), intent(in) :: co2mult(mxlay)
      real(kind=jprb), intent(in) :: colo3(mxlay)
      real(kind=jprb), intent(in) :: colo2(mxlay)
      real(kind=jprb), intent(in) :: coln2o(mxlay)
      real(kind=jprb), intent(in) :: forfac(mxlay)
      real(kind=jprb), intent(in) :: forfrac(mxlay)
      real(kind=jprb), intent(in) :: selffac(mxlay)
      real(kind=jprb), intent(in) :: selffrac(mxlay)
      real(kind=jprb), intent(in) :: fac00(mxlay)
      real(kind=jprb), intent(in) :: fac01(mxlay)
      real(kind=jprb), intent(in) :: fac10(mxlay)
      real(kind=jprb), intent(in) :: fac11(mxlay)

! Output
      real(kind=jprb), intent(out) :: pbbcd(mxlay+1)
      real(kind=jprb), intent(out) :: pbbcu(mxlay+1)
      real(kind=jprb), intent(out) :: pbbfd(mxlay+1)
      real(kind=jprb), intent(out) :: pbbfu(mxlay+1)
      real(kind=jprb), intent(out) :: pbbfddir(mxlay+1)
      real(kind=jprb), intent(out) :: pbbcddir(mxlay+1)

      real(kind=jprb), intent(out) :: puvcd(mxlay+1)
      real(kind=jprb), intent(out) :: puvfd(mxlay+1)
      real(kind=jprb), intent(out) :: puvcddir(mxlay+1)
      real(kind=jprb), intent(out) :: puvfddir(mxlay+1)

      real(kind=jprb), intent(out) :: pnicd(mxlay+1)
      real(kind=jprb), intent(out) :: pnifd(mxlay+1)
      real(kind=jprb), intent(out) :: pnicddir(mxlay+1)
      real(kind=jprb), intent(out) :: pnifddir(mxlay+1)

! Output - inactive
!      real(kind=jprb), intent(out) :: puvcu(mxlay+1)
!      real(kind=jprb), intent(out) :: puvfu(mxlay+1)
!      real(kind=jprb), intent(out) :: pnicu(mxlay+1)
!      real(kind=jprb), intent(out) :: pnifu(mxlay+1)
!      real(kind=jprb), intent(out) :: pvscd(mxlay+1)
!      real(kind=jprb), intent(out) :: pvscu(mxlay+1)
!      real(kind=jprb), intent(out) :: pvsfd(mxlay+1)
!      real(kind=jprb), intent(out) :: pvsfu(mxlay+1)


! Local

      logical :: lrtchkclr(mxlay),lrtchkcld(mxlay)

      integer(kind=jpim)  :: klev
      integer(kind=jpim) :: ib1, ib2, ibm, igt, ikl, ikp, ikx
      integer(kind=jpim) :: iw, jb, jg, jl, jk
!      integer(kind=jpim), parameter :: nuv = ?? 
!      integer(kind=jpim), parameter :: nvs = ?? 


      real(kind=jprb) :: zclear, zcloud
      real(kind=jprb) :: zdbt(mxlay+1), zdbt_nodel(mxlay+1)
      real(kind=jprb) :: zgc(mxlay), zgcc(mxlay), zgco(mxlay)
      real(kind=jprb) :: zomc(mxlay), zomcc(mxlay), zomco(mxlay)
      real(kind=jprb) :: zrdnd(mxlay+1), zrdndc(mxlay+1)
      real(kind=jprb) :: zref(mxlay+1), zrefc(mxlay+1), zrefo(mxlay+1)
      real(kind=jprb) :: zrefd(mxlay+1), zrefdc(mxlay+1), zrefdo(mxlay+1)
      real(kind=jprb) :: zrup(mxlay+1), zrupd(mxlay+1)
      real(kind=jprb) :: zrupc(mxlay+1), zrupdc(mxlay+1)
      real(kind=jprb) :: zs1(mxlay+1)
      real(kind=jprb) :: ztauc(mxlay), ztauo(mxlay)
      real(kind=jprb) :: ztdn(mxlay+1), ztdnd(mxlay+1), ztdbt(mxlay+1)
      real(kind=jprb) :: ztoc(mxlay), ztor(mxlay)
      real(kind=jprb) :: ztra(mxlay+1), ztrac(mxlay+1), ztrao(mxlay+1)
      real(kind=jprb) :: ztrad(mxlay+1), ztradc(mxlay+1), ztrado(mxlay+1)
      real(kind=jprb) :: zdbtc(mxlay+1), ztdbtc(mxlay+1)
      real(kind=jprb) :: zincflx(ngpt), zdbtc_nodel(mxlay+1) 
      real(kind=jprb) :: ztdbt_nodel(mxlay+1), ztdbtc_nodel(mxlay+1)

      real(kind=jprb) :: zdbtmc, zdbtmo, zf, zgw, zreflect
      real(kind=jprb) :: zwf, tauorig, repclc
!     real(kind=jprb) :: zincflux                                   ! inactive

! Arrays from rrtmg_sw_taumoln routines

!      real(kind=jprb) :: ztaug(mxlay,16), ztaur(mxlay,16)
!      real(kind=jprb) :: zsflxzen(16)
      real(kind=jprb) :: ztaug(mxlay,ngpt), ztaur(mxlay,ngpt)
      real(kind=jprb) :: zsflxzen(ngpt)

! Arrays from rrtmg_sw_vrtqdr routine

      real(kind=jprb) :: zcd(mxlay+1,ngpt), zcu(mxlay+1,ngpt)
      real(kind=jprb) :: zfd(mxlay+1,ngpt), zfu(mxlay+1,ngpt)

! Inactive arrays
!     real(kind=jprb) :: zbbcd(mxlay+1), zbbcu(mxlay+1)
!     real(kind=jprb) :: zbbfd(mxlay+1), zbbfu(mxlay+1)
!     real(kind=jprb) :: zbbfddir(mxlay+1), zbbcddir(mxlay+1)

! ------------------------------------------------------------------

! Initializations

      ib1 = istart
      ib2 = iend
      klev = nlayers
      iw = 0
      repclc = 1.e-12_jprb
!      zincflux = 0.0_jprb

      do jk=1,klev+1
         pbbcd(jk)=0._jprb
         pbbcu(jk)=0._jprb
         pbbfd(jk)=0._jprb
         pbbfu(jk)=0._jprb
         pbbcddir(jk)=0._jprb
         pbbfddir(jk)=0._jprb
         puvcd(jk)=0._jprb
         puvfd(jk)=0._jprb
         puvcddir(jk)=0._jprb
         puvfddir(jk)=0._jprb
         pnicd(jk)=0._jprb
         pnifd(jk)=0._jprb
         pnicddir(jk)=0._jprb
         pnifddir(jk)=0._jprb
      enddo


! Calculate the optical depths for gaseous absorption and Rayleigh scattering

      call taumol_sw(klev, &
                     colh2o, colco2, colch4, colo2, colo3, colmol, &
                     laytrop, jp, jt, jt1, &
                     fac00, fac01, fac10, fac11, &
                     selffac, selffrac, indself, forfac, forfrac, indfor, &
                     zsflxzen, ztaug, ztaur)


! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

      jb = ib1-1                  ! ???
      do jb = ib1, ib2
         ibm = jb-15
         igt = ngc(ibm)

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.ibm.ge.2) iw = ngs(ibm-1)

!        do jk=1,klev+1
!           zbbcd(jk)=0.0_jprb
!           zbbcu(jk)=0.0_jprb
!           zbbfd(jk)=0.0_jprb
!           zbbfu(jk)=0.0_jprb
!        enddo

! Top of g-point interval loop within each band (iw is cumulative counter) 
         do jg = 1,igt
            iw = iw+1

! Add adjustment for correct Earth/Sun distance
            zincflx(iw) = adjflux(jb) * zsflxzen(iw) * prmu0
!             zincflux = zincflux + adjflux(jb) * zsflxzen(iw) * prmu0           ! inactive

! Compute layer reflectances and transmittances for direct and diffuse sources, 
! first clear then cloudy

! zrefc(jk)  direct albedo for clear
! zrefo(jk)  direct albedo for cloud
! zrefdc(jk) diffuse albedo for clear
! zrefdo(jk) diffuse albedo for cloud
! ztrac(jk)  direct transmittance for clear
! ztrao(jk)  direct transmittance for cloudy
! ztradc(jk) diffuse transmittance for clear
! ztrado(jk) diffuse transmittance for cloudy
!  
! zref(jk)   direct reflectance
! zrefd(jk)  diffuse reflectance
! ztra(jk)   direct transmittance
! ztrad(jk)  diffuse transmittance
!
! zdbtc(jk)  clear direct beam transmittance
! zdbto(jk)  cloudy direct beam transmittance
! zdbt(jk)   layer mean direct beam transmittance
! ztdbt(jk)  total direct beam transmittance at levels

! Clear-sky    
!   TOA direct beam    
            ztdbtc(1)=1.0_jprb
            ztdbtc_nodel(1)=1.0_jprb
!   Surface values
            zdbtc(klev+1) =0.0_jprb
            ztrac(klev+1) =0.0_jprb
            ztradc(klev+1)=0.0_jprb
            zrefc(klev+1) =palbp(ibm)
            zrefdc(klev+1)=palbd(ibm)
            zrupc(klev+1) =palbp(ibm)
            zrupdc(klev+1)=palbd(ibm)
           
! Total sky    
!   TOA direct beam    
            ztdbt(1)=1.0_jprb
            ztdbt_nodel(1)=1.0_jprb
!   Surface values
            zdbt(klev+1) =0.0_jprb
            ztra(klev+1) =0.0_jprb
            ztrad(klev+1)=0.0_jprb
            zref(klev+1) =palbp(ibm)
            zrefd(klev+1)=palbd(ibm)
            zrup(klev+1) =palbp(ibm)
            zrupd(klev+1)=palbd(ibm)
    
    
! Top of layer loop
            do jk=1,klev

! Note: two-stream calculations proceed from top to bottom; 
!   RRTMG_SW quantities are given bottom to top and are reversed here

               ikl=klev+1-jk

! Set logical flag to do REFTRA calculation
!   Do REFTRA for all clear layers
               lrtchkclr(jk)=.true.

!   Do REFTRA only for cloudy layers in profile, since already done for clear layers
               lrtchkcld(jk)=.false.
               lrtchkcld(jk)=(pclfr(ikl) > repclc)

! Clear-sky optical parameters - this section inactive     
!   Original
!               ztauc(jk) = ztaur(ikl,iw) + ztaug(ikl,iw)
!               zomcc(jk) = ztaur(ikl,iw) / ztauc(jk)
!               zgcc(jk) = 0.0001_jprb
!   Total sky optical parameters        
!               ztauo(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptauc(ikl,ibm)
!               zomco(jk) = ptauc(ikl,ibm) * pomgc(ikl,ibm) + ztaur(ikl,iw)
!               zgco (jk) = (ptauc(ikl,ibm) * pomgc(ikl,ibm) * pasyc(ikl,ibm) + &
!                           ztaur(ikl,iw) * 0.0001_jprb) / zomco(jk)
!               zomco(jk) = zomco(jk) / ztauo(jk)

! Clear-sky optical parameters including aerosols
               ztauc(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptaua(ikl,ibm)
               zomcc(jk) = ztaur(ikl,iw) * 1.0_jprb + ptaua(ikl,ibm) * pomga(ikl,ibm)
               zgcc(jk) = pasya(ikl,ibm) * pomga(ikl,ibm) * ptaua(ikl,ibm) / zomcc(jk)
               zomcc(jk) = zomcc(jk) / ztauc(jk)


! Pre-delta-scaling clear and cloudy direct beam transmittance (must use 'orig', unscaled cloud OD)       
!   \/\/\/ This block of code is only needed for direct beam calculation
!     
               zclear = 1.0_jprb - pclfr(ikl)
               zcloud = pclfr(ikl)
               zdbtmc = exp(-ztauc(jk) / prmu0)
               tauorig = ztauc(jk) + ptaucorig(ikl,ibm)
               zdbtmo = exp(-tauorig / prmu0)
               zdbt_nodel(jk) = zclear * zdbtmc + zcloud * zdbtmo
               ztdbt_nodel(jk+1) = zdbt_nodel(jk) * ztdbt_nodel(jk)
!   Clear-sky 
               zdbtc_nodel(jk) = zdbtmc
               ztdbtc_nodel(jk+1) = zdbtc_nodel(jk) * ztdbtc_nodel(jk)
!   /\/\/\ Above code only needed for direct beam calculation


! Delta scaling - clear   
               zf = zgcc(jk) * zgcc(jk)
               zwf = zomcc(jk) * zf
               ztauc(jk) = (1.0_jprb - zwf) * ztauc(jk)
               zomcc(jk) = (zomcc(jk) - zwf) / (1.0_jprb - zwf)
               zgcc (jk) = (zgcc(jk) - zf) / (1.0_jprb - zf)

! Total sky optical parameters (cloud properties already delta-scaled)
!   Use this code if cloud properties are derived in rrtmg_sw_cldprop       
               if (icpr .ge. 1) then
                  ztauo(jk) = ztauc(jk) + ptauc(ikl,ibm)
                  zomco(jk) = ztauc(jk) * zomcc(jk) + ptauc(ikl,ibm) * pomgc(ikl,ibm) 
                  zgco (jk) = (ptauc(ikl,ibm) * pomgc(ikl,ibm) * pasyc(ikl,ibm) + &
                              ztauc(jk) * zomcc(jk) * zgcc(jk)) / zomco(jk)
                  zomco(jk) = zomco(jk) / ztauo(jk)

! Total sky optical parameters (if cloud properties not delta scaled)
!   Use this code if cloud properties are not derived in rrtmg_sw_cldprop       
               elseif (icpr .eq. 0) then
                  ztauo(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptaua(ikl,ibm) + ptauc(ikl,ibm)
                  zomco(jk) = ptaua(ikl,ibm) * pomga(ikl,ibm) + ptauc(ikl,ibm) * pomgc(ikl,ibm) + &
                              ztaur(ikl,iw) * 1.0_jprb
                  zgco (jk) = (ptauc(ikl,ibm) * pomgc(ikl,ibm) * pasyc(ikl,ibm) + &
                              ptaua(ikl,ibm)*pomga(ikl,ibm)*pasya(ikl,ibm)) / zomco(jk)
                  zomco(jk) = zomco(jk) / ztauo(jk)

! Delta scaling - clouds 
!   Use only if subroutine rrtmg_sw_cldprop is not used to get cloud properties and to apply delta scaling
                  zf = zgco(jk) * zgco(jk)
                  zwf = zomco(jk) * zf
                  ztauo(jk) = (1._jprb - zwf) * ztauo(jk)
                  zomco(jk) = (zomco(jk) - zwf) / (1.0_jprb - zwf)
                  zgco (jk) = (zgco(jk) - zf) / (1.0_jprb - zf)
               endif 

! End of layer loop
            enddo    


! Clear sky reflectivities
            call reftra_sw (klev, &
                            lrtchkclr, zgcc, prmu0, ztauc, zomcc, &
                            zrefc, zrefdc, ztrac, ztradc)

! Total sky reflectivities      
            call reftra_sw (klev, &
                            lrtchkcld, zgco, prmu0, ztauo, zomco, &
                            zrefo, zrefdo, ztrao, ztrado)


            do jk=1,klev

! Combine clear and cloudy contributions for total sky
               ikl = klev+1-jk 
               zclear = 1.0_jprb - pclfr(ikl)
               zcloud = pclfr(ikl)

               zref(jk) = zclear*zrefc(jk) + zcloud*zrefo(jk)
               zrefd(jk)= zclear*zrefdc(jk) + zcloud*zrefdo(jk)
               ztra(jk) = zclear*ztrac(jk) + zcloud*ztrao(jk)
               ztrad(jk)= zclear*ztradc(jk) + zcloud*ztrado(jk)

! Direct beam transmittance        
               zdbtmc = exp(-ztauc(jk)/prmu0)
               zdbtmo = exp(-ztauo(jk)/prmu0)
               zdbt(jk) = zclear*zdbtmc + zcloud*zdbtmo
               ztdbt(jk+1) = zdbt(jk)*ztdbt(jk)
        
! Clear-sky        
               zdbtc(jk) = zdbtmc
               ztdbtc(jk+1) = zdbtc(jk)*ztdbtc(jk)

            enddo           
                 
! Vertical quadrature for clear-sky fluxes

            call vrtqdr_sw (klev, iw, &
                            zrefc, zrefdc, ztrac, ztradc, &
                            zdbtc, zrdndc, zrupc, zrupdc, ztdbtc, &
                            zcd, zcu)
      
! Vertical quadrature for cloudy fluxes

            call vrtqdr_sw (klev, iw, &
                            zref, zrefd, ztra, ztrad, &
                            zdbt, zrdnd, zrup, zrupd, ztdbt, &
                            zfd, zfu)

! Upwelling and downwelling fluxes at levels
!   Two-stream calculations go from top to bottom; 
!   layer indexing is reversed to go bottom to top for output arrays

            do jk=1,klev+1
               ikl=klev+2-jk

! Accumulate spectral fluxes over bands - inactive
!               zbbfu(ikl) = zbbfu(ikl) + zincflx(iw)*zfu(jk,iw)  
!               zbbfd(ikl) = zbbfd(ikl) + zincflx(iw)*zfd(jk,iw)
!               zbbcu(ikl) = zbbcu(ikl) + zincflx(iw)*zcu(jk,iw)
!               zbbcd(ikl) = zbbcd(ikl) + zincflx(iw)*zcd(jk,iw)
!               zbbfddir(ikl) = zbbfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
!               zbbcddir(ikl) = zbbcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)

! Accumulate spectral fluxes over whole spectrum  
               pbbfu(ikl) = pbbfu(ikl) + zincflx(iw)*zfu(jk,iw)
               pbbfd(ikl) = pbbfd(ikl) + zincflx(iw)*zfd(jk,iw)
               pbbcu(ikl) = pbbcu(ikl) + zincflx(iw)*zcu(jk,iw)
               pbbcd(ikl) = pbbcd(ikl) + zincflx(iw)*zcd(jk,iw)
               pbbfddir(ikl) = pbbfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
               pbbcddir(ikl) = pbbcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)

! Accumulate direct fluxes for UV/visible bands
               if (ibm >= 10 .and. ibm <= 13) then
                  puvcd(ikl) = puvcd(ikl) + zincflx(iw)*zcd(jk,iw)
                  puvfd(ikl) = puvfd(ikl) + zincflx(iw)*zfd(jk,iw)
                  puvcddir(ikl) = puvcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
                  puvfddir(ikl) = puvfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
! Accumulate direct fluxes for near-IR bands
               else if (ibm == 14 .or. ibm <= 9) then  
                  pnicd(ikl) = pnicd(ikl) + zincflx(iw)*zcd(jk,iw)
                  pnifd(ikl) = pnifd(ikl) + zincflx(iw)*zfd(jk,iw)
                  pnicddir(ikl) = pnicddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
                  pnifddir(ikl) = pnifddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
               endif

            enddo

! End loop on jg, g-point interval
         enddo             

! End loop on jb, spectral band
      enddo                    

      return
      end

