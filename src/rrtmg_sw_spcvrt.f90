!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

      module rrtmg_sw_spcvrt

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use parkind, only : jpim, jprb
      use parrrsw, only : nbndsw, ngptsw, mxmol, jpband
      use rrsw_tbl, only : tblint, bpade, od_lo, exp_tbl
      use rrsw_vsn, only : hvrspv, hnamspv
      use rrsw_wvn, only : ngc, ngs
      use rrtmg_sw_reftra, only: reftra_sw
      use rrtmg_sw_taumol, only: taumol_sw
      use rrtmg_sw_vrtqdr, only: vrtqdr_sw

      implicit none

      contains

! ---------------------------------------------------------------------------
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
! ---------------------------------------------------------------------------
!
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
!
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
! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
!           Aug 2007 
!
! ------------------------------------------------------------------

! ------- Declarations ------

! -------- Input -------

      integer(kind=jpim), intent(in) :: nlayers
      integer(kind=jpim), intent(in) :: istart
      integer(kind=jpim), intent(in) :: iend
      integer(kind=jpim), intent(in) :: icpr
      integer(kind=jpim), intent(in) :: iout
      integer(kind=jpim), intent(in) :: laytrop
      integer(kind=jpim), intent(in) :: layswtch
      integer(kind=jpim), intent(in) :: laylow

      integer(kind=jpim), intent(in) :: indfor(:)
                                                                 !   Dimensions: (nlayers)
      integer(kind=jpim), intent(in) :: indself(:)
                                                                 !   Dimensions: (nlayers)
      integer(kind=jpim), intent(in) :: jp(:)
                                                                 !   Dimensions: (nlayers)
      integer(kind=jpim), intent(in) :: jt(:)
                                                                 !   Dimensions: (nlayers)
      integer(kind=jpim), intent(in) :: jt1(:)
                                                                 !   Dimensions: (nlayers)

      real(kind=jprb), intent(in) :: pavel(:)                    ! layer pressure (hPa, mb) 
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: tavel(:)                    ! layer temperature (K)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: pz(0:)                      ! level (interface) pressure (hPa, mb)
                                                                 !   Dimensions: (0:nlayers)
      real(kind=jprb), intent(in) :: tz(0:)                      ! level temperatures (hPa, mb)
                                                                 !   Dimensions: (0:nlayers)
      real(kind=jprb), intent(in) :: tbound                      ! surface temperature (K)
      real(kind=jprb), intent(in) :: wkl(:,:)                    ! molecular amounts (mol/cm2) 
                                                                 !   Dimensions: (mxmol,nlayers)
      real(kind=jprb), intent(in) :: coldry(:)                   ! dry air column density (mol/cm2)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: colmol(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: adjflux(:)                  ! Earth/Sun distance adjustment
                                                                 !   Dimensions: (jpband)

      real(kind=jprb), intent(in) :: palbd(:)                    ! surface albedo (diffuse)
                                                                 !   Dimensions: (nbndsw)
      real(kind=jprb), intent(in) :: palbp(:)                    ! surface albedo (direct)
                                                                 !   Dimensions: (nbndsw)
      real(kind=jprb), intent(in) :: prmu0                       ! cosine of solar zenith angle
      real(kind=jprb), intent(in) :: pclfr(:)                    ! cloud fraction
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: ptauc(:,:)                  ! cloud optical depth
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=jprb), intent(in) :: pasyc(:,:)                  ! cloud asymmetry parameter
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=jprb), intent(in) :: pomgc(:,:)                  ! cloud single scattering albedo
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=jprb), intent(in) :: ptaucorig(:,:)              ! cloud optical depth, non-delta scaled
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=jprb), intent(in) :: ptaua(:,:)                  ! aerosol optical depth
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=jprb), intent(in) :: pasya(:,:)                  ! aerosol asymmetry parameter
                                                                 !   Dimensions: (nlayers,nbndsw)
      real(kind=jprb), intent(in) :: pomga(:,:)                  ! aerosol single scattering albedo
                                                                 !   Dimensions: (nlayers,nbndsw)

      real(kind=jprb), intent(in) :: colh2o(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: colco2(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: colch4(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: co2mult(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: colo3(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: colo2(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: coln2o(:)
                                                                 !   Dimensions: (nlayers)

      real(kind=jprb), intent(in) :: forfac(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: forfrac(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: selffac(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: selffrac(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: fac00(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: fac01(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: fac10(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: fac11(:)
                                                                 !   Dimensions: (nlayers)

! ------- Output -------
                                                                 !   All Dimensions: (nlayers+1)
      real(kind=jprb), intent(out) :: pbbcd(:)
      real(kind=jprb), intent(out) :: pbbcu(:)
      real(kind=jprb), intent(out) :: pbbfd(:)
      real(kind=jprb), intent(out) :: pbbfu(:)
      real(kind=jprb), intent(out) :: pbbfddir(:)
      real(kind=jprb), intent(out) :: pbbcddir(:)

      real(kind=jprb), intent(out) :: puvcd(:)
      real(kind=jprb), intent(out) :: puvfd(:)
      real(kind=jprb), intent(out) :: puvcddir(:)
      real(kind=jprb), intent(out) :: puvfddir(:)

      real(kind=jprb), intent(out) :: pnicd(:)
      real(kind=jprb), intent(out) :: pnifd(:)
      real(kind=jprb), intent(out) :: pnicddir(:)
      real(kind=jprb), intent(out) :: pnifddir(:)

! Output - inactive                                              !   All Dimensions: (nlayers+1)
!      real(kind=jprb), intent(out) :: puvcu(:)
!      real(kind=jprb), intent(out) :: puvfu(:)
!      real(kind=jprb), intent(out) :: pnicu(:)
!      real(kind=jprb), intent(out) :: pnifu(:)
!      real(kind=jprb), intent(out) :: pvscd(:)
!      real(kind=jprb), intent(out) :: pvscu(:)
!      real(kind=jprb), intent(out) :: pvsfd(:)
!      real(kind=jprb), intent(out) :: pvsfu(:)


! ------- Local -------

      logical :: lrtchkclr(nlayers),lrtchkcld(nlayers)

      integer(kind=jpim)  :: klev
      integer(kind=jpim) :: ib1, ib2, ibm, igt, ikl, ikp, ikx
      integer(kind=jpim) :: iw, jb, jg, jl, jk
!      integer(kind=jpim), parameter :: nuv = ?? 
!      integer(kind=jpim), parameter :: nvs = ?? 
      integer(kind=jpim) :: itind

      real(kind=jprb) :: tblind, ze1
      real(kind=jprb) :: zclear, zcloud
      real(kind=jprb) :: zdbt(nlayers+1), zdbt_nodel(nlayers+1)
      real(kind=jprb) :: zgc(nlayers), zgcc(nlayers), zgco(nlayers)
      real(kind=jprb) :: zomc(nlayers), zomcc(nlayers), zomco(nlayers)
      real(kind=jprb) :: zrdnd(nlayers+1), zrdndc(nlayers+1)
      real(kind=jprb) :: zref(nlayers+1), zrefc(nlayers+1), zrefo(nlayers+1)
      real(kind=jprb) :: zrefd(nlayers+1), zrefdc(nlayers+1), zrefdo(nlayers+1)
      real(kind=jprb) :: zrup(nlayers+1), zrupd(nlayers+1)
      real(kind=jprb) :: zrupc(nlayers+1), zrupdc(nlayers+1)
      real(kind=jprb) :: zs1(nlayers+1)
      real(kind=jprb) :: ztauc(nlayers), ztauo(nlayers)
      real(kind=jprb) :: ztdn(nlayers+1), ztdnd(nlayers+1), ztdbt(nlayers+1)
      real(kind=jprb) :: ztoc(nlayers), ztor(nlayers)
      real(kind=jprb) :: ztra(nlayers+1), ztrac(nlayers+1), ztrao(nlayers+1)
      real(kind=jprb) :: ztrad(nlayers+1), ztradc(nlayers+1), ztrado(nlayers+1)
      real(kind=jprb) :: zdbtc(nlayers+1), ztdbtc(nlayers+1)
      real(kind=jprb) :: zincflx(ngptsw), zdbtc_nodel(nlayers+1) 
      real(kind=jprb) :: ztdbt_nodel(nlayers+1), ztdbtc_nodel(nlayers+1)

      real(kind=jprb) :: zdbtmc, zdbtmo, zf, zgw, zreflect
      real(kind=jprb) :: zwf, tauorig, repclc
!     real(kind=jprb) :: zincflux                                   ! inactive

! Arrays from rrtmg_sw_taumoln routines

!      real(kind=jprb) :: ztaug(nlayers,16), ztaur(nlayers,16)
!      real(kind=jprb) :: zsflxzen(16)
      real(kind=jprb) :: ztaug(nlayers,ngptsw), ztaur(nlayers,ngptsw)
      real(kind=jprb) :: zsflxzen(ngptsw)

! Arrays from rrtmg_sw_vrtqdr routine

      real(kind=jprb) :: zcd(nlayers+1,ngptsw), zcu(nlayers+1,ngptsw)
      real(kind=jprb) :: zfd(nlayers+1,ngptsw), zfu(nlayers+1,ngptsw)

! Inactive arrays
!     real(kind=jprb) :: zbbcd(nlayers+1), zbbcu(nlayers+1)
!     real(kind=jprb) :: zbbfd(nlayers+1), zbbfu(nlayers+1)
!     real(kind=jprb) :: zbbfddir(nlayers+1), zbbcddir(nlayers+1)

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

! Apply adjustments for correct Earth/Sun distance and zenith angle to incoming solar flux
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

! Clear
!                zdbtmc = exp(-ztauc(jk) / prmu0)
 
! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = ztauc(jk) / prmu0
               if (ze1 .le. od_lo) then
                  zdbtmc = 1._jprb - ze1 + 0.5_jprb * ze1 * ze1
               else 
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_jprb
                  zdbtmc = exp_tbl(itind)
               endif

               zdbtc_nodel(jk) = zdbtmc
               ztdbtc_nodel(jk+1) = zdbtc_nodel(jk) * ztdbtc_nodel(jk)

! Clear + Cloud
               tauorig = ztauc(jk) + ptaucorig(ikl,ibm)
!                zdbtmo = exp(-tauorig / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = tauorig / prmu0
               if (ze1 .le. od_lo) then
                  zdbtmo = 1._jprb - ze1 + 0.5_jprb * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_jprb
                  zdbtmo = exp_tbl(itind)
               endif

               zdbt_nodel(jk) = zclear * zdbtmc + zcloud * zdbtmo
               ztdbt_nodel(jk+1) = zdbt_nodel(jk) * ztdbt_nodel(jk)
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

! Clear
!                zdbtmc = exp(-ztauc(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = ztauc(jk) / prmu0
               if (ze1 .le. od_lo) then
                  zdbtmc = 1._jprb - ze1 + 0.5_jprb * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_jprb
                  zdbtmc = exp_tbl(itind)
               endif

               zdbtc(jk) = zdbtmc
               ztdbtc(jk+1) = zdbtc(jk)*ztdbtc(jk)

! Clear + Cloud
!                zdbtmo = exp(-ztauo(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = ztauo(jk) / prmu0
               if (ze1 .le. od_lo) then
                  zdbtmo = 1._jprb - ze1 + 0.5_jprb * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_jprb
                  zdbtmo = exp_tbl(itind)
               endif

               zdbt(jk) = zclear*zdbtmc + zcloud*zdbtmo
               ztdbt(jk+1) = zdbt(jk)*ztdbt(jk)
        
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

      end subroutine spcvrt_sw

      end module rrtmg_sw_spcvrt


