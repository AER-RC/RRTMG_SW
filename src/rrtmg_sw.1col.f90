!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!

      program rrtmg_sw

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                             RRTMG_SW                                     *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                 a rapid radiative transfer model                         *
! *                  for the solar spectral region                           *
! *           for application to general circulation models                  *
! *                                                                          *
! *                                                                          *
! *           Atmospheric and Environmental Research, Inc.                   *
! *                       131 Hartwell Avenue                                *
! *                       Lexington, MA 02421                                *
! *                                                                          *
! *                                                                          *
! *                          Eli J. Mlawer                                   *
! *                       Jennifer S. Delamere                               *
! *                        Michael J. Iacono                                 *
! *                        Shepard A. Clough                                 *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                      email:  miacono@aer.com                             *
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people: Steven J. Taubman, Patrick D. Brown,             *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

! ------- Description -------

! This program is the driver for RRTMG_SW, the AER SW radiation model for 
!  application to GCMs, that has been adapted from RRTM_SW for improved
!  efficiency and to provide fractional cloudiness and cloud overlap
!  capability using McICA.

! This routine
!    a) calls RRTMG_SW_INI to initialize data and to perform
!       the g-point interval reduction from 224 to 112
!    b) calls READPROF to read in the atmospheric profile;
!       all layering in RRTMG is ordered from surface to toa. 
!    c) calls CLDPROP to set cloud optical depth based on input
!       cloud properties, or CLDPRMC to set cloud optical depth
!       for McICA
!    d) calls SETCOEF to calculate various quantities needed for 
!       the radiative transfer algorithm
!    e) calls SPCVRT to call the two-stream model that in turn 
!       calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands and to perform the radiative transfer
!       with or without McICA, the Monte-Carlo Independent Column
!       Approximation to represent sub-grid scale cloud variability
!    f) writes out the calculated fluxes and cooling rates
!
! Two modes of operation are possible:
!     The mode is chosen by setting flag imca below.  
!
!    1) Standard, single forward model calculation (imca = 0); this is 
!       valid only for clear sky or fully overcast clouds
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!       For single column calculations, this method also requires setting flag
!       nmca below to the sample size of the Monte Carlo calculation; 
!       (nmca = 200 is recommended). This is method is valid for clear sky
!       or partial cloud conditions
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflag, iceflag and liqflag; see text file rrtmg_sw_instructions
!     and subroutine rrtmg_sw_cldprop.f90 for further details):
!
!    1) Input cloud fraction and cloud optical depth directly (inflgsw = 0)
!    2) Input cloud fraction and cloud physical properties (inflgsw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of iceflgsw and liqflgsw
!
! Two methods of aerosol property input are possible:
!     Aerosol properties can be input in one of two ways (controlled by input 
!     flag iaer, see text file rrtmg_sw_instructions for further details):
!
!    1) Input aerosol optical depth, single scattering albedo and asymmetry
!       parameter directly by layer and spectral band (iaer=10)
!    2) Input aerosol optical depth and 0.55 micron directly by layer and use
!       one or more of six ECMWF aerosol types (iaer=6)
!
!
! ------- Modifications -------
!
! This version of RRTMG_SW has been modified from RRTM_SW to use a reduced
! set of g-point intervals and a two-stream model for application to GCMs. 
!
!-- Original version (derived from RRTM_SW)
!     2002: AER. Inc.
!-- Conversion to F90 formatting; addition of 2-stream radiative transfer
!     Feb 2003: J.-J. Morcrette, ECMWF
!-- Additional modifications for GCM application
!     Aug 2003: M. J. Iacono, AER Inc.
!-- Total number of g-points reduced from 224 to 112.  Original
!   set of 224 can be restored by exchanging code in module parrrsw.f90 
!   and in file rrtmg_sw_init.f90.
!     Apr 2004: M. J. Iacono, AER, Inc.
!-- Modifications to include output for direct and diffuse 
!   downward fluxes.  There are output as "true" fluxes without
!   any delta scaling applied.  Code can be commented to exclude
!   this calculation in source file rrtmg_sw_spcvrt.f90.
!     Jan 2005: E. J. Mlawer, M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability.
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Reformatted for consistency with rrtmg_lw.
!     Feb 2007: M. J. Iacono, AER, Inc.

! --------- Modules ---------

      use parkind, only : jpim, jprb
      use parrrsw, only : mxlay, nbndsw, ngptsw, naerec, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi
      use rrsw_wvn, only : wavenum1, wavenum2
      use rrsw_vsn
      use mcica_subcol_gen_sw, only: mcica_subcol_sw
      use rrtmg_sw_cldprop, only: cldprop_sw
      use rrtmg_sw_cldprmc, only: cldprmc_sw
      use rrtmg_sw_init, only: rrtmg_sw_ini
      use rrtmg_sw_setcoef, only: setcoef_sw
      use rrtmg_sw_spcvrt, only: spcvrt_sw
      use rrtmg_sw_spcvmc, only: spcvmc_sw

      implicit none

! ------- Declarations

! ----- Local -----

! Control
      integer(kind=jpim) :: nlayers             ! total number of layers
      integer(kind=jpim) :: istart              ! beginning band of calculation
      integer(kind=jpim) :: iend                ! ending band of calculation
      integer(kind=jpim) :: icld                ! clear/cloud and cloud overlap flag
      integer(kind=jpim) :: icpr                ! cldprop/cldprmc use flag
      integer(kind=jpim) :: iflag               ! control flag
      integer(kind=jpim) :: iout                ! output option flag
      integer(kind=jpim) :: iaer                ! aerosol option flag
      integer(kind=jpim) :: idelm               ! delta-m scaling flag
      integer(kind=jpim) :: isccos              ! instrumental cosine response flag
      integer(kind=jpim) :: ird                 ! input unit
      integer(kind=jpim) :: iwr                 ! output unit
      integer(kind=jpim) :: i                   ! layer loop index                      ! jk
      integer(kind=jpim) :: ib                  ! band loop index                       ! jsw
      integer(kind=jpim) :: ia, ig              ! indices
      integer(kind=jpim) :: iplon               ! column loop index                     ! jl
      integer(kind=jpim) :: ims                 ! mcica statistical loop index
      integer(kind=jpim) :: imca                ! flag for mcica [0=off, 1=on]
      integer(kind=jpim) :: nmca                ! number of mcica samples (mcica mode)
      integer(kind=jpim), parameter :: ncol = 1 ! total number of columns

      integer(kind=jpim) :: iout1, iout2        ! output control flags
      integer(kind=jpim) :: indform             ! output control flag
      character page 
      character*50 outform(7)

      real(kind=jprb) :: zepsec, zepzen         ! epsilon
      real(kind=jprb) :: zdpgcp                 ! flux to heating conversion ratio


! Atmosphere
      real(kind=jprb) :: pavel(mxlay)           ! layer pressures (mb) 
      real(kind=jprb) :: tavel(mxlay)           ! layer temperatures (K)
      real(kind=jprb) :: pz(0:mxlay)            ! level (interface) pressures (hPa, mb)
      real(kind=jprb) :: tz(0:mxlay)            ! level (interface) temperatures (K)
      real(kind=jprb) :: tbound                 ! surface temperature (K)
      real(kind=jprb) :: pdp(mxlay)             ! layer pressure thickness (hPa, mb)
      real(kind=jprb) :: coldry(mxlay)          ! 
      real(kind=jprb) :: wbrodl(mxlay)          !
      real(kind=jprb) :: wkl(mxmol,mxlay)       ! molecular amounts (mol/cm-2)

      real(kind=jprb) :: cossza, zenith         ! cosine of solar zenith angle 
!      real(kind=jprb) :: earth_sun              ! function for Earth/Sun distance factor
      real(kind=jprb) :: adjflux(jpband)        ! adjustment for current Earth/Sun distance
      real(kind=jprb) :: solvar(jpband)         ! solar constant scaling factor from rrtmg_sw
                                                !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=jprb) :: semiss(jpband)         ! surface emissivity
      real(kind=jprb) :: albdir(nbndsw)         ! surface albedo, direct          ! zalbp
      real(kind=jprb) :: albdif(nbndsw)         ! surface albedo, diffuse         ! zalbd

      real(kind=jprb) :: tauaer(mxlay,jpband)   ! aerosol optical depth (iaer=10 only)
                                                ! (non-delta scaled)      
      real(kind=jprb) :: ssaaer(mxlay,jpband)   ! aerosol single scattering albedo (iaer=10 only)
                                                ! (non-delta scaled)      
      real(kind=jprb) :: asmaer(mxlay,jpband)   ! aerosol asymmetry parameter (iaer=10 only)
                                                ! (non-delta scaled)      
                                                !   first moment of input phase function
      real(kind=jprb) :: ecaer(mxlay,naerec)    ! aerosol optical thickness at 0.55 micron (iaer=6 only)
                                                ! (non-delta scaled)      

! Atmosphere - setcoef
      integer(kind=jpim) :: laytrop             ! tropopause layer index
      integer(kind=jpim) :: layswtch            ! tropopause layer index
      integer(kind=jpim) :: laylow              ! tropopause layer index
      integer(kind=jpim) :: jp(mxlay)           ! 
      integer(kind=jpim) :: jt(mxlay)           !
      integer(kind=jpim) :: jt1(mxlay)          !

      real(kind=jprb) :: colh2o(mxlay)          ! column amount (h2o)
      real(kind=jprb) :: colco2(mxlay)          ! column amount (co2)
      real(kind=jprb) :: colo3(mxlay)           ! column amount (o3)
      real(kind=jprb) :: coln2o(mxlay)          ! column amount (n2o)
      real(kind=jprb) :: colch4(mxlay)          ! column amount (ch4)
      real(kind=jprb) :: colo2(mxlay)           ! column amount (o2)
      real(kind=jprb) :: colmol(mxlay)          ! column amount
      real(kind=jprb) :: co2mult(mxlay)         ! column amount 

      integer(kind=jpim) :: indself(mxlay)
      integer(kind=jpim) :: indfor(mxlay)
      real(kind=jprb) :: selffac(mxlay)
      real(kind=jprb) :: selffrac(mxlay)
      real(kind=jprb) :: forfac(mxlay)
      real(kind=jprb) :: forfrac(mxlay)

      real(kind=jprb) :: &                      !
                         fac00(mxlay), fac01(mxlay), &
                         fac10(mxlay), fac11(mxlay) 

! Atmosphere/clouds - cldprop
      integer(kind=jpim) :: ncbands             ! number of cloud spectral bands
      integer(kind=jpim) :: inflag              ! flag for cloud property method
      integer(kind=jpim) :: iceflag             ! flag for ice cloud properties
      integer(kind=jpim) :: liqflag             ! flag for liquid cloud properties

      real(kind=jprb) :: cldfrac(mxlay)         ! layer cloud fraction
      real(kind=jprb) :: tauc(nbndsw,mxlay)     ! cloud optical depth (non-delta scaled)
      real(kind=jprb) :: ssac(nbndsw,mxlay)     ! cloud single scattering albedo (non-delta scaled)
      real(kind=jprb) :: asmc(nbndsw,mxlay)     ! cloud asymmetry parameter (non-delta scaled)
      real(kind=jprb) :: ciwp(mxlay)            ! cloud ice water path
      real(kind=jprb) :: clwp(mxlay)            ! cloud liquid water path
      real(kind=jprb) :: rei(mxlay)             ! cloud ice particle size
      real(kind=jprb) :: rel(mxlay)             ! cloud liquid particle size

      real(kind=jprb) :: taucloud(mxlay,jpband) ! cloud optical depth
      real(kind=jprb) :: taucldorig(mxlay,jpband)! cloud optical depth (non-delta scaled)
      real(kind=jprb) :: ssacloud(mxlay,jpband) ! cloud single scattering albedo
      real(kind=jprb) :: asmcloud(mxlay,jpband) ! cloud asymmetry parameter

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=jprb) :: cldfmc(ngptsw,mxlay)   ! cloud fraction [mcica]
      real(kind=jprb) :: ciwpmc(ngptsw,mxlay)   ! cloud ice water path [mcica]
      real(kind=jprb) :: clwpmc(ngptsw,mxlay)   ! cloud liquid water path [mcica]
      real(kind=jprb) :: relqmc(mxlay)          ! liquid particle size (microns)
      real(kind=jprb) :: reicmc(mxlay)          ! ice partcle size (microns)
      real(kind=jprb) :: taucmc(ngptsw,mxlay)   ! cloud optical depth [mcica]
      real(kind=jprb) :: taormc(ngptsw,mxlay)   ! unscaled cloud optical depth [mcica]
      real(kind=jprb) :: ssacmc(ngptsw,mxlay)   ! cloud single scattering albedo [mcica]
      real(kind=jprb) :: asmcmc(ngptsw,mxlay)   ! cloud asymmetry parameter [mcica]

! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      real(kind=jprb) :: ztauc(mxlay,nbndsw)    ! cloud optical depth
      real(kind=jprb) :: ztaucorig(mxlay,nbndsw)! unscaled cloud optical depth
      real(kind=jprb) :: zasyc(mxlay,nbndsw)    ! cloud asymmetry parameter 
                                                !  (first moment of phase function)
      real(kind=jprb) :: zomgc(mxlay,nbndsw)    ! cloud single scattering albedo
      real(kind=jprb) :: ztaua(mxlay,nbndsw)    ! total aerosol optical depth
      real(kind=jprb) :: zasya(mxlay,nbndsw)    ! total aerosol asymmetry parameter 
      real(kind=jprb) :: zomga(mxlay,nbndsw)    ! total aerosol single scattering albedo

      real(kind=jprb) :: zcldfmc(mxlay,ngptsw)  ! cloud fraction [mcica]
      real(kind=jprb) :: ztaucmc(mxlay,ngptsw)  ! cloud optical depth [mcica]
      real(kind=jprb) :: ztaormc(mxlay,ngptsw)  ! unscaled cloud optical depth [mcica]
      real(kind=jprb) :: zasycmc(mxlay,ngptsw)  ! cloud asymmetry parameter [mcica] 
      real(kind=jprb) :: zomgcmc(mxlay,ngptsw)  ! cloud single scattering albedo [mcica]

      real(kind=jprb) :: zbbfu(mxlay+1)         ! temporary upward shortwave flux (w/m2)
      real(kind=jprb) :: zbbfd(mxlay+1)         ! temporary downward shortwave flux (w/m2)
      real(kind=jprb) :: zbbcu(mxlay+1)         ! temporary clear sky upward shortwave flux (w/m2)
      real(kind=jprb) :: zbbcd(mxlay+1)         ! temporary clear sky downward shortwave flux (w/m2)
      real(kind=jprb) :: zbbfddir(mxlay+1)      ! temporary downward direct shortwave flux (w/m2)
      real(kind=jprb) :: zbbcddir(mxlay+1)      ! temporary clear sky downward direct shortwave flux (w/m2)
      real(kind=jprb) :: zuvfd(mxlay+1)         ! temporary UV downward shortwave flux (w/m2)
      real(kind=jprb) :: zuvcd(mxlay+1)         ! temporary clear sky UV downward shortwave flux (w/m2)
      real(kind=jprb) :: zuvfddir(mxlay+1)      ! temporary UV downward direct shortwave flux (w/m2)
      real(kind=jprb) :: zuvcddir(mxlay+1)      ! temporary clear sky UV downward direct shortwave flux (w/m2)
      real(kind=jprb) :: znifd(mxlay+1)         ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=jprb) :: znicd(mxlay+1)         ! temporary clear sky near-IR downward shortwave flux (w/m2)
      real(kind=jprb) :: znifddir(mxlay+1)      ! temporary near-IR downward direct shortwave flux (w/m2)
      real(kind=jprb) :: znicddir(mxlay+1)      ! temporary clear sky near-IR downward direct shortwave flux (w/m2)

      real(kind=jprb) :: totuflux(0:mxlay)      ! upward shortwave flux (w/m2)                  ! pfup
      real(kind=jprb) :: totdflux(0:mxlay)      ! downward shortwave flux (w/m2)                ! pfdown
      real(kind=jprb) :: fnet(0:mxlay)          ! net shortwave flux (w/m2)                     ! pfls
      real(kind=jprb) :: htr(0:mxlay)           ! shortwave heating rate (k/day)                ! pheat
      real(kind=jprb) :: totuclfl(0:mxlay)      ! clear sky upward shortwave flux (w/m2)        ! pcup 
      real(kind=jprb) :: totdclfl(0:mxlay)      ! clear sky downward shortwave flux (w/m2)      ! pcdown 
      real(kind=jprb) :: fnetc(0:mxlay)         ! clear sky net shortwave flux (w/m2)           ! pfcs
      real(kind=jprb) :: htrc(0:mxlay)          ! clear sky shortwave heating rate (k/day)      ! pheac

      real(kind=jprb) :: dirdflux(0:mxlay)      ! direct downward shortwave flux (w/m2)         ! dirdownflux
      real(kind=jprb) :: difdflux(0:mxlay)      ! diffuse downward shortwave flux (w/m2)        ! difdownflux
      real(kind=jprb) :: dflxuv(0:mxlay)        ! Total sky downward shortwave flux, UV/vis     ! pfdnuv
      real(kind=jprb) :: dflxir(0:mxlay)        ! Total sky downward shortwave flux, near-IR    ! pfdnir 
      real(kind=jprb) :: dirdnuv(0:mxlay)       ! Direct downward shortwave surface flux, UV/vis
      real(kind=jprb) :: difdnuv(0:mxlay)       ! Diffuse downward shortwave surface flux, UV/vis
      real(kind=jprb) :: dirdnir(0:mxlay)       ! Direct downward shortwave surface flux, near-IR
      real(kind=jprb) :: difdnir(0:mxlay)       ! Diffuse downward shortwave surface flux, near-IR

! Output - inactive
!      real(kind=jprb) :: zuvfu(mxlay+1)         ! temporary upward UV shortwave flux (w/m2)
!      real(kind=jprb) :: zuvfd(mxlay+1)         ! temporary downward UV shortwave flux (w/m2)
!      real(kind=jprb) :: zuvcu(mxlay+1)         ! temporary clear sky upward UV shortwave flux (w/m2)
!      real(kind=jprb) :: zuvcd(mxlay+1)         ! temporary clear sky downward UV shortwave flux (w/m2)
!      real(kind=jprb) :: zvsfu(mxlay+1)         ! temporary upward visible shortwave flux (w/m2)
!      real(kind=jprb) :: zvsfd(mxlay+1)         ! temporary downward visible shortwave flux (w/m2)
!      real(kind=jprb) :: zvscu(mxlay+1)         ! temporary clear sky upward visible shortwave flux (w/m2)
!      real(kind=jprb) :: zvscd(mxlay+1)         ! temporary clear sky downward visible shortwave flux (w/m2)
!      real(kind=jprb) :: znifu(mxlay+1)         ! temporary upward near-IR shortwave flux (w/m2)
!      real(kind=jprb) :: znifd(mxlay+1)         ! temporary downward near-IR shortwave flux (w/m2)
!      real(kind=jprb) :: znicu(mxlay+1)         ! temporary clear sky upward near-IR shortwave flux (w/m2)
!      real(kind=jprb) :: znicd(mxlay+1)         ! temporary clear sky downward near-IR shortwave flux (w/m2)

! Output (mean output for McICA calculation)
      real(kind=jprb) :: uflxsum(0:mxlay)       ! upward shortwave flux (w/m2)
      real(kind=jprb) :: dflxsum(0:mxlay)       ! downward shortwave flux (w/m2)
      real(kind=jprb) :: dirdsum(0:mxlay)       ! direct downward shortwave flux (w/m2)
      real(kind=jprb) :: difdsum(0:mxlay)       ! diffuse downward shortwave flux (w/m2)
      real(kind=jprb) :: fnetsum(0:mxlay)       ! net shortwave flux (w/m2)
      real(kind=jprb) :: htrsum(0:mxlay)        ! shortwave heating rate (k/day)


! Initializations

      hvrrtm = '$Revision$'
      hvrini = 'NOT USED'
      hvrcld = 'NOT USED'
      hvrclc = 'NOT USED'
      hvrrft = 'NOT USED'
      hvrspv = 'NOT USED'
      hvrspc = 'NOT USED'
      hvrset = 'NOT USED'
      hvrtau = 'NOT USED'
      hvrvqd = 'NOT USED'
      hvrkg  = '$Revision$'
      hvratm = 'NOT USED'
      hvrutl = 'NOT USED'
      hvrext = 'NOT USED'

      hnamrtm = '  rrtmg_sw.1col.f90:'
      hnamini = '  rrtmg_sw_init.f90:'
      hnamcld = 'rrtmgsw_cldprop.f90:'
      hnamclc = 'rrtmgsw_cldprmc.f90:'
      hnamrft = 'rrtmg_sw_reftra.f90:'
      hnamspv = 'rrtmg_sw_spcvrt.f90:'
      hnamspc = 'rrtmg_sw_spcvmc.f90:'
      hnamset = 'rrtmgsw_setcoef.f90:'
      hnamtau = 'rrtmg_sw_taumol.f90:'
      hnamvqd = 'rrtmg_sw_vrtqdr.f90:'
      hnamkg  = '   rrtmg_sw_k_g.f90:'
      hnamatm = '           rrtatm.f:'
      hnamutl = '         util_xxx.f:'
      hnamext = '            extra.f:'

      zepsec = 1.e-06_jprb
      zepzen = 1.e-10_jprb
      oneminus = 1.0_jprb - zepsec
      pi = 2._jprb * asin(1._jprb)

      ird = 9
      iwr = 10
      icpr = 0
      page = char(12)

      data outform &
         /'(1x,i3,3x,f7.6,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.5,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.4,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.3,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.2,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.1,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.1,4x,4(f10.4,4x),f11.6,4x,f10.5)'/

      uflxsum(0:) = 0._jprb
      dflxsum(0:) = 0._jprb
      dirdsum(0:) = 0._jprb
      difdsum(0:) = 0._jprb
      fnetsum(0:) = 0._jprb
      htrsum(0:) = 0._jprb

! Set imca to select calculation type:
!  (read by subroutine readprof from input file INPUT_RRTM):  
! imca = 0, use standard forward model calculation (clear and overcast only)
! imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability
!           (clear, overcast or partial cloud conditions)

! Set icld to select of clear or cloud calculation and cloud overlap method
!  (read by subroutine readprof from input file INPUT_RRTM):  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap (McICA only)
! icld = 2, with clouds using maximum/random cloud overlap (McICA only)
! icld = 3, with clouds using maximum cloud overlap (McICA only)

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 224 to 112 for input absorption
! coefficient data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  

      call rrtmg_sw_ini

! Open the input set of atmospheres
      open (ird,file='INPUT_RRTM',form='formatted')
! Open the output file
      open (iwr,file='OUTPUT_RRTM',form='formatted')
      

! This is the main longitude/column loop within rrtmg.

      do iplon = 1, ncol

! Input atmospheric profile from INPUT_RRTM. 
        call readprof(ird, nlayers, iout, imca, icld, iaer, isccos, idelm, pdp, &
                       pavel, tavel, pz, tz, tbound, semiss, zenith, adjflux, &
                       coldry, wkl, inflag, iceflag, liqflag, &
                       cldfrac, tauc, ssac, asmc, ciwp, clwp, rei, rel, &
                       tauaer, ssaaer, asmaer)

         istart = jpb1                 ! jpb1 = 16
         iend = jpb2                   ! jpb2 = 29
         iflag = iout

! Set nmca to sample size for Monte Carlo calculation
         if (imca.eq.0) nmca = 1
         if (imca.eq.1) nmca = 200

! Return here for multiple band output
 1000    continue
         if (iflag .gt. 0 .and. iflag .le. jpb2) then
            istart = iflag
            iend = iflag
         endif

! This is the statistical sampling loop for McICA

         do ims = 1, nmca

! Call sub-colum cloud generator for McICA calculations.
! Output will be averaged over all nmca samples.  The code can be modified to
! write output for each individual sample (this will be excessive if output
! is also requested for each spectral band).  

            if (imca.eq.1) then
               call mcica_subcol_sw(iplon, nlayers, icld, ims, pavel, &
                          cldfrac, ciwp, clwp, rei, rel, tauc, ssac, asmc, &
                          cldfmc, ciwpmc, clwpmc, reicmc, relqmc, taucmc, &
                          ssacmc, asmcmc)
            endif

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed in cldprop.  Cloud fraction and cloud
!  optical properties are transferred to rrtmg_sw arrays in cldprop.  
!  Note: Model will be stopped if partial cloud present without McICA.

!  If McICA is requested use cloud fraction and cloud physical properties 
!  generated by sub-column cloud generator above. 

            if (imca.eq.0) then
               do i = 1, nlayers
                  if (cldfrac(i).ne.0._jprb .and. cldfrac(i).ne.1._jprb) then
                     stop 'PARTIAL CLOUD NOT ALLOWED FOR IMCA=0'
                  endif
               enddo
               call cldprop_sw(nlayers, inflag, iceflag, liqflag, cldfrac, &
                               tauc, ssac, asmc, ciwp, clwp, rei, rel, &
                               taucldorig, taucloud, ssacloud, asmcloud)
               icpr = 1
            else
               call cldprmc_sw(nlayers, inflag, iceflag, liqflag, cldfmc, &
                               ciwpmc, clwpmc, reicmc, relqmc, &
                               taormc, taucmc, ssacmc, asmcmc)
               icpr = 1
            endif


! Calculate coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients by interpolating data from stored
! reference atmospheres.

            call setcoef_sw(nlayers, pavel, tavel, pz, tz, tbound, coldry, wkl, &
                            laytrop, layswtch, laylow, jp, jt, jt1, &
                            co2mult, colch4, colco2, colh2o, colmol, coln2o, &
                            colo2, colo3, fac00, fac01, fac10, fac11, &
                            selffac, selffrac, indself, forfac, forfrac, indfor)
  

! Cosine of the solar zenith angle 
!  Prevent using value of zero;

            cossza = zenith
            if (cossza.eq.0._jprb) cossza = zepzen

! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer 
  
! Surface albedo
            do ib=1,nbndsw
               albdif(ib) = 1._jprb - semiss(jpb1-1+ib)
               albdir(ib) = 1._jprb - semiss(jpb1-1+ib)
            enddo

! Clouds
            if (icld.eq.0) then

               ztauc(:,:) = 0._jprb
               ztaucorig(:,:) = 0._jprb
               zasyc(:,:) = 0._jprb
               zomgc(:,:) = 0._jprb
               zcldfmc(:,:) = 0._jprb
               ztaucmc(:,:) = 0._jprb
               ztaormc(:,:) = 0._jprb
               zasycmc(:,:) = 0._jprb
               zomgcmc(:,:) = 0._jprb

            elseif (icld.ge.1) then
               if (imca.eq.0) then
                  do i=1,nlayers
                     do ib=1,nbndsw
                        if (cldfrac(i) .ge. zepsec) then
                           ztauc(i,ib) = taucloud(i,jpb1-1+ib)
                           ztaucorig(i,ib) = taucldorig(i,jpb1-1+ib)
                           zasyc(i,ib) = asmcloud(i,jpb1-1+ib)
                           zomgc(i,ib) = ssacloud(i,jpb1-1+ib)
                        endif
                     enddo
                  enddo
               else
                  do i=1,nlayers
                     do ig=1,ngptsw
                        zcldfmc(i,ig) = cldfmc(ig,i)
                        ztaucmc(i,ig) = taucmc(ig,i)
                        ztaormc(i,ig) = taormc(ig,i)
                        zasycmc(i,ig) = asmcmc(ig,i)
                        zomgcmc(i,ig) = ssacmc(ig,i)
                     enddo
                  enddo
               endif   
            endif   

! Aerosol
! IAER = 0: no aerosols
            if (iaer.eq.0) then

               ztaua(:,:) = 0._jprb
               zasya(:,:) = 0._jprb
               zomga(:,:) = 0._jprb

! IAER = 6: Use ECMWF six aerosol types. See rrsw_aer.f90 for details.
! Input aerosol optical thickness at 0.55 micron for each aerosol type (ecaer)
! or set manually here for each aerosol and layer. 
            elseif (iaer.eq.6) then

               do i = 1, nlayers
                  do ia = 1, naerec
                     ecaer(i,ia) = 1.0e-15_jprb
!                      ecaer(i,ia) = 1.0e-04_jprb
                  enddo
               enddo

               do i = 1, nlayers
                  do ib = 1, nbndsw
                     ztaua(i,ib) = 0._jprb
                     zasya(i,ib) = 0._jprb
                     zomga(i,ib) = 0._jprb
                     do ia = 1, naerec
                        ztaua(i,ib) = ztaua(i,ib) + rsrtaua(ib,ia) * ecaer(i,ia)
                        zomga(i,ib) = zomga(i,ib) + rsrtaua(ib,ia) * ecaer(i,ia) * &
                                      rsrpiza(ib,ia)
                        zasya(i,ib) = zasya(i,ib) + rsrtaua(ib,ia) * ecaer(i,ia) * &
                                      rsrpiza(ib,ia) * rsrasya(ib,ia)
                     enddo
                     if (zomga(i,ib) /= 0._jprb) then
                        zasya(i,ib) = zasya(i,ib) / zomga(i,ib)
                     endif
                     if (ztaua(i,ib) /= 0._jprb) then
                        zomga(i,ib) = zomga(i,ib) / ztaua(i,ib)
                     endif
                  enddo
               enddo

! IAER=10: Direct specification of aerosol properties from IN_AER_RRTM.
            elseif (iaer.eq.10) then

               do i = 1 ,nlayers
                  do ib = 1 ,nbndsw
                     ztaua(i,ib) = tauaer(i,jpb1-1+ib)
                     zasya(i,ib) = asmaer(i,jpb1-1+ib)
                     zomga(i,ib) = ssaaer(i,jpb1-1+ib)
                  enddo
               enddo

            endif

! Call the 2-stream radiation transfer model

            do i=1,nlayers+1
               zbbcu(i) = 0._jprb
               zbbcd(i) = 0._jprb
               zbbfu(i) = 0._jprb
               zbbfd(i) = 0._jprb
               zbbcddir(i) = 0._jprb
               zbbfddir(i) = 0._jprb
               zuvcd(i) = 0._jprb
               zuvfd(i) = 0._jprb
               zuvcddir(i) = 0._jprb
               zuvfddir(i) = 0._jprb
               znicd(i) = 0._jprb
               znifd(i) = 0._jprb
               znicddir(i) = 0._jprb
               znifddir(i) = 0._jprb
            enddo

            if (imca.eq.0) then
               call spcvrt_sw &
                   (nlayers, istart, iend, icpr, iout, &
                    pavel, tavel, pz, tz, tbound, albdif, albdir, &
                    cldfrac, ztauc, zasyc, zomgc, ztaucorig, &
                    ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, &	 
                    laytrop, layswtch, laylow, jp, jt, jt1, &
                    co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
                    fac00, fac01, fac10, fac11, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, zuvcd, znifd, znicd, &
                    zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir)
            else
               call spcvmc_sw &
                   (nlayers, istart, iend, icpr, iout, &
                    pavel, tavel, pz, tz, tbound, albdif, albdir, &
                    zcldfmc, ztaucmc, zasycmc, zomgcmc, ztaormc, &
                    ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, &	 
                    laytrop, layswtch, laylow, jp, jt, jt1, &
                    co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
                    fac00, fac01, fac10, fac11, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, zuvcd, znifd, znicd, &
                    zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir)
            endif


! Prepare output up and down, clear and total flux output
            do i = 1, nlayers+1
               totuclfl(i-1) = zbbcu(i)
               totdclfl(i-1) = zbbcd(i)
               totuflux(i-1) = zbbfu(i)
               totdflux(i-1) = zbbfd(i)
! Prepare direct/diffuse flux output
               dirdflux(i-1) = zbbfddir(i)
               difdflux(i-1) = totdflux(i-1) - dirdflux(i-1)
            enddo

! Prepare net clear and total flux output
            do i = 1, nlayers+1
               fnetc(i-1) = totdclfl(i-1) - totuclfl(i-1)
               fnet(i-1) = totdflux(i-1) - totuflux(i-1)
            enddo

! Output clear and total heating rates
            do i = 1, nlayers
               zdpgcp = heatfac / pdp(i)
               htrc(i-1) = (fnetc(i) - fnetc(i-1)) * zdpgcp
               htr(i-1) = (fnet(i) - fnet(i-1)) * zdpgcp
            enddo
            htr(nlayers) = 0._jprb
            htrc(nlayers) = 0._jprb

! Process output.
            if (iout .lt. 0) goto 2000

            if (imca .eq. 0) then

               if (istart .ne. iend) then
                  iout1 = jpb2
                  iout2 = iend-1
               else
                  iout1 = istart
                  iout2 = iend
               endif

               if (isccos .eq. 1) then 
                  write(iwr,9880) 
               elseif (isccos .eq. 2) then
                  write(iwr,9881)
               else
                  write(iwr,9879)
               endif

               if (idelm .eq. 0) then
                  write(iwr,9883)
               else
                  write(iwr,9882)
               endif

               write(iwr,9899) wavenum1(iout1),wavenum2(iout2),iplon
               write(iwr,9900)
               write(iwr,9901)

               do i = nlayers, 0, -1
                  if (pz(i) .lt. 1.e-2_jprb) then
                     indform = 1
                  elseif (pz(i) .lt. 1.e-1_jprb) then
                     indform = 2
                  elseif (pz(i) .lt. 1._jprb) then
                     indform = 3
                  elseif (pz(i) .lt. 10._jprb) then
                     indform = 4
                  elseif (pz(i) .lt. 100._jprb) then
                     indform = 5
                  elseif (pz(i) .lt. 1000._jprb) then
                     indform = 6
                  else
                     indform = 7
                  endif
                  write(iwr,outform(indform)) i, pz(i), totuflux(i), &
                       difdflux(i), dirdflux(i), totdflux(i), fnet(i), htr(i)
               enddo
               write(iwr,9903)page

            elseif (imca .eq. 1) then

               do i = nlayers, 0, -1
                  uflxsum(i) = uflxsum(i) + totuflux(i)
                  dflxsum(i) = dflxsum(i) + totdflux(i)
                  dirdsum(i) = dirdsum(i) + dirdflux(i)
                  difdsum(i) = difdsum(i) + difdflux(i)
                  fnetsum(i) = fnetsum(i) + fnet(i)
                  htrsum(i) = htrsum(i) + htr(i)
               enddo

! Output average over samples when last sample reached. Comment this if-check 
! and swap flux output write statement below to output all McICA samples. 
               if (ims .eq. nmca) then

                  do i = nlayers, 0, -1
                     uflxsum(i) = uflxsum(i)/nmca
                     dflxsum(i) = dflxsum(i)/nmca
                     dirdsum(i) = dirdsum(i)/nmca
                     difdsum(i) = difdsum(i)/nmca
                     fnetsum(i) = fnetsum(i)/nmca
                     htrsum(i) = htrsum(i)/nmca
                  enddo

                  if (istart .ne. iend) then
                     iout1 = jpb2
                     iout2 = iend-1
                  else
                     iout1 = istart
                     iout2 = iend
                  endif

                  if (isccos .eq. 1) then 
                     write(iwr,9880) 
                  elseif (isccos .eq. 2) then
                     write(iwr,9881)
                  else
                     write(iwr,9879)
                  endif

                  if (idelm .eq. 0) then
                     write(iwr,9883)
                  else
                     write(iwr,9882)
                  endif

                  write(iwr,9899) wavenum1(iout1),wavenum2(iout2),iplon
                  write(iwr,9900)
                  write(iwr,9901)

                  do i = nlayers, 0, -1
                     if (pz(i) .lt. 1.e-2_jprb) then
                        indform = 1
                     elseif (pz(i) .lt. 1.e-1_jprb) then
                        indform = 2
                     elseif (pz(i) .lt. 1._jprb) then
                        indform = 3
                     elseif (pz(i) .lt. 10._jprb) then
                        indform = 4
                     elseif (pz(i) .lt. 100._jprb) then
                        indform = 5
                     elseif (pz(i) .lt. 1000._jprb) then
                        indform = 6
                     else
                        indform = 7
                     endif
                     write(iwr,outform(indform)) i, pz(i), uflxsum(i), &
                          difdsum(i), dirdsum(i), dflxsum(i), fnetsum(i), htrsum(i)
!                     write(iwr,outform(indform)) i, pz(i), totuflux(i), &
!                          difdflux(i), dirdflux(i), totdflux(i), fnet(i), htr(i)
                  enddo
                  write(iwr,9903)page

               endif

            endif

 2000       continue
! End statistical loop for McICA
         enddo

         if (iout .ge. 0 .and. iout .le. jpb2) goto 3500
         if (iflag .eq. 98) then
            iflag = jpb1
         elseif (iflag .gt. 0 .and. iflag .lt. jpb2) then
            iflag = iflag + 1
         else
            goto 3500
         endif
            goto 1000
 3500    continue

! Output module version numbers

         write(iwr,9910) hnamrtm,hvrrtm,hnamini,hvrini,hnamcld,hvrcld, &
           hnamclc,hvrclc,hnamrft,hvrrft,hnamspv,hvrspv,hnamspc,hvrspc, &
           hnamset,hvrset,hnamtau,hvrtau,hnamvqd,hvrvqd,hnamatm,hvratm, &
           hnamutl,hvrutl,hnamext,hvrext,hnamkg,hvrkg

 4000    continue

! End longitude/column loop
      enddo

      close(ird)
      close(iwr)

 9879 format(1x)
 9880 format(1x,'All output fluxes have been adjusted to account for ins&
     &trumental cosine response.') 
 9881 format(1x,'The output diffuse fluxes have been adjusted to account&
     & for instrumental cosine response.') 

 9882 format(1x,'The downwelling direct and diffuse fluxes have been com&
     &puted using the delta-M scaling approximation.') 
 9883 format(1x)

 9899 format(1x,'Wavenumbers: ',f6.0,' - ',f6.0,' cm-1, ATM ',i6)
 9900 format(1x,'LEVEL PRESSURE   UPWARD FLUX   DIFDOWN FLUX  DIRDOWN FL&
     &UX  DOWNWARD FLUX   NET FLUX    HEATING RATE')
 9901 format(1x,'         mb          W/m2          W/m2          W/m2&
     &        W/m2          W/m2       degree/day')
 9902 format(1x,i3,3x,f11.6,4x,1p,2(g12.6,2x),g13.6,3x,g16.9,0p)
 9903 format(a)
 9910 format('  Modules and versions used in this calculation:',/,/, &
              7(5x,a20,2x,a18,10x,a20,2x,a18,/))

      contains

!*************************************************************************
      real(kind=jprb) function earth_sun(idn)
!*************************************************************************
!
!  Purpose: Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  idn        : Day of the year
!  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

! ------- Modules -------

      use parkind, only : jpim, jprb
      use rrsw_con, only : pi

      implicit none

      integer(kind=jpim), intent(in) :: idn

      real(kind=jprb) :: gamma

      gamma = 2._jprb*pi*(idn-1)/365._jprb

! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110_jprb + .034221_jprb * cos(gamma) + .001289_jprb * sin(gamma) + &
                   .000719_jprb * cos(2._jprb*gamma) + .000077_jprb * sin(2._jprb*gamma)

      end function earth_sun

!************************************************************************
      subroutine readprof(ird_in, nlayers_out, iout_out, imca, icld_out, &
           iaer_out, isccos_out, idelm_out, pdp, &
           pavel_out, tavel_out, pz_out, tz_out, tbound_out, semiss_out, &
           zenith_out, adjflux_out, &
           coldry_out, wkl_out, inflag_out, iceflag_out,liqflag_out, &
           cldfrac_out, tauc, ssac, asmc, ciwp, clwp, rei, rel, &
           tauaer_out, ssaaer_out, asmaer_out)
!************************************************************************

! -------- Modules --------

      use parkind, only : jpim, jprb
      use rrsw_con, only: pi,planck,boltz,clight,avogad,alosmt,gascon, &
                          radcn1,radcn2
      use rrsw_wvn, only: wavenum1, wavenum2, delwave

! Note: COMMON blocks are left in this routine for array passing with
!       rrtatm.f for reading input profiles in single-column mode.
!       Scalars and arrays in subroutine call are renamed to avoid conflict
!       with COMMON blocks.
!                                                                         
! Purpose: Read in atmospheric profile.

      implicit integer(i-n), real(a-h,o-z)

! ------- Parameters -------
      parameter (mxlay = 203)
      parameter (nbndsw = 14)
      parameter (jpbands = 29)
      parameter (ib1 = 16, ib2 = 29)
      parameter (mg = 16)
      parameter (mxstr = 16)
      parameter (mcmu = 32)
      parameter (mxmol = 38)
      parameter (maxinpx = 35)
      parameter (maxxsec = 4)

      dimension altz(0:mxlay),ixtrans(14)
      dimension solvar(jpbands)

      common /control/  iaer, nstr, iout, istart, iend, icld, idelm, isccos
      common /constants/fluxfac,heatfac
      common /consts/   pic,planckc,boltzc,clightc,avogadc,alosmtc,gasconc, &
                        radcn1c,radcn2c
      common /swprop/   zenith, albedo, adjflux(jpbands)
      common /surface/  ireflect,semiss(jpbands)
      common /profile/  nlayers,pavel(mxlay),tavel(mxlay),pz(0:mxlay),tz(0:mxlay),tbound
      common /species/  coldry(mxlay),wkl(mxmol,mxlay),wbrodl(mxlay),colmol(mxlay),nmol
      common /ifil/     ird,ipr,ipu,idum(15)
      common /xsecctrl/ nxmol,ixindx(maxinpx)
      common /xsec/     wx(maxxsec,mxlay)
      common /pathx/    ixmax,nxmol0,ixindx0(maxinpx),wx0(maxinpx,mxlay)    
      common /xrrtatm/  ixsect

      common /cloudin/   inflag,clddat1(mxlay),clddat2(mxlay), &
                         iceflag,liqflag,clddat3(mxlay),clddat4(mxlay), &
                         clddatmom(0:16,mxlay)
      common /clouddat/  ncbands,cldfrac(mxlay), &
                         taucloud(mxlay,jpbands),ssacloud(mxlay,jpbands), &
                         xmom(0:16,mxlay,jpbands)
      common /aerdat/    ssaaer(mxlay,jpbands), phase(mcmu,mxlay,jpbands), &
                         tauaer(mxlay,jpbands)

      character*80 form1(0:1),form2(0:1),form3(0:1)
      character*1 ctest, cdollar, cdum

! Dimensions for transfer to rrtmg
      integer(kind=jpim), intent(in) :: ird_in                 ! input file unit
      integer(kind=jpim), intent(out) :: nlayers_out           ! total number of layers
      integer(kind=jpim), intent(out) :: imca                  ! McICA on/off flag (1 = use McICA)
      integer(kind=jpim), intent(out) :: icld_out              ! clear/cloud/overlap flag
      integer(kind=jpim), intent(out) :: iout_out              ! output option flag
      integer(kind=jpim), intent(out) :: iaer_out              ! aerosol flag
      integer(kind=jpim), intent(out) :: isccos_out            ! aerosol flag
      integer(kind=jpim), intent(out) :: idelm_out             ! aerosol flag

      real(kind=jprb), intent(out) :: pavel_out(mxlay)         ! layer pressures (mb) 
      real(kind=jprb), intent(out) :: tavel_out(mxlay)         ! layer temperatures (K)
      real(kind=jprb), intent(out) :: pz_out(0:mxlay)          ! level (interface) pressures (hPa, mb)
      real(kind=jprb), intent(out) :: tz_out(0:mxlay)          ! level (interface) temperatures (K)
      real(kind=jprb), intent(out) :: tbound_out               ! surface temperature (K)
      real(kind=jprb), intent(out) :: pdp(mxlay)               ! layer pressure thickness (hPa, mb)
      real(kind=jprb), intent(out) :: coldry_out(mxlay)        ! dry air molecular amount
      real(kind=jprb), intent(out) :: wkl_out(mxmol,mxlay)     ! molecular amounts (mol/cm-2)
      real(kind=jprb), intent(out) :: semiss_out(jpbands)       ! surface emissivity
      real(kind=jprb), intent(out) :: zenith_out               ! cos solar zenith angle
      real(kind=jprb), intent(out) :: adjflux_out(jpbands)      ! adjustment for current Earth/Sun distance

      integer(kind=jpim), intent(out) :: inflag_out              ! cloud property option flag
      integer(kind=jpim), intent(out) :: iceflag_out             ! ice cloud property flag
      integer(kind=jpim), intent(out) :: liqflag_out             ! liquid cloud property flag

      real(kind=jprb), intent(out) :: cldfrac_out(mxlay)         ! cloud fraction
      real(kind=jprb), intent(out) :: tauc(nbndsw,mxlay)         ! cloud optical depth (non-delta scaled)
      real(kind=jprb), intent(out) :: ssac(nbndsw,mxlay)         ! cloud single scattering albedo (non-delta scaled)
      real(kind=jprb), intent(out) :: asmc(nbndsw,mxlay)         ! cloud asymmetry parameter (non-delta scaled)
      real(kind=jprb), intent(out) :: ciwp(mxlay)                ! cloud ice water path
      real(kind=jprb), intent(out) :: clwp(mxlay)                ! cloud liquid water path
      real(kind=jprb), intent(out) :: rei(mxlay)                 ! cloud ice particle size
      real(kind=jprb), intent(out) :: rel(mxlay)                 ! cloud liquid particle size
      real(kind=jprb), intent(out) :: tauaer_out(mxlay,jpbands)   ! aerosol optical depth
      real(kind=jprb), intent(out) :: ssaaer_out(mxlay,jpbands)   ! aerosol single scattering albedo
      real(kind=jprb), intent(out) :: asmaer_out(mxlay,jpbands)   ! aerosol asymmetry parameter
                                                                 !   first momemnt of input phase function
!
! Local
      integer(kind=jpim) :: juldat                               ! day of year
      real(kind=jprb) :: fice(mxlay)                             ! cloud ice fraction

! Initializations

      data cdollar /'$'/
      data ixtrans /0,0,0,1,2,3,0,0,0,0,0,4,0,0/

      form1(0) = '(3f10.4,a3,i2,1x,2(f7.2,f8.3,f7.2))'
      form2(0) = '(3f10.4,a3,i2,23x,(f7.2,f8.3,f7.2))'
      form3(0) = '(8e10.3)'
      form1(1) = '(g15.7,g10.4,g10.4,a3,i2,1x,2(g7.2,g8.3,g7.2))'
      form2(1) = '(g15.7,g10.4,g10.4,a3,i2,23x,(g7.2,g8.3,g7.2))'
      form3(1) = '(8g15.7)'

! Pass constants to common block names for rrtatm
      pic = pi
      planckc = planck
      boltzc = boltz
      clightc = clight
      avogadc = avogad
      alosmtc = alosmt
      gasconc = gascon
      radcn1c = radcn1
      radcn2c = radcn2
 
      ixmax = maxinpx

      ird = ird_in

      do ilay = 1,mxlay
         do isp = 1,mxmol
            wkl(isp,ilay) = 0.0
         enddo
         do isp = 1,maxxsec
            wx(isp,ilay) = 0.0
         enddo
      enddo

! Top of read input loop
 1000 continue
      read (ird,9009,end=8800) ctest
      if (ctest .ne. cdollar) goto 1000

      read (ird,9011) iaer, iatm, iscat, istrm, iout, imca, icld, idelm, icos

      if (idelm.gt.0 .or. idelm.lt.0 .or. icos.gt.0 .or. icos.lt.0) then
         print *,'INVALID MEASUREMENT COMPARISON FLAG'
         stop
      endif
      isccos = icos

! No cross-sections implemented in shortwave.
      ixsect = 0

! Only 2-stream scattering implemented in rrtmg
      if (iscat .ne. 1) then
         print *,'INVALID SCATTERING OPTION CHOSEN'
         stop
      endif

      if (istrm .eq. 0) then 
         nstr = 2
      else 
         print *, 'INVALID VALUE FOR ISTRM'
         stop
      endif

      read (ird,9020) juldat, sza, isolvar, (solvar(ib),ib=ib1,ib2)

      zenith = cos(sza * pi / 180.)
      if (juldat .eq. 0) then
         adjflux_jd = 1.
      else
         adjflux_jd = earth_sun (juldat)
      endif

! If clouds are present, read in appropriate input file, IN_CLD_RRTM.
      if (icld .ge. 1) call readcld

! If aerosols are present, read in appropriate input from file, IN_AER_RRTM. 
      if (iaer .eq. 10) call readaer


      if (isolvar .eq. 0) then
         do ib = ib1,ib2
            adjflux(ib) = adjflux_jd
         enddo
      elseif (isolvar .eq. 1) then
         do ib=ib1,ib2
            adjflux(ib) = adjflux_jd * solvar(ib1)
         enddo
      elseif (isolvar .eq. 2) then
         do ib=ib1,ib2
            adjflux(ib) = adjflux_jd * solvar(ib)
         enddo
      else
         print *, 'ISOLVAR = ', isolvar, ' NOT A VALID INPUT VALUE'
         stop
      endif

      read (ird,9012) iemis, ireflect, (semiss(ib),ib=ib1,ib2)
      if (iemis .eq. 0) then
         do ib = ib1, ib2
            semiss(ib) = 1.
         enddo
      elseif (iemis .eq. 1) then
         do ib = ib1, ib2
            semiss(ib) = semiss(ib1)
         enddo
      elseif (iemis .eq. 2) then
!          print *, 'THESE ARE THE INPUT EMISSIVITY VALUES'
!          print *, semiss(ib1:ib2)
      else
          print *, 'IEMIS = ', iemis, ' NOT A VALID INPUT VALUE'
          stop
      endif
     
      if (iatm .eq. 0) then
         read (ird,9013) iform,nlayers,nmol
         if (nmol.eq.0) nmol = 7                                    
         read (ird,form1(iform)) pavel(1),tavel(1),secntk,cinp, &
              ipthak,altz(0),pz(0),tz(0),altz(1),pz(1),tz(1)
         read (ird,form3(iform)) (wkl(m,1),m=1,7), wbrodl(1)
         if(nmol .gt. 7) read (ird,form3(iform)) (wkl(m,1),m=8,nmol)

         do l = 2, nlayers
            read (ird,form2(iform)) pavel(l),tavel(l),secntk,cinp, &
                 ipthrk,altz(l),pz(l),tz(l)
            read (ird,form3(iform)) (wkl(m,l),m=1,7), wbrodl(l)
            if(nmol .gt. 7) read (ird,form3(iform)) (wkl(m,l),m=8,nmol)
         enddo
           
         if (ixsect .eq. 1) then                                 
            read (ird,9300) nxmol0
            nxmol = nxmol0
            call xsident(ird)
            read (ird,9301) iformx
     
            do l = 1, nlayers       
               read (ird,9010) cdum
               read (ird, form3(iformx)) (wx0(m,l),m=1,7),wbrodx    
               if (nxmol0 .gt. 7) read (ird,form3(iformx)) &
                  (wx0(m,l),m=8,nxmol0)
            enddo
         endif
      else
         ipu = 7
         ipr = 66
         open(unit=ipr,file='tape6',status='unknown')
         call rrtatm
         if (ixsect .eq. 1) then
            do mx = 1, nxmol0
               ixindx(mx) = ixtrans(ixindx0(mx))
            enddo
         endif
      endif

! Test for mixing ratio input.
      imix = 1
      do m = 1, nmol
         if (wkl(m,1) .gt. 1.0) then
            imix = 0
            goto 3600
         endif
      enddo
 3600 continue

      if (ixsect .eq. 1) then
         imixx = 0
         if (wx0(1,1) .le. 1.0) imixx = 1
      endif
      do l = 1, nlayers
         summol = 0.0
         do imol = 2, nmol
            summol = summol + wkl(imol,l)
         enddo
         if (imix .eq. 1) then
            coldry(l) = wbrodl(l) / (1. - summol)
            do imol = 1, nmol
               wkl(imol,l) = coldry(l) * wkl(imol,l)
            enddo
         else
            coldry(l) = wbrodl(l) + summol
         endif
         if (ixsect .eq. 1) then
            do ix = 1, nxmol0
               if (ixindx(ix) .ne. 0) then
                  if (imixx .eq. 1) then
                     wx(ixindx(ix),l) = coldry(l) * wx0(ix,l) * 1.e-20
                  else
                     wx(ixindx(ix),l) = wx0(ix,l) * 1.e-20
                  endif
               endif
            enddo 
         endif
      enddo

! Pass output arrays to new variables for transfer to rrtmg thorugh subroutine call.
      nlayers_out = nlayers
      iout_out = iout
      icld_out = icld
      iaer_out = iaer
      isccos_out = isccos
      idelm_out = idelm
      inflag_out = inflag
      iceflag_out = iceflag
      liqflag_out = liqflag

      pz_out(0) = pz(0)
      tz_out(0) = tz(0)
      tbound_out = tz(0)
      do l = 1, mxlay
         pavel_out(l) = pavel(l)
         tavel_out(l) = tavel(l)
         pz_out(l) = pz(l)
         tz_out(l) = tz(l)
         pdp(l) = (pz(l-1) - pz(l))
         coldry_out(l) = coldry(l)
         cldfrac_out(l) = cldfrac(l)
         do imol = 1,nmol
            wkl_out(imol,l) = wkl(imol,l)
         enddo
      enddo

      do l = 1, mxlay
         if (inflag.eq.0) then
            do n = 1, nbndsw
               tauc(n,l) = clddat1(l)
               ssac(n,l) = clddat2(l)
               asmc(n,l) = clddatmom(1,l)
            enddo
            ciwp(l) = 0._jprb
            clwp(l) = 0._jprb
            fice(l) = 0._jprb
            rei(l) = 0._jprb
            rel(l) = 0._jprb
         else
            do n = 1, nbndsw
               tauc(n,l) = 0._jprb
               ssac(n,l) = 1._jprb
               asmc(n,l) = 1._jprb
            enddo
            cwp = clddat1(l)
            fice(l) = clddat2(l)
            ciwp(l) = cwp * fice(l)
            clwp(l) = cwp * (1._jprb - fice(l))
            rei(l) = clddat3(l)
            rel(l) = clddat4(l)
         endif 
      enddo

      do l = 1, mxlay
         do nb = ib1,ib2
            tauaer_out(l,nb) = tauaer(l,nb)
            ssaaer_out(l,nb) = ssaaer(l,nb)
            asmaer_out(l,nb) = phase(1,l,nb)
         enddo
      enddo
      zenith_out = zenith
      do nb = ib1,ib2
        semiss_out(nb) = semiss(nb)
        adjflux_out(nb) = adjflux(nb)
      enddo
      
      goto 9000

 8800 continue
      stop ' INVALID INPUT_RRTM '

 9000 continue

 9009 format (a1,1x,i2,i2,i2)
 9010 format (a1)
 9011 format (18x,i2,29x,i1,32x,i1,1x,i1,2x,i3,3x,i1,i1,3x,i1,i1)
 9012 format (11x,i1,2x,i1,14f5.3)
 9013 format (1x,i1,i3,i5)                                     
 9020 format (12x, i3, 3x, f7.4, 4x, i1,14f5.3)
 9300 format (i5)
 9301 format (1x,i1)

      end subroutine readprof

!***************************************************************************
      subroutine readcld
!***************************************************************************

      use parkind, only : jpim, jprb

! Purpose:  To read in IN_CLD_RRTM_SW, the file that contains input 
!           cloud properties.

      implicit integer(i-n), real(a-h,o-z)

! ------- Parameters ------- 
      parameter (mxlay = 203, jpbands = 29)
      parameter (ib1 = 16, ib2 = 29)
      parameter (mg = 16)
      parameter (mxstr = 16)

      common /control/   iaer, nstr, iout, istart, iend, icld, idelm, isccos
      common /profile/   nlayers,pavel(mxlay),tavel(mxlay),pz(0:mxlay),tz(0:mxlay),tbound
      common /cloudin/   inflag,clddat1(mxlay),clddat2(mxlay), &
                         iceflag,liqflag,clddat3(mxlay),clddat4(mxlay), &
                         clddatmom(0:16,mxlay)
      common /clouddat/  ncbands,cldfrac(mxlay), &
                         taucloud(mxlay,jpbands),ssacloud(mxlay,jpbands), &
                         xmom(0:16,mxlay,jpbands)

      character*1 ctest, cpercent

      data cpercent /'%'/
      irdcld = 11

      open(irdcld,file='IN_CLD_RRTM',form='formatted')

! Read in cloud input option.  
      read(irdcld,9050) inflag, iceflag, liqflag

      do lay = 1, nlayers
         cldfrac(lay) = 0._jprb
      enddo

      if (inflag .eq. 0) then
 950     continue
!  For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
!  DAT1 are pertinent.  If CTEST = '%', then there are no more 
!  cloudy layers to process.
         read (irdcld,9099,end=8950) ctest,lay,frac,dat1,dat2,(clddatmom(istr,lay),istr=0,nstr)
         if (ctest .eq. cpercent) goto 8950
         cldfrac(lay) = frac
         clddat1(lay) = dat1
         clddat2(lay) = dat2
         goto 950
 8950    continue

      else
 1000    continue
! For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
! DAT1 are pertinent.  If CTEST = '%', then there are no more 
! cloudy layers to process.
         read (irdcld,9100,end=9000) ctest,lay,frac,dat1,dat2,dat3,dat4
         if (ctest .eq. cpercent) goto 9000
         cldfrac(lay) = frac
         clddat1(lay) = dat1
         clddat2(lay) = dat2
         clddat3(lay) = dat3
         clddat4(lay) = dat4
         goto 1000
 9000    continue
      endif

      close(irdcld)

 9050 format (3x,i2,4x,i1,4x,i1)
 9099 format (a1,1x,i3,19e10.5)
 9100 format (a1,1x,i3,5e10.5)

      end subroutine readcld

!***************************************************************************
      subroutine readaer
!***************************************************************************

      use parkind, only : jpim, jprb

! Purpose:  To read in IN_AER_RRTM, the file that contains input
!           aerosol properties.

! -------- Modules --------

      use rrsw_wvn, only : wavenum1, wavenum2

      implicit integer(i-n), real(a-h,o-z)

! ------- Parameters -------
      parameter (mxlay = 203, jpbands = 29)
      parameter (ib1 = 16, ib2 = 29)
      parameter (mg = 16)
      parameter (mxstr = 16)
      parameter (mcmu = 32)

      real aerpar(3), ssa(jpbands), asym(jpbands), aod(mxlay),aod1(jpbands)
      real rlambda(jpbands), specfac(jpbands)
      real rnu0(16:29),rnu1(23:26)
      real f1(23:26),od0(23:26),od1(23:26)
      integer lay(mxlay),ivec(mxlay)

      common /control/ iaer, nstr, iout, istart, iend, icld, idelm, isccos
      common /profile/ nlayers,pavel(mxlay),tavel(mxlay), pz(0:mxlay),tz(0:mxlay),tbound
      common /swprop/  zenith, albedo, adjflux(jpbands)
      common /aerdat/  ssaaer(mxlay,jpbands), phase(mcmu,mxlay,jpbands), tauaer(mxlay,jpbands)

      character*1 ctest, cpercent

      data cpercent /'%'/

      data rnu0 /2903._jprb,3601._jprb,4310._jprb,4892._jprb,5623._jprb,6872._jprb, &
                 7872._jprb,10590._jprb,14420._jprb,18970._jprb,25015._jprb,30390._jprb, & 
                 43507._jprb,1412._jprb/
      data rnu1 /10530.7_jprb,14293.3_jprb,18678.0_jprb,24475.1_jprb/
 
      data f1  /0.9929_jprb,0.9883_jprb,0.978_jprb,0.9696_jprb/
      data od0 /0.1084_jprb,0.167_jprb,0.245_jprb,0.3611_jprb/
      data od1 /0.3144_jprb,0.4822_jprb,0.7013_jprb,1.0239_jprb/

      eps = 1.e-10_jprb
      irdaer = 12
      open(irdaer,file='IN_AER_RRTM',form='formatted')

      aod(:) = 0.0_jprb
      tauaer(:,ib1:ib2) = 0.0_jprb

! Read in number of different aerosol models option.
      read (irdaer, 9010) naer
!       if (naer .gt. 4) then
!          print *, 'NAER (= ', naer, ') IS GREATER THAN 4'
!          stop
!       endif
        
! For each aerosol read in optical properties and layer aerosol 
! optical depths.
      do ia = 1, naer
	 read (irdaer, 9011) nlay, iaod, issa, iasym, (aerpar(i),i=1,3)

         if (iaod .eq. 0) then
! Set defaults to get standard Angstrom relation.
            if (aerpar(2) .lt. eps) aerpar(2) = 1._jprb

            do ib = ib1, ib2
  	       rlambda(ib) = 10000._jprb/rnu0(ib)
               specfac(ib) = (aerpar(2) + aerpar(3) * rlambda(ib)) / &
                            ((aerpar(2) + aerpar(3) - 1._jprb) + &
                              rlambda(ib)**aerpar(1))
            enddo
         endif

! For this aerosol, read in layers and optical depth information.
! Store a nonzero optical depth in aod to check for double specification.
         do il = 1, nlay
            read(irdaer, 9012) lay(il), (aod1(ib), ib = ib1,ib2)
            if (aod(lay(il)) .lt. eps) then
               if (iaod .eq. 0) then
                  aod(lay(il)) = aod1(ib1)
                  do ib = ib1, ib2
                     tauaer(lay(il),ib) = aod(lay(il)) * specfac(ib)
                  enddo
               else
                  do ib = ib1, ib2
                     aod(lay(il)) = max(aod(lay(il)),aod1(ib))
                     tauaer(lay(il),ib) = aod1(ib)
                  enddo
               endif
            else
               print *,'LAYER ',lay(il),' HAS MORE THAN ONE AEROSOL TYPE'
               stop
            endif
         enddo

! Build vector of aerosol layer indices 

         do il=1,nlay
            ivec(il) = lay(il) 
         end do

! Correct bands 23 through 26 for sza effect (negligible for others)
         do ib=23,26
            if (iaod.eq.0) then
                od = sum(tauaer(ivec(1:nlay),ib))/zenith
                rnu = rnu0(ib) + &
                     (rnu1(ib)-rnu0(ib))*(od-od0(ib))/(od1(ib)-od0(ib))
               rlambda_new = 10000._jprb/rnu
               specfac_new = (aerpar(2)+aerpar(3)*rlambda_new) / &
                  ((aerpar(2)+aerpar(3)- 1.)+rlambda_new**aerpar(1))
               do il=1,nlay
                  tauaer(lay(il),ib) = tauaer(lay(il),ib) * &
                         specfac_new/specfac(ib)
               end do
            endif
         end do

! For this aerosol, read and store optical properties
         read (irdaer, 9013) (ssa(ib), ib = ib1,ib2)

         do ib = ib1, ib2
            do il = 1, nlay
               if (issa .eq. 0) then 
                  ssaaer(lay(il),ib) = ssa(ib1)
               else
                  ssaaer(lay(il),ib) = ssa(ib)
               endif
            enddo
         enddo

         if (iasym .lt. 2) then
            read (irdaer, 9013) (asym(ib), ib = ib1,ib2)

            do ib = ib1, ib2
               do il = 1, nlay
                  do istr = 1,  nstr
                     if (iasym .eq. 0) then 
                        phase(istr,lay(il),ib) = asym(ib1)**istr
                     elseif (iasym .eq. 1) then
                        phase(istr,lay(il),ib) = asym(ib)**istr
                     endif
                  enddo
               enddo
            enddo
         else
            do il = 1, nlay
               do istr = 1, nstr
                  read (irdaer, 9013) (phase(istr,lay(il),ib),ib = ib1,ib2)
               enddo
            enddo
         endif

! End of naer loop
      enddo

 9000 continue
      close(irdaer)

 9010 format (3x, i2)
 9011 format (2x, i3, 4x, i1, 4x, i1, 4x, i1, 3f8.3)
 9012 format (2x, i3, 14f7.4)
 9013 format (14f5.3)

      end subroutine readaer

!**********************************************************************
      subroutine xsident(ird)
!**********************************************************************

! Purpose:  This subroutine identifies which cross-sections are to be used.

      implicit integer(i-n), real(a-h,o-z)

! ------- Parameters -------
      parameter (maxinpx=35)
      parameter (maxxsec=4)

      common /xsecctrl/ nxmol,ixindx(maxinpx)
                                                                         
!     nxmol     - number of cross-sections input by user
!     ixindx(i) - index of cross-section molecule corresponding to Ith
!                 cross-section specified by user
!                 = 0 -- not allowed in rrtm
!                 = 1 -- ccl4
!                 = 2 -- cfc11
!                 = 3 -- cfc12
!                 = 4 -- cfc22
!                                                                         
!     xsname=names, alias=aliases of the cross-section molecules          
!                                                                         
      character*10 xsname(maxinpx),alias(maxxsec,4),blank               

      data (alias(1,i),i=1,4)/ &
         'CCL4      ', 'CCL3F     ', 'CCL2F2    ', 'CHCLF2    '/ 
      data (alias(2,i),i=1,4)/ &
         ' ZZZZZZZZ ', 'CFCL3     ', 'CF2CL2    ', 'CHF2CL    '/         
      data (alias(3,i),i=1,4)/ &
         ' ZZZZZZZZ ', 'CFC11     ', 'CFC12     ', 'CFC22     '/         
      data (alias(4,i),i=1,4)/ &
         ' ZZZZZZZZ ', 'F11       ', 'F12       ', 'F22       '/        

      data blank / '          '/

      do i = 1, nxmol
         xsname(i) = blank
      enddo

! Read in the names of the molecules

      if (nxmol.gt.7) then                                               
         read (ird,'(7a10)') (xsname(i),i=1,7)                            
         read (ird,'(8a10)') (xsname(i),i=8,nxmol)                       
      else                                                                
         read (ird,'(7a10)') (xsname(i),i=1,nxmol)                       
      endif                                                               
                                                                         
!  Match the names read in against the names stored in alias           
!  and determine the index value.  
      ixmax = 4
      do i = 1, nxmol
!  Left-justify all inputed names.                                      
         call cljust (xsname(i),10)
         ixindx(i) = 0
         do j = 1, ixmax
            if ((xsname(i).eq.alias(1,j)) .or. &
                (xsname(i).eq.alias(2,j)) .or. &
                (xsname(i).eq.alias(3,j)) .or. &                           
                (xsname(i).eq.alias(4,j))) then                           
               ixindx(i) = j                                              
            endif                                                         
         enddo
      enddo   

      end subroutine xsident

      end program rrtmg_sw


