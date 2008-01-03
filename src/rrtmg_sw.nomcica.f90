!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!

       module rrtmg_sw_rad

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
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                      email:  miacono@aer.com                             *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people:  Steven J. Taubman, Patrick D. Brown,            *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

! --------- Modules ---------
      use parkind, only : jpim, jprb
      use rrsw_vsn
      use mcica_subcol_gen_sw, only: mcica_subcol_sw
      use rrtmg_sw_cldprop, only: cldprop_sw
      use rrtmg_sw_cldprmc, only: cldprmc_sw
! Move call to rrtmg_sw_ini and following use association to 
! GCM initialization area
!      use rrtmg_sw_init, only: rrtmg_sw_ini
      use rrtmg_sw_setcoef, only: setcoef_sw
      use rrtmg_sw_spcvrt, only: spcvrt_sw

      implicit none

! public interfaces/functions/subroutines
      public :: rrtmg_sw, inatm_sw, earth_sun

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine rrtmg_sw &
            (ncol    ,nlay    ,icld    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr , &
             o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  , &
             asdir   ,asdif   ,aldir   ,aldif   , &
             coszen  ,adjes   ,dyofyr  ,scon    , &
             inflgsw ,iceflgsw,liqflgsw,cldfr   , &
             taucld  ,ssacld  ,asmcld  , &
             cicewp  ,cliqwp  ,reice   ,reliq   , &
             tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc)

! ------- Description -------

! This program is the driver for RRTMG_SW, the AER SW radiation model for 
!  application to GCMs, that has been adapted from RRTM_SW for improved
!  efficiency and to provide fractional cloudiness and cloud overlap
!  capability using McICA.
!
! Note: The call to RRTMG_SW_INI should be moved to the GCM initialization 
!  area, since this has to be called only once. 
!
! This routine
!    b) calls INATM_SW to read in the atmospheric profile from GCM;
!       all layering in RRTMG is ordered from surface to toa. 
!    c) calls CLDPROP_SW to set cloud optical depth based on input
!       cloud properties
!    d) calls SETCOEF_SW to calculate various quantities needed for 
!       the radiative transfer algorithm
!    e) calls SPCVRT to call the two-stream model that in turn 
!       calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands and to perform the radiative transfer;
!    f) passes the calculated fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_sw.nomcica.f90 (to not use
!     McICA) or rrtmg_sw.f90 (to use McICA) to interface with a GCM.
!
!    1) Standard, single forward model calculation (imca = 0); this is 
!       valid only for clear sky or fully overcast clouds
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!       This method is valid for clear sky and full or partial cloud conditions.
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflag, iceflag and liqflag; see text file rrtmg_sw_instructions
!     and subroutine rrtmg_sw_cldprop.f90 for further details):
!
!    1) Input cloud fraction, cloud optical depth, single scattering albedo 
!       and asymmetry parameter directly (inflgsw = 0)
!    2) Input cloud fraction and cloud physical properties: ice fracion,
!       ice and liquid particle sizes (inflgsw = 1 or 2);  
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
!-- Reformatted for consistency with rrtmg_lw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays. 
!     Aug 2007: M. J. Iacono, AER, Inc.

! --------- Modules ---------

      use parrrsw, only : nbndsw, ngptsw, naerec, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi
      use rrsw_wvn, only : wavenum1, wavenum2

! ------- Declarations

! ----- Input -----
      integer(kind=jpim), intent(in) :: ncol            ! Number of horizontal columns     
      integer(kind=jpim), intent(in) :: nlay            ! Number of model layers
      integer(kind=jpim), intent(inout) :: icld         ! Cloud overlap method
                                                        !    0: Clear only
                                                        !    1: Random
                                                        !    2: Maximum/random
                                                        !    3: Maximum

      real(kind=jprb), intent(in) :: play(:,:)          ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: plev(:,:)          ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=jprb), intent(in) :: tlay(:,:)          ! Layer temperatures (K)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: tlev(:,:)          ! Interface temperatures (K)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=jprb), intent(in) :: tsfc(:)            ! Surface temperature (K)
                                                        !    Dimensions: (ncol)
      real(kind=jprb), intent(in) :: h2ovmr(:,:)        ! H2O volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: o3vmr(:,:)         ! O3 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: co2vmr(:,:)        ! CO2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: ch4vmr(:,:)        ! Methane volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: n2ovmr(:,:)        ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: asdir(:)           ! UV/vis surface albedo direct rad
                                                        !    Dimensions: (ncol)
      real(kind=jprb), intent(in) :: aldir(:)           ! Near-IR surface albedo direct rad
                                                        !    Dimensions: (ncol)
      real(kind=jprb), intent(in) :: asdif(:)           ! UV/vis surface albedo: diffuse rad
                                                        !    Dimensions: (ncol)
      real(kind=jprb), intent(in) :: aldif(:)           ! Near-IR surface albedo: diffuse rad
                                                        !    Dimensions: (ncol)

      integer(kind=jpim), intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
                                                        !  distance if adjflx not provided)
      real(kind=jprb), intent(in) :: adjes              ! Flux adjustment for Earth/Sun distance
      real(kind=jprb), intent(in) :: coszen(:)          ! Cosine of solar zenith angle
                                                        !    Dimensions: (ncol)
      real(kind=jprb), intent(in) :: scon               ! Solar constant (W/m2)

      integer(kind=jpim), intent(in) :: inflgsw         ! Flag for cloud optical properties
      integer(kind=jpim), intent(in) :: iceflgsw        ! Flag for ice particle specification
      integer(kind=jpim), intent(in) :: liqflgsw        ! Flag for liquid droplet specification

      real(kind=jprb), intent(in) :: cldfr(:,:)         ! Cloud fraction
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: taucld(:,:,:)      ! Cloud optical depth
                                                        !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=jprb), intent(in) :: ssacld(:,:,:)      ! Cloud single scattering albedo
                                                        !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=jprb), intent(in) :: asmcld(:,:,:)      ! Cloud asymmetry parameter
                                                        !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=jprb), intent(in) :: cicewp(:,:)        ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: cliqwp(:,:)        ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: reice(:,:)         ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: reliq(:,:)         ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: tauaer(:,:,:)      ! Aerosol optical depth (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
      real(kind=jprb), intent(in) :: ssaaer(:,:,:)      ! Aerosol single scattering albedo (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
      real(kind=jprb), intent(in) :: asmaer(:,:,:)      ! Aerosol asymmetry parameter (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
      real(kind=jprb), intent(in) :: ecaer(:,:,:)       ! Aerosol optical depth at 0.55 micron (iaer=6 only)
                                                        !    Dimensions: (ncol,nlay,naerec)
                                                        ! (non-delta scaled)      

! ----- Output -----

      real(kind=jprb), intent(out) :: swuflx(:,:)       ! Total sky shortwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=jprb), intent(out) :: swdflx(:,:)       ! Total sky shortwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=jprb), intent(out) :: swhr(:,:)         ! Total sky shortwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=jprb), intent(out) :: swuflxc(:,:)      ! Clear sky shortwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=jprb), intent(out) :: swdflxc(:,:)      ! Clear sky shortwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=jprb), intent(out) :: swhrc(:,:)        ! Clear sky shortwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)

! ----- Local -----

! Control
      integer(kind=jpim) :: nlayers             ! total number of layers
      integer(kind=jpim) :: istart              ! beginning band of calculation
      integer(kind=jpim) :: iend                ! ending band of calculation
      integer(kind=jpim) :: icpr                ! cldprop/cldprmc use flag
      integer(kind=jpim) :: iout                ! output option flag (inactive)
      integer(kind=jpim) :: iaer                ! aerosol option flag
      integer(kind=jpim) :: idelm               ! delta-m scaling flag (inactive)
      integer(kind=jpim) :: isccos              ! instrumental cosine response flag (inactive)
      integer(kind=jpim) :: iplon               ! column loop index
      integer(kind=jpim) :: i                   ! layer loop index                       ! jk
      integer(kind=jpim) :: ib                  ! band loop index                        ! jsw
      integer(kind=jpim) :: ia, ig              ! indices
      integer(kind=jpim) :: k                   ! layer loop index
      integer(kind=jpim) :: ims                 ! value for changing mcica permute seed
      integer(kind=jpim) :: imca                ! flag for mcica [0=off, 1=on]

      real(kind=jprb) :: zepsec, zepzen         ! epsilon
      real(kind=jprb) :: zdpgcp                 ! flux to heating conversion ratio

! Atmosphere
      real(kind=jprb) :: pavel(nlay+1)          ! layer pressures (mb) 
      real(kind=jprb) :: tavel(nlay+1)          ! layer temperatures (K)
      real(kind=jprb) :: pz(0:nlay+1)           ! level (interface) pressures (hPa, mb)
      real(kind=jprb) :: tz(0:nlay+1)           ! level (interface) temperatures (K)
      real(kind=jprb) :: tbound                 ! surface temperature (K)
      real(kind=jprb) :: pdp(nlay+1)            ! layer pressure thickness (hPa, mb)
      real(kind=jprb) :: coldry(nlay+1)         ! dry air column amount
      real(kind=jprb) :: wkl(mxmol,nlay+1)      ! molecular amounts (mol/cm-2)

!      real(kind=jprb) :: earth_sun             ! function for Earth/Sun distance factor
      real(kind=jprb) :: cossza                 ! Cosine of solar zenith angle
      real(kind=jprb) :: adjflux(jpband)        ! adjustment for current Earth/Sun distance
      real(kind=jprb) :: solvar(jpband)         ! solar constant scaling factor from rrtmg_sw
                                                !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=jprb) :: albdir(nbndsw)         ! surface albedo, direct          ! zalbp
      real(kind=jprb) :: albdif(nbndsw)         ! surface albedo, diffuse         ! zalbd

      real(kind=jprb) :: taua(nlay+1,nbndsw)    ! Aerosol optical depth
      real(kind=jprb) :: ssaa(nlay+1,nbndsw)    ! Aerosol single scattering albedo
      real(kind=jprb) :: asma(nlay+1,nbndsw)    ! Aerosol asymmetry parameter

! Atmosphere - setcoef
      integer(kind=jpim) :: laytrop             ! tropopause layer index
      integer(kind=jpim) :: layswtch            ! tropopause layer index
      integer(kind=jpim) :: laylow              ! tropopause layer index
      integer(kind=jpim) :: jp(nlay+1)          ! 
      integer(kind=jpim) :: jt(nlay+1)          !
      integer(kind=jpim) :: jt1(nlay+1)         !

      real(kind=jprb) :: colh2o(nlay+1)         ! column amount (h2o)
      real(kind=jprb) :: colco2(nlay+1)         ! column amount (co2)
      real(kind=jprb) :: colo3(nlay+1)          ! column amount (o3)
      real(kind=jprb) :: coln2o(nlay+1)         ! column amount (n2o)
      real(kind=jprb) :: colch4(nlay+1)         ! column amount (ch4)
      real(kind=jprb) :: colo2(nlay+1)          ! column amount (o2)
      real(kind=jprb) :: colmol(nlay+1)         ! column amount
      real(kind=jprb) :: co2mult(nlay+1)        ! column amount 

      integer(kind=jpim) :: indself(nlay+1)
      integer(kind=jpim) :: indfor(nlay+1)
      real(kind=jprb) :: selffac(nlay+1)
      real(kind=jprb) :: selffrac(nlay+1)
      real(kind=jprb) :: forfac(nlay+1)
      real(kind=jprb) :: forfrac(nlay+1)

      real(kind=jprb) :: &                      !
                         fac00(nlay+1), fac01(nlay+1), &
                         fac10(nlay+1), fac11(nlay+1) 

! Atmosphere/clouds - cldprop
      integer(kind=jpim) :: ncbands             ! number of cloud spectral bands
      integer(kind=jpim) :: inflag              ! flag for cloud property method
      integer(kind=jpim) :: iceflag             ! flag for ice cloud properties
      integer(kind=jpim) :: liqflag             ! flag for liquid cloud properties

      real(kind=jprb) :: cldfrac(nlay+1)        ! layer cloud fraction
      real(kind=jprb) :: tauc(nbndsw,nlay+1)    ! cloud optical depth (non-delta scaled)
      real(kind=jprb) :: ssac(nbndsw,nlay+1)    ! cloud single scattering albedo (non-delta scaled)
      real(kind=jprb) :: asmc(nbndsw,nlay+1)    ! cloud asymmetry parameter (non-delta scaled)
      real(kind=jprb) :: ciwp(nlay+1)           ! cloud ice water path
      real(kind=jprb) :: clwp(nlay+1)           ! cloud liquid water path
      real(kind=jprb) :: rel(nlay+1)            ! cloud liquid particle effective radius (microns)
      real(kind=jprb) :: rei(nlay+1)            ! cloud ice particle effective radius (microns)
      real(kind=jprb) :: dge(nlay+1)            ! cloud ice particle generalized effective size (microns)

      real(kind=jprb) :: taucloud(nlay+1,jpband)  ! cloud optical depth
      real(kind=jprb) :: taucldorig(nlay+1,jpband)! cloud optical depth (non-delta scaled)
      real(kind=jprb) :: ssacloud(nlay+1,jpband)  ! cloud single scattering albedo
      real(kind=jprb) :: asmcloud(nlay+1,jpband)  ! cloud asymmetry parameter

! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      real(kind=jprb) :: ztauc(nlay+1,nbndsw)     ! cloud optical depth
      real(kind=jprb) :: ztaucorig(nlay+1,nbndsw) ! unscaled cloud optical depth
      real(kind=jprb) :: zasyc(nlay+1,nbndsw)     ! cloud asymmetry parameter 
                                                  !  (first moment of phase function)
      real(kind=jprb) :: zomgc(nlay+1,nbndsw)     ! cloud single scattering albedo
      real(kind=jprb) :: ztaua(nlay+1,nbndsw)     ! total aerosol optical depth
      real(kind=jprb) :: zasya(nlay+1,nbndsw)     ! total aerosol asymmetry parameter 
      real(kind=jprb) :: zomga(nlay+1,nbndsw)     ! total aerosol single scattering albedo

      real(kind=jprb) :: zbbfu(nlay+2)          ! temporary upward shortwave flux (w/m2)
      real(kind=jprb) :: zbbfd(nlay+2)          ! temporary downward shortwave flux (w/m2)
      real(kind=jprb) :: zbbcu(nlay+2)          ! temporary clear sky upward shortwave flux (w/m2)
      real(kind=jprb) :: zbbcd(nlay+2)          ! temporary clear sky downward shortwave flux (w/m2)
      real(kind=jprb) :: zbbfddir(nlay+2)       ! temporary downward direct shortwave flux (w/m2)
      real(kind=jprb) :: zbbcddir(nlay+2)       ! temporary clear sky downward direct shortwave flux (w/m2)
      real(kind=jprb) :: zuvfd(nlay+2)          ! temporary UV downward shortwave flux (w/m2)
      real(kind=jprb) :: zuvcd(nlay+2)          ! temporary clear sky UV downward shortwave flux (w/m2)
      real(kind=jprb) :: zuvfddir(nlay+2)       ! temporary UV downward direct shortwave flux (w/m2)
      real(kind=jprb) :: zuvcddir(nlay+2)       ! temporary clear sky UV downward direct shortwave flux (w/m2)
      real(kind=jprb) :: znifd(nlay+2)          ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=jprb) :: znicd(nlay+2)          ! temporary clear sky near-IR downward shortwave flux (w/m2)
      real(kind=jprb) :: znifddir(nlay+2)       ! temporary near-IR downward direct shortwave flux (w/m2)
      real(kind=jprb) :: znicddir(nlay+2)       ! temporary clear sky near-IR downward direct shortwave flux (w/m2)

! Optional output fields 
      real(kind=jprb) :: swnflx(nlay+2)         ! Total sky shortwave net flux (W/m2)
      real(kind=jprb) :: swnflxc(nlay+2)        ! Clear sky shortwave net flux (W/m2)
      real(kind=jprb) :: dirdflux(nlay+2)       ! Direct downward shortwave surface flux
      real(kind=jprb) :: difdflux(nlay+2)       ! Diffuse downward shortwave surface flux
      real(kind=jprb) :: uvdflx(nlay+2)         ! Total sky downward shortwave flux, UV/vis  
      real(kind=jprb) :: nidflx(nlay+2)         ! Total sky downward shortwave flux, near-IR 
      real(kind=jprb) :: dirdnuv(nlay+2)        ! Direct downward shortwave flux, UV/vis
      real(kind=jprb) :: difdnuv(nlay+2)        ! Diffuse downward shortwave flux, UV/vis
      real(kind=jprb) :: dirdnir(nlay+2)        ! Direct downward shortwave flux, near-IR
      real(kind=jprb) :: difdnir(nlay+2)        ! Diffuse downward shortwave flux, near-IR

! Output - inactive
!      real(kind=jprb) :: zuvfu(nlay+2)         ! temporary upward UV shortwave flux (w/m2)
!      real(kind=jprb) :: zuvfd(nlay+2)         ! temporary downward UV shortwave flux (w/m2)
!      real(kind=jprb) :: zuvcu(nlay+2)         ! temporary clear sky upward UV shortwave flux (w/m2)
!      real(kind=jprb) :: zuvcd(nlay+2)         ! temporary clear sky downward UV shortwave flux (w/m2)
!      real(kind=jprb) :: zvsfu(nlay+2)         ! temporary upward visible shortwave flux (w/m2)
!      real(kind=jprb) :: zvsfd(nlay+2)         ! temporary downward visible shortwave flux (w/m2)
!      real(kind=jprb) :: zvscu(nlay+2)         ! temporary clear sky upward visible shortwave flux (w/m2)
!      real(kind=jprb) :: zvscd(nlay+2)         ! temporary clear sky downward visible shortwave flux (w/m2)
!      real(kind=jprb) :: znifu(nlay+2)         ! temporary upward near-IR shortwave flux (w/m2)
!      real(kind=jprb) :: znifd(nlay+2)         ! temporary downward near-IR shortwave flux (w/m2)
!      real(kind=jprb) :: znicu(nlay+2)         ! temporary clear sky upward near-IR shortwave flux (w/m2)
!      real(kind=jprb) :: znicd(nlay+2)         ! temporary clear sky downward near-IR shortwave flux (w/m2)


! Initializations

      zepsec = 1.e-06_jprb
      zepzen = 1.e-10_jprb
      oneminus = 1.0_jprb - zepsec
      pi = 2._jprb * asin(1._jprb)

      istart = jpb1
      iend = jpb2
      icpr = 0

! In a GCM with or without McICA, set nlon to the longitude dimension
!
! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation (clear and overcast only)
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability
!            (clear, overcast or partial cloud conditions)

! *** This version does not use McICA (imca = 0) ***

! Set icld to select of clear or cloud calculation and cloud 
! overlap method (read by subroutine readprof from input file INPUT_RRTM):  
! Without McICA, SW calculation is limited to clear or fully overcast conditions. 
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap (McICA only)
! icld = 2, with clouds using maximum/random cloud overlap (McICA only)
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      if (icld.lt.0.or.icld.gt.3) icld = 2

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! iaer = 6, use six ECMWF aerosol types
!           input aerosol optical depth at 0.55 microns for each aerosol type (ecaer)
! iaer = 10, input total aerosol optical depth, single scattering albedo 
!            and asymmetry parameter (tauaer, ssaaer, asmaer) directly
      iaer = 0

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 224 to 112 for input absorption
! coefficient data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  
!      call rrtmg_sw_ini

! This is the main longitude/column loop in RRTMG.
! Modify to loop over all columns (nlon) or over daylight columns

      do iplon = 1, ncol

! Prepare atmosphere profile from GCM for use in RRTMG, and define
! other input parameters

         call inatm_sw (iplon, nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, &
              h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, adjes, dyofyr, scon, &
              inflgsw, iceflgsw, liqflgsw, &
              cldfr, taucld, ssacld, asmcld, cicewp, cliqwp, &
              reice, reliq, tauaer, ssaaer, asmaer, &
              nlayers, pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
              adjflux, solvar, inflag, iceflag, liqflag, cldfrac, tauc, &
              ssac, asmc, ciwp, clwp, rei, dge, rel, taua, ssaa, asma)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed in cldprop.  Cloud fraction and cloud
!  optical properties are transferred to rrtmg_sw arrays in cldprop.  

!  Without McICA, SW calculation is limited to clear or fully overcast conditions. 
!  Stop model if partial cloudiness is present.  

         do i = 1, nlayers
            if (cldfrac(i).gt.zepsec .and. cldfrac(i).lt.oneminus) then
               stop 'PARTIAL CLOUD NOT ALLOWED'
            endif
         enddo
         call cldprop_sw(nlayers, inflag, iceflag, liqflag, cldfrac, &
                         tauc, ssac, asmc, ciwp, clwp, rei, dge, rel, &
                         taucldorig, taucloud, ssacloud, asmcloud)
         icpr = 1

! Calculate coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients by interpolating data from stored
! reference atmospheres.

         call setcoef_sw(nlayers, pavel, tavel, pz, tz, tbound, coldry, wkl, &
                         laytrop, layswtch, laylow, jp, jt, jt1, &
                         co2mult, colch4, colco2, colh2o, colmol, coln2o, &
                         colo2, colo3, fac00, fac01, fac10, fac11, &
                         selffac, selffrac, indself, forfac, forfrac, indfor)


! Cosine of the solar zenith angle 
!  Prevent using value of zero; ideally, SW model is not called from host model when sun 
!  is below horizon

         cossza = coszen(iplon)
         if (cossza .eq. 0._jprb) cossza = zepzen


! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer 

! Surface albedo
!  Near-IR bands 16-24 and 29 (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
         do ib=1,9
            albdir(ib) = aldir(iplon)
            albdif(ib) = aldif(iplon)
         enddo
         albdir(nbndsw) = aldir(iplon)
         albdif(nbndsw) = aldif(iplon)
!  UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
         do ib=10,13
            albdir(ib) = asdir(iplon)
            albdif(ib) = asdif(iplon)
         enddo


! Clouds
         if (icld.eq.0) then

            ztauc(:,:) = 0._jprb
            ztaucorig(:,:) = 0._jprb
            zasyc(:,:) = 0._jprb
            zomgc(:,:) = 1._jprb

         elseif (icld.ge.1) then
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

         endif   

! Aerosol
! IAER = 0: no aerosols
         if (iaer.eq.0) then

            ztaua(:,:) = 0._jprb
            zasya(:,:) = 0._jprb
            zomga(:,:) = 1._jprb

! IAER = 6: Use ECMWF six aerosol types. See rrsw_aer.f90 for details.
! Input aerosol optical thickness at 0.55 micron for each aerosol type (ecaer), 
! or set manually here for each aerosol and layer.
         elseif (iaer.eq.6) then

!            do i = 1, nlayers
!               do ia = 1, naerec
!                  ecaer(iplon,i,ia) = 1.0e-15_jprb
!               enddo
!            enddo

            do i = 1, nlayers
               do ib = 1, nbndsw
                  ztaua(i,ib) = 0._jprb
                  zasya(i,ib) = 0._jprb
                  zomga(i,ib) = 1._jprb
                  do ia = 1, naerec
                     ztaua(i,ib) = ztaua(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia)
                     zomga(i,ib) = zomga(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia) * &
                                   rsrpiza(ib,ia)
                     zasya(i,ib) = zasya(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia) * &
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

! IAER=10: Direct specification of aerosol optical properties from GCM
         elseif (iaer.eq.10) then

            do i = 1 ,nlayers
               do ib = 1 ,nbndsw
                  ztaua(i,ib) = taua(i,ib)
                  zasya(i,ib) = asma(i,ib)
                  zomga(i,ib) = ssaa(i,ib)
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

! Transfer up and down, clear and total sky fluxes to output arrays.
! Vertical indexing goes from bottom to top; reverse here for GCM if necessary.

         do i = 1, nlayers+1
            swuflxc(iplon,i) = zbbcu(i)
            swdflxc(iplon,i) = zbbcd(i)
            swuflx(iplon,i) = zbbfu(i)
            swdflx(iplon,i) = zbbfd(i)
            uvdflx(i) = zuvfd(i)
            nidflx(i) = znifd(i)
!  Direct/diffuse fluxes
            dirdflux(i) = zbbfddir(i)
            difdflux(i) = swdflx(iplon,i) - dirdflux(i)
!  UV/visible direct/diffuse fluxes
            dirdnuv(i) = zuvfddir(i)
            difdnuv(i) = zuvfd(i) - dirdnuv(i)
!  Near-IR direct/diffuse fluxes
            dirdnir(i) = znifddir(i)
            difdnir(i) = znifd(i) - dirdnir(i)
         enddo

!  Total and clear sky net fluxes
         do i = 1, nlayers+1
            swnflxc(i) = swdflxc(iplon,i) - swuflxc(iplon,i)
            swnflx(i) = swdflx(iplon,i) - swuflx(iplon,i)
         enddo

!  Total and clear sky heating rates
         do i = 1, nlayers
            zdpgcp = heatfac / pdp(i)
            swhrc(iplon,i) = (swnflxc(i+1) - swnflxc(i)) * zdpgcp
            swhr(iplon,i) = (swnflx(i+1) - swnflx(i)) * zdpgcp
         enddo
         swhrc(iplon,nlayers) = 0._jprb
         swhr(iplon,nlayers) = 0._jprb

! End longitude loop
      enddo

      end subroutine rrtmg_sw

!*************************************************************************
      real(kind=jprb) function earth_sun(idn)
!*************************************************************************
!
!  Purpose: Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  idn        : Day of the year
!  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

! ------- Modules -------

      use rrsw_con, only : pi

      integer(kind=jpim), intent(in) :: idn

      real(kind=jprb) :: gamma

      gamma = 2._jprb*pi*(idn-1)/365._jprb

! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110_jprb + .034221_jprb * cos(gamma) + .001289_jprb * sin(gamma) + &
                   .000719_jprb * cos(2._jprb*gamma) + .000077_jprb * sin(2._jprb*gamma)

      end function earth_sun

!***************************************************************************
      subroutine inatm_sw (iplon, nlay, icld, iaer, &
            play, plev, tlay, tlev, tsfc, &
            h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, adjes, dyofyr, scon, &
            inflgsw, iceflgsw, liqflgsw, &
            cldfr, taucld, ssacld, asmcld, cicewp, cliqwp, &
            reice, reliq, tauaer, ssaaer, asmaer, &
            nlayers, pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
            adjflux, solvar, inflag, iceflag, liqflag, cldfrac, tauc, &
            ssac, asmc, ciwp, clwp, rei, dge, rel, taua, ssaa, asma)
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_SW.
!  Set other RRTMG_SW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parrrsw, only : nbndsw, ngptsw, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2, rrsw_scon
      use rrsw_con, only : fluxfac, heatfac, oneminus, pi, grav, avogad
      use rrsw_wvn, only : ng, nspa, nspb, wavenum1, wavenum2, delwave

! ------- Declarations -------

! ----- Input -----
      integer(kind=jpim), intent(in) :: iplon           ! column loop index
      integer(kind=jpim), intent(in) :: nlay            ! number of model layers
      integer(kind=jpim), intent(in) :: icld            ! clear/cloud flag
      integer(kind=jpim), intent(in) :: iaer            ! aerosol option flag

      real(kind=jprb), intent(in) :: play(:,:)          ! Layer pressures (hPa, mb)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: plev(:,:)          ! Interface pressures (hPa, mb)
                                                        ! Dimensions: (ncol,nlay+1)
      real(kind=jprb), intent(in) :: tlay(:,:)          ! Layer temperatures (K)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: tlev(:,:)          ! Interface temperatures (K)
                                                        ! Dimensions: (ncol,nlay+1)
      real(kind=jprb), intent(in) :: tsfc(:)            ! Surface temperature (K)
                                                        ! Dimensions: (ncol)
      real(kind=jprb), intent(in) :: h2ovmr(:,:)        ! H2O volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: o3vmr(:,:)         ! O3 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: co2vmr(:,:)        ! CO2 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: ch4vmr(:,:)        ! Methane volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: n2ovmr(:,:)        ! Nitrous oxide volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)

      integer(kind=jpim), intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
                                                        !  distance if adjflx not provided)
      real(kind=jprb), intent(in) :: adjes              ! Flux adjustment for Earth/Sun distance
      real(kind=jprb), intent(in) :: scon               ! Solar constant (W/m2)

      integer(kind=jpim), intent(in) :: inflgsw         ! Flag for cloud optical properties
      integer(kind=jpim), intent(in) :: iceflgsw        ! Flag for ice particle specification
      integer(kind=jpim), intent(in) :: liqflgsw        ! Flag for liquid droplet specification

      real(kind=jprb), intent(in) :: cldfr(:,:)         ! Cloud fraction
                                                        ! Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: taucld(:,:,:)      ! Cloud optical depth (optional)
                                                        ! Dimensions: (nbndsw,ncol,nlay)
      real(kind=jprb), intent(in) :: ssacld(:,:,:)      ! Cloud single scattering albedo
                                                        ! Dimensions: (nbndsw,ncol,nlay)
      real(kind=jprb), intent(in) :: asmcld(:,:,:)      ! Cloud asymmetry parameter
                                                        ! Dimensions: (nbndsw,ncol,nlay)
      real(kind=jprb), intent(in) :: cicewp(:,:)        ! Cloud ice water path (g/m2)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: cliqwp(:,:)        ! Cloud liquid water path (g/m2)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: reice(:,:)         ! Cloud ice effective radius (microns)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=jprb), intent(in) :: reliq(:,:)         ! Cloud water drop effective radius (microns)
                                                        ! Dimensions: (ncol,nlay)

      real(kind=jprb), intent(in) :: tauaer(:,:,:)      ! Aerosol optical depth
                                                        ! Dimensions: (ncol,nlay,nbndsw)
      real(kind=jprb), intent(in) :: ssaaer(:,:,:)      ! Aerosol single scattering albedo
                                                        ! Dimensions: (ncol,nlay,nbndsw)
      real(kind=jprb), intent(in) :: asmaer(:,:,:)      ! Aerosol asymmetry parameter
                                                        ! Dimensions: (ncol,nlay,nbndsw)

! Atmosphere
      integer(kind=jpim), intent(out) :: nlayers        ! number of layers

      real(kind=jprb), intent(out) :: pavel(:)          ! layer pressures (mb) 
                                                        ! Dimensions: (nlay)
      real(kind=jprb), intent(out) :: tavel(:)          ! layer temperatures (K)
                                                        ! Dimensions: (nlay)
      real(kind=jprb), intent(out) :: pz(0:)            ! level (interface) pressures (hPa, mb)
                                                        ! Dimensions: (0:nlay)
      real(kind=jprb), intent(out) :: tz(0:)            ! level (interface) temperatures (K)
                                                        ! Dimensions: (0:nlay)
      real(kind=jprb), intent(out) :: tbound            ! surface temperature (K)
      real(kind=jprb), intent(out) :: pdp(:)            ! layer pressure thickness (hPa, mb)
                                                        ! Dimensions: (nlay)
      real(kind=jprb), intent(out) :: coldry(:)         ! dry air column density (mol/cm2)
                                                        ! Dimensions: (nlay)
      real(kind=jprb), intent(out) :: wkl(:,:)          ! molecular amounts (mol/cm-2)
                                                        ! Dimensions: (mxmol,nlay)

      real(kind=jprb), intent(out) :: adjflux(:)        ! adjustment for current Earth/Sun distance
                                                        ! Dimensions: (jpband)
      real(kind=jprb), intent(out) :: solvar(:)         ! solar constant scaling factor from rrtmg_sw
                                                        ! Dimensions: (jpband)
                                                        !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=jprb), intent(out) :: taua(:,:)         ! Aerosol optical depth
                                                        ! Dimensions: (nlay,nbndsw)
      real(kind=jprb), intent(out) :: ssaa(:,:)         ! Aerosol single scattering albedo
                                                        ! Dimensions: (nlay,nbndsw)
      real(kind=jprb), intent(out) :: asma(:,:)         ! Aerosol asymmetry parameter
                                                        ! Dimensions: (nlay,nbndsw)

! Atmosphere/clouds - cldprop
      integer(kind=jpim), intent(out) :: inflag         ! flag for cloud property method
      integer(kind=jpim), intent(out) :: iceflag        ! flag for ice cloud properties
      integer(kind=jpim), intent(out) :: liqflag        ! flag for liquid cloud properties

      real(kind=jprb), intent(out) :: cldfrac(:)        ! layer cloud fraction
                                                        ! Dimensions: (nlay)
      real(kind=jprb), intent(out) :: tauc(:,:)         ! cloud optical depth (non-delta scaled)
                                                        ! Dimensions: (nbndsw,nlay)
      real(kind=jprb), intent(out) :: ssac(:,:)         ! cloud single scattering albedo (non-delta-scaled)
                                                        ! Dimensions: (nbndsw,nlay)
      real(kind=jprb), intent(out) :: asmc(:,:)         ! cloud asymmetry parameter (non-delta scaled)
                                                        ! Dimensions: (nbndsw,nlay)
      real(kind=jprb), intent(out) :: ciwp(:)           ! cloud ice water path
                                                        ! Dimensions: (nlay)
      real(kind=jprb), intent(out) :: clwp(:)           ! cloud liquid water path
                                                        ! Dimensions: (nlay)
      real(kind=jprb), intent(out) :: rel(:)            ! cloud liquid particle effective radius (microns)
                                                        ! Dimensions: (nlay)
      real(kind=jprb), intent(out) :: rei(:)            ! cloud ice particle effective radius (microns)
                                                        ! Dimensions: (nlay)
      real(kind=jprb), intent(out) :: dge(:)            ! cloud ice particle generalized effective size (microns)
                                                        ! Dimensions: (nlay)

! ----- Local -----
      real(kind=jprb), parameter :: amd = 28.9660_jprb   ! Effective molecular weight of dry air (g/mol)
      real(kind=jprb), parameter :: amw = 18.0160_jprb   ! Molecular weight of water vapor (g/mol)
!      real(kind=jprb), parameter :: amc = 44.0098_jprb   ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=jprb), parameter :: amo = 47.9998_jprb   ! Molecular weight of ozone (g/mol)
!      real(kind=jprb), parameter :: amo2 = 31.9999_jprb  ! Molecular weight of oxygen (g/mol)
!      real(kind=jprb), parameter :: amch4 = 16.0430_jprb ! Molecular weight of methane (g/mol)
!      real(kind=jprb), parameter :: amn2o = 44.0128_jprb ! Molecular weight of nitrous oxide (g/mol)

! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=jprb), parameter :: amdw = 1.607793_jprb  ! Molecular weight of dry air / water vapor
      real(kind=jprb), parameter :: amdc = 0.658114_jprb  ! Molecular weight of dry air / carbon dioxide
      real(kind=jprb), parameter :: amdo = 0.603428_jprb  ! Molecular weight of dry air / ozone
      real(kind=jprb), parameter :: amdm = 1.805423_jprb  ! Molecular weight of dry air / methane
      real(kind=jprb), parameter :: amdn = 0.658090_jprb  ! Molecular weight of dry air / nitrous oxide
      real(kind=jprb), parameter :: amdo2 = 0.905140_jprb ! Molecular weight of dry air / oxygen

      real(kind=jprb), parameter :: sbc = 5.67e-08_jprb   ! Stefan-Boltzmann constant (W/m2K4)
      real(kind=jprb), parameter :: o2mmr = 0.23143_jprb  ! o2 mass mixing ratio

      integer(kind=jpim) :: isp, l, ix, n, imol, ib       ! Loop indices
      real(kind=jprb) :: amm, summol                      ! 
      real(kind=jprb) :: adjflx                           ! flux adjustment for Earth/Sun distance
      real(kind=jprb) :: earth_sun                        ! function for Earth/Sun distance adjustment

! Add one to nlayers here to include extra model layer at top of atmosphere
      nlayers = nlay

!  Initialize all molecular amounts to zero here, then pass input amounts
!  into RRTM array WKL below.

      wkl(:,:) = 0.0_jprb
      cldfrac(:) = 0.0_jprb
      tauc(:,:) = 0.0_jprb
      ssac(:,:) = 1.0_jprb
      asmc(:,:) = 0.0_jprb
      ciwp(:) = 0.0_jprb
      clwp(:) = 0.0_jprb
      rei(:) = 0.0_jprb
      dge(:) = 0.0_jprb
      rel(:) = 0.0_jprb
      taua(:,:) = 0.0_jprb
      ssaa(:,:) = 1.0_jprb
      asma(:,:) = 0.0_jprb
 
! Set flux adjustment for current Earth/Sun distance (two options).
! 1) Use Earth/Sun distance flux adjustment provided by GCM (input as adjes);
      adjflx = adjes
!
! 2) Calculate Earth/Sun distance from DYOFYR, the cumulative day of the year.
!    (Set adjflx to 1. to use constant Earth/Sun distance of 1 AU). 
      if (dyofyr .gt. 0) then
         adjflx = earth_sun(dyofyr)
      endif

! Set incoming solar flux adjustment to include adjustment for
! current Earth/Sun distance (ADJFLX) and scaling of default internal
! solar constant (rrsw_scon = 1368.22 Wm-2) by band (SOLVAR).  SOLVAR can be set 
! to a single scaling factor as needed, or to a different value in each 
! band, which may be necessary for paleoclimate simulations. 
! 
      do ib = jpb1,jpb2
!         solvar(ib) = 1._jprb
         solvar(ib) = scon / rrsw_scon
         adjflux(ib) = adjflx * solvar(ib)
      enddo

!  Set surface temperature.
      tbound = tsfc(iplon)

!  Install input GCM arrays into RRTMG_SW arrays for pressure, temperature,
!  and molecular amounts.  
!  Pressures are input in mb, or are converted to mb here.
!  Molecular amounts are input in volume mixing ratio, or are converted from 
!  mass mixing ratio (or specific humidity for h2o) to volume mixing ratio
!  here. These are then converted to molecular amount (molec/cm2) below.  
!  The dry air column COLDRY (in molec/cm2) is calculated from the level 
!  pressures, pz (in mb), based on the hydrostatic equation and includes a 
!  correction to account for h2o in the layer.  The molecular weight of moist 
!  air (amm) is calculated for each layer.  
!  Note: In RRTMG, layer indexing goes from bottom to top, and coding below
!  assumes GCM input fields are also bottom to top. Input layer indexing
!  from GCM fields should be reversed here if necessary.

      pz(0) = plev(iplon,1)
      tz(0) = tlev(iplon,1)
      do l = 1, nlayers
         pavel(l) = play(iplon,l)
         tavel(l) = tlay(iplon,l)
         pz(l) = plev(iplon,l+1)
         tz(l) = tlev(iplon,l+1)
         pdp(l) = pz(l-1) - pz(l)
! For h2o input in vmr:
         wkl(1,l) = h2ovmr(iplon,l)
! For h2o input in mmr:
!         wkl(1,l) = h2o(iplon,l)*amdw
! For h2o input in specific humidity;
!         wkl(1,l) = (h2o(iplon,l)/(1._jprb - h2o(iplon,l)))*amdw
         wkl(2,l) = co2vmr(iplon,l)
         wkl(3,l) = o3vmr(iplon,l)
         wkl(4,l) = n2ovmr(iplon,l)
         wkl(6,l) = ch4vmr(iplon,l)
         wkl(7,l) = o2mmr*amdo2
         amm = (1._jprb - wkl(1,l)) * amd + wkl(1,l) * amw            
         coldry(l) = (pz(l-1)-pz(l)) * 1.e3_jprb * avogad / &
                     (1.e2_jprb * grav * amm * (1._jprb + wkl(1,l)))
      enddo

! The following section can be used to set values for an additional layer (from
! the GCM top level to 1.e-4 mb) for improved calculation of TOA fluxes. 
! Temperature and molecular amounts in the extra model layer are set to 
! their values in the top GCM model layer, though these can be modified
! here if necessary. 
! If this feature is utilized, increase nlayers by one above, limit the two
! loops above to (nlayers-1), and set the top most (nlayers) layer values here. 

!      pavel(nlayers) = 0.5_jprb * pz(nlayers-1)
!      tavel(nlayers) = tavel(nlayers-1)
!      pz(nlayers) = 1.e-4jprb
!      tz(nlayers-1) = 0.5_jprb * (tavel(nlayers)+tavel(nlayers-1))
!      tz(nlayers) = tz(nlayers-1)
!      pdp(nlayers) = pz(nlayers-1) - pz(nlayers)
!      wkl(1,nlayers) = wkl(1,nlayers-1)
!      wkl(2,nlayers) = wkl(2,nlayers-1)
!      wkl(3,nlayers) = wkl(3,nlayers-1)
!      wkl(4,nlayers) = wkl(4,nlayers-1)
!      wkl(6,nlayers) = wkl(6,nlayers-1)
!      wkl(7,nlayers) = wkl(7,nlayers-1)
!      amm = (1._jprb - wkl(1,nlayers-1)) * amd + wkl(1,nlayers-1) * amw
!      coldry(nlayers) = (pz(nlayers-1)) * 1.e3_jprb * avogad / &
!                        (1.e2_jprb * grav * amm * (1._jprb + wkl(1,nlayers-1)))

! At this point all molecular amounts in wkl are in volume mixing ratio; 
! convert to molec/cm2 based on coldry for use in rrtm.  

      do l = 1, nlayers
         do imol = 1, nmol
            wkl(imol,l) = coldry(l) * wkl(imol,l)
         enddo
      enddo

! Transfer aerosol optical properties to RRTM variables;
! modify to reverse layer indexing here if necessary.

      if (iaer .ge. 1) then 
         do l = 1, nlayers
            do ib = 1, nbndsw
               taua(l,ib) = tauaer(iplon,l,ib)
               ssaa(l,ib) = ssaaer(iplon,l,ib)
               asma(l,ib) = asmaer(iplon,l,ib)
            enddo
         enddo
      endif

! Transfer cloud fraction and cloud optical properties to RRTM variables;
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then 
         inflag = inflgsw
         iceflag = iceflgsw
         liqflag = liqflgsw

! Move incoming GCM cloud arrays to RRTMG cloud arrays.
! For GCM input, incoming reice is in effective radius; for Fu parameterization (iceflag = 3)
! convert effective radius to generalized effective size using method of Mitchell, JAS, 2002:

         do l = 1, nlayers
            cldfrac(l) = cldfr(iplon,l)
            ciwp(l) = cicewp(iplon,l)
            clwp(l) = cliqwp(iplon,l)
            rei(l) = reice(iplon,l)
            if (iceflag .eq. 3) then
               dge(l) = 1.5396_jprb * reice(iplon,l)
            endif
            rel(l) = reliq(iplon,l)
            do n = 1,nbndsw
               tauc(n,l) = taucld(n,iplon,l)
               ssac(n,l) = ssacld(n,iplon,l)
               asmc(n,l) = asmcld(n,iplon,l)
            enddo
         enddo

! If an extra layer is being used in RRTMG, set all cloud properties to zero in the extra layer.

!         cldfrac(nlayers) = 0.0_jprb
!         tauc(:,nlayers) = 0.0_jprb
!         ssac(:,nlayers) = 1.0_jprb
!         asmc(:,nlayers) = 0.0_jprb
!         ciwp(nlayers) = 0.0_jprb
!         clwp(nlayers) = 0.0_jprb
!         rei(nlayers) = 0.0_jprb
!         dge(nlayers) = 0.0_jprb
!         rel(nlayers) = 0.0_jprb
      
      endif

      end subroutine inatm_sw

      end module rrtmg_sw_rad


