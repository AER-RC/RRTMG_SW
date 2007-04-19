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

      subroutine rrtmg_sw &
            (play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr , &
             o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  , &
             asdir   ,asdif   ,aldir   ,aldif   , &
             coszen  ,adjes   ,dyofyr  , &
             inflgsw ,iceflgsw,liqflgsw,cldfmcl , &
             taucmcl ,ssacmcl ,asmcmcl , &
             ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc)


! ------- Description -------

! This program is the driver for RRTMG_SW, the AER SW radiation model for 
!  application to GCMs, that has been adapted from RRTM_SW for improved
!  efficiency and to provide fractional cloudiness and cloud overlap
!  capability using McICA.
!
! Note: The call to RRTMG_SW_INIT should be moved to the GCM initialization 
!  area, since this has to be called only once. 
!
! This routine
!    b) calls INATM_SW to read in the atmospheric profile;
!       all layering in RRTMG is ordered from surface to toa. 
!    c) calls CLDPRMC_SW to set cloud optical depth for McICA based
!       on input cloud properties
!    d) calls SETCOEF_SW to calculate various quantities needed for 
!       the radiative transfer algorithm
!    e) calls SPCVMC to call the two-stream model that in turn 
!       calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands and to perform the radiative transfer
!       using McICA, the Monte-Carlo Independent Column Approximation,
!       to represent sub-grid scale cloud variability
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
!       This method is valid for clear sky or partial cloud conditions.
!
! This call to RRTMG_SW must be preceeded by a call to the module
!     mcica_subcol_gen_sw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngpt) dimension.
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
!-- Revised to add McICA capability.
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Reformatted for consistency with rrtmg_lw.
!     Feb 2007: M. J. Iacono, AER, Inc.

! --------- Modules ---------

      use parkind, only : jpim, jprb
      use parrrsw, only : nlon, mxlay, naer, nbndsw, nstr, nmol, ngpt, mxmol, &
                          jpband, jpb1, jpb2
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi
      use rrsw_wvn, only : wavenum1, wavenum2
      use rrsw_vsn

      implicit none

! ------- Declarations

! ----- Input -----
! Add necessary GCM include and parameter statements here (e.g. longitude, layer dimensions)

! Input arrays from GCM
      real(kind=jprb), intent(in) :: play(nlon,mxlay)           ! Layer pressures (hPa, mb)
      real(kind=jprb), intent(in) :: plev(nlon,mxlay+1)         ! Interface pressures (hPa, mb)
      real(kind=jprb), intent(in) :: tlay(nlon,mxlay)           ! Layer temperatures (K)
      real(kind=jprb), intent(in) :: tlev(nlon,mxlay+1)         ! Interface temperatures (K)
      real(kind=jprb), intent(in) :: tsfc(nlon)                 ! Surface temperature (K)
      real(kind=jprb), intent(in) :: h2ovmr(nlon,mxlay)         ! H2O volume mixing ratio
      real(kind=jprb), intent(in) :: o3vmr(nlon,mxlay)          ! O3 volume mixing ratio
      real(kind=jprb), intent(in) :: co2vmr                     ! CO2 volume mixing ratio
      real(kind=jprb), intent(in) :: ch4vmr(nlon,mxlay)         ! Methane volume mixing ratio
      real(kind=jprb), intent(in) :: n2ovmr(nlon,mxlay)         ! Nitrous oxide volume mixing ratio
      real(kind=jprb), intent(in) :: asdir(nlon)                ! UV/vis surface albedo direct rad
      real(kind=jprb), intent(in) :: aldir(nlon)                ! Near-IR surface albedo direct rad
      real(kind=jprb), intent(in) :: asdif(nlon)                ! UV/vis surface albedo: diffuse rad
      real(kind=jprb), intent(in) :: aldif(nlon)                ! Near-IR surface albedo: diffuse rad

      integer(kind=jpim), intent(in) :: dyofyr                  ! Day of the year (used to get Earth/Sun
                                                                !  distance if adjflx not provided)
      real(kind=jprb), intent(in) :: adjes                      ! Flux adjustment for Earth/Sun distance
      real(kind=jprb), intent(in) :: coszen                     ! Cosine of solar zenith angle

      integer(kind=jpim), intent(in) :: inflgsw                 ! Flag for cloud optical properties
      integer(kind=jpim), intent(in) :: iceflgsw                ! Flag for ice particle specification
      integer(kind=jpim), intent(in) :: liqflgsw                ! Flag for liquid droplet specification

      real(kind=jprb), intent(in) :: cldfmcl(ngpt,nlon,mxlay)   ! Cloud fraction
      real(kind=jprb), intent(in) :: taucmcl(ngpt,nlon,mxlay)   ! Cloud optical depth
      real(kind=jprb), intent(in) :: ssacmcl(ngpt,nlon,mxlay)   ! Cloud single scattering albedo
      real(kind=jprb), intent(in) :: asmcmcl(ngpt,nlon,mxlay)   ! Cloud asymmetry parameter
      real(kind=jprb), intent(in) :: ciwpmcl(ngpt,nlon,mxlay)   ! Cloud ice water path (g/m2)
      real(kind=jprb), intent(in) :: clwpmcl(ngpt,nlon,mxlay)   ! Cloud liquid water path (g/m2)
      real(kind=jprb), intent(in) :: reicmcl(nlon,mxlay)        ! Cloud ice effective radius (microns)
      real(kind=jprb), intent(in) :: relqmcl(nlon,mxlay)        ! Cloud water drop effective radius (microns)
      real(kind=jprb), intent(in) :: tauaer(nlon,mxlay,nbndsw)  ! Aerosol optical depth (iaer=10 only)
                                                                ! (non-delta scaled)      
      real(kind=jprb), intent(in) :: ssaaer(nlon,mxlay,nbndsw)  ! Aerosol single scattering albedo (iaer=10 only)
                                                                ! (non-delta scaled)      
      real(kind=jprb), intent(in) :: asmaer(nlon,mxlay,nbndsw)  ! Aerosol asymmetry parameter (iaer=10 only)
                                                                ! (non-delta scaled)      
      real(kind=jprb), intent(in) :: ecaer(nlon,mxlay,naer)     ! Aerosol optical depth at 0.55 micron (iaer=6 only)
                                                                ! (non-delta scaled)      

! ----- Output -----

      real(kind=jprb), intent(out) :: swuflx(nlon,mxlay+1)      ! Total sky shortwave upward flux (W/m2)
      real(kind=jprb), intent(out) :: swdflx(nlon,mxlay+1)      ! Total sky shortwave downward flux (W/m2)
      real(kind=jprb), intent(out) :: swhr(nlon,mxlay)          ! Total sky shortwave radiative heating rate (K/d)
      real(kind=jprb), intent(out) :: swuflxc(nlon,mxlay+1)     ! Clear sky shortwave upward flux (W/m2)
      real(kind=jprb), intent(out) :: swdflxc(nlon,mxlay+1)     ! Clear sky shortwave downward flux (W/m2)
      real(kind=jprb), intent(out) :: swhrc(nlon,mxlay)         ! Clear sky shortwave radiative heating rate (K/d)

! ----- Local -----

! Control
      integer(kind=jpim) :: nlayers             ! total number of layers
      integer(kind=jpim) :: istart              ! beginning band of calculation
      integer(kind=jpim) :: iend                ! ending band of calculation
      integer(kind=jpim) :: icld                ! clear/cloud flag
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
      real(kind=jprb) :: pavel(mxlay)           ! layer pressures (mb) 
      real(kind=jprb) :: tavel(mxlay)           ! layer temperatures (K)
      real(kind=jprb) :: pz(0:mxlay)            ! level (interface) pressures (hPa, mb)
      real(kind=jprb) :: tz(0:mxlay)            ! level (interface) temperatures (K)
      real(kind=jprb) :: tbound                 ! surface temperature (K)
      real(kind=jprb) :: pdp(mxlay)             ! layer pressure thickness (hPa, mb)
      real(kind=jprb) :: coldry(mxlay)          ! dry air column amount
      real(kind=jprb) :: wkl(mxmol,mxlay)       ! molecular amounts (mol/cm-2)

      real(kind=jprb) :: earth_sun              ! function for Earth/Sun distance factor
      real(kind=jprb) :: cossza                 ! Cosine of solar zenith angle
      real(kind=jprb) :: adjflux(jpband)        ! adjustment for current Earth/Sun distance
      real(kind=jprb) :: solvar(jpband)         ! solar constant scaling factor from rrtmg_sw
                                                !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=jprb) :: albdir(nbndsw)         ! surface albedo, direct          ! zalbp
      real(kind=jprb) :: albdif(nbndsw)         ! surface albedo, diffuse         ! zalbd

      real(kind=jprb) :: taua(mxlay,nbndsw)     ! Aerosol optical depth
      real(kind=jprb) :: ssaa(mxlay,nbndsw)     ! Aerosol single scattering albedo
      real(kind=jprb) :: asma(mxlay,nbndsw)     ! Aerosol asymmetry parameter

! Atmosphere - setcoef
      integer(kind=jpim) :: laytrop            ! tropopause layer index
      integer(kind=jpim) :: layswtch           ! tropopause layer index
      integer(kind=jpim) :: laylow             ! tropopause layer index
      integer(kind=jpim) :: jp(mxlay)          ! 
      integer(kind=jpim) :: jt(mxlay)          !
      integer(kind=jpim) :: jt1(mxlay)         !

      real(kind=jprb) :: colh2o(mxlay)         ! column amount (h2o)
      real(kind=jprb) :: colco2(mxlay)         ! column amount (co2)
      real(kind=jprb) :: colo3(mxlay)          ! column amount (o3)
      real(kind=jprb) :: coln2o(mxlay)         ! column amount (n2o)
      real(kind=jprb) :: colch4(mxlay)         ! column amount (ch4)
      real(kind=jprb) :: colo2(mxlay)          ! column amount (o2)
      real(kind=jprb) :: colmol(mxlay)         ! column amount
      real(kind=jprb) :: co2mult(mxlay)        ! column amount 

      integer(kind=jpim) :: indself(mxlay)
      integer(kind=jpim) :: indfor(mxlay)
      real(kind=jprb) :: selffac(mxlay)
      real(kind=jprb) :: selffrac(mxlay)
      real(kind=jprb) :: forfac(mxlay)
      real(kind=jprb) :: forfrac(mxlay)

      real(kind=jprb) :: &                     !
                         fac00(mxlay), fac01(mxlay), &
                         fac10(mxlay), fac11(mxlay) 

! Atmosphere/clouds - cldprop
      integer(kind=jpim) :: ncbands             ! number of cloud spectral bands
      integer(kind=jpim) :: inflag              ! flag for cloud property method
      integer(kind=jpim) :: iceflag             ! flag for ice cloud properties
      integer(kind=jpim) :: liqflag             ! flag for liquid cloud properties

      real(kind=jprb) :: cldfrac(mxlay)           ! layer cloud fraction
      real(kind=jprb) :: tauc(mxlay)              ! cloud optical depth (non-delta scaled)
      real(kind=jprb) :: ssac(mxlay)              ! cloud single scattering albedo (non-delta scaled)
      real(kind=jprb) :: asmc(mxlay)              ! cloud asymmetry parameter (non-delta scaled)
      real(kind=jprb) :: ciwp(mxlay)              ! cloud ice water path
      real(kind=jprb) :: clwp(mxlay)              ! cloud liquid water path
      real(kind=jprb) :: rei(mxlay)               ! cloud ice particle size
      real(kind=jprb) :: rel(mxlay)               ! cloud liquid particle size

      real(kind=jprb) :: taucloud(mxlay,jpband)   ! cloud optical depth
      real(kind=jprb) :: taucldorig(mxlay,jpband) ! cloud optical depth (non-delta scaled)
      real(kind=jprb) :: ssacloud(mxlay,jpband)   ! cloud single scattering albedo
      real(kind=jprb) :: asmcloud(mxlay,jpband)   ! cloud asymmetry parameter

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=jprb) :: cldfmc(ngpt,mxlay)     ! cloud fraction [mcica]
      real(kind=jprb) :: ciwpmc(ngpt,mxlay)     ! cloud ice water path [mcica]
      real(kind=jprb) :: clwpmc(ngpt,mxlay)     ! cloud liquid water path [mcica]
      real(kind=jprb) :: relqmc(mxlay)          ! liquid particle size (microns)
      real(kind=jprb) :: reicmc(mxlay)          ! ice partcle size (microns)
      real(kind=jprb) :: taucmc(ngpt,mxlay)     ! cloud optical depth [mcica]
      real(kind=jprb) :: taormc(ngpt,mxlay)     ! unscaled cloud optical depth [mcica]
      real(kind=jprb) :: ssacmc(ngpt,mxlay)     ! cloud single scattering albedo [mcica]
      real(kind=jprb) :: asmcmc(ngpt,mxlay)     ! cloud asymmetry parameter [mcica]

! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      real(kind=jprb) :: ztauc(mxlay,nbndsw)    ! cloud optical depth
      real(kind=jprb) :: ztaucorig(mxlay,nbndsw)! unscaled cloud optical depth
      real(kind=jprb) :: zasyc(mxlay,nbndsw)    ! cloud asymmetry parameter 
                                                !  (first moment of phase function)
      real(kind=jprb) :: zomgc(mxlay,nbndsw)    ! cloud single scattering albedo
      real(kind=jprb) :: ztaua(mxlay,nbndsw)    ! total aerosol optical depth
      real(kind=jprb) :: zasya(mxlay,nbndsw)    ! total aerosol asymmetry parameter 
      real(kind=jprb) :: zomga(mxlay,nbndsw)    ! total aerosol single scattering albedo

      real(kind=jprb) :: zcldfmc(mxlay,ngpt)    ! cloud fraction [mcica]
      real(kind=jprb) :: ztaucmc(mxlay,ngpt)    ! cloud optical depth [mcica]
      real(kind=jprb) :: ztaormc(mxlay,ngpt)    ! unscaled cloud optical depth [mcica]
      real(kind=jprb) :: zasycmc(mxlay,ngpt)    ! cloud asymmetry parameter [mcica] 
      real(kind=jprb) :: zomgcmc(mxlay,ngpt)    ! cloud single scattering albedo [mcica]

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

! Optional output fields 
      real(kind=jprb) :: swnflx(nlon,mxlay+1)        ! Total sky shortwave net flux (W/m2)
      real(kind=jprb) :: swnflxc(nlon,mxlay+1)       ! Clear sky shortwave net flux (W/m2)
      real(kind=jprb) :: dirdflux(nlon,mxlay+1)      ! Direct downward shortwave surface flux
      real(kind=jprb) :: difdflux(nlon,mxlay+1)      ! Diffuse downward shortwave surface flux
      real(kind=jprb) :: uvdflx(nlon,mxlay+1)        ! Total sky downward shortwave flux, UV/vis     ! pfdnuv
      real(kind=jprb) :: nidflx(nlon,mxlay+1)        ! Total sky downward shortwave flux, near-IR    ! pfdnir 
      real(kind=jprb) :: dirdnuv(nlon,mxlay+1)       ! Direct downward shortwave flux, UV/vis
      real(kind=jprb) :: difdnuv(nlon,mxlay+1)       ! Diffuse downward shortwave flux, UV/vis
      real(kind=jprb) :: dirdnir(nlon,mxlay+1)       ! Direct downward shortwave flux, near-IR
      real(kind=jprb) :: difdnir(nlon,mxlay+1)       ! Diffuse downward shortwave flux, near-IR

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


! Initializations

      zepsec = 1.e-06_jprb
      zepzen = 1.e-10_jprb
      oneminus = 1.0_jprb - zepsec
      pi = 2._jprb * sin(1._jprb)

      istart = jpb1
      iend = jpb2
      icpr = 0
      ims = 2

! In a GCM with or without McICA, set nlon to the longitude dimension
!
! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation (clear and overcast only)
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability
!            (clear, overcast or partial cloud conditions)

! *** This version uses McICA (imca = 1) ***

! Set icld to select of clear or cloud calculation and cloud 
! overlap method (read by subroutine readprof from input file INPUT_RRTM):  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap (McICA only)
! icld = 2, with clouds using maximum/random cloud overlap (McICA only)
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      icld = 2

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

      call rrtmg_sw_init

! This is the main longitude/column loop in RRTMG.
! Modify to loop over all columns (nlon) or over daylight columns

      do iplon = 1, nlon

! Prepare atmosphere profile from GCM for use in RRTMG, and define
! other input parameters

         call inatm_sw (iplon, icld, iaer, play, plev, tlay, tlev, tsfc, &
              h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, adjes, dyofyr, &
              inflgsw, iceflgsw, liqflgsw, &
              cldfmcl, taucmcl, ssacmcl, asmcmcl, ciwpmcl, clwpmcl, &
              reicmcl, relqmcl, tauaer, ssaaer, asmaer, &
              nlayers, pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
              adjflux, solvar, inflag, iceflag, liqflag, cldfmc, taucmc, &
              ssacmc, asmcmc, ciwpmc, clwpmc, reicmc, relqmc, &
              taua, ssaa, asma)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed in cldprop.  Cloud fraction and cloud
!  optical properties are transferred to rrtmg_sw arrays in cldprop.  

         call cldprmc_sw(nlayers, inflag, iceflag, liqflag, cldfmc, &
                         ciwpmc, clwpmc, reicmc, relqmc, &
                         taormc, taucmc, ssacmc, asmcmc)
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

         cossza = coszen
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

            zcldfmc(:mxlay,:ngpt) = 0._jprb
            ztaucmc(:mxlay,:ngpt) = 0._jprb
            ztaormc(:mxlay,:ngpt) = 0._jprb
            zasycmc(:mxlay,:ngpt) = 0._jprb
            zomgcmc(:mxlay,:ngpt) = 0._jprb

         elseif (icld.ge.1) then
            do i=1,nlayers
               do ig=1,ngpt
                  zcldfmc(i,ig) = cldfmc(ig,i)
                  ztaucmc(i,ig) = taucmc(ig,i)
                  ztaormc(i,ig) = taormc(ig,i)
                  zasycmc(i,ig) = asmcmc(ig,i)
                  zomgcmc(i,ig) = ssacmc(ig,i)
               enddo
            enddo

         endif   

! Aerosol
! IAER = 0: no aerosols
         if (iaer.eq.0) then

            ztaua(:mxlay,:nbndsw) = 0._jprb
            zasya(:mxlay,:nbndsw) = 0._jprb
            zomga(:mxlay,:nbndsw) = 0._jprb

! IAER = 6: Use ECMWF six aerosol types. See rrsw_aer.f90 for details.
! Input aerosol optical thickness at 0.55 micron for each aerosol type (ecaer), 
! or set manually here for each aerosol and layer.
         elseif (iaer.eq.6) then

!            do i = 1, nlayers
!               do ia = 1, naer
!                  ecaer(iplon,i,ia) = 1.0e-15_jprb
!               enddo
!            enddo

            do i = 1, nlayers
               do ib = 1, nbndsw
                  ztaua(i,ib) = 0._jprb
                  zasya(i,ib) = 0._jprb
                  zomga(i,ib) = 0._jprb
                  do ia = 1, naer
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

! Transfer up and down, clear and total sky fluxes to output arrays.
! Vertical indexing goes from bottom to top; reverse here for GCM if necessary.

         do i = 1, nlayers+1
            swuflxc(iplon,i) = zbbcu(i)
            swdflxc(iplon,i) = zbbcd(i)
            swuflx(iplon,i) = zbbfu(i)
            swdflx(iplon,i) = zbbfd(i)
            uvdflx(iplon,i) = zuvfd(i)
            nidflx(iplon,i) = znifd(i)
!  Direct/diffuse fluxes
            dirdflux(iplon,i) = zbbfddir(i)
            difdflux(iplon,i) = swdflx(iplon,i) - dirdflux(iplon,i)
!  UV/visible direct/diffuse fluxes
            dirdnuv(iplon,i) = zuvfddir(i)
            difdnuv(iplon,i) = zuvfd(i) - dirdnuv(iplon,i)
!  Near-IR direct/diffuse fluxes
            dirdnir(iplon,i) = znifddir(i)
            difdnir(iplon,i) = znifd(i) - dirdnir(iplon,i)
         enddo

!  Total and clear sky net fluxes
         do i = 1, nlayers+1
            swnflxc(iplon,i) = swdflxc(iplon,i) - swuflxc(iplon,i)
            swnflx(iplon,i) = swdflx(iplon,i) - swuflx(iplon,i)
         enddo

!  Total and clear sky heating rates
         do i = 1, nlayers
            zdpgcp = heatfac / pdp(i)
            swhrc(iplon,i) = (swnflxc(iplon,i+1) - swnflxc(iplon,i)) * zdpgcp
            swhr(iplon,i) = (swnflx(iplon,i+1) - swnflx(iplon,i)) * zdpgcp
         enddo
         swhrc(iplon,nlayers) = 0._jprb
         swhr(iplon,nlayers) = 0._jprb

! End longitude loop
      enddo

      return
      end subroutine rrtmg_sw

!*************************************************************************
      real function earth_sun(idn)
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

      return
      end function earth_sun

!***************************************************************************
      subroutine swdatinit
!***************************************************************************

! --------- Modules ----------

      use parkind, only : jpim, jprb 
      use rrsw_con, only: heatfac, grav, planck, boltz, &
                          clight, avogad, alosmt, gascon, radcn1, radcn2 
      use rrsw_wvn, only: ng, nspa, nspb, wavenum1, wavenum2, delwave
      use rrsw_vsn

      implicit none
      save 
 
! Shortwave spectral band limits (wavenumbers)
      wavenum1(:) = (/2600._jprb, 3250._jprb, 4000._jprb, 4650._jprb, 5150._jprb, 6150._jprb, 7700._jprb, &
                      8050._jprb,12850._jprb,16000._jprb,22650._jprb,29000._jprb,38000._jprb,  820._jprb/)
      wavenum2(:) = (/3250._jprb, 4000._jprb, 4650._jprb, 5150._jprb, 6150._jprb, 7700._jprb, 8050._jprb, &
                     12850._jprb,16000._jprb,22650._jprb,29000._jprb,38000._jprb,50000._jprb, 2600._jprb/)
      delwave(:) =  (/ 650._jprb,  750._jprb,  650._jprb,  500._jprb, 1000._jprb, 1550._jprb,  350._jprb, &
                      4800._jprb, 3150._jprb, 6650._jprb, 6350._jprb, 9000._jprb,12000._jprb, 1780._jprb/)

! Spectral band information
      ng(:) = (/16,16,16,16,16,16,16,16,16,16,16,16,16,16/)
      nspa(:) = (/9,9,9,9,1,9,9,1,9,1,0,1,9,1/)
      nspb(:) = (/1,5,1,1,1,5,1,0,1,0,0,1,5,1/)

!     Heatfac is the factor by which one must multiply delta-flux/ 
!     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
!     the heating rate in units of degrees/day.  It is equal to 
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(86400)(1e-5)/(1.004)
      heatfac = 8.4391_jprb

!     Modified values for consistency with CAM3:
!        =  (9.80616)(86400)(1e-5)/(1.00464)
!      heatfac = 8.43339130434_jprb

!    Constants from NIST 01/11/2002

      grav = 9.8066_jprb
      planck = 6.62606876e-27_jprb
      boltz = 1.3806503e-16_jprb
      clight = 2.99792458e+10_jprb
      avogad = 6.02214199e+23_jprb
      alosmt = 2.6867775e+19_jprb
      gascon = 8.31447200e+07_jprb
      radcn1 = 1.191042722e-12_jprb
      radcn2 = 1.4387752_jprb
!
!     units are generally cgs
!
!     The first and second radiation constants are taken from NIST.
!     They were previously obtained from the relations:
!          radcn1 = 2.*planck*clight*clight*1.e-07
!          radcn2 = planck*clight/boltz

      return
      end

!***************************************************************************
      subroutine inatm_sw (iplon, icld, iaer, play, plev, tlay, tlev, tsfc, &
            h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, adjes, dyofyr, &
            inflgsw, iceflgsw, liqflgsw, &
            cldfmcl, taucmcl, ssacmcl, asmcmcl, ciwpmcl, clwpmcl, &
            reicmcl, relqmcl, tauaer, ssaaer, asmaer, &
            nlayers, pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
            adjflux, solvar, inflag, iceflag, liqflag, cldfmc, taucmc, &
            ssacmc, asmcmc, ciwpmc, clwpmc, reicmc, relqmc, &
            taua, ssaa, asma)
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_SW.
!  Set other RRTMG_SW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parkind, only : jpim, jprb 
      use parrrsw, only : nlon, mxlay, nbndsw, nstr, nmol, ngpt, mxmol, &
                          jpband, jpb1, jpb2
      use rrsw_con, only : heatfac, oneminus, pi, grav, avogad
      use rrsw_wvn, only : ng, nspa, nspb, wavenum1, wavenum2, delwave

      implicit none

! ------- Declarations -------

! ----- Input -----
! Add necessary GCM include and parameter statements here (e.g. longitude, layer dimensions)

! Input arrays from GCM
      integer(kind=jpim), intent(in) :: iplon                   ! column loop index
      integer(kind=jpim), intent(in) :: icld                    ! clear/cloud flag
      integer(kind=jpim), intent(in) :: iaer                    ! aerosol option flag

      real(kind=jprb), intent(in) :: play(nlon,mxlay)           ! Layer pressures (hPa, mb)
      real(kind=jprb), intent(in) :: plev(nlon,mxlay+1)         ! Interface pressures (hPa, mb)
      real(kind=jprb), intent(in) :: tlay(nlon,mxlay)           ! Layer temperatures (K)
      real(kind=jprb), intent(in) :: tlev(nlon,mxlay+1)         ! Interface temperatures (K)
      real(kind=jprb), intent(in) :: tsfc(nlon)                 ! Surface temperature (K)
      real(kind=jprb), intent(in) :: h2ovmr(nlon,mxlay)         ! H2O volume mixing ratio
      real(kind=jprb), intent(in) :: o3vmr(nlon,mxlay)          ! O3 volume mixing ratio
      real(kind=jprb), intent(in) :: co2vmr                     ! CO2 volume mixing ratio
      real(kind=jprb), intent(in) :: ch4vmr(nlon,mxlay)         ! Methane volume mixing ratio
      real(kind=jprb), intent(in) :: n2ovmr(nlon,mxlay)         ! Nitrous oxide volume mixing ratio

      integer(kind=jpim), intent(in) :: dyofyr                  ! Day of the year (used to get Earth/Sun
                                                                !  distance if adjflx not provided)
      real(kind=jprb), intent(in) :: adjes                      ! Flux adjustment for Earth/Sun distance

      integer(kind=jpim), intent(in) :: inflgsw                 ! Flag for cloud optical properties
      integer(kind=jpim), intent(in) :: iceflgsw                ! Flag for ice particle specification
      integer(kind=jpim), intent(in) :: liqflgsw                ! Flag for liquid droplet specification

      real(kind=jprb), intent(in) :: cldfmcl(ngpt,nlon,mxlay)   ! Cloud fraction
      real(kind=jprb), intent(in) :: taucmcl(ngpt,nlon,mxlay)   ! Cloud optical depth (optional)
      real(kind=jprb), intent(in) :: ssacmcl(ngpt,nlon,mxlay)   ! Cloud single scattering albedo
      real(kind=jprb), intent(in) :: asmcmcl(ngpt,nlon,mxlay)   ! Cloud asymmetry parameter
      real(kind=jprb), intent(in) :: ciwpmcl(ngpt,nlon,mxlay)   ! Cloud ice water path (g/m2)
      real(kind=jprb), intent(in) :: clwpmcl(ngpt,nlon,mxlay)   ! Cloud liquid water path (g/m2)
      real(kind=jprb), intent(in) :: reicmcl(nlon,mxlay)        ! Cloud ice effective radius (microns)
      real(kind=jprb), intent(in) :: relqmcl(nlon,mxlay)        ! Cloud water drop effective radius (microns)

      real(kind=jprb), intent(in) :: tauaer(nlon,mxlay,nbndsw)  ! Aerosol optical depth
      real(kind=jprb), intent(in) :: ssaaer(nlon,mxlay,nbndsw)  ! Aerosol single scattering albedo
      real(kind=jprb), intent(in) :: asmaer(nlon,mxlay,nbndsw)  ! Aerosol asymmetry parameter

! Atmosphere
      integer(kind=jpim), intent(out) :: nlayers             ! number of layers

      real(kind=jprb), intent(out) :: pavel(mxlay)           ! layer pressures (mb) 
      real(kind=jprb), intent(out) :: tavel(mxlay)           ! layer temperatures (K)
      real(kind=jprb), intent(out) :: pz(0:mxlay)            ! level (interface) pressures (hPa, mb)
      real(kind=jprb), intent(out) :: tz(0:mxlay)            ! level (interface) temperatures (K)
      real(kind=jprb), intent(out) :: tbound                 ! surface temperature (K)
      real(kind=jprb), intent(out) :: pdp(mxlay)             ! layer pressure thickness (hPa, mb)
      real(kind=jprb), intent(out) :: coldry(mxlay)          ! 
      real(kind=jprb), intent(out) :: wkl(mxmol,mxlay)       ! molecular amounts (mol/cm-2)

      real(kind=jprb), intent(out) :: adjflux(jpband)        ! adjustment for current Earth/Sun distance
      real(kind=jprb), intent(out) :: solvar(jpband)         ! solar constant scaling factor from rrtmg_sw
                                                             !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=jprb), intent(out) :: taua(mxlay,nbndsw)     ! Aerosol optical depth
      real(kind=jprb), intent(out) :: ssaa(mxlay,nbndsw)     ! Aerosol single scattering albedo
      real(kind=jprb), intent(out) :: asma(mxlay,nbndsw)     ! Aerosol asymmetry parameter

! Atmosphere/clouds - cldprop
      integer(kind=jpim), intent(out) :: inflag              ! flag for cloud property method
      integer(kind=jpim), intent(out) :: iceflag             ! flag for ice cloud properties
      integer(kind=jpim), intent(out) :: liqflag             ! flag for liquid cloud properties

      real(kind=jprb), intent(out) :: cldfmc(ngpt,mxlay)     ! layer cloud fraction
      real(kind=jprb), intent(out) :: taucmc(ngpt,mxlay)     ! cloud optical depth (non-delta scaled)
      real(kind=jprb), intent(out) :: ssacmc(ngpt,mxlay)     ! cloud single scattering albedo (non-delta-scaled)
      real(kind=jprb), intent(out) :: asmcmc(ngpt,mxlay)     ! cloud asymmetry parameter (non-delta scaled)
      real(kind=jprb), intent(out) :: ciwpmc(ngpt,mxlay)     ! cloud ice water path
      real(kind=jprb), intent(out) :: clwpmc(ngpt,mxlay)     ! cloud liquid water path
      real(kind=jprb), intent(out) :: reicmc(mxlay)          ! cloud ice particle size
      real(kind=jprb), intent(out) :: relqmc(mxlay)          ! cloud liquid particle size


! ----- Local -----
      real(kind=jprb), parameter :: amd = 28.9660_jprb       ! Effective molecular weight of dry air (g/mol)
      real(kind=jprb), parameter :: amw = 18.0160_jprb       ! Molecular weight of water vapor (g/mol)
!      real(kind=jprb), parameter :: amc = 44.0098_jprb       ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=jprb), parameter :: amo = 47.9998_jprb       ! Molecular weight of ozone (g/mol)
!      real(kind=jprb), parameter :: amo2 = 31.9999_jprb      ! Molecular weight of oxygen (g/mol)
!      real(kind=jprb), parameter :: amch4 = 16.0430_jprb     ! Molecular weight of methane (g/mol)
!      real(kind=jprb), parameter :: amn2o = 44.0128_jprb     ! Molecular weight of nitrous oxide (g/mol)
!      real(kind=jprb), parameter :: amc11 = 137.3684_jprb    ! Molecular weight of CFC11 (g/mol) - CFCL3
!      real(kind=jprb), parameter :: amc12 = 120.9138_jprb    ! Molecular weight of CFC12 (g/mol) - CF2CL2
! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=jprb), parameter :: amdw = 1.607793_jprb      ! Molecular weight of dry air / water vapor
      real(kind=jprb), parameter :: amdc = 0.658114_jprb      ! Molecular weight of dry air / carbon dioxide
      real(kind=jprb), parameter :: amdo = 0.603428_jprb      ! Molecular weight of dry air / ozone
      real(kind=jprb), parameter :: amdm = 1.805423_jprb      ! Molecular weight of dry air / methane
      real(kind=jprb), parameter :: amdn = 0.658090_jprb      ! Molecular weight of dry air / nitrous oxide
      real(kind=jprb), parameter :: amdo2 = 0.905140_jprb     ! Molecular weight of dry air / oxygen
      real(kind=jprb), parameter :: amdc1 = 0.210852_jprb     ! Molecular weight of dry air / CFC11
      real(kind=jprb), parameter :: amdc2 = 0.239546_jprb     ! Molecular weight of dry air / CFC12

      real(kind=jprb), parameter :: sbc = 5.67e-08_jprb       ! Stefan-Boltzmann constant (W/m2K4)
      real(kind=jprb), parameter :: o2mmr = 0.23143_jprb      ! o2 mass mixing ratio

      integer(kind=jpim) :: isp, l, ix, n, imol, ib, ig       ! Loop indices
      real(kind=jprb) :: amm, summol                          ! 
      real(kind=jprb) :: adjflx                               ! flux adjustment for Earth/Sun distance
      real(kind=jprb) :: earth_sun                            ! function for Earth/Sun distance adjustment

      nlayers = mxlay

!  Initialize all molecular amounts to zero here, then pass input amounts
!  into RRTM array WKL below.

       wkl(:,:) = 0.0_jprb
       cldfmc(:,:) = 0.0_jprb
       taucmc(:,:) = 0.0_jprb
       ssacmc(:,:) = 0.0_jprb
       asmcmc(:,:) = 0.0_jprb
       ciwpmc(:,:) = 0.0_jprb
       clwpmc(:,:) = 0.0_jprb
       reicmc(:) = 0.0_jprb
       relqmc(:) = 0.0_jprb
       taua(:,:) = 0.0_jprb
       ssaa(:,:) = 0.0_jprb
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
! solar constant (1368.22 Wm-2) by band (SOLVAR).  SOLVAR can be set 
! to a single scaling factor as needed, or to a different value in each 
! band, which may be necessary for paleoclimate simulations. 
! 
      do ib = jpb1,jpb2
         solvar(ib) = 1._jprb
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
! For h2o input in vmr:
         wkl(1,l) = h2ovmr(iplon,l)
! For h2o input in mmr:
!         wkl(1,l) = h2o(iplon,l)*amdw
! For h2o input in specific humidity;
!         wkl(1,l) = (h2o(iplon,l)/(1._jprb - h2o(iplon,l)))*amdw
         wkl(2,l) = co2vmr
         wkl(3,l) = o3vmr(iplon,l)
         wkl(4,l) = n2ovmr(iplon,l)
         wkl(6,l) = ch4vmr(iplon,l)
         wkl(7,l) = o2mmr*amdo2
         amm = (1._jprb - wkl(1,l)) * amd + wkl(1,l) * amw            
         coldry(l) = (pz(l-1)-pz(l)) * 1.e3_jprb * avogad / &
                     (1.e2_jprb * grav * amm * (1._jprb + wkl(1,l)))
      enddo

! The following section can be used to set values for an additional layer (from
! the GCM top level to 0. mb) for improved calculation of TOA fluxes. 
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
!      wkl(1,nlayers) = wkl(1,nlayers-1)
!      wkl(2,nlayers) = wkl(2,nlayers-1)
!      wkl(3,nlayers) = wkl(3,nlayers-1)
!      wkl(4,nlayers) = wkl(4,nlayers-1)
!      wkl(6,nlayers) = wkl(6,nlayers-1)
!      wkl(7,nlayers) = wkl(7,nlayers-1)
!      amm = (1._jprb - wkl(1,nlayers-1)) * amd + wkl(1,nlayers-1) * amw
!      coldry(nlayers) = (pz(nlayers-1)) * 1.e3_jprb * avogad / &
!                        (1.e2_jprb * grav * amm * (1._jprb + wkl(1,nlayers-1)))

! At this point all moleculular amounts in wkl are in volume mixing ratio; 
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
            do ig = 1, ngpt
               cldfmc(ig,l) = cldfmcl(ig,iplon,l)
               taucmc(ig,l) = taucmcl(ig,iplon,l)
               ssacmc(ig,l) = ssacmcl(ig,iplon,l)
               asmcmc(ig,l) = asmcmcl(ig,iplon,l)
               ciwpmc(ig,l) = ciwpmcl(ig,iplon,l)
               clwpmc(ig,l) = clwpmcl(ig,iplon,l)
            enddo
            reicmc(l) = reicmcl(iplon,l)
            if (iceflag .eq. 3) then
               reicmc(l) = 1.5396 * reicmcl(iplon,l)
            endif
            relqmc(l) = relqmcl(iplon,l)
         enddo

! If an extra layer is being used in RRTMG, set all cloud properties to zero in the extra layer.

!         cldfmc(:,nlayers) = 0.0_jprb
!         taucmc(:,nlayers) = 0.0_jprb
!         ssacmc(:,nlayers) = 0.0_jprb
!         asmcmc(:,nlayers) = 0.0_jprb
!         ciwpmc(:,nlayers) = 0.0_jprb
!         clwpmc(:,nlayers) = 0.0_jprb
!         reicmc(nlayers) = 0.0_jprb
!         relqmc(nlayers) = 0.0_jprb
      
      endif

      return
      end 

