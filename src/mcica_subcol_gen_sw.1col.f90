
      module mcica_subcol_gen_sw

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! Purpose: Create McICA stochastic arrays for cloud physical or optical properties.   
! Two options are possible:
! 1) Input cloud physical properties: cloud fraction, ice and liquid water
!    paths, ice fraction, and particle sizes.  Output will be stochastic
!    arrays of these variables.  (inflag = 1)
! 2) Input cloud optical properties directly: cloud optical depth, single
!    scattering albedo and asymmetry parameter.  Output will be stochastic
!    arrays of these variables.  (inflag = 0)

! --------- Modules ----------

      use parkind, only : jpim, jprb 
      use parrrsw, only : mxlay, nbndsw, ngpt, nlon
      use rrsw_con, only: grav
      use rrsw_wvn, only: ngb
      use rrsw_vsn

      implicit none

      private

! public interfaces/functions/subroutines
      public mcica_subcol_sw, generate_stochastic_clouds

      contains

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine mcica_subcol_sw(iplon, nlayers, icld, ims, pavel, &
                       cldfrac, ciwp, clwp, rei, rel, tauc, ssac, asmc, &
                       cldfmc, ciwpmc, clwpmc, reicmc, relqmc, taucmc, &
                       ssacmc, asmcmc)

! Control
      integer(kind=jpim), intent(in) :: nlayers             ! number of model layers
      integer(kind=jpim), intent(in) :: icld                ! clear/cloud, cloud overlap flag
      integer(kind=jpim), intent(in) :: ims                 ! mcica statistical loop index; also
                                                            ! value for changing mcica permute seed
      integer(kind=jpim), intent(in) :: iplon               ! column/longitude dimension

! Atmosphere
      real(kind=jprb), intent(in) :: pavel(mxlay)           ! layer pressures (mb) 
!      real(kind=jprb), intent(in) :: pdp(mxlay)             ! layer pressure thickness (mb) 

! Atmosphere/clouds - cldprop
      real(kind=jprb), intent(in) :: cldfrac(mxlay)         ! layer cloud fraction
      real(kind=jprb), intent(in) :: tauc(nbndsw,mxlay)     ! cloud optical depth
      real(kind=jprb), intent(in) :: ssac(nbndsw,mxlay)     ! cloud single scattering albedo (non-delta scaled)
      real(kind=jprb), intent(in) :: asmc(nbndsw,mxlay)     ! cloud asymmetry parameter (non-delta scaled)
      real(kind=jprb), intent(in) :: ciwp(mxlay)            ! cloud ice water path
      real(kind=jprb), intent(in) :: clwp(mxlay)            ! cloud liquid water path
      real(kind=jprb), intent(in) :: rei(mxlay)             ! cloud ice particle size
      real(kind=jprb), intent(in) :: rel(mxlay)             ! cloud liquid particle size

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=jprb), intent(out) :: cldfmc(ngpt,mxlay)       ! cloud fraction [mcica]
      real(kind=jprb), intent(out) :: ciwpmc(ngpt,mxlay)       ! cloud ice water path [mcica]
      real(kind=jprb), intent(out) :: clwpmc(ngpt,mxlay)       ! cloud liquid water path [mcica]
      real(kind=jprb), intent(out) :: relqmc(mxlay)            ! liquid particle size (microns)
      real(kind=jprb), intent(out) :: reicmc(mxlay)            ! ice partcle size (microns)
      real(kind=jprb), intent(out) :: taucmc(ngpt,mxlay)       ! cloud optical depth [mcica]
      real(kind=jprb), intent(out) :: ssacmc(ngpt,mxlay)       ! cloud single scattering albedo [mcica]
      real(kind=jprb), intent(out) :: asmcmc(ngpt,mxlay)       ! cloud asymmetry parameter [mcica]

! Local

! Stochastic cloud generator variables [mcica]
      integer(kind=jpim), parameter :: nsubcsw = ngpt  ! number of sub-columns (g-point intervals)
      integer(kind=jpim) :: permuteseed                ! if the cloud generator is called multiple times, 
                                                       ! permute the seed between each call.    
      integer(kind=jpim) :: km, im, nm                 ! loop indices

      real(kind=jprb) :: pmid(mxlay)                          ! layer pressures (Pa) 
!      real(kind=jprb) :: pdel(mxlay)                          ! layer pressure thickness (Pa) 
!      real(kind=jprb) :: qi(mxlay)                            ! ice water (specific humidity)
!      real(kind=jprb) :: ql(mxlay)                            ! liq water (specific humidity)


! Return if clear sky; or stop if icld out of range
      if (icld.eq.0) return
      if (icld.lt.0.or.icld.gt.3) then 
         stop 'RRTMG_SW: INVALID ICLD'
      endif 

! Set permute seed for random number generator
! For single column mode, permuteseed must be different for each (ims) sample performed

      permuteseed = (ims+1)*nsubcsw

! Pass particle sizes to new arrays, no subcolumns for these properties yet
! Convert pressures from mb to Pa

      reicmc(:mxlay) = rei(:mxlay)
      relqmc(:mxlay) = rel(:mxlay)
      pmid(:mxlay) = pavel(:mxlay)*1.e2_jprb
!      pdel(:mxlay) = pdp(:mxlay)*1.e2_jprb

! Convert input ice and liquid cloud water paths to specific humidity ice and liquid components 

!      cwp =  (q * pdel * 1000.) / gravit)
!           = (kg/kg * kg m-1 s-2 *1000.) / m s-2
!           = (g m-2)
!
!      q  = (cwp * gravit) / (pdel *1000.)
!         = (g m-2 * m s-2) / (kg m-1 s-2 * 1000.)
!         =  kg/kg

!      do km = 1, nlayers
!         qi(km) = (ciwp(km) * grav) / (pdel(km) * 1000._jprb)
!         ql(km) = (clwp(km) * grav) / (pdel(km) * 1000._jprb)
!      enddo

!  Generate the stochastic subcolumns and their optical properties for the longwave;
      call generate_stochastic_clouds (nlayers, icld, pmid, cldfrac, clwp, ciwp, tauc, &
                               ssac, asmc, cldfmc, clwpmc, ciwpmc, taucmc, &
                               ssacmc, asmcmc, permuteseed)

      end subroutine mcica_subcol_sw


!-------------------------------------------------------------------------------------------------
      subroutine generate_stochastic_clouds(nlayers, icld, pmid, cld, clwp, ciwp, tauc, & 
                                   ssac, asmc, cld_stoch, clwp_stoch, ciwp_stoch, tauc_stoch, &
                                   ssac_stoch, asmc_stoch, changeSeed) 
!-------------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------------------------------------
  ! ---------------------
  ! Contact: Cecile Hannay (hannay@ucar.edu)
  ! 
  ! Modifications: Generalized for use with RRTMG and added Mersenne Twister as the default
  !   random number generator, which can be changed to the optional kissvec random number generator
  !   with flag 'irnd' below . Some extra functionality has been commented or removed.  
  !   Michael J. Iacono, AER, Inc., February 2007
  !
  ! Given a profile of cloud fraction, cloud water and cloud ice, we produce a set of subcolumns.
  ! Each layer within each subcolumn is homogeneous, with cloud fraction equal to zero or one 
  ! and uniform cloud liquid and cloud ice concentration.
  ! The ensemble as a whole reproduces the probability function of cloud liquid and ice within each layer 
  ! and obeys an overlap assumption in the vertical.   
  ! 
  ! Overlap assumption:
  !  The cloud are consistent with 4 overlap assumptions: random, maximum, maximum-random and exponential. 
  !  The default option is maximum-random (option 3)
  !  The options are: 1=random overlap, 2=max/random, 3=maximum overlap, 4=exponential overlap
  !  This is set with the variable "overlap" 
  !mji - Exponential overlap option (overlap=4) has been deactivated in this version
  !  The exponential overlap uses also a length scale, Zo. (real,    parameter  :: Zo = 2500. ) 
  ! 
  ! Seed:
  !  If the stochastic cloud generator is called several times during the same timestep, 
  !  one should change the seed between the call to insure that the subcolumns are different.
  !  This is done by changing the argument 'changeSeed'
  !  For example, if one wants to create a set of columns for the shortwave and another set for the longwave ,
  !  use 'changeSeed = 1' for the first call and'changeSeed = 2' for the second call 
  !
  ! PDF assumption:
  !  We can use arbitrary complicated PDFS. 
  !  In the present version, we produce homogeneuous clouds (the simplest case).  
  !  Future developments include using the PDF scheme of Ben Johnson. 
  !
  ! History file:
  !  Option to add diagnostics variables in the history file. (using FINCL in the namelist)
  !  nsubcol = number of subcolumns
  !  overlap = overlap type (1-3)
  !  Zo = length scale 
  !  CLOUD_S = mean of the subcolumn cloud fraction ('_S" means Stochastic)
  !  CLDLIQ_S = mean of the subcolumn cloud water
  !  CLDICE_S = mean of the subcolumn cloud ice 
  !
  ! Note:
  !   Here: we force that the cloud condensate to be consistent with the cloud fraction 
  !   i.e we only have cloud condensate when the cell is cloudy. 
  !   In CAM: The cloud condensate and the cloud fraction are obtained from 2 different equations 
  !   and the 2 quantities can be inconsistent (i.e. CAM can produce cloud fraction 
  !   without cloud condensate or the opposite).
  !---------------------------------------------------------------------------------------------------------------

      use mcica_random_numbers
! The Mersenne Twister random number engine
      use MersenneTwister, only: randomNumberSequence, &   
                                 new_RandomNumberSequence, getRandomReal

      implicit none

      type(randomNumberSequence) :: randomNumbers

! -- Arguments

      integer(kind=jpim), intent(in) :: nlayers                  ! number of layers
      integer(kind=jpim), intent(in) :: icld                     ! clear/cloud, cloud overlap flag
      integer(kind=jpim), optional, intent(in) :: changeSeed     ! allows permuting seed

! Column state (cloud fraction, cloud water, cloud ice) + variables needed to read physics state 
      real(kind=jprb), intent(in) :: pmid(mxlay)                 ! layer pressure (Pa)
      real(kind=jprb), intent(in) :: cld(mxlay)                  ! cloud fraction 
      real(kind=jprb), intent(in) :: clwp(mxlay)                 ! cloud liquid water path (g/m2)
      real(kind=jprb), intent(in) :: ciwp(mxlay)                 ! cloud ice water path (g/m2)
      real(kind=jprb), intent(in) :: tauc(nbndsw,mxlay)          ! cloud optical depth (non-delta scaled)
      real(kind=jprb), intent(in) :: ssac(nbndsw,mxlay)          ! cloud single scattering albedo (non-delta scaled)
      real(kind=jprb), intent(in) :: asmc(nbndsw,mxlay)          ! cloud asymmetry parameter (non-delta scaled)

      real(kind=jprb), intent(out) :: cld_stoch(ngpt, mxlay)  ! subcolumn cloud fraction 
      real(kind=jprb), intent(out) :: clwp_stoch(ngpt, mxlay) ! subcolumn cloud liquid water path
      real(kind=jprb), intent(out) :: ciwp_stoch(ngpt, mxlay) ! subcolumn cloud ice water path
      real(kind=jprb), intent(out) :: tauc_stoch(ngpt, mxlay) ! subcolumn cloud optical depth
      real(kind=jprb), intent(out) :: ssac_stoch(ngpt, mxlay) ! subcolumn cloud single scattering albedo
      real(kind=jprb), intent(out) :: asmc_stoch(ngpt, mxlay) ! subcolumn cloud asymmetry parameter

! -- Local variables
      integer(kind=jpim), parameter :: nsubcol = ngpt            ! number of sub-columns (g-point intervals)
      real(kind=jprb) :: cldf(mxlay)                             ! cloud fraction 
    
! Mean over the subcolumns (cloud fraction, cloud water , cloud ice) - inactive
!      real(kind=jprb) :: mean_cld_stoch(mxlay)                   ! cloud fraction 
!      real(kind=jprb) :: mean_clwp_stoch(mxlay)                  ! cloud water
!      real(kind=jprb) :: mean_ciwp_stoch(mxlay)                  ! cloud ice
!      real(kind=jprb) :: mean_tauc_stoch(mxlay)                  ! cloud optical depth
!      real(kind=jprb) :: mean_ssac_stoch(mxlay)                  ! cloud single scattering albedo
!      real(kind=jprb) :: mean_asmc_stoch(mxlay)                  ! cloud asymmetry parameter

! Set overlap
      integer(kind=jpim) :: overlap                              ! 1 = random overlap, 2 = maximum/random,
                                                                 ! 3 = maximum overlap, 
!      real(kind=jprb), parameter  :: Zo = 2500.                 ! length scale (m) 
!      real(kind=jprb) :: zm(mxlay)                              ! Height of midpoints (above surface)
!      real(kind=jprb), dimension(mxlay) :: alpha=0.0            ! overlap parameter  

! Constants (min value for cloud fraction and cloud water and ice)
      real(kind=jprb), parameter :: cldmin = 1.0e-2              ! min cloud fraction
!      real(kind=jprb), parameter :: qmin   = 1.0e-10             ! min cloud water and cloud ice (not used)

! Variables related to random number and seed 
      integer(kind=jpim) :: irnd                                 ! flag for random number generator
                                                                 !  0 = kissvec
                                                                 !  1 = Mersenne Twister

      real(kind=jprb), dimension(nsubcol, mxlay) :: CDF, CDF2    ! random numbers
      integer(kind=jpim) :: seed1, seed2, seed3, seed4           ! seed to create random number (kissvec)
      real(kind=jprb) :: rand_num                                ! random number (kissvec)
      integer(kind=jpim) :: iseed                                ! seed to create random number (Mersenne Twister)
      real(kind=jprb) :: rand_num_mt                             ! random number (Mersenne Twister)

! Flag to identify cloud fraction in subcolumns
      logical,  dimension(nsubcol, mxlay) :: iscloudy            ! flag that says whether a gridbox is cloudy

! Indices
      integer(kind=jpim) :: ilev, isubcol, i, n                  ! indices

!------------------------------------------------------------------------------------------ 
! Set randum number generator to use (0 = kissvec; 1 = mersennetwister)
!      irnd = 0
      irnd = 1

! Pass input cloud overlap setting to local variable
    overlap = icld

! ensure that cloud fractions are in bounds 
! to avoid to get ql_stoch getting to big (ql_stoch = ql/cld).
    cldf(:) = cld(:)
    where (cldf(:) < cldmin)
        cldf(:) = 0._jprb
    end where

! ----- Create seed  --------
   
! Advance randum number generator by changeseed values
      if (irnd.eq.0) then   
! For kissvec, create a seed that depends on the state of the columns. Maybe not the best way, but it works.  
! Must use pmid from bottom four layers. 
         seed1 = (pmid(1) - int(pmid(1)))  * 1000000000_jpim
         seed2 = (pmid(2) - int(pmid(2)))  * 1000000000_jpim
         seed3 = (pmid(3) - int(pmid(3)))  * 1000000000_jpim
         seed4 = (pmid(4) - int(pmid(4)))  * 1000000000_jpim
         do i=1,changeSeed
            call kissvec(seed1, seed2, seed3, seed4, rand_num)
         enddo
      elseif (irnd.eq.1) then
         randomNumbers = new_RandomNumberSequence(seed = changeSeed)
      endif 


! ------ Apply overlap assumption --------

! generate the random numbers  

      select case (overlap)

      case(1) 
! Random overlap
! i) pick a random value at every level
  
         if (irnd.eq.0) then 
            do isubcol = 1,nsubcol
               do ilev = 1,nlayers
                  call kissvec(seed1, seed2, seed3, seed4, rand_num)  ! we get different random number for each level
                  CDF(isubcol, ilev) = rand_num
               enddo
            enddo
         elseif (irnd.eq.1) then
            do isubcol = 1, nsubcol
               do ilev = 1, nlayers
                  rand_num_mt = getRandomReal(randomNumbers)
                  CDF(isubcol,ilev) = rand_num_mt
               enddo
             enddo
         endif

      case(2) 
! Maximum-Random overlap
! i) pick  a random number for top layer.
! ii) walk down the column: 
!    - if the layer above is cloudy, we use the same random number than in the layer above
!    - if the layer above is clear, we use a new random number 

         if (irnd.eq.0) then 
            do isubcol = 1,nsubcol
               do ilev = 1,nlayers
                  call kissvec(seed1, seed2, seed3, seed4, rand_num) 
                  CDF(isubcol, ilev) = rand_num
               enddo
            enddo
         elseif (irnd.eq.1) then
            do isubcol = 1, nsubcol
               do ilev = 1, nlayers
                  rand_num_mt = getRandomReal(randomNumbers)
                  CDF(isubcol,ilev) = rand_num_mt
               enddo
             enddo
         endif

         do ilev = 2,nlayers
            where (CDF(:, ilev-1) > spread(1._jprb - cldf(ilev-1), dim=1, nCopies=nsubcol) )
               CDF(:,ilev) = CDF(:,ilev-1) 
            elsewhere
               CDF(:,ilev) = CDF(:,ilev) *  spread(1._jprb - cldf(ilev-1), dim=1, nCopies=nsubcol) 
            end where
         enddo

    case(3) 
! Maximum overlap
! i) pick same random numebr at every level  

         if (irnd.eq.0) then 
            do isubcol = 1,nsubcol
               call kissvec(seed1, seed2, seed3, seed4, rand_num) 
               do ilev = 1,nlayers
                  CDF(isubcol, ilev) = rand_num
               enddo
            enddo
         elseif (irnd.eq.1) then
            do isubcol = 1, nsubcol
                  rand_num_mt = getRandomReal(randomNumbers)
               do ilev = 1, nlayers
                  CDF(isubcol,ilev) = rand_num_mt
               enddo
             enddo
         endif

!    case(4) - inactive
!       ! Exponential overlap: weighting between maximum and random overlap increases with the distance. 
!       ! The random numbers for exponential overlap verify:
!       ! j=1   RAN(j)=RND1
!       ! j>1   if RND1 < alpha(j,j-1) => RAN(j) = RAN(j-1)
!       !                                 RAN(j) = RND2
!       ! alpha is obtained from the equation
!       ! alpha = exp(- (Zi-Zj-1)/Zo) where Zo is a characteristic length scale    


!       ! compute alpha
!       zm    = state%zm     
!       alpha(:, 1) = 0.
!       do ilev = 2,nlayers
!          alpha(:, ilev) = exp( -( zm (:, ilev-1) -  zm (:, ilev)) / Zo)
!       end do
       
!       ! generate 2 streams of random numbers
!       do isubcol = 1,nsubcol
!          do ilev = 1,nlayers
!             call kissvec(seed1, seed2, seed3, seed4, rand_num)
!             CDF(isubcol, ilev) = rand_num
!             call kissvec(seed1, seed2, seed3, seed4, rand_num)
!             CDF2(isubcol, ilev) = rand_num
!          end do
!       end do

!       ! generate random numbers
!       do ilev = 2,nlayers
!          where (CDF2(:, ilev) < spread(alpha (:,ilev), dim=1, nCopies=nsubcol) )
!             CDF(:,ilev) = CDF(:,ilev-1) 
!          end where
!       end do

      end select

 
! -- generate subcolumns for homogeneous clouds -----

      do ilev = 1, nlayers
         iscloudy(:,ilev) = (CDF(:,ilev) >= 1._jprb - spread(cldf(ilev), dim=1, nCopies=nsubcol) )
      enddo

! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0

      do ilev = 1, nlayers
         where (iscloudy(:,ilev) )
            cld_stoch(:,ilev) = 1._jprb
         elsewhere (.not. iscloudy(:,ilev) )
            cld_stoch(:,ilev) = 0._jprb
         end where
      enddo

! where there is a cloud, set the subcolumn cloud properties;
! cloud-averaged and not grid-averaged cloud properties are assumed for single-column input;
! do NOT divide by cldf here to convert grid-averaged to cloud-averaged quantities

      do ilev = 1, nlayers
         where ( iscloudy(:,ilev) .and. (spread(cldf(ilev), dim=1, nCopies=nsubcol) > 0._jprb) )
            clwp_stoch(:,ilev) = spread(clwp(ilev), dim=1, nCopies=nsubcol)
            ciwp_stoch(:,ilev) = spread(ciwp(ilev), dim=1, nCopies=nsubcol)
         elsewhere
            clwp_stoch(:,ilev) = 0._jprb
            ciwp_stoch(:,ilev) = 0._jprb
         end where
      enddo
      do ilev = 1, nlayers
         do isubcol = 1, ngpt
            if ( iscloudy(isubcol,ilev) .and. (cldf(ilev) > 0._jprb) ) then
               n = ngb(isubcol)  
               tauc_stoch(isubcol,ilev) = tauc(n,ilev)
               ssac_stoch(isubcol,ilev) = ssac(n,ilev)
               asmc_stoch(isubcol,ilev) = asmc(n,ilev)
            else
               tauc_stoch(isubcol,ilev) = 0._jprb
               ssac_stoch(isubcol,ilev) = 1._jprb
               asmc_stoch(isubcol,ilev) = 1._jprb
            endif
         enddo
      enddo

! -- compute the means of the subcolumns ---
!      mean_cld_stoch(:) = 0._jprb
!      mean_clwp_stoch(:) = 0._jprb
!      mean_ciwp_stoch(:) = 0._jprb
!      mean_tauc_stoch(:) = 0._jprb
!      mean_ssac_stoch(:) = 0._jprb
!      mean_asmc_stoch(:) = 0._jprb
!      do i = 1, nsubcol
!         mean_cld_stoch(:) =  cld_stoch(i,:) + mean_cld_stoch(:) 
!         mean_clwp_stoch( :) =  clwp_stoch( i,:) + mean_clwp_stoch( :) 
!         mean_ciwp_stoch( :) =  ciwp_stoch( i,:) + mean_ciwp_stoch( :) 
!         mean_tauc_stoch( :) =  tauc_stoch( i,:) + mean_tauc_stoch( :) 
!         mean_ssac_stoch( :) =  ssac_stoch( i,:) + mean_ssac_stoch( :) 
!         mean_asmc_stoch( :) =  asmc_stoch( i,:) + mean_asmc_stoch( :) 
!      end do
!      mean_cld_stoch(:) = mean_cld_stoch(:) / nsubcol
!      mean_clwp_stoch( :) = mean_clwp_stoch( :) / nsubcol
!      mean_ciwp_stoch( :) = mean_ciwp_stoch( :) / nsubcol
!      mean_tauc_stoch( :) = mean_tauc_stoch( :) / nsubcol
!      mean_ssac_stoch( :) = mean_ssac_stoch( :) / nsubcol
!      mean_asmc_stoch( :) = mean_asmc_stoch( :) / nsubcol

      end subroutine generate_stochastic_clouds


!-------------------------------------------------------------------------------------------------- 
      subroutine kissvec(seed1,seed2,seed3,seed4,ran_arr)
!-------------------------------------------------------------------------------------------------- 

! public domain code
! made available from http://www.fortran.com/
! downloaded by pjr on 03/16/04
! converted to vector form, functions inlined by pjr,mvr on 05/10/2004

! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123; 
!
      implicit none

      real(kind=jprb), intent(inout)  :: ran_arr
      integer(kind=jpim), intent(inout) :: seed1,seed2,seed3,seed4
!      integer(kind=jpim) :: i,sz,kiss
      integer(kind=jpim) :: kiss
      integer(kind=jpim) :: m, k, n

! inline function 
      m(k, n) = ieor (k, ishft (k, n) )

!      sz = size(ran_arr)
!      do i = 1, sz
         seed1 = 69069_jpim * seed1 + 1327217885_jpim
         seed2 = m (m (m (seed2, 13_jpim), - 17_jpim), 5_jpim)
         seed3 = 18000_jpim * iand (seed3, 65535_jpim) + ishft (seed3, - 16_jpim)
         seed4 = 30903_jpim * iand (seed4, 65535_jpim) + ishft (seed4, - 16_jpim)
         kiss = seed1 + seed2 + ishft (seed3, 16_jpim) + seed4
         ran_arr = kiss*2.328306e-10_jprb + 0.5_jprb
!     end do
    
      end subroutine kissvec

      end module mcica_subcol_gen_sw


