C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

C************************  SUBROUTINE READPROF  *****************************C

      SUBROUTINE READPROF(NLAYERSX,IOUTX,ICLDX,IAERX,ISCCOSX,IDELMX,
     &      PAVELX,TAVELX,PZX,TZX,ZENITH,ADJFLUX,SEMISS,
     &      WKLX,COLDRYX,INFLAGX,ICEFLAGX,LIQFLAGX,
     &      CLDFRACX,CLDDAT1X,CLDDAT2X,CLDDAT3X,CLDDAT4X,CLDDATMOMX,
     &      TAUAERX,SSAAERX,PHASEX)

C     Read in atmospheric profile.

      IMPLICIT DOUBLE PRECISION (V)                                      
                                                                         
      PARAMETER (MAXINPX=35)
      PARAMETER (MAXXSEC=4)
C      PARAMETER (MAXPROD = MXLAY*MAXXSEC)

      PARAMETER (MCMU = 32)

      INCLUDE 'param.f'
c

      DIMENSION ALTZ(0:MXLAY),IXTRANS(14)

      COMMON /CONTROL/  IAER, NSTR, IOUT, ISTART, IEND, ICLD,
     &                  idelm, isccos
      COMMON /CONSTANTS/FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
c      COMMON /SWPROP/   ZENITH, ALBEDO, ADJFLUX
c      COMMON /SURFACE/  IREFLECT,SEMISS(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBRODL(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /IFIL/     IRD,IPR,IPU,IDUM(15)
      COMMON /XSECCTRL/ NXMOL,IXINDX(MAXINPX)
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /PATHX/    IXMAX,NXMOL0,IXINDX0(MAXINPX),WX0(MAXINPX,MXLAY)    
      COMMON /XRRTATM/  IXSECT

      COMMON /CLOUDIN/   INFLAG,CLDDAT1(MXLAY),CLDDAT2(MXLAY),
     &                   ICEFLAG,LIQFLAG,CLDDAT3(MXLAY),CLDDAT4(MXLAY),
     &                   CLDDATMOM(0:16,MXLAY)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),
     &     TAUCLOUD(MXLAY,NBANDS),SSACLOUD(MXLAY,NBANDS),
     &     xmom(0:16,MXLAY,NBANDS)
      common /AERDAT/    ssaaer(mxlay,nbands), phase(mcmu,mxlay,nbands), 
     &                   tauaer(mxlay,nbands)

      REAL*8 PAVELX(MXLAY),TAVELX(MXLAY),PZX(0:MXLAY),TZX(0:MXLAY)
      REAL*8 COLDRYX(MXLAY),WKLX(35,MXLAY)
      REAL*8 CLDFRACX(MXLAY),CLDDAT1X(MXLAY),CLDDAT2X(MXLAY),
     &          CLDDAT3X(MXLAY),CLDDAT4X(MXLAY),
     &          CLDDATMOMX(0:16,MXLAY)
      REAL*8 SEMISS(NBANDS),ZENITH,ADJFLUX
      REAL*8 ssaaerx(mxlay,nbands), phasex(mcmu,mxlay,nbands), 
     &       tauaerx(mxlay,nbands)

      CHARACTER*80 FORM1(0:1),FORM2(0:1),FORM3(0:1)
      CHARACTER*1 CTEST, CDOLLAR, CDUM

      DATA CDOLLAR /'$'/
      DATA IXTRANS /0,0,0,1,2,3,0,0,0,0,0,4,0,0/
C      DATA WX /MAXPROD*0.0/

      FORM1(0) = '(3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))'
      FORM2(0) = '(3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))'
      FORM3(0) = '(8E10.3)'
      FORM1(1) = '(G15.7,G10.4,G10.4,A3,I2,1X,2(G7.2,G8.3,G7.2))'
      FORM2(1) = '(G15.7,G10.4,G10.4,A3,I2,23X,(G7.2,G8.3,G7.2))'
      FORM3(1) = '(8G15.7)'

C  Initialize molecular amount and cross section arrays to zero here.
      
      DO 1200 ILAY = 1,MXLAY
         DO 1100 ISP = 1,35
 1100       WKL(ISP,ILAY) = 0.0
         DO 1150 ISP = 1,MAXXSEC
 1150       WX(ISP,ILAY) = 0.0
 1200 CONTINUE

      IXMAX = MAXINPX
      IRD = 9
      OPEN (IRD,FILE='INPUT_RRTM',FORM='FORMATTED')

 1000 CONTINUE
      READ (IRD,9009,END=8800) CTEST
      IF (CTEST .NE. CDOLLAR) GO TO 1000
      READ (IRD,9011) IAER, IATM, ISCAT, ISTRM, IOUT, ICLD, IDELM, ICOS

      if (idelm.gt.1 .or. idelm.lt.0 .or. icos.gt.2 .or. icos.lt.0) then
         print *,'INVALID MEASUREMENT COMPARISON FLAG'
         stop
      endif
      isccos = icos

C     No cross-sections implemented in shortwave.
      IXSECT = 0

      IF (ISCAT .NE. 1) THEN
         PRINT *,' INVALID SCATTERING OPTION CHOSEN'
         STOP
      ENDIF

c mji - 2 stream version only
      IF (ISTRM .EQ. 0) THEN 
         NSTR = 2
      ELSE 
         PRINT *, 'INVALID VALUE FOR ISTRM'
         STOP
      ENDIF

      IF (IAER.EQ.10) CALL READAER

C     If clouds are present, read in appropriate input file, IN_CLD_RRTM.
      IF (ICLD .EQ. 1) CALL READCLD


      READ (IRD,9020) JULDAT, SZA
      ZENITH = COS(SZA * PI / 180.)
      IF (JULDAT .EQ. 0) THEN
         ADJFLUX = 1.
      ELSE
         ADJFLUX = EARTH_SUN (JULDAT)
      ENDIF
      READ (IRD,9012) IEMIS, IREFLECT, (SEMISS(IB),IB=IB1,IB2)
      IF (IEMIS .EQ. 0) THEN
         DO 1500 IB = IB1, IB2
            SEMISS(IB) = 1.
 1500    CONTINUE
      ELSEIF (IEMIS .EQ. 1) THEN
         DO 1600 IB = IB1, IB2
            SEMISS(IB) = SEMISS(IB1)
 1600    CONTINUE
      ELSEIF (IEMIS .EQ. 2) THEN
C          PRINT *, 'THESE ARE THE INPUT EMISSIVITY VALUES'
C          PRINT *, SEMISS(IB1:IB2)
      ELSE
          PRINT *, 'IEMIS = ', IEMIS, ' NOT A VALID INPUT VALUE'
          STOP
      ENDIF
     
      IF (IATM .EQ. 0) THEN
         READ (IRD,9013) IFORM,NLAYERS,NMOL
         IF (NMOL.EQ.0) NMOL = 7                                    
         READ (IRD,FORM1(IFORM)) PAVEL(1),TAVEL(1),SECNTK,CINP,
     &        IPTHAK,ALTZ(0),PZ(0),TZ(0),ALTZ(1),PZ(1),TZ(1)
         READ (IRD,FORM3(IFORM)) (WKL(M,1),M=1,7), WBRODL(1)
         IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,1),M=8,NMOL)

         DO 2000 L = 2, NLAYERS
            READ (IRD,FORM2(IFORM)) PAVEL(L),TAVEL(L),SECNTK,CINP,
     &           IPTHRK,ALTZ(L),PZ(L),TZ(L)
            READ (IRD,FORM3(IFORM)) (WKL(M,L),M=1,7), WBRODL(L)
            IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,L),M=8,NMOL)
 2000    CONTINUE                                                            
           
         IF (IXSECT .EQ. 1) THEN                                 
            READ (IRD,9300) NXMOL0
            NXMOL = NXMOL0
            CALL XSIDENT(IRD)
            READ (IRD,9301) IFORMX
C     
            DO 3000 L = 1, NLAYERS       
               READ (IRD,9010) CDUM
               READ (IRD, FORM3(IFORMX)) (WX0(M,L),M=1,7),WBRODX    
               IF (NXMOL0 .GT. 7) READ (IRD,FORM3(IFORMX)) 
     &              (WX0(M,L),M=8,NXMOL0)
 3000       CONTINUE
         ENDIF
      ELSE
         IPU = 7
         IPR = 66
         OPEN(UNIT=IPR,FILE='TAPE6',STATUS='UNKNOWN')
         CALL RRTATM
         IF (IXSECT .EQ. 1) THEN
            DO 3300 MX = 1, NXMOL0
               IXINDX(MX) = IXTRANS(IXINDX0(MX))
 3300       CONTINUE
         ENDIF
      ENDIF

C     Test for mixing ratio input.
      IMIX = 1
      DO 3500 M = 1, NMOL
         IF (WKL(M,1) .GT. 1.0) THEN
            IMIX = 0
            GO TO 3600
         ENDIF
 3500 CONTINUE
 3600 CONTINUE

      IF (IXSECT .EQ. 1) THEN
         IMIXX = 0
         IF (WX0(1,1) .LE. 1.0) IMIXX = 1
      ENDIF
      DO 5000 L = 1, NLAYERS
         SUMMOL = 0.0
         DO 4100 IMOL = 2, NMOL
            SUMMOL = SUMMOL + WKL(IMOL,L)
 4100    CONTINUE
         IF (IMIX .EQ. 1) THEN
            COLDRY(L) = WBRODL(L) / (1. - SUMMOL)
            DO 4200 IMOL = 1, NMOL
               WKL(IMOL,L) = COLDRY(L) * WKL(IMOL,L)
 4200       CONTINUE
         ELSE
            COLDRY(L) = WBRODL(L) + SUMMOL
         ENDIF
         IF (IXSECT .EQ. 1) THEN
            DO 4400 IX = 1, NXMOL0
               IF (IXINDX(IX) .NE. 0) THEN
                  IF (IMIXX .EQ. 1) THEN
                     WX(IXINDX(IX),L) = COLDRY(L) * WX0(IX,L) * 1.E-20
                  ELSE
                     WX(IXINDX(IX),L) = WX0(IX,L) * 1.E-20
                  ENDIF
               ENDIF
 4400       CONTINUE
         ENDIF
 5000 CONTINUE

c mji out
c      DO L = 1, NLAYERS
c        print*, L,PAVEL(L),TAVEL(L),
c     &              (WKL(JA,L),JA=1,4),WKL(6,L),COLDRY(L)
cc9200  format(1x,'READPROF ',I3,2F7.1,6E13.5)
c      END DO
c!

      CLOSE(IRD)

C Pass RRTM COMMON arrays into dummy arrays 
      NLAYERSX = NLAYERS
      IOUTX = IOUT
      ICLDX = ICLD
      IAERX = IAER
      ISCCOSX = ISCCOS
      IDELMX = IDELM
      INFLAGX = INFLAG
      ICEFLAGX = ICEFLAG
      LIQFLAGX = LIQFLAG
      PZX(0) = PZ(0)
      TZX(0) = TZ(0)
      DO 7000 L = 1, NLAYERS
         PAVELX(L) = PAVEL(L)
         TAVELX(L) = TAVEL(L)
         PZX(L) = PZ(L)
         TZX(L) = TZ(L)
         COLDRYX(L) = COLDRY(L)
         CLDFRACX(L) = CLDFRAC(L)
         CLDDAT1X(L) = CLDDAT1(L)
         CLDDAT2X(L) = CLDDAT2(L)
         CLDDAT3X(L) = CLDDAT3(L)
         CLDDAT4X(L) = CLDDAT4(L)
         DO 7100 IMOL = 1,NMOL
            WKLX(IMOL,L) = WKL(IMOL,L)
 7100    CONTINUE
         DO 7200 NS = 0,NSTR 
           CLDDATMOMX(NS,L) = CLDDATMOM(NS,L)
 7200    CONTINUE
 7000 CONTINUE
      DO 7500 L = 1, NLAYERS
         DO 7400 NB = IB1,IB2
           SSAAERX(L,NB) = SSAAER(L,NB)
           TAUAERX(L,NB) = TAUAER(L,NB)
         DO 7300 NS = 1,NSTR 
           PHASEX(NS,L,NB) = PHASE(NS,L,NB)
 7300    CONTINUE
 7400    CONTINUE
 7500 CONTINUE

c mji out
c      DO JB=16,29
c        DO JK=1,NLAYERS
c         print*, JB,JK,TAUAERX(JK,JB),SSAAERX(JK,JB),
c     &      PHASEX(1,JK,JB)
c 9102    format(1x,'RDPROF Aerx :',2I3,7E12.5)
c        END DO
c      END DO
c      DO L = 1, NLAYERS
c        print*, L,PAVELX(L),TAVELX(L),
c     &              (WKLX(JA,L),JA=1,4),WKLX(6,L),COLDRYX(L)
c9210  format(1x,'READPROF(X) ',I3,2F7.1,6E13.5)
c      END DO
c      DO L = 1, NLAYERS
c        print*, L,CLDFRACX(L),CLDDAT1X(L),CLDDAT2X(L),CLDDAT3X(L),
c     &              CLDDAT4X(L),(CLDDATMOMX(NS,L),NS=0,NSTR)
c9220  format(1x,'READPROF(CLDX) ',I3,22E13.5)
c      END DO
c!

      GO TO 9000

 8800 CONTINUE
      STOP ' INVALID INPUT_RRTM '

 9000 CONTINUE

 9009 FORMAT (A1,1X,I2,I2,I2)
 9010 FORMAT (A1)
 9011 FORMAT (18X,I2,29X,I1,32X,I1,1X,I1,2X,I3,4X,I1,3x,i1,i1)
 9012 FORMAT (11X,I1,2X,I1,14F5.3)
 9013 FORMAT (1X,I1,I3,I5)                                     
 9020 format (12X, I3, 3X, F7.4)
 9300 FORMAT (I5)
 9301 FORMAT (1X,I1)

      RETURN
      END 

C************************  SUBROUTINE READCLD  *****************************C

      SUBROUTINE READCLD

C     Purpose:  To read in IN_CLD_RRTM_SW, the file that contains input 
C               cloud properties.

      INCLUDE 'param.f'
c
      COMMON /CONTROL/   IAER, NSTR, IOUT, ISTART, IEND, ICLD,
     &                   idelm, isccos
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /CLOUDIN/   INFLAG,CLDDAT1(MXLAY),CLDDAT2(MXLAY),
     &                   ICEFLAG,LIQFLAG,CLDDAT3(MXLAY),CLDDAT4(MXLAY),
     &                   CLDDATMOM(0:16,MXLAY)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),
     &     TAUCLOUD(MXLAY,NBANDS),SSACLOUD(MXLAY,NBANDS),
     &     xmom(0:16,MXLAY,NBANDS)
      CHARACTER*1 CTEST, CPERCENT

      DATA CPERCENT /'%'/
      IRDCLD = 11

      OPEN(IRDCLD,FILE='IN_CLD_RRTM',FORM='FORMATTED')

C     Read in cloud input option.  
      READ(IRDCLD,9050) INFLAG, ICEFLAG, LIQFLAG

      DO 500 LAY = 1, NLAYERS
         CLDFRAC(LAY) = 0.
 500  CONTINUE

      IF (INFLAG .EQ. 0) THEN
 950     CONTINUE
c     For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
c     DAT1 are pertinent.  If CTEST = '%', then there are no more 
C     cloudy layers to process.
         READ (IRDCLD,9099,END=8950) CTEST,LAY,FRAC,
     &        DAT1,DAT2,(CLDDATMOM(ISTR,LAY),ISTR=0,NSTR)
         IF (CTEST .EQ. CPERCENT) GO TO 8950
         CLDFRAC(LAY) = FRAC
         CLDDAT1(LAY) = DAT1
         CLDDAT2(LAY) = DAT2
         GO TO 950
 8950    CONTINUE
      ELSE
 1000    CONTINUE
c     For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
c     DAT1 are pertinent.  If CTEST = '%', then there are no more 
C     cloudy layers to process.
         READ (IRDCLD,9100,END=9000) CTEST,LAY,FRAC,
     &        DAT1,DAT2,DAT3,DAT4
         IF (CTEST .EQ. CPERCENT) GO TO 9000
         CLDFRAC(LAY) = FRAC
         CLDDAT1(LAY) = DAT1
         CLDDAT2(LAY) = DAT2
         CLDDAT3(LAY) = DAT3
         CLDDAT4(LAY) = DAT4
         GO TO 1000
 9000    CONTINUE
      ENDIF

c mji out
c      DO L = 1, NLAYERS
c        print 9210, L,CLDFRAC(L),CLDDAT1(L),CLDDAT2(L),CLDDAT3(L),
c     &              CLDDAT4(L),(CLDDATMOM(NS,L),NS=0,NSTR)
c9210  format(1x,'READPROF(CLD) ',I3,22E13.5)
c      END DO


      CLOSE(IRDCLD)

 9050 FORMAT (3X,I2,4X,I1,4X,I1)
 9099 FORMAT (A1,1X,I3,19E10.5)
 9100 FORMAT (A1,1X,I3,5E10.5)
      RETURN
      END

C************************  SUBROUTINE READAER  *****************************C

      SUBROUTINE READAER

C     Purpose:  To read in IN_AER_RRTM, the file that contains input
C               aerosol properties.

      INCLUDE 'param.f'
c

      PARAMETER (MCMU = 32)
      real aerpar(3), ssa(nbands), asym(nbands), aod(mxlay),aod1(nbands)
      real rl1(nbands), rl2(nbands), rlambda(nbands), specfac(nbands)

      COMMON /CONTROL/   IAER, NSTR, IOUT, ISTART, IEND, ICLD,
     &                   idelm, isccos
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      common /AERDAT/    ssaaer(mxlay,nbands), phase(mcmu,mxlay,nbands), 
     &                   tauaer(mxlay,nbands)
      CHARACTER*1 CTEST, CPERCENT

      DATA CPERCENT /'%'/

      integer lay(mxlay)

      eps = 1.e-10
      IRDAER = 12
      OPEN(IRDAER,FILE='IN_AER_RRTM',FORM='FORMATTED')

      do ib = ib1, ib2
         DO 500 ILAY = 1, MXLAY
            AOD(ILAY) = 0.
            tauaer(ILAY,ib) = 0.
 500     CONTINUE
         rl1(ib) = 10000. / wavenum1(ib)	
         rl2(ib) = 10000. / wavenum2(ib)	
      enddo

C     Read in number of different aerosol models option.
      read (irdaer, 9010) naer
c      if (naer .gt. 4) then
c         print *, 'NAER (= ', naer, ') IS GREATER THAN 4'
c         stop
c      endif
        
c     For each aerosol read in optical properties and layer aerosol 
c     optical depths.
      do ia = 1, naer
	 read (irdaer, 9011) nlay, iaod, issa, iasym, (aerpar(i),i=1,3)

         if (iaod .eq. 0) then
c           Set defaults to get standard Angstrom relation.
            if (aerpar(2) .lt. eps) aerpar(2) = 1.

            omaer  = 1. - aerpar(1)
            do ib = ib1, ib2
               if (omaer .ne. 0.) then
                  factor = 1. / omaer
                  aodbar = factor *(rl2(ib)**omaer - rl1(ib)**omaer)/ 
     &                 (rl2(ib) - rl1(ib))
               else
                  aodbar = alog(rl2(ib)/rl1(ib))/(rl2(ib) - rl1(ib))
               endif
               if (aerpar(1) .ne. 0.) then
                  rlambda(ib) = (1. / aodbar) ** (1./ aerpar(1))
               else
                  rlambda(ib) = 1.0
               endif
               specfac(ib) = (aerpar(2) + aerpar(3) * rlambda(ib)) /
     &              ((aerpar(2) + aerpar(3) - 1.) + 
     &              rlambda(ib)**aerpar(1))
            enddo
         endif
C        For this aerosol, read in layers and optical depth information.
C        Store a nonzero optical depth in aod to check for double
C        specification.
         do il = 1, nlay
            read(irdaer, 9012) lay(il), (aod1(ib), ib = ib1,ib2)

            if (aod(lay(il)) .lt. eps) then
               if (iaod .eq. 0) then
                  aod(lay(il)) = aod1(ib1)
                  do ib = ib1, ib2
                     tauaer(lay(il),ib) = aod1(ib1) * specfac(ib)
                  enddo
               else
                  do ib = ib1, ib2
                     aod(lay(il)) = max(aod(lay(il)),aod1(ib))
                     tauaer(lay(il),ib) = aod1(ib)
                  enddo
               endif
            else
               print *,'LAYER ',lay(il),' HAS MORE THAN 
     &              ONE AEROSOL TYPE'
               stop
            endif
         enddo

c        For this aerosol, read and store optical properties
         read (irdaer, 9013) (ssa(ib), ib = ib1,ib2)

         DO 2000 IB = IB1, IB2
            do il = 1, nlay
               if (issa .eq. 0) then 
                  ssaaer(lay(il),IB) = ssa(ib1)
               else
                  ssaaer(lay(il),IB) = ssa(IB)
               endif
            enddo
 2000    CONTINUE

         if (iasym .lt. 2) then
            read (irdaer, 9013) (asym(ib), ib = ib1,ib2)

            DO 3000 IB = IB1, IB2
               do il = 1, nlay
                  do istr = 1,  nstr
                     if (iasym .eq. 0) then 
                        phase(istr,lay(il),IB) = asym(ib1)**istr
                     elseif (iasym .eq. 1) then
                        phase(istr,lay(il),IB) = asym(IB)**istr
                     endif
                  enddo
               enddo
 3000       CONTINUE
         else
            do il = 1, nlay
               do istr = 1, nstr
                  read (irdaer, 9013) (phase(istr,lay(il),ib), 
     &                 ib = ib1,ib2)
               enddo
            enddo
         endif

      enddo


 9000 CONTINUE
      CLOSE(IRDAER)

 9010 format (4x, i1)
 9011 format (2x, i3, 4x, i1, 4x, i1, 4x, i1, 3f8.3)
 9012 format (2x, i3, 14f7.4)
 9013 format (14f5.3)

      RETURN
      END

C************************  SUBROUTINE XSIDENT  *****************************C

      SUBROUTINE XSIDENT(IRD)
C                                                                         
C     This subroutine identifies which cross-sections are to be used.

      PARAMETER (MAXINPX=35)
      PARAMETER (MAXXSEC=4)

      IMPLICIT DOUBLE PRECISION (V)                                     ! 
C                                                                         
      COMMON /XSECCTRL/ NXMOL,IXINDX(MAXINPX)
C                                                                         
C     NXMOL     - number of cross-sections input by user
C     IXINDX(I) - index of cross-section molecule corresponding to Ith
C                 cross-section specified by user
C                 = 0 -- not allowed in RRTM
C                 = 1 -- CCL4
C                 = 2 -- CFC11
C                 = 3 -- CFC12
C                 = 4 -- CFC22
C                                                                         
C     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          
C                                                                         
      CHARACTER*10 XSNAME(MAXINPX),ALIAS(MAXXSEC,4),BLANK               
C                                                                         
      DATA (ALIAS(1,I),I=1,4)/                                           
     *    'CCL4      ', 'CCL3F     ', 'CCL2F2    ', 'CHCLF2    '/ 
      DATA (ALIAS(2,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'CFCL3     ', 'CF2CL2    ', 'CHF2CL    '/         
      DATA (ALIAS(3,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'CFC11     ', 'CFC12     ', 'CFC22     '/         
      DATA (ALIAS(4,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'F11       ', 'F12       ', 'F22       '/        

      DATA BLANK / '          '/                                          
C                                                                         
      DO 10 I = 1, NXMOL                                                 
         XSNAME(I) = BLANK                                                
   10 CONTINUE                                                            
C                                                                         
C     READ IN THE NAMES OF THE MOLECULES                                  
C                                                                         
      IF (NXMOL.GT.7) THEN                                               
         READ (IRD,'(7A10)') (XSNAME(I),I=1,7)                            
         READ (IRD,'(8A10)') (XSNAME(I),I=8,NXMOL)                       
      ELSE                                                                
         READ (IRD,'(7A10)') (XSNAME(I),I=1,NXMOL)                       
      ENDIF                                                               
C                                                                         
C     MATCH THE NAMES READ IN AGAINST THE NAMES STORED IN ALIAS           
C     AND DETERMINE THE INDEX VALUE.  
      IXMAX = 4                                                          
      DO 40 I = 1, NXMOL                                                 
C        Left-justify all inputed names.                                      
         CALL CLJUST (XSNAME(I),10)
         IXINDX(I) = 0
         DO 20 J = 1, IXMAX
            IF ((XSNAME(I).EQ.ALIAS(1,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(2,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(3,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(4,J))) THEN                           
               IXINDX(I) = J                                              
            ENDIF                                                         
   20    CONTINUE
   40 CONTINUE                                                            

      RETURN
      END

*****************************************************************
      BLOCK DATA

      INCLUDE 'param.f'

      COMMON /HVERSN/ HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *                HVDUM1(4),HVRUTL,HVREXT,
     *                HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
     *                HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR

      COMMON /HVRSNB/    HVRKG(16:29)

      CHARACTER*15 HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVDUM1,HVRUTL,HVREXT,
     *            HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
     *            HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR


      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /FEATURES/  NG(IB1:IB2),NSPA(IB1:IB2),NSPB(IB1:IB2)

      DATA HVRRTM / 'NOT USED' /,  
     *     HVRRTR / 'NOT USED' /,   HVRATM / 'NOT USED' /,
     *     HVRSET / 'NOT USED' /,   HVRTAU / 'NOT USED' /,
     *     HVDUM1 / 4*'NOT USED' /, HVRUTL / 'NOT USED' /,
     *     HVREXT / 'NOT USED' /, 
     *     HVRD1M / 'NOT USED' /,   HVRR1M / 'NOT USED' /,
     *     HVREPK / 'NOT USED' /,   HVRLPK / 'NOT USED' /,
     *     HVRAER / 'NOT USED' /,
     *     HVRCLD / 'NOT USED' /,
     *     HVRDIS / 'NOT USED' /,   HVRLAM / 'NOT USED' /,
     *     HVRPAR / 'NOT USED' /


      DATA WAVENUM1(16) /2600./,WAVENUM2(16) /3250./,DELWAVE(16) /650./
      DATA WAVENUM1(17) /3250./,WAVENUM2(17) /4000./,DELWAVE(17) /750./
      DATA WAVENUM1(18) /4000./,WAVENUM2(18) /4650./,DELWAVE(18) /650./
      DATA WAVENUM1(19) /4650./,WAVENUM2(19) /5150./,DELWAVE(19) /500./
      DATA WAVENUM1(20) /5150./,WAVENUM2(20) /6150./,DELWAVE(20) /1000./
      DATA WAVENUM1(21) /6150./,WAVENUM2(21) /7700./,DELWAVE(21) /1550./
      DATA WAVENUM1(22) /7700./,WAVENUM2(22) /8050./,DELWAVE(22) /350./
      DATA WAVENUM1(23) /8050./,WAVENUM2(23)/12850./,DELWAVE(23) /4800./
      DATA WAVENUM1(24)/12850./,WAVENUM2(24)/16000./,DELWAVE(24) /3150./
      DATA WAVENUM1(25)/16000./,WAVENUM2(25)/22650./,DELWAVE(25) /6650./
      DATA WAVENUM1(26)/22650./,WAVENUM2(26)/29000./,DELWAVE(26) /6350./
      DATA WAVENUM1(27)/29000./,WAVENUM2(27)/38000./,DELWAVE(27) /9000./
      DATA WAVENUM1(28)/38000./,WAVENUM2(28)/50000./,DELWAVE(28)/12000./
      DATA WAVENUM1(29)/820./,  WAVENUM2(29)/2600./, DELWAVE(29)/1780./

      DATA NG  /16,16,16,16,16,16,16,16,16,16,16,16,16,16/
      DATA NSPA /9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 0, 1, 9, 1/
      DATA NSPB /1, 5, 1, 1, 1, 5, 1, 0, 1, 0, 0, 1, 5, 1/

C     HEATFAC is the factor by which one must multiply delta-flux/ 
C     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
C     the heating rate in units of degrees/day.  It is equal to 
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(3600)(1e-5)/(1.004)
      DATA HEATFAC /8.4391/


      END
c**********************************************************************
      Block Data phys_consts
c
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
c
      DATA PI / 3.1415927410125732 /
c
c    Constants from NIST 01/11/2002
c
      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
     *     CLIGHT / 2.99792458E+10 /, 
     *     AVOGAD / 6.02214199E+23 /, ALOSMT / 2.6867775E+19 /,
     *     GASCON / 8.314472  E+07 /
     *     RADCN1 / 1.191042722E-12 /, RADCN2 / 1.4387752    /
c
c     Pi was obtained from   PI = 2.*ASIN(1.)                             A03980
c
c     units are genrally cgs
c
c     The first and second radiation constants are taken from NIST.
c     They were previously obtained from the relations:
c                            RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07      A03990
c                            RADCN2 = PLANCK*CLIGHT/BOLTZ                 A04000
      end
c

