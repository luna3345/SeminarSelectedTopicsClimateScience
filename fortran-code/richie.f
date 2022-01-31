                           PROGRAM RICHIE
C
C     ******************************************************************
C     *                                                                *
C     *        FINE-MESH LIMITED AREA PRIMITIVE EQUATION MODEL.        *
C     *                                                                *
C     *                                                                *
C     *        THIS MODEL IS A REALIZATION OF THE NUMERICAL            *
C     *            SCHEME DESCRIBED IN RICHARDSON'S BOOK               *
C     *                                                                *
C     *           Weather Prediction by Numerical Process              *
C     *                                                                *
C     ******************************************************************
C
C     **************************************************************
C                                                                 
C         Coding started on 29th Jan, 1992, IMS, Dublin. 
C		  *Edited by Luna Lehmann in 2022, KFU, Graz.
C
C     **************************************************************
C                                                                
C        Units here are SI. LFR's Momenta are ten times bigger. 
C        LFR's pressures are in microbars or deci-Pascals.     
C                                                             
C        Trick: prognostic variables on E-grid are interpolated 
C               to a finer A-grid. The forecast step is done
C               on the A-grid (with due care in choosing
C               finite differences). Then, the interpolation
C               is repeated at the beginning of each step.
C               The result is to double the computation, but
C               the program is MUCH easier to code on the A-grid.
C                                                             
C     **************************************************************
C
c          Parameters to determine size of grid
           PARAMETER
     *                ( IM= 9, JM=5, LM=5,
     *                  IMM1=IM-1, IMM2=IM-2, 
     *                  JMM1=JM-1, JMM2=JM-2,
     *                  LMM1=LM-1, LMM2=LM-2, 
     *                  KB=IM+2*JM-3          )  
C
c=================================================================
c
c      2-D Work fields common to all modules.
       COMMON /WORK/ WORK1(IM,JM), WORK2(IM,JM), WORK3(IM,JM)
c
c=================================================================
C
          DIMENSION
     *    COSPHI(JM,2), FCOR(JM,2), H1(JM,2), RH1(JM,2)

C         Divergence and Vorticity
          DIMENSION DIV (IM,JM,LM)
CCC       DIMENSION VORT(IM,JM,LM)
C
c=================================================================
c
c       space for the basic INPUT variables 
        REAL PS(IM,JM), Z(IM,JM,LMM1)
        REAL U(IM,JM,LM),V(IM,JM,LM), DDDFF(IM,JM,LM)
        REAL OROG(IM,2*JM)
        COMMON / WINPUT / PS, Z, U, V, DDDFF, OROG
c
c       space for the basic PROGNOSTIC variables
        real P(IM,JM,LM),MU(IM,JM,LM),MV(IM,JM,LM),TTOP(IM,JM)
        COMMON / WPROG / P, MU, MV, TTOP
c
c    -     -    -    -    -    -    -    -    -    -    -    -
c
c       size of fine A-grid
        PARAMETER (JM2=10,JM2M1=JM2-1,JM2M2=JM2-2,JM2M3=JM2-3)
C
c       space for the basic PROGNOSTIC variables on FINE grid
        real P2(IM,JM2,LM),MU2(IM,JM2,LM),
     X       MV2(IM,JM2,LM),TTOP2(IM,JM2)
        COMMON / WPROG2 / P2, MU2, MV2, TTOP2
C
c       space for various DIAGNOSTIC variables on FINE grid
        real PP2(IM,JM2,LM),U2(IM,JM2,LM),
     X       V2(IM,JM2,LM),W2(IM,JM2,LM),
     X       RR2(IM,JM2,LM),PMID2(IM,JM2,LM),
     X       RHO2(IM,JM2,LM),T2(IM,JM2,LM),
     X       RHO2I(IM,JM2,LM),T2I(IM,JM2,LM)
        COMMON / WDIAG2 / PP2, U2, V2, W2, 
     X       RR2, PMID2, RHO2, T2, RHO2I,T2I
c
c       space for the basic HISTORICAL variables on FINE grid
        real P2OLD(IM,JM2,LM),MU2OLD(IM,JM2,LM),
     X       MV2OLD(IM,JM2,LM),TT2OLD(IM,JM2)
        real P2NEW(IM,JM2,LM),MU2NEW(IM,JM2,LM),
     X       MV2NEW(IM,JM2,LM),TT2NEW(IM,JM2)
        COMMON / WOLDNEW / P2OLD, MU2OLD, MV2OLD, TT2OLD,
     X                     P2NEW, MU2NEW, MV2NEW, TT2NEW
C
c       Space for various TENDENCY variables on FINE grid
        real DP2DT (IM,JM2,LM),DRR2DT(IM,JM2,LM),DTT2DT(IM,JM2),
     X       DMU2DT(IM,JM2,LM),DMV2DT(IM,JM2,LM)
        COMMON /TENDY / DP2DT, DRR2DT, DTT2DT, DMU2DT, DMV2DT 
c
        COMMON / WMISC / BRAKET(IM,JM2,LM), DWDZ(IM,JM2,LM)
c
CCC     COMMON / PGRAD / DPDX(IM,JM2,LM), DPDY(IM,JM2,LM)
c
        REAL COSPH2(JM2), FCOR2(JM2), TAN2(JM2)
c
c       Work fields for plotting (on fine A-grid)
CCC     real ZPLOT(IM,2*JM)
        real PSEA2(IM,JM2),Z2(IM,JM2,LM)
c
c      2-D Work fields on A-grid.
       COMMON /WORKK/ WORKK1(IM,JM2), WORKK2(IM,JM2), 
     +                WORKK3(IM,JM2), WORKK4(IM,JM2)
C
C      Divergence of momentum and velocity.
       REAL  DIV2M(IM,JM2,LM), DIV2V(IM,JM2,LM)
c
c      1-D Work fields for vertical integrals, etc. 
       COMMON /WVERT / VERT1(LM), VERT2(LM), 
     +                 VERT3(LM), VERT4(LM),
     +                 VERT5(LM), VERT6(LM) 
C
c       space for the basic INITIALIZED variables on FINE grid
        real P2INI(IM,JM2,LM),MU2INI(IM,JM2,LM),
     X       MV2INI(IM,JM2,LM),TT2INI(IM,JM2)
        COMMON / WINIT2 / P2INI, MU2INI, MV2INI, TT2INI
C
c=================================================================
C
c       Conventional interface levels (km).
        real ZBAR(LM), PBAR(LM),TBAR(LM)
c
c      Work fields for other vertical quantities (SI units).
       COMMON /VERTIC/ ZLEV(LM), DELZ(LM)
C
        REAL KAPPA , GAMMA
C
        REAL NOISE1, NOISE2, BRUS
        REAL N1MID , N2MID , BRMID
c
        CHARACTER*10 STRING
        INTEGER YEAR, MONTH, DAY, HOUR
c
        LOGICAL OLDDAT, LFRDAT 
        LOGICAL TIMEFILT, MASSFILT, WINDFILT, DIVDAMP
c
        REAL WHOUR(10)
        CHARACTER*6 WNAME(10), WFILE
        CHARACTER*10 WIFILE, OLDFIL
c
c   	*the X indicates a line break in fortran since its on the 6th column
        PARAMETER ( NSMAX=1000)
        REAL 
     X       XVAR(0:NSMAX),YVAR(0:NSMAX),
     X       UVAR(0:NSMAX),VVAR(0:NSMAX),DVAR(0:NSMAX),
     X       KVAR(0:NSMAX),AVAR(0:NSMAX),PVAR(0:NSMAX),
     X       MVAR(0:NSMAX),
     X       N1VAR(0:NSMAX),N2VAR(0:NSMAX),BRVAR(0:NSMAX),
     X       N1PVAR(0:NSMAX),N2PVAR(0:NSMAX),BRPVAR(0:NSMAX)
C
C-----------------------------------------------------------------------
c
c       Parameters relating to initialization.
        LOGICAL INIT
        REAL HH(0:500), HH2(0:1000)
        DATA TAUC /6./, IWINDOW  /0/, WPARAM /0./
c
C       Control for geostrophic momentum initial conditions.
        LOGICAL ICGEOS
        DATA ICGEOS /.FALSE./
C
C-----------------------------------------------------------------------
C
        DATA ZBAR / 11.8,   7.2,   4.2 ,   2.0 ,   0.0 /
        DATA PBAR / 200.,  400.,   600.,   800.,  1013./
        DATA TBAR / -55.,  -32.,   -12.,    +2.,   +15./
c
                             DATA
     X   A/6376000./,  OMEGA2/.00014584/, GRAV/9.80/,
     X   PI/3.141592654/,  RGAS/287./, CP/1004./
c
         DATA NDISC /88/, TSTART /0./, TEND /24./
c
         DATA IMID,JMID / 7 , 4 /
         DATA JMID2Z, JMID2V  / 7 , 8 /
c
         DATA YEAR /1910/, MONTH /05/, DAY /20/, HOUR / 07/
c
         DATA DT / 600./
         DATA EpsTime / 0./, EpSpace / 1./
c
CCC      DATA WHOUR / 00 , 0.25, 06, 12, 24, 48, 4*0 /
         DATA WHOUR / 00 ,   01, 06, 12, 24, 48, 4*0 /
         DATA WNAME /'FCST00','FCST01','FCST06',
     X               'FCST12','FCST24','FCST48',
     X              4*' ' /
c
C***********************************************************************
c
c       OPEN CONTROL CARD FILE and READ PARAMETERS.
C		*cds file contains parameters for simulation like start-time, end-time etc.
        OPEN(UNIT=5,FILE='richie.cds',FORM='FORMATTED')
        READ(5,*) string,  YEAR, MONTH, DAY, HOUR
        READ(5,*) string,  OLDDAT , OLDFIL
        READ(5,*) string,  LFRDAT
        READ(5,*) string,  DT 
        READ(5,*) string,  TSTART 
        READ(5,*) string,  TEND 
        READ(5,*) string,  TFORWARD
        READ(5,*) string,  TIMEFILT,EpsTime
        READ(5,*) string,  MASSFILT, WINDFILT, EpSpace,EPSSDT
        READ(5,*) string,  DIVDAMP, CDDAMP 
        READ(5,*) string,  SURFAC
        READ(5,*) string,  INIT
        READ(5,*) string,  TAUC 
        READ(5,*) string,  IWINDOW , WPARAM
        READ(5,*) string,  ICGEOS
C
        print*, ' Date:   ',  YEAR, MONTH, DAY, HOUR
        print*, ' OLDdat  ',  OLDDAT, OLDFIL
        print*, ' LFRdat  ',  LFRDAT
        print*, ' Delta t ',  DT 
        print*, ' tStart  ',  TSTART
        print*, ' tEnd    ',  TEND   
        print*, ' tForward',  TFORWARD   
        print*, ' EpsTime ',  TIMEFILT, EpsTime
        print*, ' EpSpace ',  MASSFILT, WINDFILT, EpSpace,EPSSDT
        print*, ' DivDamp ',  DIVDAMP,  CDDAMP
        print*, ' INIT    ',  INIT
        print*, ' TauC    ',  TAUC 
        print*, ' Window  ',  IWINDOW , WPARAM
        print*, ' ICGEOS  ',  ICGEOS
c
C***********************************************************************
c
C-----------------------------------------------------------------------
C
C       SET THE AREA AND OTHER PARAMETERS.
c
C       Specify the grid for plotting.
C 		*DLAMD is the longitude and DPHID the latitude of the grid
        WBLMD = -7.0
        PHISB = +37.8
        DLAMD = 3.00
        DPHID = 1.50
C
        GLM0D = 0.00
        ETH0D = 0.00
        IMCOPY = 13
        JMCOPY = 14
        LMCOPY = 1
C
C***********************************************************************
c
      IP=IM
      JP=JM
      LP=LM
c
      print  2200, IP,JP,LP,DLAMD,DPHID,WBLMD,PHISB
 2200 FORMAT(' IM=',I3,' JM=',I2,' LM=',I2/
     +       ' DLAMD=',F4.2,' DPHID=',F4.2,
     +       ' WBLMD=',F3.0,' PHISB=',F3.0)
C
      print  2240, NDISC,TSTART,TEND
 2240 FORMAT(' NDISC=',I2,' TSTART=',F4.0,' TEND=',F4.0)
c
C--------------RUN CONTROL CONSTANTS------------------------------------
C
      NTSPH=3600./ABS(DT)
C
      NSTART=ABS(TSTART)*NTSPH   +.5
      NSTEPS=ABS(TEND  )*NTSPH   +.5
c
      IF (ABS(DT).gt.3600.) THEN
         XNTSPH=3600./ABS(DT)
         NSTEPS=NINT(ABS(TEND)*XNTSPH)
      ENDIF
C
      print  2250, DT, NSTART, NSTEPS 
 2250 FORMAT(' DELTAT=',F8.0,'   NSTART=',I6,'   NSTEPS=',I6)
c
      IF ( NSTEPS.GT.NSMAX ) STOP ' Increase NSMAX '
C
C--------------DERIVED SPACE-TIME GRID AND OTHER CONSTANTS--------------
C
C 	  *actual grid parameters used in calculations are in radiant
      DRAD=PI/180.
      DLAM=DLAMD*DRAD
      DPHI=DPHID*DRAD
C
      DO 102 J=1,JM
      DO 102 K=1,2
        PHI=(PHISB+(2*J+K-3)*DPHID)*DRAD
        COSPHI(J,K)=COS(PHI)
        H1 (J,K)=A*COSPHI(J,K)
        RH1(J,K)=1./H1(J,K)
        FCOR(J,K)=OMEGA2*SIN(PHI)
 102  CONTINUE
C
      DO 602 J=1,JM2
        PHI=(PHISB+(J-1)*DPHID)*DRAD
        COSPH2(J)=COS(PHI)
        FCOR2(J)=OMEGA2*SIN(PHI)
        TAN2(J)=SIN(PHI)/COS(PHI)
 602  CONTINUE
C
C=======================================================================
c
C 		*ZLEV gives elevation in meters, ZBAR in km
        DO L=1,LM
          ZLEV(L) = ZBAR(L) * 1000.
        ENDDO
C		
C		*calculate deltaZ for all levels	
        DO L=2,LM
          DELZ(L) = ZLEV(L-1) - ZLEV(L)
        ENDDO
C
C-----------------------------------------------------------------------
c
        KAPPA = RGAS/CP
        GAMMA = 1/(1-KAPPA)
c       print*, ' KAPPA, GAMMA ', KAPPA, GAMMA
c
C=======================================================================
c
C     Calculate the Filter Weights for DFI.
C 	  *note: in richie.cds INIT is set to false
      IF ( INIT ) THEN 
C
        IF      (IWINDOW.eq.1) THEN
          CALL HFILT1 (NSTEPS,DT,TAUC,IWINDOW,WPARAM,HH)
          print*,' HFILT1: ', (HH(NNN),NNN=0,NSTEPS)
        ELSE IF (IWINDOW.eq.-1) THEN
          NCHOICE = 1
          TAUS = TAUC
          NSH = NSTEPS/2
          CALL  DOLPH(ABS(DT), NCHOICE, TauS, ripple, NSTEPS, HH2)
          DO NNN=0,NSTEPS
            HH(NNN) = HH2(NNN+NSTEPS)
          ENDDO
          print*,' DOLPH : ', (HH(NNN),NNN=0,NSTEPS)
        ELSE 
          STOP ' IWINDOW not recognized '
        ENDIF
C
      ENDIF
c
C=======================================================================
C
C       READ IN THE INITIAL DATA.
CCC     OPEN(UNIT=NDISC,FILE='RICHIE.DAT',FORM='FORMATTED')
CCC     print*, ' READ in the INITIAL fields '
CCC     READ(NDISC,777) P,MU,MV,TTOP,OROG
c 777   format( 5E20.8 )
C
C***********************************************************************
c
c       get the initial fields
C		*note: in richie.cds OLDDAT is set to false
        IF (OLDDAT ) THEN
c
C           READ IN THE INITIAL DATA FOR RE_RUNS.
            OPEN(UNIT=60,FILE=OLDFIL,FORM='UNFORMATTED')
            print*, ' Read in the INITIAL fields (B1)'
            READ (60) PSEA2,Z2,U2,V2,TTOP2,OROG
            print*, ' Read in the INITIAL fields (B2)'
            READ (60) P2,MU2,MV2
            CLOSE(UNIT=60)
            GO TO 1000
c        
        ELSE
c
            print*,' @ @    Entry to GETDAT  @ @  '
            CALL GETDAT
            print*,' @@@@  Back from GETDAT  @@@@ '
c
c           make alterations in data to agree with
c           Richardson's values (Table, LFR, p 185).
            CALL TABLOT(P,MU,MV,TTOP,OROG,IM,JM,LM,'   BEFOR   ')
            IF ( LFRDAT ) THEN
               CALL RICDAT
               CALL TABLOT(P,MU,MV,TTOP,OROG,IM,JM,LM,'   AFTER   ')
            ENDIF
c
        ENDIF
c
C***********************************************************************
c
C-----------   COMPUTATION OF DIVERGENCE on E-GRID   -------------------
C-----------   (Old code, for checking purposes only)-------------------
c
      DO 106 J=1,JM
      DO 106 I=1,IM
         WORK1(I,J)=0.
 106  CONTINUE
 
      DO 402 L=1,LM
c
      ncount = 0
      sumdiv = 0.
      absdiv = 0.
      DO 401 I=3,IMM2
      MI2=MOD(I,2)
      KV=1+MI2
      KH=2-MI2
      JF=JMM2+MI2
      DO 401 J=2,JF
      JM1=J-MI2
      JP1=JM1+1
      DIV(I,J,L)= (MU(I+1,J,L)-MU(I-1,J,L))/(H1(J,KH)*2*DLAM)
     +           +( MV(I,JP1,L)*COSPHI(JP1,KV)
     +             -MV(I,JM1,L)*COSPHI(JM1,KV) )/(H1(J,KH)*2*DPHI)
            ncount = ncount + 1
            sumdiv = sumdiv + DIV(I,J,L)
            absdiv = absdiv + ABS(DIV(I,J,L))
  401 CONTINUE
      sumdiv = sumdiv / ncount
      absdiv = absdiv / ncount
c
        DMID = DIV(IMID,JMID,L)
        print9441, L, DMID, sumdiv, absdiv, ncount
 9441   FORMAT(' Div(7,4), lev:',I2,1PE12.3,
     +         ' sumdiv,absdiv :', 1P2E12.3,    i3)
c
  402 CONTINUE
      DSAVE = DIV(IMID,JMID,4)
      print*, ' DSAVE: ', DSAVE
C
C---------VERTICAL INTEGRATION TO GET SURF. PRESS. TENDENCY---------
c
      DO 202 L=1,LM
      DO 202 I=3,IMM2
      JF=JMM2+MOD(I,2)
      DO 202 J=2,JF
      WORK1(I,J)=WORK1(I,J)-DIV(I,J,L)*GRAV
 202  CONTINUE
c
      DPSDT = WORK1(IMID,JMID)
      print9442, DPSDT
 9442 FORMAT(' Pressure Tendency (7,4)', 1PE12.3,' Pa/s')
      DPSDT2 = DPSDT*6*3600/100
      print9443, DPSDT2
 9443 FORMAT(' Pressure Tendency (7,4)', F12.1,' mb/6h')
C
C-----------------------------------------------------------------------
C
c       Interpolate PROGNOSTIC Variables from E-grid to A-grid.
        DO L=1,LM
          CALL FILLZ(P (1,1,L),P2 (1,1,L),IM,JM)
          CALL FILLV(MU(1,1,L),MU2(1,1,L),IM,JM)
          CALL FILLV(MV(1,1,L),MV2(1,1,L),IM,JM)
        ENDDO
        CALL FILLZ(TTOP,TTOP2,IM,JM)
C
C***********************************************************************
C
C     Define geostrophic initial momenta if desired
      IF (ICGEOS) THEN

      DO I=1,IM
      DO J=1,JM2 
         PMID2(I,J,1) = P2(I,J,1) * EXP(-1.) 
         DO L=2,LMM1
           PMID2(I,J,L)  = ( P2(I,J,L-1) + P2(I,J,L) ) / 2.
         ENDDO
         PMID2(I,J,LM) = ( P2(I,J,L-1) + 100000. ) / 2.
         IF(I.EQ.IMID .and. J.EQ.JMID2Z )
     +   print*, ' PMID0 (7,7) ', (PMID2(I,J,L),l=1,LM)
      ENDDO
      ENDDO
c
         FCOR0 = FCOR2(JMID2Z)
         DO L=1,LM
         DO I=2,IMM1
         DO J=1,JM2
           DDXP=(PMID2(I+1,J,L)-PMID2(I-1,J,L))/(A*COSPH2(J)*2*DLAM)
           MV2(I,J,L) = (20000/grav) * DDXP/FCOR0
           IF(I.EQ.IMID .and. J.EQ.JMID2V) THEN
              print*, ' MV2GEOS ', MV2(I,J,L)
           ENDIF
         ENDDO
         ENDDO
         ENDDO
c
         DO L=1,LM
         DO I=1,IM
         DO J=2,JM2M1
           DDYP = ( PMID2(I,J+1,L)-PMID2(I,J-1,L) ) / (A*2*DPHI)
           MU2(I,J,L) = - (20000/grav) * DDYP/FCOR0
           IF(I.EQ.IMID .and. J.EQ.JMID2V) THEN
              print*, ' MU2GEOS ', MU2(I,J,L)
           ENDIF
         ENDDO
         ENDDO
         ENDDO

C        CALL MMTOUV(P2,MU2,MV2, U2,V2,IM,JM2,LM)
      ELSE
C  HACK of indices: set LM to LM -1
         DO L=1,LM -1
            print*, ' MV2ORIG ', MV2(IMID,JMID2V,L)
         ENDDO
         print*, LM, IMID, JMIDV, shape(MU2)
         DO L=1,LM -1
            print*, ' MU2ORIG ', MU2(IMID,JMIDV,L)
         ENDDO

       ENDIF
C
C***********************************************************************
C
 1000   CONTINUE
c
c       Copy the initial values into the HISTORICAL fields. 
        DO J=1,JM2
        DO I=1,IM
          DO L=1,LM
            P2OLD (I,J,L) = P2 (I,J,L)
            MU2OLD(I,J,L) = MU2(I,J,L)
            MV2OLD(I,J,L) = MV2(I,J,L)
            P2NEW (I,J,L) = P2 (I,J,L)
            MU2NEW(I,J,L) = MU2(I,J,L)
            MV2NEW(I,J,L) = MV2(I,J,L)
          ENDDO
        TT2OLD(I,J) = TTOP2(I,J)
        TT2NEW(I,J) = TTOP2(I,J)
        ENDDO
        ENDDO
c
         NSTEP = 0
         print9091, 
     X    ( NSTEP, L, P2NEW (IMID,JMID2Z,L), MU2NEW(IMID,JMID2V,L),
     X                MV2NEW(IMID,JMID2V,L) ,  L=1,LM )
 9091    FORMAT(' NsteP=',I6,'   LeveL=',I1,'   P2/MU2/MV2/ ',
     X                 F12.1,2x,2f9.1) 
         print9092, NSTEP, TT2NEW(IMID,JMID2Z)
 9092    FORMAT(' NsteP=',I6,'   TTOP2/',F12.1) 
C
c     Put the initial values into the INITIALIZED arrays. 
C     (First pass is with DT < 0; second is with DT > 0).
      IF ( INIT ) THEN
         IF ( DT . GT. 0 ) THEN
C           READ IN THE semi-INITIALIZED DATA
            OPEN(UNIT=55,FILE='FINI1',FORM='UNFORMATTED')
            print*, ' Read in the FINI1 fields '
            READ (55) P2INI,MU2INI,MV2INI, TT2INI
            CLOSE(UNIT=55)
            DO J=1,JM2
            DO I=1,IM
             DO L=1,LM
               P2INI (I,J,L) = P2INI (I,J,L)+0.5*hh(nstep)*P2 (I,J,L)
               MU2INI(I,J,L) = MU2INI(I,J,L)+0.5*hh(nstep)*MU2(I,J,L)
               MV2INI(I,J,L) = MV2INI(I,J,L)+0.5*hh(nstep)*MV2(I,J,L)
             ENDDO
             TT2INI(I,J) = TT2INI(I,J)+0.5*hh(nstep)*TTOP2(I,J)
            ENDDO
            ENDDO
         ELSE
            DO J=1,JM2
            DO I=1,IM
               DO L=1,LM
                  P2INI (I,J,L) = 0.5*hh(nstep)*P2 (I,J,L)
                  MU2INI(I,J,L) = 0.5*hh(nstep)*MU2(I,J,L)
                  MV2INI(I,J,L) = 0.5*hh(nstep)*MV2(I,J,L)
               ENDDO
               TT2INI(I,J) = 0.5*hh(nstep)*TTOP2(I,J)
            ENDDO
            ENDDO
         ENDIF
      ENDIF
c
c--------------------------------------------------------------
c
C       Calculate the total Mass or areal pressure integral.
        ZAREA = 0.
        ZMASS = 0.
        DO J=1,JM2
        DO I=1,IM
          IF(((I+J)/2*2).EQ.(I+J)) THEN 
            DAREA = 2*(A*COSPH2(J)*DLAM)*(A*DPHI)
            ZAREA = ZAREA + DAREA
            DMASS = P2 (I,J,LM)*DAREA
            ZMASS = ZMASS + DMASS
          ENDIF
        ENDDO
        ENDDO
        print9095, NSTEP, ZAREA, ZMASS, ZMASS/ZAREA 
 9095   FORMAT(' NsteP=',I6,' Area, MASS, pMEAN:',1P2E12.3,-2PF12.1) 
c
C-----------------------------------------------------------------------
c
C      Convert fields to PSEA and Z and U and V for OUTPUT.
       CALL PTOZ(P2,TTOP2,OROG, Z2,PSEA2,IM,JM2,LM)
       CALL MMTOUV(P2,MU2,MV2, U2,V2,IM,JM2,LM)
C
C
C       Write out the Initial Data.
        OPEN(UNIT=50,FILE=WNAME(1),FORM='UNFORMATTED')
        print*, ' Write out the INITIAL fields '
        WRITE(50) PSEA2,Z2,U2,V2,TTOP2,OROG
        WRITE(50) P2,MU2,MV2
        CLOSE(UNIT=50)
C
C***********************************************************************
C
C       Calculate Energy and other diagnostic quantities.
        NTSTEP = 0
C       CALL ENERGY(NTSTEP, FI,U,V,FI0,IM,JM, DX,DY,
C    X                    EK,EA,EP,AK,PK,EM, .TRUE.)
C
  		XVAR(NTSTEP)=NTSTEP
  		YVAR(NTSTEP)= P2NEW (IMID,JMID2Z,LM)
		UVAR(NTSTEP)=USAVE
		VVAR(NTSTEP)=VSAVE
 		DVAR(NTSTEP)=0.
                KVAR(NTSTEP)=EK
                AVAR(NTSTEP)=EA
		PVAR(NTSTEP)=EP
                MVAR(NTSTEP) =  ZMASS/ZAREA
		N1VAR(NTSTEP)=  0.
		N2VAR(NTSTEP)=  0.
                BRVAR(NTSTEP) = 0
		N1PVAR(NTSTEP)=  N1MID
		N2PVAR(NTSTEP)=  N2MID
                BRPVAR(NTSTEP) = BRMID
C
C***********************************************************************
C
C                            MAIN LOOP
c
C-----------------------------------------------------------------------
c
      DO 5000 NSTEP=1,NSTEPS
c
        print*,' >>>>> Start of Main Loop, NSTEP = ', NSTEP
c
c       First step is forward, others are centered.
        IF ( NSTEP.eq.1 ) THEN
          DELTAT = DT
        ELSE
          DELTAT = 2.*DT
        ENDIF
c
        TIME = NSTEP*DT
        HR= TIME/(3600.)
        print*,' Time Step ending at', HR, '   Hours '
c
c       Take a forward step now and then.
        ITF = INT(TFORWARD*3600.)
        ITIME = INT(TIME)
        If ( ITIME/ITF*ITF .EQ. ITIME ) THEN
           print*,' Forward Step at ', HR,'   Hours '
           DELTAT = DT
c          Copy the MIDDLE values into the HISTORICAL fields. 
           DO L=1,LM
             CALL SMEARZ(P2 (1,1,L),P2OLD (1,1,L),IM,JM)
             CALL SMEARV(MU2(1,1,L),MU2OLD(1,1,L),IM,JM)
             CALL SMEARV(MV2(1,1,L),MV2OLD(1,1,L),IM,JM)
           ENDDO
           CALL SMEARZ(TTOP2,TT2OLD,IM,JM)
        ENDIF
C
C************************************************************
C*********************    ( i )  ****************************
C************************************************************
c
C------  CENTRAL PRESSURE (LAYER INTEGRAL OF P) -------------
C------  LAYER INTEGRAL OF DENSITY ( RR ) -------------------
C
      DO 1011 I=1,IM
      DO 1011 J=1,JM2 
         PP2(I,J,1) = ( RGAS*TTOP2(I,J)/GRAV ) * P2(I,J,1)
         DO 1010 L=2,LM
         IF (L.eq.LM) THEN
            DELTAZ = ZLEV(L-1)-OROG(I,J)
         ELSE
            DELTAZ = DELZ(L)
         ENDIF
         PP2(I,J,L) = DELTAZ*(P2(I,J,L)-P2(I,J,L-1))
     +              / ( LOG(P2(I,J,L)) - LOG(P2(I,J,L-1)) )
 1010    CONTINUE
         IF(I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 )
     +   print*, ' PP2(7,7) ', (PP2(I,J,L),l=1,LM)
 1011 CONTINUE
c
C
      DO 1015 L=1,LM
      DO 1014 I=1,IM
      DO 1014 J=1,JM2
        IF ( L.eq.1 ) THEN
           RR2(I,J,L) = ( P2(I,J,1)            ) / GRAV
        ELSE
           RR2(I,J,L) = ( P2(I,J,L)-P2(I,J,L-1) ) / GRAV
        ENDIF
           IF(I.eq.IMID .and. J.eq.JMID2Z .and. NSTEP.EQ.1) THEN
              RMID=RR2(I,J,L)
              print*,' CHECK: RR2: ', L,RMID
           ENDIF
 1014 CONTINUE
 1015 CONTINUE
c
C************************************************************
C*********************    ( ii )  ***************************
C************************************************************
c
C------  CENTRAL VALUE OF PRESSURE ( PMID ) -----------------
C------  CENTRAL VALUE OF DENSITY ( RHO2 ) ------------------
c
      DO 1017 I=3,IMM2
      DO 1017 J=3,JM2M3 
CCC      PMID2(I,J,1) = P2(I,J,1) + ( RGAS*TTOP2(I,J)/GRAV )
         PMID2(I,J,1) = P2(I,J,1) * EXP (-1.) 
         DO 1016 L=2,LM
            IF (L.eq.LM) THEN
               DELTAZ = ZLEV(L-1)-OROG(I,J)
            ELSE
               DELTAZ = DELZ(L)
            ENDIF
            PMID2(I,J,L) = PP2(I,J,L)/DELTAZ
 1016    CONTINUE
         IF(I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 )
     +   print*, ' PMID2(7,7) ', (PMID2(I,J,L),l=1,LM)
 1017 CONTINUE
c
      DO 1019 I=3,IMM2
      DO 1019 J=3,JM2M3 
         DO 1018 L=2,LM
            IF (L.eq.LM) THEN
               DELTAZ = ZLEV(L-1)-OROG(I,J)
            ELSE
               DELTAZ = DELZ(L)
            ENDIF
            RHO = RR2(I,J,L) / DELTAZ
            TEMP = PMID2(I,J,L) / ( RGAS*RHO )
            RHO2(I,J,L) = RHO 
            T2(I,J,L) = TEMP 
 1018    CONTINUE
      IF ( I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 ) THEN
         print*, ' RHO   ', (RHO2(I,J,L),L=2,LM)
         print*, ' TEMP  ', (T2 (I,J,L),L=2,LM)
      ENDIF
 1019 CONTINUE
c
C************************************************************
C*********************    (iii)  ****************************
C************************************************************
c
C--------------COMPUTATION OF DIVERGENCE (of Momentum) -----------------
c
      DO 1004 L=1,LM
c
      ncount = 0
      sumdiv = 0.
      absdiv = 0.
CCC   DO 1002 I=3,IMM2
CCC   DO 1002 J=3,JM2M3
                            DO 1002 I=2,IMM1
                            DO 1002 J=2,JM2M2
      DIV2M(I,J,L)= (MU2(I+1,J,L)-MU2(I-1,J,L))/(A*COSPH2(J)*2*DLAM)
     +            +(MV2(I,J+1,L)*COSPH2(J+1)
     +             -MV2(I,J-1,L)*COSPH2(J-1) )/(A*COSPH2(J)*2*DPHI)
      if(((I+J)/2*2).EQ.(I+J)) ncount = ncount + 1
      if(((I+J)/2*2).EQ.(I+J)) sumdiv = sumdiv + DIV2M(I,J,L)
      if(((I+J)/2*2).EQ.(I+J)) absdiv = absdiv + ABS(DIV2M(I,J,L))
 1002 CONTINUE
      sumdiv = sumdiv / ncount
      absdiv = absdiv / ncount
c

        DMID = DIV2M(IMID,JMID2Z,L)
        print9444, L, DMID, sumdiv, absdiv, ncount
 9444   FORMAT(' DivM(7,7),lev:',I2,1PE12.3,
     +         ' sumdiv,absdiv :', 1P2E12.3,  i3)
c
 1004 CONTINUE
      DSAVE = DIV2M(IMID,JMID2Z,4)
      print*, ' DSAVE: ', DSAVE
C
C----   Divergence Analysis Check.
      DIVXAS = 0
      DIVYAS = 0
      DIVAS = 0
      IF ( NSTEP.EQ.1 ) THEN
        DO 10049 L=1,LM
        I=IMID
        J=JMID2Z
        DIVX2M =
     +   (MU2(I+1,J,L)-MU2(I-1,J,L))/(A*COSPH2(J)*2*DLAM)
        DIVY2M =
     +   (MV2(I,J+1,L)*COSPH2(J+1)-MV2(I,J-1,L)*COSPH2(J-1))
     +    /  (A*COSPH2(J)*2*DPHI)
        DMIDXY = DIVX2M + DIVY2M 
        print 94449, L, DIVX2M, DIVY2M, DMIDXY
94449   FORMAT(' Div. Anal., Lev:',I2,' DIVX2M DIVY2M',
     +           2F10.4,' DMIDXY ', F10.4)
        DIVXAS = DIVXAS + ABS(DIVX2M)
        DIVYAS = DIVYAS + ABS(DIVY2M)
        DIVAS = DIVAS + ABS(DMIDXY)
10049   CONTINUE
        print 94459, DIVXAS/LM, DIVYAS/LM, DIVAS/LM
94459   FORMAT(' Div. Mean Abs Vals DIVXAS DIVYAS DIVAS',
     +           3F10.4)
      ENDIF
C
C************************************************************
C*********************    (iv)  *****************************
C************************************************************
c
C---------VERTICAL INTEGRATION TO GET SURF. PRESS. TENDENCY---------
c
      DO 1001 J=1,JM2
      DO 1001 I=1,IM
         WORKK2(I,J)=0.
 1001 CONTINUE
c
      DO 1006 L=1,LM
        DO 1006 I=3,IMM2
        DO 1006 J=3,JM2M3
          WORKK2(I,J)=WORKK2(I,J)-DIV2M(I,J,L)*GRAV
 1006  CONTINUE
c
      DPSDT = WORKK2(IMID,JMID2Z)
      print9445, DPSDT
 9445 FORMAT(' Pressure Tendency (7,7)', 1PE12.3,' Pa/s')
      DPSDT2 = DPSDT*6*3600/100
      print9446, DPSDT2
 9446 FORMAT(' Pressure Tendency (7,7)', F12.1,' mb/6h')
C
C--------
C
C     GET the NOISE Parameter based on mean modulus of dps/dt.
      ncount = 0
      PNOIS1 = 0.
      PNOIS2 = 0.
      DO I=3,IMM2
      DO J=3,JM2M3
        IF ( ((I+J)/2*2).EQ.(I+J) ) THEN 
          ncount = ncount + 1
          PNOIS1 = PNOIS1 + ABS(WORKK2(I,J))
        ENDIF
        IF(I.eq.IMID .and. J.eq.JMID2Z) THEN
          PNOIS2 = PNOIS2 + WORKK2(I,J)
        ENDIF
      ENDDO
      ENDDO
      PNOIS1 = PNOIS1 / ncount
      PNOIS2 = ABS(PNOIS2)
c
      print94581, PNOIS1, NSTEP
94581 FORMAT(' PNOIS1 Param, Mean dps/dt',1PE12.3,' NSTEP:',I4)
      PNOISE = PNOIS1*6*3600/100
      print94582, PNOISE
94582 FORMAT(' PNOISE:  Mean abs(dps/dt)',1PE12.3,' hPa/6h ')
      print94583, PNOIS2*6*3600/100
94583 FORMAT(' PNOIS2: dps/dt at mid-pt ',1PE12.3,' hPa/6h ')
C
C************************************************************
C*********************    ( v )  ****************************
C************************************************************
c
C--------------COMPUTATION OF DIVERGENCE (of Velocity) -----------------
c
      DO 1020 J=1,JM2
      DO 1020 I=1,IM
         WORKK2(I,J)=0.
 1020 CONTINUE
C
C     First, get the Velocities.
      DO 1022 L=1,LM
c
      DO 1021 I=1,IM
      DO 1021 J=1,JM2
        RR = RR2(I,J,L)
        U2(I,J,L) = MU2(I,J,L)/RR
        V2(I,J,L) = MV2(I,J,L)/RR
 1021 CONTINUE
CCC     print9455, L, RMID
 9455   FORMAT(' RR(7,7), lev:',I2,1PE12.3)
c
 1022 CONTINUE
c
C     Now, get the Divergences
      DO 1024 L=1,LM
c
      ncount = 0
      sumdiv = 0.
      absdiv = 0.
      DO 1023 I=3,IMM2
      DO 1023 J=3,JM2M3
      DIV2V(I,J,L)= (U2(I+1,J,L)- U2(I-1,J,L))/(A*COSPH2(J)*2*DLAM)
     +             +(V2(I,J+1,L)*COSPH2(J+1)
     +              -V2(I,J-1,L)*COSPH2(J-1) )/(A*COSPH2(J)*2*DPHI)
C = = = Approximation to obtain closer agreement with LFR = = = =
          IF(NSTEP.EQ.1) DIV2V(I,J,L)= DIV2M(I,J,L)/RR2(I,J,L)
C = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      if(((I+J)/2*2).EQ.(I+J)) ncount = ncount + 1
      if(((I+J)/2*2).EQ.(I+J)) sumdiv = sumdiv + DIV2V(I,J,L)
      if(((I+J)/2*2).EQ.(I+J)) absdiv = absdiv + ABS(DIV2V(I,J,L))
c
         IF ( I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 ) THEN
            VERT1(L) = U2(I-1,J,L)
            VERT2(L) = U2(I+1,J,L)
            VERT3(L) = V2(I,J-1,L)
            VERT4(L) = V2(I,J+1,L)
         ENDIF
c
 1023 CONTINUE
      sumdiv = sumdiv / ncount
      absdiv = absdiv / ncount
            DMID = DIV2V(IMID,JMID2Z,L)
            IF(NSTEP.EQ.1)
     X      print9454, L, DMID, sumdiv, absdiv, ncount
 9454       FORMAT(' DivV(7,7),lev:',I2,1PE12.3,
     +             ' sumdiv,absdiv :', 1P2E12.3,  i3)
c
 1024 CONTINUE
c
C     print9453, (L,VERT1(L),VERT2(L),VERT3(L),VERT4(L),M)
 9453 FORMAT(' UUVV(7,7),lev:',I2,2F10.2,3x,2F10.2)
C
C     GET the NOISE Parameters of the Divergences
c
      ncount = 0
      NOISE1 = 0.
      NOISE2 = 0.
         N1MID = 0.
         N2MID = 0.
      DO 1028 I=3,IMM2
      DO 1028 J=3,JM2M3
        sum1 = 0
        sum2 = 0
         DO 1027 L=1,LM
            IF ( ((I+J)/2*2).EQ.(I+J) ) THEN 
               ncount = ncount+1
               sum1 = sum1 + DIV2M(I,J,L) 
               sum2 = sum2 + ABS( DIV2M(I,J,L) )
            ENDIF
            IF(I.eq.IMID .and. J.eq.JMID2Z) THEN
               N1MID = N1MID + DIV2M(I,J,L) 
               N2MID = N2MID + ABS( DIV2M(I,J,L) )
            ENDIF
 1027    CONTINUE
         NOISE1 = NOISE1 + ABS(sum1)
         NOISE2 = NOISE2 + sum2
 1028 CONTINUE
      NOISE1 = NOISE1 / ncount
      NOISE2 = NOISE2 / ncount
C     Convert to hPa per hour.
      NOISE1 = (3600/100) * GRAV * NOISE1 
      NOISE2 = (3600/100) * GRAV * NOISE2
      BRUS = ( NOISE1 / NOISE2 ) * 100.
      N1MID = (3600/100) * GRAV * ABS(N1MID)
      N2MID = (3600/100) * GRAV * N2MID
      BRMID = ( N1MID / N2MID ) * 100.
c
        print9458, NOISE1, NOISE2, BRUS, NSTEP
 9458   FORMAT(' Noise Params (hPa/h), N1 N2 Br : ',1P2E12.3,0PF8.1,
     X         '%   NSTEP:',I4)
        print9459, N1MID, N2MID, BRMID
 9459   FORMAT(' N1 N2 Br at mid-point  : ',1P2E12.3,0PF8.1,'%')
C
      USAVE = U2(IMID,JMID2V,4)
      VSAVE = V2(IMID,JMID2V,4)
c     
        PZERO = 10**5
        DPZERO = 10**5 / LM
        print*,'  PNOIS1/(p0*N1)  ', PNOIS1/(NOISE1*PZERO)
        print*,' DPSDT/(DP*N1MID) ', DPSDT /(N1MID*DPZERO)
        print*,'PNOIS2/(DP*N1MID) ', PNOIS2/(N1MID*DPZERO)
C
C************************************************************
C*********************    ( vi )  ***************************
C************************************************************
c
C--------------THE VERTICAL VELOCITY GRADIENT -------------------------
c
      DO 1032 I=3,IMM2
      DO 1032 J=3,JM2M3 
         PMID = PMID2(I,J,1)
         VERT1(1) = DIV2V(I,J,1) * PMID
         VHALF1   = DIV2V(I,J,1) * ( P2(I,J,1)-PMID )
         DWDZ(I,J,1) = - DIV2V(I,J,1) + VERT1(1)/(GAMMA*PMID)
         DO 1031 L=2,LM
           PMID = PMID2(I,J,L)
           VHALF2 = DIV2V(I,J,L  )*( PMID-P2(I,J,L-1) )
           VERT1(L) = VERT1(L-1) + VHALF1 + VHALF2
           VHALF1 = DIV2V(I,J,L  )*( P2(I,J,L)-PMID )
           DWDZ(I,J,L) = - DIV2V(I,J,L) + VERT1(L)/(GAMMA*PMID)
 1031    CONTINUE
         IF(I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 )
     +   print91031, (DWDZ(I,J,L),l=1,LM)
91031    FORMAT(' DWDZ:', 1P5E15.6)
 1032 CONTINUE
c
C     Correction for small integral term.
      DO 1034 I=3,IMM2
      DO 1034 J=3,JM2M3 
         DUDP = (U2(I,J,2)-U2(I,J,1))/(PMID2(I,J,2)-PMID2(I,J,1))
         DVDP = (V2(I,J,2)-V2(I,J,1))/(PMID2(I,J,2)-PMID2(I,J,1))
         DPDX = (P2(I+1,J,1)-P2(I-1,J,1))/(A*COSPH2(J)*2*DLAM)
         DPDY = (P2(I,J+1,1)-P2(I,J-1,1))/(A*          2*DPHI)
         DVGP = DUDP*DPDX + DVDP*DPDY
         PMID = PMID2(I,J,1)
         VSTEP = DVGP * PMID
         VERT1(1) =  VSTEP
         DWDZ(I,J,1) = DWDZ(I,J,1) - VERT1(1)/(GAMMA*PMID)
        DO 1033 L=2,LM
         DUDP = (U2(I,J,L)-U2(I,J,L-1))/(PMID2(I,J,L)-PMID2(I,J,L-1))
         DVDP = (V2(I,J,L)-V2(I,J,L-1))/(PMID2(I,J,L)-PMID2(I,J,L-1))
         DPDX = (P2(I+1,J,L-1)-P2(I-1,J,L-1))/(A*COSPH2(J)*2*DLAM)
         DPDY = (P2(I,J+1,L-1)-P2(I,J-1,L-1))/(A*          2*DPHI)
         DVGP = DUDP*DPDX + DVDP*DPDY
         DELP = PMID2(I,J,L)-PMID2(I,J,L-1)
           VSTEP = DVGP * DELP
           PMID = PMID2(I,J,L)
           VERT1(L) = VERT1(L-1) + VSTEP
           DWDZ(I,J,L) = DWDZ(I,J,L) - VERT1(L)/(GAMMA*PMID)
 1033    CONTINUE
         IF(I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 )
     +   print91033, (DWDZ(I,J,L),l=1,LM)
91033    FORMAT(' DwDz:', 1P5E15.6)
 1034 CONTINUE
c
C--------------THE VERTICAL VELOCITY ----------------------------------
c
C     Scale Factor for surface wind 
c     (as multiple of bottom layer wind)
CCC   SURFAC = 0.2  !!!   Now read in on a card.
c
      DO 1038 I=3,IMM2
      DO 1038 J=3,JM2M3 
c        First, get the value of w at the ground.
         DHDX = ( OROG(I+1,J)-OROG(I-1,J) ) / (A*COSPH2(J)*2*DLAM)
         DHDY = ( OROG(I,J+1)-OROG(I,J-1) ) / (A*2*DPHI)
         WLAYER = U2(I,J,LM)*DHDX + V2(I,J,LM)*DHDY
         W2(I,J,LM) = SURFAC * WLAYER
C        Now, work up from the bottom, using the gradient. 
         DO 1035 L=LM,2,-1
            IF (L.eq.LM) THEN
               DELTAZ = ZLEV(L-1)-OROG(I,J)
            ELSE
               DELTAZ = DELZ(L)
            ENDIF
            W2(I,J,L-1) = W2(I,J,L) + DELTAZ*DWDZ(I,J,L)
 1035    CONTINUE
         IF(I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 )
     +   print91035, (W2(I,J,L),l=1,LM)
91035    FORMAT('  W2 :', 1P5E15.6)
 1038 CONTINUE
C
C************************************************************
C*********************    ( vii )  **************************
C************************************************************
c
C--------------THE STRATOSPHERIC TEMPERATURE TENDENCY -----------------
c
      DO 1048 I=3,IMM2
      DO 1048 J=3,JM2M3 
         DTT2DT(I,J) = TTOP2(I,J) * DWDZ(I,J,1)
 1048 CONTINUE
c
      DTDT = DTT2DT(IMID,JMID2Z)
C     print9434, TTOP2(IMID,JMID2Z)
 9434 FORMAT(' Temperature (Strato) (7,7)', F12.1,' K')
      print9435, DTDT
 9435 FORMAT(' Temperature Tendency (7,7)', 1PE12.3,' K/s')
      DTDT2 = DTDT*6*3600
      print9436, DTDT2
 9436 FORMAT(' Temperature Tendency (7,7)', F12.1,' K/6h')
C
C************************************************************
C*********************    ( viii )  *************************
C************************************************************
C
C     Get temperature and, thence, density at the interfaces.
      DO 1052 I=3,IMM2
      DO 1052 J=3,JM2M3 
        DO 1051 l=1,LM
          IF ( L.eq.1 ) then
            T2I(I,J,L)  = ( TTOP2(I,J) + T2(I,J,2) ) / 2.
          ELSE IF ( L.eq.LM ) THEN
            T2I(I,J,L)  = ( 3.*T2(I,J,LM)-T2(I,J,LMM1) ) / 2.
          ELSE
            T2I(I,J,L)  = ( T2(I,J,L) + T2(I,J,L+1) ) / 2.
          ENDIF
          RHO2I(I,J,L) = P2(I,J,L) / ( RGAS*T2I(I,J,L) )
 1051   CONTINUE
              IF(I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 ) THEN
                print*, ' T2I  (7,7) ', (T2I  (I,J,L),l=1,LM)
                print*, ' RHO2I(7,7) ', (RHO2I(I,J,L),l=1,LM)
              ENDIF
 1052 CONTINUE
c
c     Get the vertical momentum at the interfaces
c     and the layer boundary term [mh].
      DO 1055 I=3,IMM2
      DO 1055 J=3,JM2M3 
        DO 1053 L=1,LM
          VERT1(L) = RHO2I(I,J,L)*W2(I,J,L)  
 1053   CONTINUE
              IF(I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 )
     +        print*, ' MH(7,7) ', (VERT1(L),l=1,LM)
        DO 1054 L=1,LM
          IF ( L.eq.1 ) then
            VERT2(L) =            - VERT1(L) 
          ELSE IF ( L.eq.LM ) THEN
            VERT2(L) = VERT1(L-1)
          ELSE
            VERT2(L) = VERT1(L-1) - VERT1(L) 
          ENDIF
          BRAKET(I,J,L) = VERT2(L)
 1054   CONTINUE
          IF(I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 ) THEN
            print 90155, (BRAKET(I,J,L),l=1,LM)
90155       FORMAT( ' BRAKET  ', 5F10.4) 
            BRAKAS = 0
            DO L=1,LM
              BRAKAS = BRAKAS + ABS(BRAKET(I,J,L))
            ENDDO
          ENDIF
 1055 CONTINUE
          print 90255, BRAKAS
90255     FORMAT( ' BRAKAS  ', F10.4) 
c
C************************************************************
C*********************    ( ix )  ***************************
C************************************************************
C
C--------------THE CONTINUITY EQUATION ---------------------------------
C
c     Get the tendency of layer 'density' RR: DRR2DT.
c
      SUMM1 = 0
      SUMM2 = 0
      SUMM3 = 0
      DO 1057 I=3,IMM2
      DO 1057 J=3,JM2M3 
        DO 1056 L=1,LM
          DRR2DT(I,J,L) = - ( DIV2M(I,J,L) + BRAKET(I,J,L) )
          DRR2DTX = DRR2DT(I,J,L)*(6*3600)
           IF(I.EQ.IMID .and. J.EQ.JMID2Z .and. NSTEP.EQ.1 ) THEN
              print*, ' DRR2DT ', L, DRR2DT(I,J,L), DRR2DTX
              VERT1(L) =   DRR2DT(I,J,L)*(6*3600)*(grav/100.)
              VERT2(L) = - DIV2M(I,J,L) *(6*3600)*(grav/100.)
              VERT3(L) = - BRAKET(I,J,L)*(6*3600)*(grav/100.)
              SUMM1 = SUMM1 + VERT1(L)
              SUMM2 = SUMM2 + VERT2(L)
              SUMM3 = SUMM3 + VERT3(L)
           ENDIF
 1056   CONTINUE
 1057 CONTINUE
        IF(NSTEP.EQ.1 ) THEN
           print90901, ( VERT1(L), VERT2(L),VERT3(L),L=1,LM )
90901      FORMAT( 'PPPPP: DDP, HORIZ, VERT ', 3F8.1)
           print90902,   SUMM1, SUMM2, SUMM3           
90902      FORMAT( ':::::::  SUMS        ', 3F8.1)
        ENDIF
C
C************************************************************
C*********************    ( x )  ****************************
C************************************************************
C
c     Get the tendency of interface pressures DP2DT.
c
      DO 1059 I=3,IMM2
      DO 1059 J=3,JM2M3 
        DP2DT(I,J,1) = GRAV*DRR2DT(I,J,1)
        DO 1058 L=2,LM
          DP2DT(I,J,L) = DP2DT(I,J,L-1) + GRAV*DRR2DT(I,J,L)
 1058   CONTINUE
 1059 CONTINUE
      IF(NSTEP.EQ.1)
     X    print9641, (L, (6*3600/100)*DP2DT(IMID,JMID2Z,L),L=1,LM)
 9641 FORMAT(' DP2DT, Lev',I2,F12.1,' mb/6h')
      DPSDT = DP2DT(IMID,JMID2Z,LM)
      print9642, DPSDT
 9642 FORMAT(' Pressure Tendency (7,7)', 1PE12.3,' Pa/s')
      DPSDT = DPSDT*6*3600/100
      print9643, DPSDT,NSTEP
 9643 FORMAT(' Pressure Tendency (7,7)', F12.1,' mb/6h ### NSTEP:',I4)
C
C-----------------------------------------------------------------------
c
c       Interpolate Fields from scalar to vector points
        DO L=1,LM
          CALL SMEARZ(RR2(1,1,L),RR2(1,1,L),IM,JM)
          CALL SMEARZ(RHO2I(1,1,L),RHO2I(1,1,L),IM,JM)
          CALL SMEARZ(W2(1,1,L),W2(1,1,L),IM,JM)
        ENDDO
C
C***********************************************************************
C***********************************************************************
C
C***********************************************************************
C***********************************************************************
c
C       MOMENTUM EQUATIONS.
C
C-----------------------------------------------------------------------
      print*, ' MMMMM:  Solution of Momentum Equations '
c
      DO 2000 L=1,LM
c
          IF(NSTEP.EQ.1) print9921,L
 9921     FORMAT('MMM: Momentum Tendency, Level:',I2)
c
C************************************************************
C*********************    ( xi )  ***************************
C************************************************************
C
C--------------THE PRESSURE GRADIENT ----------------------------------
c
      IF ( L.LT.LM ) THEN
c
        DO 1102 I=3,IMM2
        DO 1102 J=3,JM2M3 
          DDXP = ( PP2(I+1,J,L)-PP2(I-1,J,L) ) / (A*COSPH2(J)*2*DLAM)
          DDYP = ( PP2(I,J+1,L)-PP2(I,J-1,L) ) / (A*2*DPHI)
              IF(I.EQ.IMID .and. J.EQ.JMID2V .and. NSTEP.EQ.1) THEN
                 print*, ' DDXP, DDYP ', DDXP, DDYP
                 VERT1(L) = DDXP*6*3600
                 VERT4(L) = DDYP*6*3600
              ENDIF
         DMU2DT(I,J,L) = DDXP
         DMV2DT(I,J,L) = DDYP
 1102   CONTINUE
c
      ELSE
C
        DO 1104 I=3,IMM2
        DO 1104 J=3,JM2M3 
          TERM1X = (PP2(I+1,J,L)-PP2(I-1,J,L))
          TERM2X = (OROG(I+1,J)-OROG(I-1,J))*
     X             (P2(I+1,J,L)-P2(I-1,J,L))/
     X             (LOG(P2(I+1,J,L))-LOG(P2(I-1,J,L)))
          DDXP = ( TERM1X + TERM2X ) / (A*COSPH2(J)*2*DLAM)
          TERM1Y = (PP2(I,J+1,L)-PP2(I,J-1,L))
          TERM2Y = (OROG(I,J+1)-OROG(I,J-1))*
     X             (P2(I,J+1,L)-P2(I,J-1,L))/
     X             (LOG(P2(I,J+1,L))-LOG(P2(I,J-1,L)))
          DDYP = ( TERM1Y + TERM2Y ) / (A*2*DPHI)
              IF (I.EQ.IMID .and. J.EQ.JMID2V .and.  NSTEP.EQ.1 ) THEN
                 print*, ' DDXP, DDYP ', DDXP, DDYP
                 print*, ' TM1X, TM1Y ',
     X           TERM1X/(A*COSPH2(J)*2*DLAM),TERM1Y/(A*2*DPHI)
                 print*, ' TERM2X, TERM2Y ',
     X           TERM2X/(A*COSPH2(J)*2*DLAM),TERM2Y/(A*2*DPHI)
                 VERT1(L) = DDXP*6*3600
                 VERT4(L) = DDYP*6*3600
              ENDIF
         DMU2DT(I,J,L) = DDXP
         DMV2DT(I,J,L) = DDYP
 1104   CONTINUE
c
      ENDIF
C
C************************************************************
C*********************    ( xii )  **************************
C************************************************************
C
C--------------THE CORIOLIS FORCE (and Spherical Terms) ---------------
c
        DO 1105 I=3,IMM2
        DO 1105 J=3,JM2M3 
          CORX = - MV2(I,J,L) * FCOR2(J)
          CORY = + MU2(I,J,L) * FCOR2(J)
              IF(I.EQ.IMID .and. J.EQ.JMID2V .and. NSTEP.EQ.1) THEN
                 print*, ' CORX, CORY ', CORX, CORY
                 VERT2(L) = CORX*6*3600
                 VERT5(L) = CORY*6*3600
              ENDIF
          DMU2DT(I,J,L) =  DMU2DT(I,J,L) + CORX
          DMV2DT(I,J,L) =  DMV2DT(I,J,L) + CORY
c
          SPHERU = -2.*MU2(I,J,L)*MV2(I,J,L)*TAN2(J)/(A*RR2(I,J,L)) 
          SPHERV = (MU2(I,J,L)**2-MV2(I,J,L)**2)*TAN2(J)/(A*RR2(I,J,L)) 
              IF(I.EQ.IMID .and. J.EQ.JMID2V .and. NSTEP.EQ.1 )
     +        print*, ' SPHERIC    ', SPHERU, SPHERV
          DMU2DT(I,J,L) =  DMU2DT(I,J,L) + SPHERU
          DMV2DT(I,J,L) =  DMV2DT(I,J,L) + SPHERV
 1105   CONTINUE
c
C************************************************************
C*********************    ( xiii )  *************************
C************************************************************
C
C--------- Horizontal Advection or flux terms --------------------------
C
        DO 1106 I=3,IMM2
        DO 1106 J=3,JM2M3 
          XPLUS = ( MU2(I+2,J,L)**2 )/RR2(I+2,J,L)
          XMINS = ( MU2(I-2,J,L)**2 )/RR2(I-2,J,L)
          FLUXUX = ( XPLUS-XMINS ) / (A*COSPH2(J)*4*DLAM)
          YPLUS = ( MU2(I,J+2,L)*MV2(I,J+2,L) )/RR2(I,J+2,L)
          YMINS = ( MU2(I,J-2,L)*MV2(I,J-2,L) )/RR2(I,J-2,L)
          FLUXUY = ( YPLUS-YMINS ) / (A*4*DPHI)
c
          XPLUS = ( MU2(I+2,J,L)*MV2(I+2,J,L) )/RR2(I+2,J,L)
          XMINS = ( MU2(I-2,J,L)*MV2(I-2,J,L) )/RR2(I-2,J,L)
          FLUXVX = ( XPLUS-XMINS ) / (A*COSPH2(J)*4*DLAM)
          YPLUS = ( MV2(I,J+2,L)**2 )/RR2(I,J+2,L)
          YMINS = ( MV2(I,J-2,L)**2 )/RR2(I,J-2,L)
          FLUXVY = ( YPLUS-YMINS ) / (A*4*DPHI)
              IF(I.EQ.IMID .and. J.EQ.JMID2V .and. NSTEP.EQ.1 ) THEN
                 print*, ' FLUX(east)', FLUXUX, FLUXVX
                 print*, ' FLUX(nrth)', FLUXUY, FLUXVY
              ENDIF
         DMU2DT(I,J,L) = DMU2DT(I,J,L) + (FLUXUX+FLUXUY)
         DMV2DT(I,J,L) = DMV2DT(I,J,L) + (FLUXVX+FLUXVY)
 1106   CONTINUE
c
C************************************************************
C*********************    ( xiv )  **************************
C************************************************************
C
C--------- Vertical Advection or flux terms ----------------------------
c
C     (Note: Perhaps w Should be Smoothed from M to P points.)
      DO 1109 I=3,IMM2
      DO 1109 J=3,JM2M3 
          IF ( L.eq.1 ) then
            UHI = 0.
            VHI = 0.
            WHI = 0.
            ULO = ( U2(I,J,L)+U2(I,J,L+1) ) / 2.
            VLO = ( V2(I,J,L)+V2(I,J,L+1) ) / 2.
            WLO = W2(I,J,L)*RHO2I(I,J,L)  
ccccc             WLO = W2(I,J-1,L)*RHO2I(I,J-1,L)  
          ELSE IF ( L.eq.LM ) THEN
            UHI = ( U2(I,J,L)+U2(I,J,L-1) ) / 2.
            VHI = ( V2(I,J,L)+V2(I,J,L-1) ) / 2.
            WHI = W2(I,J,L-1)*RHO2I(I,J,L-1)  
ccccc             WHI = W2(I,J-1,L-1)*RHO2I(I,J-1,L)  
            ULO = 0.
            VLO = 0.
            WLO = 0.
          ELSE
            UHI = ( U2(I,J,L)+U2(I,J,L-1) ) / 2.
            VHI = ( V2(I,J,L)+V2(I,J,L-1) ) / 2.
            WHI = W2(I,J,L-1)*RHO2I(I,J,L-1)  
ccccc             WHI = W2(I,J-1,L-1)*RHO2I(I,J,L-1)  
            ULO = ( U2(I,J,L)+U2(I,J,L+1) ) / 2.
            VLO = ( V2(I,J,L)+V2(I,J,L+1) ) / 2.
            WLO = W2(I,J,L)*RHO2I(I,J,L)  
ccccc             WLO = W2(I,J-1,L)*RHO2I(I,J-1,L)  
          ENDIF
          VFLXU1 =   UHI*WHI 
          VFLXU2 = - ULO*WLO
          VFLXV1 =   VHI*WHI 
          VFLXV2 = - VLO*WLO
              IF(I.EQ.IMID .and. J.EQ.JMID2V .and. NSTEP.EQ.1 ) THEN
                 print*, ' WFLUX UPR ', VFLXU1, VFLXV1
                 print*, ' WFLUX LWR ', VFLXU2, VFLXV2
              ENDIF
         DMU2DT(I,J,L) = DMU2DT(I,J,L) + ( VFLXU1 + VFLXU2 )
         DMV2DT(I,J,L) = DMV2DT(I,J,L) + ( VFLXV1 + VFLXV2 ) 
 1109 CONTINUE
c
C-----------------------------------------------------------------------
c
c       Divergence Damping
        IF ( DIVDAMP ) THEN
           DO 1112 I=3,IMM2
           DO 1112 J=3,JM2M3 
             DDXDM=(DIV2M(I+1,J,L)-DIV2M(I-1,J,L))/(A*COSPH2(J)*2*DLAM)
             DDYDM=(DIV2M(I,J+1,L)-DIV2M(I,J-1,L))/(A*2*DPHI)
                 IF(I.EQ.IMID .and. J.EQ.JMID2V .and. NSTEP.EQ.1 )
     +           print*, ' DDXDM, DDYDM ', DDXDM, DDYDM
             DMU2DT(I,J,L) = DMU2DT(I,J,L) - CDDAMP*DDXDM
             DMV2DT(I,J,L) = DMV2DT(I,J,L) - CDDAMP*DDYDM 
 1112      CONTINUE
        ENDIF
c
C-----------------------------------------------------------------------
c
 2000 CONTINUE
C
C************************************************************
C*********************    ( xv )  ***************************
C************************************************************
C
C----------------- Tendency of Momenta ---------------------------------
c
        DO 1400 L=1,LM
c
        DO 1200 I=3,IMM2
        DO 1200 J=3,JM2M3 
         DMU2DT(I,J,L) = - DMU2DT(I,J,L)
         DMV2DT(I,J,L) = - DMV2DT(I,J,L)
 1200   CONTINUE
c
      DDTMU = DMU2DT(IMID,JMID2V,L)
      DDTMV = DMV2DT(IMID,JMID2V,L)
      DDTMU2 = DDTMU*6*3600
      DDTMV2 = DDTMV*6*3600
      print9742, DDTMU, DDTMV, DDTMU2, DDTMV2
 9742 FORMAT('Momentum Tendencies:(', 2F7.2,')kg/m/s; (',
     X                             -3P2F6.1,')T/m/6h')
      VERT3(L) = DDTMU2 + ( VERT1(L)+VERT2(L) )
      VERT6(L) = DDTMV2 + ( VERT4(L)+VERT5(L) )
c
        IF(NSTEP.EQ.1 ) THEN
           print90911, DDTMU2,-VERT1(L),-VERT2(L), VERT3(L)
           print90912, DDTMV2,-VERT4(L),-VERT5(L), VERT6(L)
90911      FORMAT( 'UUUUU:  DU, PGF, COR, RES ', -3P4F8.1)
90912      FORMAT( 'VVVVV:  DV, PGF, COR, RES ', -3P4F8.1)
        ENDIF
c
 1400   CONTINUE
C
C***********************************************************************
C***********************************************************************
C
C***********************************************************************
C***********************************************************************
c
C       CALCULATE THE FIELDS AT THE NEW TIME.
        DO 1600 I=3,IMM2
        DO 1600 J=3,JM2M3 
           DO L=1,LM
             P2NEW (I,J,L) = P2OLD (I,J,L) + DELTAT*DP2DT (I,J,L)
             MU2NEW(I,J,L) = MU2OLD(I,J,L) + DELTAT*DMU2DT(I,J,L)
             MV2NEW(I,J,L) = MV2OLD(I,J,L) + DELTAT*DMV2DT(I,J,L)
           ENDDO
           TT2NEW(I,J) = TT2OLD(I,J) + DELTAT*DTT2DT(I,J)
 1600   CONTINUE
c
C------- Filtering and Smoothing ---------------------------------------
c
c     Apply the Robert-Asselin time filter.
      IF ( TIMEFILT ) THEN
          ALPHA = EpsTime/2.
          DO 1700 I=3,IMM2
          DO 1700 J=3,JM2M3 
             DO L=1,LM
               P2 (I,J,L) = P2 (I,J,L) + 
     X           ALPHA*(P2OLD (I,J,L)-2*P2 (I,J,L)+P2NEW (I,J,L))
               MU2(I,J,L) = MU2(I,J,L) + 
     X           ALPHA*(MU2OLD(I,J,L)-2*MU2(I,J,L)+MU2NEW(I,J,L))
               MV2(I,J,L) = MV2(I,J,L) + 
     X           ALPHA*(MV2OLD(I,J,L)-2*MV2(I,J,L)+MV2NEW(I,J,L))
             ENDDO
             TTOP2(I,J) = TTOP2(I,J) + 
     X         ALPHA*(TT2OLD(I,J)-2*TTOP2(I,J)+TT2NEW(I,J))
 1700     CONTINUE
      ENDIF
c
C     Spatial Smoothing with Raymond Filter.
C     ISWEST = 2
C     ISEAST = IMM1
C     JSSOUTH = 2
C     JSNORTH = JM2M2
                ISWEST = 1
                ISEAST = IM
                JSSOUTH = 1
                JSNORTH = JM2M1
c
C     Allow for DAMP COEFF to depend on the time.
      IF ( EPSSDT.eq.0. ) THEN
         EPSILON = EpSpace 
      ELSE
         EPSILON = EpSpace * ((NSTEP-1)*DT)/(EPSSDT*3600.-DT)
      ENDIF
      print*, ' EPSILON  ', EPSILON
C
      IF ( MASSFILT ) THEN
c          First, average the fields. 
           DO L=1,LM
             CALL SMEARZ(P2NEW (1,1,L),P2NEW (1,1,L),IM,JM)
           ENDDO
           CALL SMEARZ(TT2NEW,TT2NEW,IM,JM)
C          Then, Apply the Implicit Filter (Surface pressure
C          requires special treatment).
           DO L=1,LMM1
             CALL IMPFIL(P2NEW (1,1,L),P2NEW (1,1,L),IM,JM2,IM,JM2,
     X            ISWEST,ISEAST,JSSOUTH,JSNORTH, EPSILON)
           ENDDO
           CALL IMPFIL(TT2NEW,TT2NEW,IM,JM2,IM,JM2,
     X            ISWEST,ISEAST,JSSOUTH,JSNORTH, EPSILON)
           CALL PSFILT(P2NEW(1,1,LM),OROG,PSEA2,IM,JM2,
     X            ISWEST,ISEAST,JSSOUTH,JSNORTH, EPSILON)
      ENDIF
c
      IF ( WINDFILT ) THEN
c          First, average the fields. 
           DO L=1,LM
             CALL SMEARV(MU2NEW(1,1,L),MU2NEW(1,1,L),IM,JM)
             CALL SMEARV(MV2NEW(1,1,L),MV2NEW(1,1,L),IM,JM)
           ENDDO
C          Then, Apply the Implicit Filter.
           DO L=1,LM
             CALL IMPFIL(MU2NEW(1,1,L),MU2NEW(1,1,L),IM,JM2,IM,JM2,
     X            ISWEST,ISEAST,JSSOUTH,JSNORTH, EPSILON)
             CALL IMPFIL(MV2NEW(1,1,L),MV2NEW(1,1,L),IM,JM2,IM,JM2,
     X            ISWEST,ISEAST,JSSOUTH,JSNORTH, EPSILON)
           ENDDO
      ENDIF
c
C-----------------------------------------------------------------------
c
c       print out a few of the new values.
         print9191, 
     X    ( NSTEP, L, P2NEW (IMID,JMID2Z,L), MU2NEW(IMID,JMID2V,L),
     X                MV2NEW(IMID,JMID2V,L) ,  L=1,LM )
 9191    FORMAT(' NsteP=',I6,'   LeveL=',I1,'   P2/MU2/MV2/ ',
     X                 F12.1,2x,2f9.1) 
         print9192, NSTEP, TT2NEW(IMID,JMID2Z)
 9192    FORMAT(' NsteP=',I6,'   TTOP2/',F12.1) 
C
C       Calculate the total Mass or areal pressure integral.
        ZAREA = 0.
        ZMASS = 0.
        DO J=1,JM2
        DO I=1,IM
            DAREA = (A*COSPH2(J)*DLAM)*(A*DPHI)
            ZAREA = ZAREA + DAREA
            DMASS = P2 (I,J,LM)*DAREA
            ZMASS = ZMASS + DMASS
        ENDDO
        ENDDO
        print9195, NSTEP, ZAREA, ZMASS, ZMASS/ZAREA 
 9195   FORMAT(' NsteP=',I6,' Area, MASS, pMEAN:',1P2E12.3,-2PF12.1) 
c
C-----------------------------------------------------------------------
C------ Prepare for the Next Time Step ---------------------------------
c
c       Copy the MIDDLE values into the HISTORICAL fields. 
        DO L=1,LM
          CALL SMEARZ(P2 (1,1,L),P2OLD (1,1,L),IM,JM)
          CALL SMEARV(MU2(1,1,L),MU2OLD(1,1,L),IM,JM)
          CALL SMEARV(MV2(1,1,L),MV2OLD(1,1,L),IM,JM)
        ENDDO
        CALL SMEARZ(TTOP2,TT2OLD,IM,JM)
C
c       Copy the NEW values into the MIDDLE fields. 
        DO L=1,LM
          CALL SMEARZ(P2NEW (1,1,L),P2 (1,1,L),IM,JM)
          CALL SMEARV(MU2NEW(1,1,L),MU2(1,1,L),IM,JM)
          CALL SMEARV(MV2NEW(1,1,L),MV2(1,1,L),IM,JM)
        ENDDO
        CALL SMEARZ(TT2NEW,TTOP2,IM,JM)
C
C-----------------------------------------------------------------------
c
C       Calculate Energy and other diagnostic quantities.
        NTSTEP = NSTEP
C       CALL ENERGY(NTSTEP, FI,U,V,FI0,IM,JM, DX,DY,
C    X                    EK,EA,EP,AK,PK,EM, .TRUE.)
C
  		XVAR(NTSTEP)=NTSTEP
  		YVAR(NTSTEP)= P2NEW (IMID,JMID2Z,LM)
		UVAR(NTSTEP)=USAVE
		VVAR(NTSTEP)=VSAVE
 		DVAR(NTSTEP)=0.
                KVAR(NTSTEP)=EK
                AVAR(NTSTEP)=EA
		PVAR(NTSTEP)=EP
                MVAR(NTSTEP) =  ZMASS/ZAREA
		N1VAR(NTSTEP-1)= NOISE1
		N2VAR(NTSTEP-1)= NOISE2
                BRVAR(NTSTEP-1)= BRUS
		N1PVAR(NTSTEP-1)=  N1MID
		N2PVAR(NTSTEP-1)=  N2MID
                BRPVAR(NTSTEP-1) = BRMID
C
c---------------------------------------------------------------
c
C      Write out Results at intermittent times.
c
      DO 3000 NWRITE=1,10
         TIME = NSTEP*DT
         HR= TIME/(3600.)
         IF ( ABS(HR).EQ.WHOUR(NWRITE) ) THEN
         print*,' Write-out at ', HR, '   Hours '
c
C           Convert fields to PSEA and Z and U and V for OUTPUT.
            CALL PTOZ(P2,TTOP2,OROG, Z2,PSEA2,IM,JM2,LM)
            CALL MMTOUV(P2,MU2,MV2, U2,V2,IM,JM2,LM)
C
            WFILE = WNAME(NWRITE)
C           Write out the Data.
            OPEN(UNIT=50,FILE=WFILE,FORM='UNFORMATTED')
            WRITE(50) PSEA2,Z2,U2,V2,TTOP2,OROG
            WRITE(50) P2,MU2,MV2
            CLOSE(UNIT=50)
         ENDIF
 3000 CONTINUE
C
C***********************************************************************
c
c     Update the values of the INITIALIZED arrays. 
      IF ( INIT ) THEN
        DO J=1,JM2
        DO I=1,IM
          DO L=1,LM
            P2INI (I,J,L) = P2INI (I,J,L) + hh(nstep)*P2 (I,J,L)
            MU2INI(I,J,L) = MU2INI(I,J,L) + hh(nstep)*MU2(I,J,L)
            MV2INI(I,J,L) = MV2INI(I,J,L) + hh(nstep)*MV2(I,J,L)
          ENDDO
        TT2INI(I,J) = TT2INI(I,J) + hh(nstep)*TTOP2(I,J)
        ENDDO
        ENDDO
      ENDIF
c
c---------------------------------------------------------------
c
      DPTRUE = (P2(IMID,JMID2Z,LM)-P2OLD(IMID,JMID2Z,LM))/DT
      print9842, DPTRUE
 9842 FORMAT(' True Pressure Tendency (7,7)', 1PE12.3,' Pa/s')
      DPTRUE = DPTRUE*6*3600/100
      print9843, DPTRUE,NSTEP
 9843 FORMAT('True Pressure Tendency (7,7)', 
     X         F12.1,' mb/6h #%# NSTEP:',I4)
c
c---------------------------------------------------------------
c
 5000 CONTINUE
c
c-----------------------------------------------------------------------
C--------------OUTPUT SECTION-------------------------------------------
C
        OPEN(UNIT=33,FILE='graf.dat',FORM='UNFORMATTED')
        print*,' Opening Graph-File'
          WRITE (33)  NTSTEP,XVAR,YVAR
          WRITE (33)  NTSTEP,XVAR,UVAR
          WRITE (33)  NTSTEP,XVAR,VVAR
          WRITE (33)  NTSTEP,XVAR,DVAR
          WRITE (33)  NTSTEP,XVAR,KVAR
          WRITE (33)  NTSTEP,XVAR,AVAR
          WRITE (33)  NTSTEP,XVAR,PVAR
          WRITE (33)  NTSTEP,XVAR,MVAR
          WRITE (33)  NTSTEP,XVAR,N1VAR
          WRITE (33)  NTSTEP,XVAR,N2VAR
          WRITE (33)  NTSTEP,XVAR,BRVAR
          WRITE (33)  NTSTEP,XVAR,N1PVAR
          WRITE (33)  NTSTEP,XVAR,N2PVAR
          WRITE (33)  NTSTEP,XVAR,BRPVAR
        CLOSE(UNIT=33)
C
C-----------------------------------------------------------------------
c
C     Write out Initialized fields .
      IF ( INIT ) THEN
        IF ( DT .LT. 0 ) THEN
           WIFILE = 'FINI1'
           OPEN(UNIT=60,FILE=WIFILE,FORM='UNFORMATTED')
           WRITE(60) P2INI,MU2INI,MV2INI, TT2INI
           CLOSE(UNIT=60)
        ELSE IF ( DT .GT. 0 ) THEN
           CALL PTOZ(P2INI,TT2INI,OROG, Z2,PSEA2,IM,JM2,LM)
           CALL MMTOUV(P2INI,MU2INI,MV2INI, U2,V2,IM,JM2,LM)
           WIFILE = 'FINIT'
           OPEN(UNIT=65,FILE=WIFILE,FORM='UNFORMATTED')
           WRITE(65) PSEA2,Z2,U2,V2,TT2INI,OROG
           WRITE(65) P2INI,MU2INI,MV2INI
           CLOSE(UNIT=65)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
c
      STOP 1000
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
  998 STOP 998  
  999 STOP 999  
c
      END
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
