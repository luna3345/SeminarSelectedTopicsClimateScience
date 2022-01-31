      PROGRAM RMODEL
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
C     ******************************************************************
C     *                                                                *
C     *      Coding started on 13th March, 1991, at MISU.              *
C     *                                                                *
C     ******************************************************************
C
                              P A R A M E T E R
     *   IM= 13, JM=07, LM=5,
     *   KB=IM+2*JM-3, 
     *   IMM1=IM-1, IMM2=IM-2, IMM3=IM-3, IMM4=IM-4, 
     *   JMM1=JM-1, JMM2=JM-2,
     *   LMP1=LM+1, LMM1=LM-1, LMM2=LM-2
C
                               L O G I C A L
     1 FIRST,DRCA,LEQ1,RESTRT,NODIS,   CBNDC,ABHAD,STPGF,CRVAD,THETA
                             D I M E N S I O N
     * DATE(4), DSG(LM), SGML(LM), PAPAI(LM), FFT(LM), FFV(LM),
     * EF4T(LM), FST(LMM1), FSV(LMM1), SG(LMP1), SDP(LMP1), STPL(7),
     * DISL(7), Z0L(7), ISHDE(10),
     * CDCN(JM,2),CDCS(JM,2),CPGFU(JM,2),CPH(JM,2),
     * CVORT(JM,2),CVPOTX(JM,2),DDAMPU(JM,2),DFCN(JM,2),DFCS(JM,2),
     * EM(JM,2),FDTHF(JM,2),H1(JM,2),RH1(JM,2)
                             D I M E N S I O N
     1            WF2(IM,JM),WF3(IM,JM),PDT(IM,JM),ALGPLO(IM,JM),
     2 ZSQLO(IM,JM),PHILO(IM,JM),DPHDZE(IM,JM),PHIM(IM,JM),ZSQM(IM,JM),
     3 CU(IM,JM),CV(IM,JM),PGFU(IM,JM),PGFV(IM,JM),UFLUX(IM,JM),
     4 VFLUX(IM,JM),ADPD(IM,JM),QVORT(IM,JM),VPOTV(IM,JM),VPOTH(IM,JM),
     5 WF4(IM,JM),WF5(IM,JM),WF6(IM,JM),WF7(IM,JM)
                             D I M E N S I O N
     1 ALPIFC(IM,JM,LMP1), PHIFC(IM,JM,LMP1), DIV(IM,JM,LM),
     2 DFI(IM,JM,LM), SD(IM,JM,LMM1)
                             D I M E N S I O N
     1 U(IM,JM,LM), V(IM,JM,LM), T(IM,JM,LM), PD(IM,JM),
     2 PDB(KB,2), TB(KB,LM,2), UB(KB,LM,2), VB(KB,LM,2), PHIS(IM,JM)
                             D I M E N S I O N
     * PDVP(IM,JM), PDU(IM,JM), PDV(IM,JM), PADU(IMM4,JMM2,LM),
     * PADV(IMM4,JMM2,LM), PADT(IMM4,JMM2,LM), CHAL1(JM,2), CHAP1(JM,2),
     * CHAL2(JM,2), CHAP2(JM,2), TUTLN(IM,JM), TUTLS(IM,JM),
     * TVTPE(IM,JM), TVTPW(IM,JM), PDSD(IM,JM), SPDSD(IM,JM,LMM1),
     * TH(IM,JM,LM), CDVL(JM,2), CDVP(JM,2), FVSQ2(IM,JM)
                           E Q U I V A L E N C E
     1 (ALPIFC(1),ALGPLO(1),DIV(1),WF2(1)), (ALPIFC(1,1,2),TVTPE(1)),
     2 (ALPIFC(1,1,3),TVTPW(1)), (ALPIFC(1,1,6),CU(1),PDT(1)),
     * (ALPIFC(1,1,4),FVSQ2(1)),
     3 (PHIFC(1),ZSQLO(1),DFI(1),TH(1)),
     4 (PHIFC(1,1,6),CV(1),PDSD(1),VPOTV(1),WF3(1)),
     5 (WF4(1),ZSQM(1),PGFV(1),VFLUX(1),SD(1),SPDSD(1),PDU(1)),
     6 (WF5(1),PHIM(1),ADPD(1),SD(1,1,2),PDV(1)),
     7 (WF6(1),PHILO(1),SD(1,1,3),TUTLN(1),QVORT(1)),
     8 (WF7(1),DPHDZE(1),PGFU(1),UFLUX(1),SD(1,1,4),TUTLS(1),VPOTH(1))
                               C O M M O N
     1 /ARRAYS/ U,V,T,PD,H1,DSG,RH1,CPH
     * ,ALPIFC,PHIFC,PHIS,SG,SGML,DATE,STPL,DISL,Z0L,WF4,WF5,WF6,WF7
     * ,ISHDE,   PDB,TB,UB,VB
     2 /CONSTS/    NTSD,          H2,DPH,DLM,RTDLM,RTDPH
     * ,DISLP,DIST,DLMD,DPHD,IHRST,                LSTPLM,NTSPH,PT,R
     * ,Z0SLP,Z0T,   IOUT,NSHDE,                  TAU
                             D A T A
     1 A/6376000./,G/9.81/,TWOM/.00014584/,PI/3.141592654/,RN/.1/
     2,R/287.04/,CP/1004./,FIRST/.TRUE./,RESTRT/.FALSE./
C
C-----------------------------------------------------------------------
 2200 FORMAT(' IM=',I3,' JM=',I2,' LM=',I2,' DLMD=',F4.2,' DPHD=',F4.2,'
     " PHISB=',F3.0)
 2210 FORMAT(' DSG=',5(F4.3,1X),'PT=',F6.0,' DT=',F7.2,' W=',F4.3)
 2215 FORMAT(' IDTCF=',I1,' VMTC1=',E8.2,' VMTC2=',E8.2,' VMTC3=',E8.2)
 2220 FORMAT(' COAC=',E7.1,' CODAMP=',F3.0,' TDDAMP=',F3.0,' CRILAT=',F3
     ".0,' CRILTV=',F3.0)
 2230 FORMAT(' NBLHDT=',I2,' NBLHDV=',I2,' NBLUD=',I2)
 2240 FORMAT(' **** NDISC=',I1,' TSTART=',F3.0,' TEND=',F3.0,' TCP=',F3.
     "0,' ****')
 2300 FORMAT(' LSTPLM=',I1,' STPL=',7(F5.0))
 2310 FORMAT(' Z0L=',7(F6.0),' DISL=',F3.0,6(1X,F3.0))
 2320 FORMAT(' Z0T=',F4.0,' DIST=',F2.0,' Z0SLP=',F5.0,' DISLP=',F3.1)
 2290 FORMAT(' ISHDE=',10I3)
 2620 FORMAT('0REFERENCE ISOLINE IS LOCATED ALONG THE BORDER OF A FIELD
     "OF ZEROS, ON THE SIDE FACING THE NEIGHBORING FIELD OF NINES')
 2000 FORMAT(' CHECK POINT AT NTSD=',I4)
C
 2390 FORMAT ('0VALUES 1,2 OF THE CONSTANT MULTIPLYING JACOBIAN1 RESULT
     *IN A SECOND,FOURTH ORDER VORT,T,Q HORIZONTAL ADVECTION')
 3300 FORMAT (' CJCB1=',F2.0)
 3500 FORMAT ('0VALUES 1 AND 4/3 OF CKEA1 RESULT IN THE SECOND AND APPRO
     *XIMATELY FOURTH ORDER KE ADVECTION TERMS')
 3510 FORMAT (' CKEA1=',F8.6)
 3520 FORMAT ('0ADAMS-BASHFORTH HAD, A,A**2 CONSERVING VAD, THETA OR T T
     *YPE SCHEMES')
 3530 FORMAT ('               ABHAD=',L1,' CRVAD=',L1,' THETA=',L1,
     * ' STPGF=',L1)
 3420 FORMAT ('0CONSTANT END CONDITIONS, OR TIME-DEPENDENT')
 3310 FORMAT ('               CBNDC=',L1)
C WE DONT WANT TO READ THE DATA FROM CARDS , HENCE
C THE COMMENTING OUT OF THE READ STATEMENTS BELOW .
C
C	THE ORIGINAL PROGRAM READ A LOT OF DATA FROM CARDS .
C	WE USE DATA STATEMENTS INSTEAD
C
	                   D A T A 
     1 DLMD/1.50/,  DPHD/1.00/,  PHISB/30./,
     2 (DSG(L), L=1,5)/.250,.250,.250,.125,.125/,  
     3 PT/20000./,  DT/450.00/,  W/.250/,
     4 IDTCF/4/,  VMTC1/0.63E 05/,  VMTC2/0.13E 06/,  VMTC3/0.18E-02/,
     5 COAC/0.1E 07/,  CODAMP/3./,  TDDAMP/9./,  CRILAT/99./, 
     6 CRILTV/60./,
     7 NBLHDT/5/,  NBLHDV/0/,  NBLUD/3/
                           D A T A
     1 NDISC/8/,  TSTART/0./,  TEND/0.375/,  TCP/36./,
     2 LSTPLM/4/,  (STPL(L), L=1,7)/300.,500.,700.,850.,0.,0.,0./,
     3 (Z0L(L), L=1,7)/8960.,5440.,2880.,1320.,0.,0.,0./,
     4 (DISL(L), L=1,7)/40.,20.,20.,20.,0.,0.,0./,
     5 Z0T/273./,  DIST/2./,  Z0SLP/1000./,  DISLP/2.5/,
     6 (ISHDE(L), L=1,10)/24,99,99,99,99,99,99,99,99,99/,
     7 CJCB1/2./,  CKEA1/1.33333/,
     8 ABHAD/.FALSE. /,  CRVAD/.FALSE./,  THETA/.FALSE./,  
     9 STPGF/.FALSE./,  CBNDC/.TRUE./
      IP=IM
      JP=JM
      LP=LM
      PRINT   2200, IP,JP,LP,DLMD,DPHD,PHISB
C     READ    2210, DSG,PT,DT,W
      PRINT   2210, DSG,PT,DT,W
C     READ    2215, IDTCF,VMTC1,VMTC2,VMTC3
      PRINT   2215, IDTCF,VMTC1,VMTC2,VMTC3
C     READ    2220, COAC,CODAMP,TDDAMP,CRILAT,CRILTV
      PRINT   2220, COAC,CODAMP,TDDAMP,CRILAT,CRILTV
C     READ    2230, NBLHDT,NBLHDV,NBLUD
      PRINT   2230, NBLHDT,NBLHDV,NBLUD
C     READ    2240, NDISC,TSTART,TEND,TCP
      PRINT   2240, NDISC,TSTART,TEND,TCP
C     READ    2300, LSTPLM,STPL
      PRINT   2300, LSTPLM,STPL
C     READ    2310, Z0L,DISL
      PRINT   2310, Z0L,DISL
C     READ    2320, Z0T,DIST,Z0SLP,DISLP
      PRINT   2320, Z0T,DIST,Z0SLP,DISLP
C     READ    2290, ISHDE
      PRINT   2290, ISHDE
      PRINT   2620
C
      PRINT   2390
C     READ    3300, CJCB1
      PRINT   3300, CJCB1
      CJ1TH=CJCB1/3.
      CJ2TH=(1.-CJCB1)/3.
      PRINT   3500
C     READ    3510, CKEA1
      PRINT   3510, CKEA1
      CKEA2=1.-CKEA1
C
      PRINT   3520
C     READ    3530, ABHAD,CRVAD,THETA,STPGF
      PRINT   3530, ABHAD,CRVAD,THETA,STPGF
      PRINT   3420
C     READ    3310, CBNDC
      PRINT   3310, CBNDC
C
      DO 107 L=1,LSTPLM
      STPL(L)=STPL(L)*100.
      Z0L(L)=Z0L(L)*9.81
 107  DISL(L)=DISL(L)*9.81
      Z0SLP=Z0SLP*100.
      DISLP=DISLP*100.
      JCRIT=(CRILAT-PHISB)*.5/DPHD
      JCRITV=(CRILTV-PHISB)*.5/DPHD
C
C--------------RUN CONTROL CONSTANTS------------------------------------
      NTSPH=3600./DT
C--------------NUMBER OF TIME STEPS BETWEEN SUCCESSIVE MAP OUTPUTS
      NSHDE=ISHDE(1)*NTSPH  +.5
      IOUT=1
C--------------NUMBER OF TIME STEPS BETWEEN SUCCESSIVE CHECK POINTS
      NCP=TCP*NTSPH         +.5
C--------------NUMBER OF TIME STEPS AT WHICH THE RUN IS TO BE TERMINATED
      NSTART=TSTART*NTSPH   +.5
      NTSTM=TEND*NTSPH      +.5
C--------------NUMBER OF TIME STEPS WITH DIVEREGENCE DAMPING
      NDDAMP=TDDAMP*NTSPH   +.5
C
C--------------DERIVED SPACE-TIME GRID AND OTHER CONSTANTS--------------
      SG(1)=0.
      DO 272 L=1,LM
 272  SG(L+1) = SG(L)+DSG(L)
      DO 453 L=1,LM
 453  SGML(L)=.5*(SG(L+1)+SG(L))
C
      JBCV=2*JM-2
      JBCH=2*JM-3
      DTR=PI/180.
      DLM=DLMD*DTR
      DPH=DPHD*DTR
      RTDLM=0.5/DLM
      RTDPH=0.5/DPH
      WDT=W*DT
      WPD=(100000.-PT)*W
      DTHF=.5*DT
      VMTC1=VMTC1*IDTCF
      VMTC2=VMTC2*IDTCF
      VMTC3=VMTC3*IDTCF
      ACDTD=DT*COAC*(A*SQRT(0.5*DLM**2+DPH**2)/300000.)**(4./3.)
      AC=DT*COAC*IDTCF
      ALGPT=ALOG(PT)
      ZSQPT=ALGPT*ALGPT
C
      DO 510 L=1,LM
      EF4T(L)=DTHF/DSG(L)/CP
      FFT(L) =DTHF/DSG(L)
 510  FFV(L) =0.25*FFT(L)
      DO 620 L=1,LMM1
      FST(L)=DTHF/(SGML(L+1)-SGML(L))
 620  FSV(L)=0.25*FST(L)
C
      DO 102 J=1,JM
      DO 102 K=1,2
      PH  =(PHISB+(2*J+K-3)*DPHD)*DTR
      CPH(J,K)=COS(PH)
      H1 (J,K)=A*CPH(J,K)
      RH1(J,K)=1./H1(J,K)
      EM(J,K)=0.5*DT/H1(J,K)/DLM
      FDTHF(J,K)=TWOM*SIN(PH)*DTHF
      IF(J.GT.JCRIT) GO TO 101
      ACDT     =AC*(SQRT((H1(J,K)*DLM)**2+(A*DPH)**2)/300000.)**(4./3.)
C
 101  CSD05 =SIN(PH)*SIN(PH+DPH)+COS(PH)*COS(PH+DPH)*COS(DLM)
      CSD08 =SIN(PH)*SIN(PH-DPH)+COS(PH)*COS(PH-DPH)*COS(DLM)
      SND05 =SQRT(1.-CSD05**2)
      SND08 =SQRT(1.-CSD08**2)
      CSB05 =( SIN(PH+DPH)-SIN(PH)*CSD05)/(COS(PH)*SND05)
      CSB08 =(-SIN(PH-DPH)+SIN(PH)*CSD08)/(COS(PH)*SND08)
      SND05H=SQRT(0.5*(1.-CSD05))
      SND08H=SQRT(0.5*(1.-CSD08))
      TND05H=SND05H/SQRT(1.-SND05H**2)
      TND08H=SND08H/SQRT(1.-SND08H**2)
      T58   =TND05H/TND08H
      T85   =TND08H/TND05H
      SNB05 =SQRT(1.-CSB05**2)
      SNB08 =SQRT(1.-CSB08**2)
      CSB   =CSB05*CSB08-SNB05*SNB08
      SNB   =SNB05*CSB08+CSB05*SNB08
      TNA   =(1.+T58*CSB)/(T58*SNB)
      TNB   =(1.+T85*CSB)/(T85*SNB)
      DL1L2 =ATAN(SND05H*TNA)+ATAN(SND05H*SNB05/CSB05)
      DL8L7 =ATAN(SND08H*TNB)+ATAN(SND08H*SNB08/CSB08)
      D05   =ACOS(CSD05)
      D08   =ACOS(CSD08)
      RAR   =1./(0.5*(DL1L2*D05+DL8L7*D08)*A*A)
      ACDTA =ACDT*RAR
      WDTPDA=0.5*(100000.-PT)*WDT*RAR
      DFCN(J,K)=ACDTA *DL1L2/D05
      DFCS(J,K)=ACDTA *DL8L7/D08
      CDCN(J,K)=WDTPDA*DL1L2/D05
      CDCS(J,K)=WDTPDA*DL8L7/D08
C
      CPGFU(J,K)=DT/H1(J,K)*RTDLM*.5
      CVORT(J,K)=.25*RH1(J,K)*DT
      CVPOTX(J,K)=.25*CPGFU(J,K)
      DDAMPU(J,K)=CODAMP*ACDTD*RTDLM/H1(J,K)
      CDVL(J,K)=1./(A*COS(PH)*2.*DLM)
      CDVP(J,K)=1./(A*COS(PH)*2.*DPH)
C
      CHAL1(J,K)=DT/(8.*A*COS(PH)*DLM) *(1.+(CJCB1-1.)/6.)
      CHAP1(J,K)=DT/(8.*A*COS(PH)*DPH) *(1.+(CJCB1-1.)/6.)
      CHAL2(J,K)=-(CJCB1-1.)*CHAL1(J,K)/7.
 102  CHAP2(J,K)=-(CJCB1-1.)*CHAP1(J,K)/7.
C
      H2=A
      CPGFV=DT/H2*RTDPH*.5
      CVPOTY=.25*CPGFV
      EN=.5*DT/H2/DPH
      GDT=G*DT
      OCC=1.1E-03/G*0.25
      FCPIN=-.25/CP*.25
      DDAMPV     =CODAMP*ACDTD*RTDPH/H2
      P0SQ=(1000.*100.)**2
      APAH=0.286/2.
      DO 632 L=1,LMP1
 632  SDP(L)=0.
C
      TAU=0.
      IF(NSTART.EQ.0) GO TO 303
 309  READ(NDISC) NTSD,DATE,IHRST
      IF(NTSD.LT.0) GO TO 307
      IF(NTSD.GT.NSTART) GO TO 307
      READ(NDISC)
      GO TO 309
C
C        RESTARTING SECTION: IF THE RESTARTING TIME IS 24 HOURS, AND
C             .NOT.CBNDC, THE PROGRAM NEEDS A NEW SET OF BOUNDARY VALUES
 307  BACKSPACE NDISC
      BACKSPACE NDISC
      BACKSPACE NDISC
      RESTRT=.TRUE.
C
 303  READ(NDISC) NTSD,DATE,IHRST,U,V
      READ(NDISC) T,PD,PDB,TB,UB,VB,PHIS
      IF (RESTRT .AND. NSTART.EQ.24*NTSPH .AND. .NOT.CBNDC)
     *                                             READ (9) PDB,TB,UB,VB
      IF (NTSD.GE.NTSTM) STOP
      IF (.NOT.RESTRT) GO TO 2900
  11  IF (NSHDE.GE.NTSD) GO TO 2900
      IOUT=IOUT+1
      NSHDE=ISHDE(IOUT)*NTSPH  +.5
      GO TO 11
 420  FIRST=.FALSE.
      RESTRT=.FALSE.
C
C--------------THE CONTINUITY EQUATION, COMPUTATION OF PDT--------------
 149  NTSD=NTSD+1
      NODIS=.TRUE.
      IF(MOD(NTSD,IDTCF).EQ.0) NODIS=.FALSE.
C
      DO 106 I=3,IMM2
      DO 106 J=2,JMM1
 106  PDT(I,J)=0.
C--------------VERTICAL INTEGRATION-------------------------------------
      DO 202 L=1,LM
      DO 202 I=3,IMM2
      JF=JMM2+MOD(I,2)
      DO 202 J=2,JF
 202  PDT(I,J)=PDT(I,J)-DIV(I,J,L)*DSG(L)
C
C--------------COMPUTATION OF SD----------------------------------------
      DO 616 L=1,LMM1
      DO 619 I=2,IMM1,2
      DO 619 J=1,JMM1,JMM2
 619  SD(I,J,L)=0.
      DO 616 I=2,IMM1,IMM3
      DO 616 J=2,JMM2
 616  SD(I,J,L)=0.
C
      DO 203 I=3,IMM2
      JF=JMM2+MOD(I,2)
      DO 203 J=2,JF
      RPDIJ=1./PD(I,J)
      PDTIJ=PDT(I,J)
      SD(I,J,LMM1)=(PDTIJ+DIV(I,J,LM))*DSG(LM)*RPDIJ
      DO 203 IVI=1,LMM2
      L=LMM1-IVI
 203  SD(I,J,L)=SD(I,J,L+1)+(PDTIJ+DIV(I,J,L+1))*DSG(L+1)*RPDIJ
C
C--------------KINETIC ENERGY GENERATION TERMS IN T EQUATION------------
      DO 340 I=3,IMM2
      JF=JMM2+MOD(I,2)
      DO 340 J=2,JF
      ALPDT=PDT(I,J)/PD(I,J)
      IF (THETA) GO TO 629
      DO 341 L=2,LM
 341  SDP(L)=SD(I,J,L-1)
 629  DO 340 L=1,LM
 340  T(I,J,L)=T(I,J,L)+EF4T(L)*DFI(I,J,L)
     *                         *(SDP(L+1)+SDP(L)+(SG(L+1)+SG(L))*ALPDT)
C
C--------------SPLIT VERTICAL ADVECTION---------------------------------
C--------------VERTICAL ADVECTION OF TEMPERATURE------------------------
      IF (.NOT.THETA) GO TO 630
      DO 631 I=1,IM
      JF=JMM1+MOD(I,2)
      DO 631 J=1,JF
      PA=PT
      DO 631 L=1,LM
      PB=PT+SG(L+1)*PD(I,J)
      TH(I,J,L)=T(I,J,L)*(P0SQ/(PA*PB))**APAH
 631  PA=PB
C
      IF (.NOT.CRVAD) GO TO 721
      DO 701 I=3,IMM2
      JF=JMM2+MOD(I,2)
      DO 701 J=2,JF
      TTB=0.
      DO 702 L=1,LMM1
      TTA=TTB
      TTB=SD(I,J,L)*(TH(I,J,L+1)-TH(I,J,L))
 702  T(I,J,L) =T(I,J,L) -FFT(L) *(T(I,J,L) /TH(I,J,L))*(TTA+TTB)
 701  T(I,J,LM)=T(I,J,LM)-FFT(LM)*(T(I,J,LM)/TH(I,J,LM))*TTB
      GO TO 622
C
 721  DO 723 I=3,IMM2
      JF=JMM2+MOD(I,2)
      DO 723 J=2,JF
      TTB=0.
      DO 724 L=1,LMM1
      TTA=TTB
      TTB=SD(I,J,L)*FST(L)*(TH(I,J,L+1)-TH(I,J,L))
 724  T(I,J,L) =T(I,J,L) -(T(I,J,L) /TH(I,J,L))*(TTA+TTB)
 723  T(I,J,LM)=T(I,J,LM)-(T(I,J,LM)/TH(I,J,LM))*TTB
      GO TO 622
C
 630  IF (.NOT.CRVAD) GO TO 621
      DO 501 I=3,IMM2
      JF=JMM2+MOD(I,2)
      DO 501 J=2,JF
      TTB=0.
      DO 502 L=1,LMM1
      TTA=TTB
      TTB=SD(I,J,L)*(T(I,J,L+1)-T(I,J,L))
 502  T(I,J,L) =T(I,J,L) -FFT(L)*(TTA+TTB)
 501  T(I,J,LM)=T(I,J,LM)-FFT(LM)*TTB
      GO TO 622
C
 621  DO 623 I=3,IMM2
      JF=JMM2+MOD(I,2)
      DO 623 J=2,JF
      TTB=0.
      DO 624 L=1,LMM1
      TTA=TTB
      TTB=SD(I,J,L)*FST(L)*(T(I,J,L+1)-T(I,J,L))
 624  T(I,J,L) =T(I,J,L) -TTA-TTB
 623  T(I,J,LM)=T(I,J,LM)-TTB
C
C--------------VERTICAL ADVECTION OF MOMENTUM---------------------------
 622  DO 617 L=1,LMM1
      DO 618 I=2,IMM1
      JS=1+MOD(I,2)
      DO 618 J=JS,JMM1
 618  PDSD(I,J)=PD(I,J)*SD(I,J,L)
C
      DO 617 I=3,IMM2
      MI2=MOD(I,2)
      JD=-1+2*MI2
      JF=JMM1-MI2
      DO 617 J=2,JF
 617  SPDSD(I,J,L)=PDSD(I+1,J)+PDSD(I,J+JD)+PDSD(I-1,J)+PDSD(I,J)
C
      IF(.NOT.CRVAD) GO TO 625
      DO 504 I=3,IMM2
      JF=JMM1-MOD(I,2)
      DO 504 J=2,JF
      TUB=0.
      TVB=0.
      DO 505 L=1,LMM1
      TUA=TUB
      TVA=TVB
      TUB=SPDSD(I,J,L)*(U(I,J,L+1)-U(I,J,L))
      TVB=SPDSD(I,J,L)*(V(I,J,L+1)-V(I,J,L))
      U(I,J,L) =U(I,J,L) -FFV(L)*(TUA+TUB)/PDVP(I,J)
 505  V(I,J,L) =V(I,J,L) -FFV(L)*(TVA+TVB)/PDVP(I,J)
      U(I,J,LM)=U(I,J,LM)-FFV(LM)*TUB     /PDVP(I,J)
 504  V(I,J,LM)=V(I,J,LM)-FFV(LM)*TVB     /PDVP(I,J)
      GO TO 626
C
 625  DO 627 I=3,IMM2
      JF=JMM1-MOD(I,2)
      DO 627 J=2,JF
      TUB=0.
      TVB=0.
      DO 628 L=1,LMM1
      TUA=TUB
      TVA=TVB
      TUB=SPDSD(I,J,L)*FSV(L)*(U(I,J,L+1)-U(I,J,L))
      TVB=SPDSD(I,J,L)*FSV(L)*(V(I,J,L+1)-V(I,J,L))
      U(I,J,L) =U(I,J,L)-(TUA+TUB)/PDVP(I,J)
 628  V(I,J,L) =V(I,J,L)-(TVA+TVB)/PDVP(I,J)
      U(I,J,LM)=U(I,J,LM)-TUB     /PDVP(I,J)
 627  V(I,J,LM)=V(I,J,LM)-TVB     /PDVP(I,J)
C
C--------------VERTICAL TURBULENT TRANSPORT OF MOMENTUM-----------------
 626  IF(NODIS) GO TO 390
C--------------( SIMPLIFIED MIXING-LENGTH THEORY )----------------------
      DO 360 I=3,IMM2
      MI2=MOD(I,2)
      JF=JMM1-MI2
      DO 360 J=2,JF
      JP1=J+MI2
      JM1=JP1-1
      DU2=U(I,J,LM)-U(I,J,LMM1)
      DU1=U(I,J,LMM1)-U(I,J,LMM2)
      DV2=V(I,J,LM)-V(I,J,LMM1)
      DV1=V(I,J,LMM1)-V(I,J,LMM2)
      US=U(I,J,LM)+.5*DU2
      VS=V(I,J,LM)+.5*DV2
      AVS=SQRT(US*US+VS*VS)
      APHIS=PHIS(I+1,J)+PHIS(I-1,J)+PHIS(I,JP1)+PHIS(I,JM1)
      F7=-VMTC3  *AVS*(1.+OCC*APHIS)
      IF(APHIS.LT.0.2) F7=0.
      TU3=F7*US
      TV3=F7*VS
      AVDP2=SQRT(DU2*DU2+DV2*DV2)
      AVDP1=SQRT(DU1*DU1+DV1*DV1)
      DP2=PDVP(I,J)*(SGML(LM)  -SGML(LMM1))
      DP1=PDVP(I,J)*(SGML(LMM1)-SGML(LMM2))
      FAC2=VMTC2 *AVDP2/(DP2*DP2)
      FAC1=VMTC1 *AVDP1/(DP1*DP1)
      TU2=FAC2*DU2
      TV2=FAC2*DV2
      TU1=FAC1*DU1
      TV1=FAC1*DV1
      FAC0=GDT/PDVP(I,J)
      FAC3=FAC0/DSG(LM)
      FAC2=FAC0/DSG(LMM1)
      FAC1=FAC0/DSG(LMM2)
      U(I,J,LM)=FAC3*(TU3-TU2)  +U(I,J,LM)
      U(I,J,LMM1)=FAC2*(TU2-TU1) +U(I,J,LMM1)
      U(I,J,LMM2)=FAC1*TU1       +U(I,J,LMM2)
      V(I,J,LM)=FAC3*(TV3-TV2)+V(I,J,LM)
      V(I,J,LMM1)=FAC2*(TV2-TV1) +V(I,J,LMM1)
 360  V(I,J,LMM2)=FAC1*TV1       +V(I,J,LMM2)
C
C-    ---------ADVECTION LOOP-------------------------------------------
 390  DO 206 L=1,LM
C
      DO 601 I=1,IM
      JF=JM-MOD(I,2)
      DO 601 J=1,JF
      PDU(I,J)=PDVP(I,J)*U(I,J,L)
 601  PDV(I,J)=PDVP(I,J)*V(I,J,L)
C
      DO 635 I=1,IMM1
      MI2=MOD(I,2)
      KV=1+MI2
      KH=2-MI2
      DO 635 J=1,JMM1
      JP1=J+1-MI2
      IF (.NOT.THETA) GO TO 636
      TUTLN(I,J)=(PDU(I+1,J)+PDU(I,JP1))*(TH(I+1,JP1,L)-TH(I,J,L))
      TVTPE(I,J)=(PDV(I+1,J)*CPH(J,KH)+PDV(I,JP1)*CPH(JP1,KV))
     *                                  *(TH(I+1,JP1,L)-TH(I,J,L))
      GO TO 635
 636  TUTLN(I,J)=(PDU(I+1,J)+PDU(I,JP1))*(T (I+1,JP1,L)-T (I,J,L))
      TVTPE(I,J)=(PDV(I+1,J)*CPH(J,KH)+PDV(I,JP1)*CPH(JP1,KV))
     *                                  *(T (I+1,JP1,L)-T (I,J,L))
 635  CONTINUE
C
      DO 637 I=2,IM
      MI2=MOD(I,2)
      KV=1+MI2
      KH=2-MI2
      DO 637 J=1,JMM1
      JP1=J+1-MI2
      IF (.NOT.THETA) GO TO 638
      TUTLS(I,J)=(PDU(I,JP1)+PDU(I-1,J))*(TH(I,J,L)-TH(I-1,JP1,L))
      TVTPW(I,J)=(PDV(I,JP1)*CPH(JP1,KV)+PDV(I-1,J)*CPH(J,KH))
     *                                  *(TH(I-1,JP1,L)-TH(I,J,L))
      GO TO 637
 638  TUTLS(I,J)=(PDU(I,JP1)+PDU(I-1,J))*(T (I,J,L)-T (I-1,JP1,L))
      TVTPW(I,J)=(PDV(I,JP1)*CPH(JP1,KV)+PDV(I-1,J)*CPH(J,KH))
     *                                  *(T (I-1,JP1,L)-T (I,J,L))
 637  CONTINUE
C
C--------------THE FIRST LAW OF THERMODYNAMICS--------------------------
      DO 605 I=3,IMM2
      MI2=MOD(I,2)
      KV=1+MI2
      KH=2-MI2
      JF=JMM2+MI2
      MI2P3=MI2+3
      IDTCB=MIN0(I-3,IMM2-I)
      JMIBC=JBCH+MI2
C
      DO 605 J=2,JF
      JM1=J-MI2
      JP1=JM1+1
      JMMI2=J-MI2
      ITWJ=J+J
      TIJL=T(I,J,L)
      NLTCB=MIN0(IDTCB,ITWJ-MI2P3,JMIBC-ITWJ)
      HDTT=0.
      IF(NODIS) GO TO 250
      IF(PHIS(I,J).GT.4905.) GO TO 250
      IF(NLTCB.GE.NBLHDT) GO TO 250
C
C        HORIZONTAL DIFFUSION
      HDTT= DFCN(J,KH)*(T(I+1,JP1,L)+T(I-1,JP1,L)-TIJL-TIJL)
     *     -DFCS(J,KH)*(TIJL+TIJL-T(I-1,JM1,L)-T(I+1,JM1,L))
C
 250  IF(NLTCB.LT.NBLUD) GO TO 640
      ADT=-(CHAL1(J,KH)*(TUTLN(I-1,JM1)+TUTLN(I,J)
     2                  +TUTLS(I,J)    +TUTLS(I+1,JM1))
     3     +CHAP1(J,KH)*(TVTPE(I-1,JM1)+TVTPE(I,J)
     4                  +TVTPW(I+1,JM1)+TVTPW(I,J))
     5     +CHAL2(J,KH)*(TUTLN(I-2,J-1)+TUTLN(I+1,JP1)
     6                  +TUTLS(I-1,JP1)+TUTLS(I+2,J-1))
     7     +CHAP2(J,KH)*(TVTPE(I-2,J-1)+TVTPE(I+1,JP1)
     8                  +TVTPW(I+2,J-1)+TVTPW(I-1,JP1)))/PD(I,J)
      GO TO 308
C
 640  UIJL=0.25*(U(I+1,J,L)+U(I,JP1,L)+U(I-1,J,L)+U(I,JM1,L))
      VIJL=0.25*(V(I+1,J,L)+V(I,JP1,L)+V(I-1,J,L)+V(I,JM1,L))
      TTA=EM(J,KH)*UIJL
      TTB=EN      *VIJL
      P=-TTA-TTB
      Q=TTA-TTB
      ISP=SIGN(1.,P)
      ISQ=SIGN(1.,Q)
      I1=I+ISP
      J1=JMMI2+(1+ISP)/2
      I2=I-ISQ
      J2=JMMI2+(1+ISQ)/2
      I3=I+ISP-ISQ
      J3=J+(ISP+ISQ)/2
      P=ABS(P)
      Q=ABS(Q)
      F3=P*Q
      F0=F3-P-Q
      F1=P-F3
      F2=Q-F3
           ADT=F0*T (I,J,L)+F1*T (I1,J1,L)+F2*T (I2,J2,L)+F3*T (I3,J3,L)
      IF (THETA)
     *     ADT=F0*TH(I,J,L)+F1*TH(I1,J1,L)+F2*TH(I2,J2,L)+F3*TH(I3,J3,L)
 308  IF (THETA) ADT=T(I,J,L)/TH(I,J,L)*ADT
      DTCH=0.
      IF (.NOT.ABHAD) GO TO 605
      PADTS = PADT(I-2,J-1,L)
      PADT(I-2,J-1,L) = ADT
      IF (NTSD-NSTART.EQ.1) GO TO 605
      ADT = 0.5*(ADT+ADT+ADT-PADTS)
 605  WF2(I,J)=TIJL+ADT+HDTT+DTCH
C
      DO 512 I=3,IMM2
      JF=JMM2+MOD(I,2)
      DO 512 J=2,JF
 512  T(I,J,L)=WF2(I,J)
C
C--------------THE EQUATION OF MOTION, HORIZONTAL ADVECTION-------------
      DO 351 I=1,IM
      MI2=MOD(I,2)
      JF=JM-MI2
      KV=1+MI2
      DO 351 J=1,JF
      PDV(I,J)=PDV(I,J)*CPH(J,KV)
 351  VPOTV(I,J)=U(I,J,L)**2+V(I,J,L)**2
C
      DO 613 I=2,IMM1
      MI2=MOD(I,2)
      KV=1+MI2
      KH=2-MI2
      JS=1+MI2
      DO 613 J=JS,JMM1
      JM1=J-MI2
      JP1=JM1+1
      QVORT(I,J)=CVORT(J,KH)*((V(I+1,J,L)-V(I-1,J,L))*RTDLM
     2           -(U(I,JP1,L)*CPH(JP1,KV)-U(I,JM1,L)*CPH(JM1,KV))*RTDPH)
     3           /PD(I,J)
 613  VPOTH(I,J)=VPOTV(I+1,J)+VPOTV(I-1,J)+VPOTV(I,JP1)+VPOTV(I,JM1)
C
      DO 641 I=3,IMM2
      MI2=MOD(I,2)
      JF=JMM2+MI2
      DO 641 J=2,JF
      JM1=J-MI2
      JP1=JM1+1
 641  FVSQ2(I,J)=0.5*
     *     (VPOTV(I+2,JP1)+VPOTV(I+1,J+1)+VPOTV(I-1,J+1)+VPOTV(I-2,JP1)
     *     +VPOTV(I-2,JM1)+VPOTV(I-1,J-1)+VPOTV(I+1,J-1)+VPOTV(I+2,JM1))
C
      DO 117 I=3,IMM2
      MI2=MOD(I,2)
      KV=1+MI2
      KH=2-MI2
      JF=JMM1-MI2
      MI2M4=MI2-4
      IDTCB=MIN0(I-3,IMM2-I)
      JMIBC=JBCV-MI2
C
      DO 117 J=2,JF
      JP1=J+MI2
      JM1=JP1-1
      JPMI2=J+MI2
      ITWJ=J+J
      UIJL=U(I,J,L)
      VIJL=V(I,J,L)
      NLTCB=MIN0(IDTCB,ITWJ+MI2M4,JMIBC-ITWJ)
      HDTU=0.
      HDTV=0.
      IF(NODIS) GO TO 252
      IF(NLTCB.GE.NBLHDV.AND.J.LT.JCRITV) GO TO 252
C
C        HORIZONTAL DIFFUSION
      HDTU= DFCN(J,KV)*(U(I+1,JP1,L)+U(I-1,JP1,L)-UIJL-UIJL)
     *     -DFCS(J,KV)*(UIJL+UIJL-U(I-1,JM1,L)-U(I+1,JM1,L))
      HDTV= DFCN(J,KV)*(V(I+1,JP1,L)+V(I-1,JP1,L)-VIJL-VIJL)
     *     -DFCS(J,KV)*(VIJL+VIJL-V(I-1,JM1,L)-V(I+1,JM1,L))
C
 252  IF(NLTCB.LT.NBLUD) GO TO 642
      HVORT=2.*(QVORT(I+1,J)+QVORT(I-1,J)+QVORT(I,JP1)+QVORT(I,JM1))
      ADU=-CVPOTX(J,KV)*(CKEA1*(VPOTH(I+1,J)-VPOTH(I-1,J))
     1                  +CKEA2*(FVSQ2(I+1,J)-FVSQ2(I-1,J)))
     2    +(CJ1TH*(PDV(I,J)*HVORT
     3            +(PDV(I+1,JP1)+PDV(I+1,JM1))*QVORT(I+1,J)
     4            +(PDV(I-1,JP1)+PDV(I-1,JM1))*QVORT(I-1,J))
     5     +CJ2TH*
     6     (PDV(I+1,JP1)*(QVORT(I+2,JP1)+QVORT(I,JP1)+QVORT(I,JM1))
     7     +PDV(I-1,JP1)*(QVORT(I-2,JP1)+QVORT(I,JP1)+QVORT(I,JM1))
     8     +PDV(I-1,JM1)*(QVORT(I-2,JM1)+QVORT(I,JM1)+QVORT(I,JP1))
     9     +PDV(I+1,JM1)*(QVORT(I+2,JM1)+QVORT(I,JM1)+QVORT(I,JP1))))
     A     /CPH(J,KV)
      ADV=-CVPOTY      *(CKEA1*(VPOTH(I,JP1)-VPOTH(I,JM1))
     1                  +CKEA2*(FVSQ2(I,JP1)-FVSQ2(I,JM1)))
     2     -CJ1TH*(PDU(I,J)*HVORT
     3            +(PDU(I+1,JP1)+PDU(I-1,JP1))*QVORT(I,JP1)
     4            +(PDU(I+1,JM1)+PDU(I-1,JM1))*QVORT(I,JM1))
     5     -CJ2TH*
     6     (PDU(I+1,JP1)*(QVORT(I+1,J+1)+QVORT(I+1,J)+QVORT(I-1,J))
     7     +PDU(I-1,JP1)*(QVORT(I-1,J+1)+QVORT(I-1,J)+QVORT(I+1,J))
     8     +PDU(I-1,JM1)*(QVORT(I-1,J-1)+QVORT(I-1,J)+QVORT(I+1,J))
     9     +PDU(I+1,JM1)*(QVORT(I+1,J-1)+QVORT(I+1,J)+QVORT(I-1,J)))
      GO TO 306
C
 642  TTA=EM(J,KV)*UIJL
      TTB=EN      *VIJL
      P=-TTA-TTB
      Q=TTA-TTB
      ISP=SIGN(1.,P)
      ISQ=SIGN(1.,Q)
      I1=I+ISP
      J1=JPMI2-(1-ISP)/2
      I2=I-ISQ
      J2=JPMI2-(1-ISQ)/2
      I3=I+ISP-ISQ
      J3=J+(ISP+ISQ)/2
      P=ABS(P)
      Q=ABS(Q)
      F3=P*Q
      F0=F3-P-Q
      F1=P-F3
      F2=Q-F3
      ADU=F0*UIJL+F1*U(I1,J1,L)+F2*U(I2,J2,L)+F3*U(I3,J3,L)
      ADV=F0*VIJL+F1*V(I1,J1,L)+F2*V(I2,J2,L)+F3*V(I3,J3,L)
 306  IF (.NOT.ABHAD) GO TO 606
      PADUS = PADU(I-2,J-1,L)
      PADVS = PADV(I-2,J-1,L)
      PADU(I-2,J-1,L) = ADU
      PADV(I-2,J-1,L) = ADV
      IF (NTSD-NSTART.EQ.1) GO TO 606
      ADU = 0.5*(ADU+ADU+ADU-PADUS)
      ADV = 0.5*(ADV+ADV+ADV-PADVS)
 606  WF2(I,J)=UIJL+ADU+HDTU
 117  WF3(I,J)=VIJL+ADV+HDTV
C
      DO 206 I=3,IMM2
      JF=JMM1-MOD(I,2)
      DO 206 J=2,JF
      U(I,J,L)=WF2(I,J)
 206  V(I,J,L)=WF3(I,J)
C
C--------------COMPUTATION OF PD----------------------------------------
      DO 211 I=3,IMM2
      DO 211 J=2,JMM1
 211  PD(I,J)=PD(I,J)+DT*PDT(I,J)
C
C-    ---------DRY CONVECTIVE ADJUSTMENT--------------------------------
      IF(MOD(NTSD,NTSPH).NE.0) GO TO 450
      DO 230 I=3,IMM2
      JF=JMM2+MOD(I,2)
      DO 230 J=2,JF
      PDIJ=PD(I,J)
      DO 451 L=1,LM
 451  PAPAI(L)=(PT+SGML(L)*PDIJ)**(-.2858964143)
 232  DRCA=.FALSE.
      PPLP1=PAPAI(1)
      DO 231 L=1,LMM1
      PPL=PPLP1
      PPLP1=PAPAI(L+1)
      THL=T(I,J,L)*PPL
      THLP1=T(I,J,L+1)*PPLP1
      DTH=THLP1-THL
      IF(DTH.LT.RN) GO TO 231
      RDSD=DSG(L)/DSG(L+1)
      DTL=DTH/(PPL+PPLP1*RDSD)
      DTLP1=-DTL*RDSD
      T(I,J,L)=T(I,J,L)+DTL
      T(I,J,L+1)=T(I,J,L+1)+DTLP1
      DRCA=.TRUE.
 231  CONTINUE
      IF(DRCA) GO TO 232
 230  CONTINUE
C
C------------- TIME INTERPOLATION OF PD AND T AT THE OUTER BOUNDARY ----
 450  IF (CBNDC) GO TO 603
      TAU=NTSD*DT
      IF (TAU.GT.86400.) TAU=TAU-86400.
      IF (NTSD.EQ.24*NTSPH+1)                   READ (9) PDB,TB,UB,VB
      CALL RTIPDT
C
C------------- SPACE INTERPOLATION OF PD AND T AT THE INNER BOUNDARY ---
 603  DO 145 I=2,IMM1,2
      PD(I,   1)=0.25*(PD(I+1, 2)+PD(I-1, 2)
     2           +PD(I-1,   1)+PD(I+1,   1))
      PD(I,JMM1)=0.25*(PD(I+1,JM)+PD(I-1,JM)
     2           +PD(I-1,JMM1)+PD(I+1,JMM1))
      DO 145 L=1,LM
      T(I,   1,L)=0.25*(T(I+1,2,L)+T(I-1, 2,L)
     2            +T(I-1,   1,L)+T(I+1,   1,L))
 145  T(I,JMM1,L)=0.25*(T(I+1,JM,L)+T(I-1,JM,L)
     2            +T(I-1,JMM1,L)+T(I+1,JMM1,L))
C
      DO 146 J=2,JMM2
      PD(   2,J)=0.25*(PD( 3,J+1)+PD(   1,J+1)
     2           +PD(   1,J)+PD( 3,J))
      PD(IMM1,J)=0.25*(PD(IM,J+1)+PD(IMM2,J+1)
     2           +PD(IMM2,J)+PD(IM,J))
      DO 146 L=1,LM
      T(   2,J,L)=0.25*(T( 3,J+1,L)+T(   1,J+1,L)
     2            +T(   1,J,L)+T( 3,J,L))
 146  T(IMM1,J,L)=0.25*(T(IM,J+1,L)+T(IMM2,J+1,L)
     2            +T(IMM2,J,L)+T(IM,J,L))
C
C--------------OUTPUT SECTION-------------------------------------------
 2900 IF(NTSD.EQ.NSHDE) CALL ROUTPT
      IF(NTSD.EQ.0)     CALL RDGNS
C
C-----------------------------------------------------------------------
C              PREPARATION FOR THE LOOP 401
      DO 602 I=2,IMM1
      MI2=MOD(I,2)
      JS=2-MI2
      DO 602 J=JS,JMM1
      JP1=J+MI2
      JM1=JP1-1
 602  PDVP(I,J)=0.25*(PD(I+1,J)+PD(I-1,J)+PD(I,JP1)+PD(I,JM1))
      DO 633 I=2,IMM1,2
      DO 633 J=1,JM,JMM1
 633  PDVP(I,J)=0.50*(PD(I+1,J)+PD(I-1,J))
      DO 634 I=1,IM,IMM1
      DO 634 J=1,JMM1
 634  PDVP(I,J)=0.50*(PD(I,J+1)+PD(I,J))
C
      DO 400 I=1,IM
      JF=JMM1+MOD(I,2)
      DO 400 J=1,JF
      ALGPL=ALOG(PT+PD(I,J))
      ALGPLO(I,J)=ALGPL
      ZSQLO(I,J)=ALGPL*ALGPL
 400  PHILO(I,J)=PHIS(I,J)
C
      DO 401 IVI=1,LM
      L=LMP1-IVI
      LEQ1=.FALSE.
      IF(L.EQ.1) LEQ1=.TRUE.
      DO 402 I=1,IM
      JF=JMM1+MOD(I,2)
      DO 402 J=1,JF
      IF(LEQ1) GO TO 410
      ALGPUP=ALOG(PT+SG(L)*PD(I,J))
      ZSQUP=ALGPUP*ALGPUP
      GO TO 411
 410  ALGPUP=ALGPT
      ZSQUP=ZSQPT
 411  ALGPL=ALGPLO(I,J)
      DPHDZE(I,J)=R*T(I,J,L)/(ALGPUP+ALGPL)
      ZSQM(I,J)=.5*(ZSQLO(I,J)+ZSQUP)
      PHILIJ=PHILO(I,J)
      DFIJL=DPHDZE(I,J)*(ZSQLO(I,J)-ZSQUP)
      PHIUP=PHILIJ+DFIJL
      DFI(I,J,L)=DFIJL
      PHIM(I,J)=PHILIJ+PHIUP
      IF(LEQ1) GO TO 402
      ALGPLO(I,J)=ALGPUP
      PHILO(I,J)=PHIUP
      ZSQLO(I,J)=ZSQUP
 402  CONTINUE
C
C--------------CALCULATE DIAGONAL PART OF CORRECTION--------------------
      DO 403 I=3,IMM2
      MI2=MOD(I,2)
      KH=2-MI2
      JF=JMM2+MI2
      DO 403 J=2,JF
      PHIJ=PHIM(I,J)
      ZSQIJ=ZSQM(I,J)
      DFDZIJ=DPHDZE(I,J)
      JM1=J-MI2
      JP1=JM1+1
 403  DIV(I,J,L)=-CDCN(J,KH)  *(PHIM(I+1,JP1)+PHIM(I-1,JP1)-PHIJ-PHIJ
     2             +(DPHDZE(I+1,JP1)+DFDZIJ)*(ZSQM(I+1,JP1)-ZSQIJ)
     3             +(DPHDZE(I-1,JP1)+DFDZIJ)*(ZSQM(I-1,JP1)-ZSQIJ))
     4           +CDCS(J,KH)  *(PHIJ+PHIJ-PHIM(I-1,JM1)-PHIM(I+1,JM1)
     5             +(DPHDZE(I-1,JM1)+DFDZIJ)*(ZSQIJ-ZSQM(I-1,JM1))
     6             +(DPHDZE(I+1,JM1)+DFDZIJ)*(ZSQIJ-ZSQM(I+1,JM1)))
C
C---------------CALCULATION OF CU AND CV--------------------------------
      DO 404 I=2,IMM1
      MI2=MOD(I,2)
      KV=1+MI2
      JS=2-MI2
      DO 404 J=JS,JMM1
      JP1=J+MI2
      JM1=JP1-1
      DPIP1=DPHDZE(I+1,J)
      ZMIP1=ZSQM(I+1,J)
      DPJP1=DPHDZE(I,JP1)
      ZMJP1=ZSQM(I,JP1)
      DPIM1=DPHDZE(I-1,J)
      ZMIM1=ZSQM(I-1,J)
      DPJM1=DPHDZE(I,JM1)
      ZMJM1=ZSQM(I,JM1)
      IF (.NOT.STPGF) GO TO 643
      CU(I,J)=-CPGFU(J,KV)*(DPIP1+DPIM1)*(ZMIP1-ZMIM1)
      CV(I,J)=-CPGFV      *(DPJP1+DPJM1)*(ZMJP1-ZMJM1)
      GO TO 404
 643  CU(I,J)=-.5*CPGFU(J,KV)*
     2        ((DPJP1+DPIM1)*(ZMJP1-ZMIM1)
     3        +(DPIP1+DPJP1)*(ZMIP1-ZMJP1)
     4        +(DPJM1+DPIM1)*(ZMJM1-ZMIM1)
     5        +(DPIP1+DPJM1)*(ZMIP1-ZMJM1))
      CV(I,J)=-.5*CPGFV*
     2        ((DPIM1+DPJM1)*(ZMIM1-ZMJM1)
     3        +(DPJP1+DPIM1)*(ZMJP1-ZMIM1)
     4        +(DPIP1+DPJM1)*(ZMIP1-ZMJM1)
     5        +(DPJP1+DPIP1)*(ZMJP1-ZMIP1))
 404  CONTINUE
C
C--------------CALCULATION OF PGFU AND PGFV-----------------------------
      DO 405 I=2,IMM1
      MI2=MOD(I,2)
      KV=1+MI2
      JS=2-MI2
      DO 405 J=JS,JMM1
      JP1=J+MI2
      JM1=JP1-1
      PGFU(I,J)=-CPGFU(J,KV)*(PHIM(I+1,J)-PHIM(I-1,J))+CU(I,J)
 405  PGFV(I,J)=-CPGFV      *(PHIM(I,JP1)-PHIM(I,JM1))+CV(I,J)
      IF(RESTRT)  GO TO 415
C
C        HORIZONTAL PART OF THE OMEGA ALPHA TERM
      IF (THETA) GO TO 412
      DO 430 I=2,IMM1
      MI2=MOD(I,2)
      JD=-1+MI2+MI2
      JS=2-MI2
      DO 430 J=JS,JMM1
      PDBAR=PD(I+1,J)+PD(I,J)+PD(I,J+JD)+PD(I-1,J)
 430  ADPD(I,J)=PDBAR*(U(I,J,L)*CU(I,J)+V(I,J,L)*CV(I,J))
      DO 431 I=3,IMM2
      MI2=MOD(I,2)
      JF=JMM2+MI2
      DO 431 J=2,JF
      JM1=J-MI2
      JP1=JM1+1
 431   T(I,J,L)=T(I,J,L)+FCPIN*(ADPD(I,JP1)+ADPD(I,JM1)+ADPD(I-1,J)
     2                    +ADPD(I+1,J))/PD(I,J)
C
C--------------UPDATE U AND V UNLESS THIS IS THE FIRST PASS-------------
C--------------CORIOLIS AND PRESSURE GRADIENT FORCE---------------------
 412  IF(FIRST) GO TO 415
      DO 215 I=3,IMM2
      MI2=MOD(I,2)
      KV=1+MI2
      JF=JMM1-MI2
C
      DO 215 J=2,JF
      UIJL=U(I,J,L)
      UHLD=UIJL
      VIJL=V(I,J,L)
      RF0SP1=1./(1.+FDTHF(J,KV)**2)
      UIJL=UIJL+PGFU(I,J)+FDTHF(J,KV)*VIJL
      VIJL=VIJL+PGFV(I,J)-FDTHF(J,KV)*UHLD
      U(I,J,L)=(UIJL+FDTHF(J,KV)*VIJL)*RF0SP1
 215  V(I,J,L)=VIJL-FDTHF(J,KV)*U(I,J,L)
C
 415  DO 103 I=2,IMM1
      JS=2-MOD(I,2)
      DO 103 J=JS,JMM1
      UFLUX(I,J)=PDVP(I,J)*U(I,J,L)-WPD*PGFU(I,J)
 103  VFLUX(I,J)=PDVP(I,J)*V(I,J,L)-WPD*PGFV(I,J)
C
C-    ---------COMPUTATION OF DIVERGENCE--------------------------------
      DO 401 I=3,IMM2
      MI2=MOD(I,2)
      KV=1+MI2
      KH=2-MI2
      JF=JMM2+MI2
      DO 401 J=2,JF
      JM1=J-MI2
      JP1=JM1+1
 401  DIV(I,J,L)=DIV(I,J,L) +(UFLUX(I+1,J)-UFLUX(I-1,J))*CDVL(J,KH)
     2         +(VFLUX(I,JP1)*CPH(JP1,KV)-VFLUX(I,JM1)*CPH(JM1,KV))*
     3          CDVP(J,KH)
      IF(FIRST) GO TO 420
C
C--------------DIVERGENCE DAMPING  -------------------------------------
      IF(NTSD.GT.NDDAMP) GO TO 471
      DO 472 L=1,LM
      DO 639 I=3,IMM2
      MI2=MOD(I,2)
      KV=1+MI2
      KH=2-MI2
      JF=JMM2+MI2
      DO 639 J=2,JF
      JM1=J-MI2
      JP1=JM1+1
 639  WF3(I,J)=CDVL(J,KH)*(U(I+1,J,L)-U(I-1,J,L))
     *       +CDVP(J,KH)*(V(I,JP1,L)*CPH(JP1,KV)-V(I,JM1,L)*CPH(JM1,KV))
C
      DO 472 I=4,IMM3
      MI2=MOD(I,2)
      KV=1+MI2
      JS=3-MI2
      DO 472 J=JS,JMM2
      JP1=J+MI2
      JM1=JP1-1
      U(I,J,L)=U(I,J,L)+DDAMPU(J,KV)*(WF3(I+1,J)-WF3(I-1,J))
 472  V(I,J,L)=V(I,J,L)+DDAMPV      *(WF3(I,JP1)-WF3(I,JM1))
C
C------------- TIME INTERPOLATION AND/OR SPACE EXTRAPOLATION OF U AND V
C              AT THE OUTER BOUNDARY
 471  CALL RTISEV
C
C------------- SPACE INTERPOLATION OF U AND V AT THE INNER BOUNDARY ----
      DO 147 I=3,IMM2,2
      DO 147 L=1,LM
      U(I,   1,L)=0.25*(U(I+1, 2,L)+U(I-1, 2,L)
     2            +U(I-1,   1,L)+U(I+1,   1,L))
      V(I,   1,L)=0.25*(V(I+1, 2,L)+V(I-1, 2,L)
     2            +V(I-1,   1,L)+V(I+1,   1,L))
      U(I,JMM1,L)=0.25*(U(I+1,JM,L)+U(I-1,JM,L)
     2            +U(I-1,JMM1,L)+U(I+1,JMM1,L))
147   V(I,JMM1,L)=0.25*(V(I+1,JM,L)+V(I-1,JM,L)

     2            +V(I-1,JMM1,L)+V(I+1,JMM1,L))
C
      DO 148 J=2,JMM1
      DO 148 L=1,LM
      U(   2,J,L)=0.25*(U( 3,J,L)+U(   1,J,L)
     2            +U(   1,J-1,L)+U( 3,J-1,L))
      V(   2,J,L)=0.25*(V( 3,J,L)+V(   1,J,L)
     2            +V(   1,J-1,L)+V( 3,J-1,L))
      U(IMM1,J,L)=0.25*(U(IM,J,L)+U(IMM2,J,L)
     2            +U(IMM2,J-1,L)+U(IM,J-1,L))
 148  V(IMM1,J,L)=0.25*(V(IM,J,L)+V(IMM2,J,L)
     2            +V(IMM2,J-1,L)+V(IM,J-1,L))
C
C--------------RUN CONTROL STATEMENTS-----------------------------------
      IF (MOD(NTSD,2*NTSPH).EQ.0) CALL RDGNS
      IF (MOD(NTSD,NCP).NE.0) GO TO 301
      NNTSD = -NTSD
      WRITE(NDISC) NTSD,DATE,IHRST,U,V
      WRITE(NDISC) T,PD,PDB,TB,UB,VB,PHIS
      WRITE(NDISC) NNTSD,DATE,IHRST
      BACKSPACE NDISC
      PRINT 2000, NTSD
C
 301  IF(NTSD.EQ.NTSTM) STOP
      GO TO 149
      END
      SUBROUTINE RDGNS
C   PARAMETER ONLY ACCEPTS CONSTANTS HENCE THE FOLLOWING:
                             P A R A M E T E R
     * IM= 51, JM=19, LM=5,
     1  KB=86, IMM1=50, IMM2=49, IMM4=47, 
     2  JMM1=18, JMM2=17, LMP1=6, LMM1=4
C    * KB=IM+2*JM-3, IMM1=IM-1, IMM2=IM-2, IMM4=IM-4, JMM1=JM-1,
C    *               JMM2=JM-2, LMP1=LM+1, LMM1=LM-1
                          D I M E N S I O N
     1 DAREA(2)
                             C O M M O N
     1 /ARRAYS/ U(IM,JM,LM), V(IM,JM,LM), T(IM,JM,LM), PD(IM,JM),
     * H1(JM,2), DSG(LM), RH1(JM,2), CPH(JM,2), ALPIFC(IM,JM,LMP1),
     * PHIFC(IM,JM,LMP1), PHIS(IM,JM), SG(LMP1), SGML(LM), DATE(4),
     * STPL(7), DISL(7), Z0L(7), WF4(IM,JM), WF5(IM,JM), WF6(IM,JM),
     * WF7(IM,JM), ISHDE(10), PDB(KB,2), TB(KB,LM,2), UB(KB,LM,2),
     * VB(KB,LM,2)
     2 /CONSTS/    NTSD,          H2,DPH,DLM,RTDLM,RTDPH
     * ,DISLP,DIST,DLMD,DPHD,IHRST,                LSTPLM,NTSPH,PT,R
     * ,Z0SLP,Z0T,   IOUT,NSHDE,                  TAU
      PRINT 201
 201  FORMAT('0NTSD    AKEU     AKEV     TAKE       APE        TE
     "APD     APENS*10**15')
      AKEU=0.
      AKEV=0.
      APE=0.
      ARKE=0.
      APENS=0.
      AVPD=0.
      ARPD=0.
      DS2=2.*H2*DPH
      FACDS1=2.*DLM
      DO 101 J=2,JMM1
      DAREA(1)=FACDS1*H1(J,1)*DS2
      DAREA(2)=FACDS1*H1(J,2)*DS2
      DO 102 I=3,IMM2
      MI2=MOD(I,2)
      JM1=J-MI2
      JP1=JM1+1
      KH=2-MI2
      KV=1+MI2
      IF(J.EQ.JMM1.AND.MI2.EQ.0) GO TO 103
      DARPD=DAREA(KH)
      IF(I.EQ.3.OR.I.EQ.IMM2) DARPD=0.5*DARPD
      IF(J.EQ.2.AND.MI2.NE.0) DARPD=0.5*DARPD
      IF(J.EQ.JMM1) DARPD=0.5*DARPD
      AVPD=AVPD+PD(I,J)*DARPD
      ARPD=ARPD+DARPD
      CPPDDA=1004.*PD(I,J)*DARPD
      SVSDS=0.
      DAP=0.
      DO 105 L=1,LM
      VORT      =RH1  (J,KH)*((V(I+1,J,L)-V(I-1,J,L))*RTDLM
     2           -(U(I,JP1,L)*CPH(JP1,KV)-U(I,JM1,L)*CPH(JM1,KV))*RTDPH)
      SVSDS=SVSDS+VORT**2*DSG(L)
 105  DAP=DAP+T(I,J,L)*DSG(L)
      APENS=APENS+SVSDS/PD(I,J)*DARPD
      APE=APE+DAP*CPPDDA
 103  IF(J.EQ.JMM1.AND.MI2.NE.0) GO TO 102
      DARKE=DAREA(KV)
      IF(I.EQ.3.OR.I.EQ.IMM2) DARKE=0.5*DARKE
      IF(J.EQ.2.AND.MI2.EQ.0) DARKE=0.5*DARKE
      IF(J.EQ.JMM1) DARKE=0.5*DARKE
      USQH=0.
      VSQH=0.
      DO 104 L=1,LM
      DSGHF=0.5*DSG(L)
      USQH=U(I,J,L)**2*DSGHF+USQH
 104  VSQH=V(I,J,L)**2*DSGHF+VSQH
      PDVP=0.25*(PD(I+1,J)+PD(I,J)+PD(I,J-1+2*MI2)+PD(I-1,J))
      DM=PDVP*DARKE
      AKEU=AKEU+USQH*DM
      AKEV=AKEV+VSQH*DM
      ARKE=ARKE+DARKE
 102  CONTINUE
 101  CONTINUE
      AVPD=AVPD/ARPD
      APENS=APENS/ARPD*10.**15   *0.5
      APE=APE/(ARPD*AVPD)
      TM=ARKE*AVPD
      AKEU=AKEU/TM
      AKEV=AKEV/TM
      TAKE=AKEU+AKEV
      TE=APE+TAKE
      PRINT 110, NTSD,AKEU,AKEV,TAKE,APE,TE,AVPD,APENS
 110  FORMAT (' ',I4,3F9.2,2F11.2,F10.0,F15.4)
      RETURN
      END
      SUBROUTINE ROUTPT
C   PARAMETER ONLY ACCEPTS CONSTANTS HENCE THE FOLLOWING:
                             P A R A M E T E R
     * IM= 51, JM=19, LM=5,
     1  KB=86, IMM1=50, IMM2=49, IMM4=47, 
     2  JMM1=18, JMM2=17, LMP1=6, LMM1=4
C    * KB=IM+2*JM-3, IMM1=IM-1, IMM2=IM-2, IMM4=IM-4, JMM1=JM-1,
C    *               JMM2=JM-2, LMP1=LM+1, LMM1=LM-1
                             C O M M O N
     1 /ARRAYS/ U(IM,JM,LM), V(IM,JM,LM), T(IM,JM,LM), PD(IM,JM),
     * H1(JM,2), DSG(LM), RH1(JM,2), CPH(JM,2), ALPIFC(IM,JM,LMP1),
     * PHIFC(IM,JM,LMP1), PHIS(IM,JM), SG(LMP1), SGML(LM), DATE(4),
     * STPL(7), DISL(7), Z0L(7), WF4(IM,JM), WF5(IM,JM), WF6(IM,JM),
     * WF7(IM,JM), ISHDE(10), PDB(KB,2), TB(KB,LM,2), UB(KB,LM,2),
     * VB(KB,LM,2)
     2 /CONSTS/    NTSD,          H2,DPH,DLM,RTDLM,RTDPH
     * ,DISLP,DIST,DLMD,DPHD,IHRST,                LSTPLM,NTSPH,PT,R
     * ,Z0SLP,Z0T,   IOUT,NSHDE,                  TAU
C
 2400 FORMAT(' AT',-2PF5.0,'MB',4A4,I2,' GMT + ',I2)
 2500 FORMAT(' PRESSURE AT THE MEAN SEA LEVEL',4A4,I2,' GMT + ',I2)
 2600 FORMAT(' TEMPERATURE AT',-2PF5.0,'MB',4A4,I2,' GMT + ',I2)
C
C--------------PREPARATION----------------------------------------------
 2901 IHR=NTSD/NTSPH
      DO 526 I=1,IM
      JF=JMM1+MOD(I,2)
      DO 526 J=1,JF
      PDIJ=PD(I,J)
      ALPIFC(I,J,LMP1)=ALOG(PT+PDIJ)
      PHIFC(I,J,LMP1)=PHIS(I,J)
      DO 526 IVI=1,LM
      L=LMP1-IVI
      ALPIJ=ALOG(PT+SG(L)*PDIJ)
      PHIFC(I,J,L)=PHIFC(I,J,L+1)+R*T(I,J,L)*(ALPIFC(I,J,L+1)-ALPIJ)
 526  ALPIFC(I,J,L)=ALPIJ
C
C--------------GEOPOTENTIAL AT STANDARD PRESSURE LEVELS-----------------
C--------------SEA LEVEL PRESSURE MAP-----------------------------------
      DO 375 I=1,IM
      JF=JMM1+MOD(I,2)
      DO 375 J=1,JF
      PDIJ=PD(I,J)
      IF(PHIS(I,J).LT.1.) GO TO 375
      TIJLM=T(I,J,LM)
      ALPLM1=ALPIFC(I,J,LMM1)
      ALPLP1=ALPIFC(I,J,LMP1)
      SLOP=2.*(TIJLM-T(I,J,LMM1))/(ALPLP1-ALPLM1)
      IF(SLOP.LT..5) GO TO 560
      TTT=TIJLM-SLOP*.5*(ALPLP1+ALPIFC(I,J,LM))
      PDIJ=(-TTT+SQRT(TTT*TTT+2.*SLOP*(PHIS(I,J)/R+(TTT+.5*SLOP*
     2        ALPLP1)*ALPLP1)))/SLOP
      GO TO 561
 560  PDIJ=ALPLP1+PHIS(I,J)/(R*TIJLM)
 561  PDIJ=EXP(PDIJ)-PT
 375  WF4 (I,J)=PDIJ+PT
      DO 528 I=1,IM
      MI2=MOD(I,2)
      JF=JMM1+MI2
      DO 528 J=1,JF
      IF(PHIS(I,J).LT.9810.) GO TO 529
      JM1=J-MI2
      JP1=JM1+1
      WF5(I,J)=.1464466*
     2 (WF4(I+1,JP1)+WF4(I-1,JP1)+WF4(I+1,JM1)+WF4(I-1,JM1))
     3        +.1035534*
     4 (WF4(I+2,J)+WF4(I-2,J)+WF4(I,J+1)+WF4(I,J-1))
      GO TO 528
 529  WF5 (I,J)=WF4 (I,J)
 528  CONTINUE
      CALL RSTMAP (WF5,  DLMD,DPHD,DISLP,Z0SLP)
      PRINT 2500, DATE,IHRST,IHR
C
      DO 520 LSTPL=1,LSTPLM
      REFALP=ALOG(STPL(LSTPL))
      TWOREF=2.*REFALP
      DO 521 I=1,IM
      JF=JMM1+MOD(I,2)
      DO 521 J=1,JF
      PDIJ=PD(I,J)
      DO 522 L=2,LMP1
      IF((REFALP-ALPIFC(I,J,L)).GT.0) GO TO 522
      NL1=L
      GO TO 523
 522  CONTINUE
      NL1=LMP1
 523  IF((TWOREF-ALPIFC(I,J,NL1)-ALPIFC(I,J,NL1-1)).LE.0.) NL1=NL1-1
      PHI1=PHIFC(I,J,NL1)
      ALP1=ALPIFC(I,J,NL1)
      IF(NL1.EQ.1) GO TO 524
      IF(NL1.EQ.LMP1) GO TO 525
      COFB=T(I,J,NL1)
      FAC=2.*ALOG(PT+SGML(NL1)*PDIJ)
      COFAHF=(COFB-T(I,J,NL1-1))/(ALPIFC(I,J,NL1+1)-ALPIFC(I,J,NL1-1))
      GO TO 527
 524  COFB=T(I,J,1)
      FAC=2.*ALOG(PT+SGML(1)*PDIJ)
      COFAHF=(T(I,J,2)-COFB)/(ALPIFC(I,J,3)-ALPIFC(I,J,1))
      GO TO 527
 525  COFB=T(I,J,LM)
      FAC=2.*ALOG(PT+SGML(LM)*PDIJ)
      COFAHF=(COFB-T(I,J,LMM1))/(ALPIFC(I,J,LMP1)-ALPIFC(I,J,LMM1))
 527  WF5(I,J)=COFB+COFAHF*(TWOREF-FAC)
      COFB=COFB-FAC*COFAHF
 521  WF4 (I,J)=PHI1-R*(REFALP-ALP1)*(COFAHF*(REFALP+ALP1)+COFB)
      DO 550 I=1,IM
      MI2=MOD(I,2)
      JF=JMM1+MI2
      DO 550 J=1,JF
      IF(PHIS(I,J).LT.9810.) GO TO 552
      JM1=J-MI2
      JP1=JM1+1
      WF6(I,J)=.125*(WF4(I+1,JP1)+WF4(I-1,JP1)+WF4(I+1,JM1)+WF4(I-1,JM1)
     2                +WF4(I,J)+WF4(I,J)+WF4(I,J)+WF4(I,J))
      WF7(I,J)=.25*(WF5(I+1,JP1)+WF5(I-1,JP1)+WF5(I+1,JM1)+WF5(I-1,JM1))
      GO TO 550
 552  WF6(I,J)=WF4(I,J)
      WF7(I,J)=WF5(I,J)
 550  CONTINUE
      STNIV=STPL(LSTPL)
      DIS=DISL(LSTPL)
      Z0=Z0L(LSTPL)
      CALL RSTMAP(WF6,DLMD,DPHD,DIS,Z0)
      PRINT 2400,STNIV,DATE,IHRST,IHR
      CALL RSTMAP(WF7,DLMD,DPHD,DIST,Z0T)
      PRINT 2600,STNIV,DATE,IHRST,IHR
 520  CONTINUE
C-----------------------------------------------------------------------
      IOUT=IOUT+1
      NSHDE=ISHDE(IOUT)*NTSPH + .5
C
      RETURN
      END
      SUBROUTINE RSTMAP (W,DLMD,DPHD,DIS,Z0)
C        THIS IS A STEREOGRAPHIC SHADING SUBROUTINE, DESIGNED FOR AN
C           8 LINES/INCH PRINTING
C        STEREOGRAPHIC MAP REFERENCE, E.G., GODSKE ET AL. 1957, P.132
C
                             P A R A M E T E R
     * IM= 51, JM=19
      REAL W(IM,JM)
      INTEGER   BLK, PRINT(141), CTBL(20)
      DATA BLK/'    '/, CTBL/'0   ','    ','1   ','    ','2   ','    ',
     "   '3   ','    ','4   ','    ','5   ','    ','6   ','    ','7   ',
     "   '    ','8   ','    ','9   ','    '/
C
C        MAP GEOMETRY CONSTANTS
C        UNIT OF DISTANCE OF THE IMMAGE SURFACE X1,Y1 AND X,Y SYSTEMS IS
C           THE SHORTER SIDE OF THE PRINTING CELL (2.54/10. CM)
C        X1P,Y1P ARE THE NORTH POLE COORDINATES
      X1P = 108.
      Y1P = 185.
C        DP30 IS THE DISTANCE, ON THE PHMI SYNOPTIC MAP, BETWEEN THE
C           NORTH POLE AND 30 DEG N LATITUDE CIRCLE, IN CM
      DP30 = 45.78
C        CENTRAL MAP MERIDIAN LONGITUDE, IN DEG
      CMLD = 2.5
C        WESTERN AND SOUTHERN BOUNDARIES OF THE MAP REGION, IN DEG
      WBD = -35.+DLMD
      SBD = 30.+DPHD
C
C        BASIC AND DERIVED CONSTANTS
      PI = 3.141593
      PIHF = PI/2.
      DTR = PI/180.
      RC = ( SIN(30.*DTR)/COS(30.*DTR))/(DP30/0.254)
      DLM = DLMD*DTR
      DPH = DPHD*DTR
      RDLM = 1./DLM
      RDPH = 1./DPH
      ALMWB = WBD*DTR
      APHSB = SBD*DTR
      ALMEB = ALMWB +(IM-3)*DLM
      APHNB = APHSB +2*(JM-2)*DPH
      AGM = -PIHF -CMLD*DTR
      ALM0 = (WBD-2.*DLMD)*DTR
      APH0 = (SBD-2.*DPHD)*DTR
C
      PRINT 199
  199 FORMAT (1H1)
      DO 101 LINE=1,200
      DO 102 K=1,141
  102 PRINT(K) = BLK
      X = LINE*1.25-X1P
C
      DO 103 K=4,141
      Y = K-Y1P
      R = SQRT(X**2+Y**2)
      A = ATAN2(Y,X)
C
C        ALM,APH IS LONGITUDE,LATITUDE OF THE PRINTING CELL, IN RADIAN
C           MEASURE
      ALM = A-AGM
      IF (ALM.LE.ALMWB) GO TO 101
      APH = PIHF -2.*ATAN(RC*R)
      IF (APH.GE.APHNB) GO TO 101
      IF (APH.GT.APHSB) GO TO 105
      IF (ALM.LT.ALMEB) GO TO 103
      GO TO 106
  105 IF (ALM.GE.ALMEB) GO TO 101
C
C        THE PRINTING CELL HAS BEEN FOUND TO BE INSIDE THE MAP REGION
C           BOUNDARIES
C        X2,Y2 IS A GRID COORDINATE SYSTEM, WITH DLM,DPH AS LENGTH UNITS
      X2 = (ALM-ALM0)*RDLM
      Y2 = (APH-APH0)*RDPH
C
C        I2,J2 ARE THE X2,Y2 COORDINATES OF THE GRID POINT AT THE LEFT
C           LOWER CORNER OF THE GRID SQUARE
      I2 = INT(X2)
      J2 = INT(Y2)
C
C        NOW CHANGE I2,J2, IF NECESSARY, TO OBTAIN THE X2,Y2 COORDINATES
C           OF THE H GRID POINT TO BE USED AS THE ORIGIN OF THE
C           INTERPOLATION P,Q COORDINATE SYSTEM
      IF (MOD(I2+J2,2).EQ.0) GO TO 109
      IF (X2-I2+Y2-J2.GT.1.) GO TO 107
      J2 = J2-1
      GO TO 108
  107 I2 = I2+1
      GO TO 108
  109 IF (Y2-J2.GT.X2-I2) GO TO 108
      I2 = I2+1
      J2 = J2-1
C
C        I00,J00 ARE THE I,J SUBSCRIPTS OF THE I2,J2 H GRID POINT
  108 I00 = I2
      MI2 = MOD(I00,2)
      J00 = (J2+MI2)/2
C
C        CALCULATE REMAINING PARAMETERS NEEDED FOR BILINEAR
C           INTERPOLATION (NOTATION IS USED FOLLOWING MESINGER, JAS
C           1965, EQ.16)
      J00P1 = J00+1
      I10 = I00+1
      J10 = J00P1-MI2
      I01 = I00-1
      J01 = J10
      I11 = I00
      J11 = J00P1
C
      X2MI2 = X2-I2
      Y2MJ2 = Y2-J2
      P = 0.5*( X2MI2+Y2MJ2)
      Q = 0.5*(-X2MI2+Y2MJ2)
      PQ = P*Q
C
C        INTERPOLATION AND PRINTING
      Z = (1.-P-Q+PQ)*W(I00,J00) +(P-PQ)*W(I10,J10) +(Q-PQ)*W(I01,J01)
     "   +PQ*W(I11,J11)
      ZR = Z-Z0
      IF (ZR.LT.0.) ZR=-ZR
      IQ = ZR/DIS
      MIQ = MOD(IQ,20)
      N = MIQ+1
      IF (Z.LT.Z0) N=20-MIQ
      PRINT(K) = CTBL(N)
  103 CONTINUE
C
  101 PRINT 104, (PRINT(K),K=4,134)
  104 FORMAT (1H ,131A1)
  106 PRINT 110
  110 FORMAT (1H ,//)
      RETURN
      END
      SUBROUTINE RTIPDT
C   PARAMETER ONLY ACCEPTS CONSTANTS HENCE THE FOLLOWING:
                             P A R A M E T E R
     * IM= 51, JM=19, LM=5,
     1  KB=86, IMM1=50, IMM2=49, IMM4=47, 
     2  JMM1=18, JMM2=17, LMP1=6, LMM1=4
C    * KB=IM+2*JM-3, IMM1=IM-1, IMM2=IM-2, IMM4=IM-4, JMM1=JM-1,
C    *               JMM2=JM-2, LMP1=LM+1, LMM1=LM-1
                             C O M M O N
     1 /ARRAYS/ U(IM,JM,LM), V(IM,JM,LM), T(IM,JM,LM), PD(IM,JM),
     * H1(JM,2), DSG(LM), RH1(JM,2), CPH(JM,2), ALPIFC(IM,JM,LMP1),
     * PHIFC(IM,JM,LMP1), PHIS(IM,JM), SG(LMP1), SGML(LM), DATE(4),
     * STPL(7), DISL(7), Z0L(7), WF4(IM,JM), WF5(IM,JM), WF6(IM,JM),
     * WF7(IM,JM), ISHDE(10), PDB(KB,2), TB(KB,LM,2), UB(KB,LM,2),
     * VB(KB,LM,2)
     2 /CONSTS/    NTSD,          H2,DPH,DLM,RTDLM,RTDPH
     * ,DISLP,DIST,DLMD,DPHD,IHRST,                LSTPLM,NTSPH,PT,R
     * ,Z0SLP,Z0T,   IOUT,NSHDE,                  TAU
C
      N = 1
      DO 327 I=1,IM,2
      PD(I,1)=PDB(N,1)+PDB(N,2)*TAU
 327  N=N+1
      DO 328 I=1,IM,2
      PD(I,JM)=PDB(N,1)+PDB(N,2)*TAU
 328  N=N+1
      DO 329 J=2,JMM1
      PD(1,J)=PDB(N,1)+PDB(N,2)*TAU
 329  N=N+1
      DO 330 J=2,JMM1
      PD(IM,J)=PDB(N,1)+PDB(N,2)*TAU
 330  N=N+1
C
      DO 144 L=1,LM
      N = 1
      DO 141 I=1,IM,2
      T(I, 1,L) = TB(N,L,1)+TB(N,L,2)*TAU
 141  N = N+1
      DO 142 I=1,IM,2
      T(I,JM,L) = TB(N,L,1)+TB(N,L,2)*TAU
 142  N = N+1
      DO 143 J=2,JMM1
      T( 1,J,L) = TB(N,L,1)+TB(N,L,2)*TAU
 143  N = N+1
      DO 144 J=2,JMM1
      T(IM,J,L) = TB(N,L,1)+TB(N,L,2)*TAU
 144  N = N+1
C
      RETURN
      END
      SUBROUTINE RTISEV
C   PARAMETER ONLY ACCEPTS CONSTANTS HENCE THE FOLLOWING:
                             P A R A M E T E R
     * IM= 51, JM=19, LM=5,
     1  KB=86, IMM1=50, IMM2=49, IMM4=47, 
     2  JMM1=18, JMM2=17, LMP1=6, LMM1=4
C    * KB=IM+2*JM-3, IMM1=IM-1, IMM2=IM-2, IMM4=IM-4, JMM1=JM-1,
C    *               JMM2=JM-2, LMP1=LM+1, LMM1=LM-1
                             C O M M O N
     1 /ARRAYS/ U(IM,JM,LM), V(IM,JM,LM), T(IM,JM,LM), PD(IM,JM),
     * H1(JM,2), DSG(LM), RH1(JM,2), CPH(JM,2), ALPIFC(IM,JM,LMP1),
     * PHIFC(IM,JM,LMP1), PHIS(IM,JM), SG(LMP1), SGML(LM), DATE(4),
     * STPL(7), DISL(7), Z0L(7), WF4(IM,JM), WF5(IM,JM), WF6(IM,JM),
     * WF7(IM,JM), ISHDE(10), PDB(KB,2), TB(KB,LM,2), UB(KB,LM,2),
     * VB(KB,LM,2)
     2 /CONSTS/    NTSD,          H2,DPH,DLM,RTDLM,RTDPH
     * ,DISLP,DIST,DLMD,DPHD,IHRST,                LSTPLM,NTSPH,PT,R
     * ,Z0SLP,Z0T,   IOUT,NSHDE,                  TAU
C
      DO 140 L=1,LM
      N = 1
      DO 137 I=2,IMM1,2
      U(I, 1,L) = UB(N,L,1)+UB(N,L,2)*TAU
      V(I, 1,L) = VB(N,L,1)+VB(N,L,2)*TAU
      IF (V(I, 1,L).LT.0.)
     "                  U(I, 1,L) =U(I,   2,L)+(U(I,   2,L)-U(I,   3,L))
 137  N = N+1
      DO 138 I=2,IMM1,2
      U(I,JM,L) = UB(N,L,1)+UB(N,L,2)*TAU
      V(I,JM,L) = VB(N,L,1)+VB(N,L,2)*TAU
      IF (V(I,JM,L).GT.0.)
     "                  U(I,JM,L) =U(I,JMM1,L)+(U(I,JMM1,L)-U(I,JMM2,L))
 138  N = N+1
      DO 139 J=1,JMM1
      U( 1,J,L) = UB(N,L,1)+UB(N,L,2)*TAU
      V( 1,J,L) = VB(N,L,1)+VB(N,L,2)*TAU
      IF (U( 1,J,L).LT.0.)
     "                  V( 1,J,L) =V(   3,J,L)+(V(   3,J,L)-V(   5,J,L))
 139  N = N+1
      DO 140 J=1,JMM1
      U(IM,J,L) = UB(N,L,1)+UB(N,L,2)*TAU
      V(IM,J,L) = VB(N,L,1)+VB(N,L,2)*TAU
      IF (U(IM,J,L).GT.0.)
     "                  V(IM,J,L) =V(IMM2,J,L)+(V(IMM2,J,L)-V(IMM4,J,L))
 140  N = N+1
C
      RETURN
      END
