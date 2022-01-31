        PROGRAM CDIFF
C
C********************************************************************
C
C       COMPARE TWO OUTPUT FILES MMDDHH.FCn AND CALCULATE THE
C       STATISTICS OF THE DIFFERENCES BETWEEN THEM.
C      
C       Write out the DIFFERENCE fields to a file.
C
C********************************************************************
C
                       P A R A M E T E R 
     X  IM=40,JM=26
     X ,IMM1=im-1,IMM2=im-2,IMM3=im-3,IMM4=im-4,
     X  JMM1=jm-1,JMM2=jm-2
C
			DIMENSION
     1  IDATE0(4),FI0(IM,JM),UZ0(IM,JM),VZ0(IM,JM),
     2  IDATE1(4),FI1(IM,JM),UZ1(IM,JM),VZ1(IM,JM),
     3  FIDIF(IM,JM),UZDIF(IM,JM),VZDIF(IM,JM)
C
        character*50 FILE1, FILE2, FILE3
C       character*10 STRING
C
        LOGICAL CORNER, PDIFF, FASCI1, FASCI2, FASCI3
        DATA CORNER /.FALSE./
        DATA PDIFF  /.FALSE./
C
C********************************************************************
C
C       OPEN THE CARD FILE.
C     OPEN (UNIT=55,form='formatted',status='old',
C    1   FILE='cdiff.cds')
C
C       OPEN THE LINEPRINTER FILE.
      OPEN (UNIT=66,form='formatted',FILE='cdiff.lpt')
C
C--------------------------------------------------------------------
C
C       READ (55,*) STRING, CORNER
C       READ (55,*) STRING, PDIFF
C
C--------------------------------------------------------------------
        TYPE 9001
 9001   FORMAT('  TYPE IN NAME OF  FIRST FILE '$)
        ACCEPT 9101, FILE1
 9101   FORMAT(A)
        WRITE (66,9701) FILE1
 9701   FORMAT(/' FILE 1: ',A)
        TYPE 90001
90001   FORMAT('  FORMAT OF  FIRST FILE ASCII (T/F) '$)
        ACCEPT 91001, FASCI1
91001   FORMAT(L1)
C
        TYPE 9002
 9002   FORMAT('  TYPE IN NAME OF SECOND FILE '$)
        ACCEPT 9101, FILE2
        WRITE (66,9702) FILE2
 9702   FORMAT(' FILE 2: ',A/)
        TYPE 90002
90002   FORMAT('  FORMAT OF SECOND FILE ASCII (T/F) '$)
        ACCEPT 91001, FASCI2
C
        TYPE 9003
 9003   FORMAT('  TYPE IN NAME OF OUTPUT FILE '$)
        ACCEPT 9101, FILE3
           IF(FILE3.eq.'') GO TO 9900
        WRITE (66,9703) FILE3
 9703   FORMAT(' FILE 3: ',A/)
        TYPE 90003
90003   FORMAT('  FORMAT OF OUTPUT FILE ASCII (T/F) '$)
        ACCEPT 91001, FASCI3
c
 9900   CONTINUE
C--------------------------------------------------------------------
C
        GRAV = 9.81
C
      IF ( FASCI1 ) THEN 
	OPEN(UNIT=24,form='formatted',status='old',FILE=FILE1)
	READ(24,*)  NTSD0,IDATE0,IHRST0,FI0,UZ0,VZ0
	CLOSE(UNIT=24)
	TYPE 30 , NTSD0,IDATE0,IHRST0
      ELSE
	OPEN(UNIT=24,form='unformatted',status='old',FILE=FILE1)
	READ(24)  NTSD0,IDATE0,IHRST0,FI0,UZ0,VZ0
	CLOSE(UNIT=24)
	TYPE 30 , NTSD0,IDATE0,IHRST0
      ENDIF
30    FORMAT( ' NTSD0=',I4,' IDATE0=',4I4,' IHRST0=',I4)
C
      IF ( FASCI2 ) THEN 
	OPEN(UNIT=25,form='formatted',status='old',FILE=FILE2)
	READ(25,*)  NTSD1,IDATE1,IHRST1,FI1,UZ1,VZ1
	CLOSE(UNIT=25)
	TYPE 40 , NTSD1,IDATE1,IHRST1
      ELSE
	OPEN(UNIT=25,form='unformatted',status='old',FILE=FILE2)
	READ(25)  NTSD1,IDATE1,IHRST1,FI1,UZ1,VZ1
	CLOSE(UNIT=25)
	TYPE 40 , NTSD1,IDATE1,IHRST1
      ENDIF
40    FORMAT( ' NTSD1=',I4,' IDATE1=',4I4,' IHRST1=',I4)
C
C--------------------------------------------------------------------
C
	NCOUNT=0
	SMUDIF=0.
	SMUSQ=0.
	SMVDIF=0.
	SMVSQ=0.
	SMFDIF=0.
	SMFSQ=0.
C
        FMAX = 0.
        UMAX = 0.
        VMAX = 0.
        FMIN = 0.
        UMIN = 0.
        VMIN = 0.
C
        NBSKIP = 0
        ISTART = 1  + NBSKIP
        ISTOP  = IM - NBSKIP
        JSTART = 1  + NBSKIP
        JSTOP  = JM - NBSKIP
C
        DO 100 I=ISTART,ISTOP
        DO 100 J=JSTART,JSTOP
            UDIF=UZ0(I,J)-UZ1(I,J)
            UZDIF(I,J)=UDIF
            VDIF=VZ0(I,J)-VZ1(I,J)
            VZDIF(I,J)=VDIF
            FDIF      =FI0(I,J)-FI1(I,J)
            FIDIF(I,J)=FDIF
            UDIFSQ=UDIF*UDIF
            VDIFSQ=VDIF*VDIF
            FDIFSQ=FDIF*FDIF
            NCOUNT=NCOUNT+1
            SMUSQ=SMUSQ+UDIFSQ
            SMVSQ=SMVSQ+VDIFSQ
            SMFSQ=SMFSQ+FDIFSQ
            SMFDIF=SMFDIF+FDIF
            SMUDIF=SMUDIF+UDIF
            SMVDIF=SMVDIF+VDIF
C
            IF((UDIF).GT.UMAX) THEN
              IUMAX = I
              JUMAX = J
              UMAX = (UDIF)
            END IF
            IF((VDIF).GT.VMAX) THEN
              IVMAX = I
              JVMAX = J
              VMAX = (VDIF)
            END IF
            IF((FDIF).GT.FMAX) THEN
              IFMAX = I
              JFMAX = J
              FMAX = (FDIF)
            END IF
C
            IF((UDIF).LT.UMIN) THEN
              IUMIN = I
              JUMIN = J
              UMIN = (UDIF)
            END IF
            IF((VDIF).LT.VMIN) THEN
              IVMIN = I
              JVMIN = J
              VMIN = (VDIF)
            END IF
            IF((FDIF).LT.FMIN) THEN
              IFMIN = I
              JFMIN = J
              FMIN = (FDIF)
            END IF
100	CONTINUE
C
	ROOTF=SQRT(SMFSQ/NCOUNT)
	ROOTu=SQRT(SMuSQ/NCOUNT)
	ROOTv=SQRT(SMvSQ/NCOUNT)
	ROOTuv=SQRT((SMuSQ+SMvSQ)/NCOUNT)
C
	STDEVU=SQRT((NCOUNT*SMUSQ-SMUDIF*SMUDIF)/(NCOUNT*NCOUNT))
	STDEVV=SQRT((NCOUNT*SMVSQ-SMVDIF*SMVDIF)/(NCOUNT*NCOUNT))
	STDEVF=SQRT((NCOUNT*SMFSQ-SMFDIF*SMFDIF)/(NCOUNT*NCOUNT))
C
  	ROOTF=ROOTF/GRAV
  	STDEVF=STDEVF/GRAV
        DFMEAN = SMFDIF/NCOUNT/GRAV
        FMAX=FMAX/GRAV
        FMIN=FMIN/GRAV
C
C--------------------------------------------------------------------
C
	TYPE     234, ROOTF,STDEVF,DFMEAN
        WRITE (66,234) ROOTF,STDEVF,DFMEAN
 234  FORMAT(' FI:  RMS DIFF. =',F8.2,' STANDARD DEV.=',
     x  F8.2,' MEAN DIFF.=',F8.2,' METRES')
C
	TYPE     235, ROOTuv,ROOTu,ROOTv
        WRITE (66,235) ROOTuv,ROOTu,ROOTv
 235  FORMAT('  V:  RMS DIFF. =',F8.2,' u,v: RMS DIFF.=',2F8.2,' M/S '/)
C
C       TYPE 237,STDEVU,STDEVV
237	FORMAT( ' U, V: STANDARD DEVS.=',2F10.2,' METRES/SEC')
C
	TYPE      2391, IFMAX,JFMAX, FMAX, IFMIN,JFMIN, FMIN
	TYPE      2392, IUMAX,JUMAX, UMAX, IUMIN,JUMIN, UMIN
	TYPE      2393, IVMAX,JVMAX, VMAX, IVMIN,JVMIN, VMIN
	WRITE (66, 2391) IFMAX,JFMAX, FMAX, IFMIN,JFMIN, FMIN
	WRITE (66, 2392) IUMAX,JUMAX, UMAX, IUMIN,JUMIN, UMIN
	WRITE (66, 2393) IVMAX,JVMAX, VMAX, IVMIN,JVMIN, VMIN
C
 2391   FORMAT( ' MAX AND MIN DELTA-FI :',2(2I4,6X,F10.4))
 2392   FORMAT( ' MAX AND MIN DELTA-U  :',2(2I4,6X,F10.4))
 2393   FORMAT( ' MAX AND MIN DELTA-V  :',2(2I4,6X,F10.4))
C
CCC            STOP ' STOP 777 '
C--------------------------------------------------------------------
C
	NTSD=48
	IHRST=0
	TYPE 55 , NTSD,IDATE0,IHRST
55	FORMAT( ' NTSD=',I4,' IDATE=',4I4,' IHRST=',I4)
C
        GO TO 5678
c       OPEN(UNIT=26,form='unformatted',status='new',
c    1       FILE='CDIFF.OUT')
c       WRITE(26)  NTSD,IDATE0,IHRST,FIDIF,UZDIF,VZDIF
c       CLOSE(UNIT=26)
 5678   CONTINUE
C
           IF(FILE3.eq.'') GO TO 9950
      IF ( FASCI3 ) THEN 
CCC     OPEN(UNIT=26,form='formatted',status='new',FILE=FILE3)
	OPEN(UNIT=26,form='formatted',             FILE=FILE3)
        WRITE (26,*)  NTSD1,IDATE1,IHRST1,FIDIF,UZDIF,VZDIF
        CLOSE(UNIT=26)
        TYPE 50 , NTSD1,IDATE1,IHRST1
      ELSE
CCC    	OPEN(UNIT=26,form='unformatted',status='new',FILE=FILE3)
	OPEN(UNIT=26,form='unformatted',             FILE=FILE3)
	WRITE (26)  NTSD1,IDATE1,IHRST1,FIDIF,UZDIF,VZDIF
	CLOSE(UNIT=26)
	TYPE 50 , NTSD1,IDATE1,IHRST1
      ENDIF
50    FORMAT( ' NTSD3=',I4,' IDATE3=',4I4,' IHRST3=',I4)
c
 9950 CONTINUE
C
C--------------------------------------------------------------------
C
      IF(CORNER) CALL KORNER(FI1  ,IM,JM, 6,'  FI REF  ')
      IF(PDIFF ) CALL KORNER(FIDIF,IM,JM, 6,'  FI DIF  ')
      IF(CORNER) CALL KORNER(UZ1  ,IM,JM, 6,'  U  REF  ')
      IF(PDIFF ) CALL KORNER(UZDIF,IM,JM, 6,'  U DIF   ')
      IF(CORNER) CALL KORNER(VZ1  ,IM,JM, 6,'  V  REF  ')
      IF(PDIFF ) CALL KORNER(VZDIF,IM,JM, 6,'  V DIF   ')
C
C--------------------------------------------------------------------
C
        STOP
	END
C--------------------------------------------------------------------
      SUBROUTINE KORNER(FIELD,KIM,KJM,KPRINT,MESAGE)
C --- [C-GRID ANALYSIS * J.HAMILTON * 21-MAY-1984*MODD PLYNCH DEC 86]
C
C *** PRINT OUT THE FIELD WITH STAGGERING AS SELECTED ******************
C
C-----------------------------------------------------------------------
      DIMENSION FIELD(KIM,KJM)
      DIMENSION XFIELD(10)
      DIMENSION XINDEX(10)
      character*10 MESAGE
C-----------------------------------------------------------------------
C
C#PRINT HEADER (IF REQUIRED)
      IF(KPRINT.LE.0)  RETURN
      WRITE(KPRINT,10) MESAGE
   10 FORMAT(1H1//////48X,' KORNER DUMP OF ',A10//)
C
C#LOOP FOR THE TOP AND BOTTOM LINES
      DO 65 JLOOP=1,2
C
C#THE TOP LINES
         IYMAX=KJM
         IYMIN=IYMAX-4
C
C#THE BOTTOM LINES
         IF(JLOOP.EQ.2) IYMIN=1
         IF(JLOOP.EQ.2) IYMAX=5
C
C#PRINT THE LINES IN REVERSE
         DO 50 JREV=IYMAX,IYMIN,-1
            JY=JREV
C#LEFT HAND SIDE
            IXMIN=1
            IXMAX=5
            INDEX=0
            DO 15 JX=IXMIN,IXMAX
               INDEX=INDEX+1
               XINDEX(INDEX)=FLOAT(JX) + 0.01*FLOAT(JY)
               XFIELD(INDEX)=FIELD(JX,JY)
   15          CONTINUE
C
C#RIGHT HAND SIDE
            IXMIN=KIM-4
            IXMAX=KIM
            DO 20 JX=IXMIN,IXMAX
               INDEX=INDEX+1
               XINDEX(INDEX)=FLOAT(JX) + 0.01*FLOAT(JY)
               XFIELD(INDEX)=FIELD(JX,JY)
   20          CONTINUE
C
C#NOW PRINT OUT THE LINES
            WRITE(KPRINT,25)XINDEX
   25       FORMAT(/4X,5F11.2,'      ',5F11.2)
            WRITE(KPRINT,30)XFIELD
   30       FORMAT(4X,1P5E11.3,'   ---',1P5E11.3)
C#END OF THE INNER LOOP
   50       CONTINUE
         IF((JLOOP.EQ.1)) WRITE(KPRINT,60)
   60    FORMAT(/4X,5('        ---'),'   ---',5('        ---'))
C
C#END OF MAIN LOOP
   65    CONTINUE
C
      RETURN
      END
C--------------------------------------------------------------------
