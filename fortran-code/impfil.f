c***************************************************************
c***************************************************************
C
        SUBROUTINE IMPFIL(FIN,FOUT,IM,JM,IMAX,JMAX,
     +                    IWEST,IEAST,JSOUTH,JNORTH,EPSILON)
C
C       IMPLICIT FILTERING OF A TWO DIMENSIONAL FIELD USING
C       THE SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER
C       DESCRIBED IN RAYMOND, M.W.R., 116, 2132-2141.
C
C     FIN (IMAX,JMAX): ----  INPUT  FIELD
C     FOUT(IMAX,JMAX): ----  OUTPUT FIELD
C
C     IM, JM:                ACTIVE DIMENSIONS OF FIN, FOUT
C     IMAX, JMAX:      ----  MAXIMUM DIMENSIONS OF FIN, FOUT
C
C     IWEST, IEAST,    ----  LIMITS OF AREA OVER WHICH THE
C     JSOUTH,JNORTH:   ----  FILTER IS TO BE APPLIED (THESE
C                            LINES ARE UNCHANGED BY THE FILTER)
C
C     EPSILON:         ----  FILTER PARAMETER TO DETERMINE CUTOFF
C                            WHERE FILTER RESPONSE IS ONE-HALF.
C
        REAL FIN(IMAX,JMAX),FOUT(IMAX,JMAX)
        DOUBLE PRECISION FDUBL(250), EDUBL
C
C       CHECK FOR ADEQUATE WORK-SPACE DIMENSIONS
        IF (IMAX.GT.250 .OR. JMAX.GT.250) THEN
           WRITE(6,9999)
           STOP
 9999      FORMAT(' DIMENSIONS TOO LARGE IN IMPFIL ')
        END IF
C
CCC     DETERMINE CUTOFF FREQUENCY/WAVENUMBER.
CCC     PI = 2*ASIN(1.)
CCC     THETAC = 2*ATAN(EPSILON**(-1./6.))
CCC     XLENC = 2*PI/THETAC 
CCC     TYPE *,' XLENC, THETAC, EPSILON ',XLENC,THETAC,EPSILON
C
C       PARAMETER IS DOUBLE PRECISION IN LOWPAS.
        EDUBL = EPSILON
C
C       FIRST, COPY THE INPUT TO THE OUTPUT. 
        DO J=1,JM
        DO I=1,IM
            FOUT(I,J) = FIN(I,J)
        ENDDO
        ENDDO
C
C       NEXT, FILTER THE ROWS IN THE X-DIRECTION.
        NLEN = IEAST-IWEST+1
        DO J=JSOUTH+1,JNORTH-1
          DO I=IWEST,IEAST
            FDUBL(I-IWEST+1) = FOUT(I,J)
          ENDDO
          CALL LOWPAS(FDUBL,NLEN,EDUBL)
          DO I=IWEST,IEAST
            FOUT(I,J) = FDUBL(I-IWEST+1)
          ENDDO
        ENDDO
C
C       FINALLY, FILTER THE ROWS IN THE Y-DIRECTION.
        NLEN = JNORTH-JSOUTH+1
        DO I=IWEST+1,IEAST-1
          DO J=JSOUTH,JNORTH
            FDUBL(J-JSOUTH+1) = FOUT(I,J)
          ENDDO
          CALL LOWPAS(FDUBL,NLEN,EDUBL)
          DO J=JSOUTH,JNORTH
            FOUT(I,J) = FDUBL(J-JSOUTH+1)
          ENDDO
        ENDDO
C
        RETURN
        END
C
c***************************************************************
c***************************************************************
C
        SUBROUTINE LOWPAS(XY,N,EPS)
C
C       SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER
C       (RAYMOND, MWR, 116, 2132-2141)
C
C*************************************************************
C***     THIS CODE IS COPIED FROM A LISTING PROVIDED       ***
C***     BY WILLIAM H RAYMOND. SOME NOTATIONAL CHANGES     ***
C***     HAVE BEEN MADE IN THE ROUTINE LOWPAS. THE         ***
C***     ROUTINE INVLOW HAS BEEN COPIED ALMOST VERBATIM.   ***
C*************************************************************
C
C       XY     UNFILTERED VALUES ON INPUT
C              FILTERED VALUES ON OUTPUT.
C       N      NUMBER OF VALUES.
C       EPS    FILTER PARAMETER 
C              (DETERMINES CUTOFF)
C
C       LOCAL VARIABLE NSAVE KEPT AND COMPARED TO
C       INCOMING VALUE N. (INITIALLY SET TO ZERO).
C       IF VALUE IS UNCHANGED, ISKIP IS SET TO .TRUE.
C
C---------------------------------------------------------------
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION XY(N),RHS(481),XDASH(481)
        LOGICAL ISKIP
        DATA NSAVE /0/
        DATA ISKIP /.FALSE./
C
        SAVE
C
C---------------------------------------------------------------
C
        IF ( N .GT. 481 ) STOP ' Increase N in LOWPAS '
C
C---------------------------------------------------------------
C
C       SKIP FIRST PART OF GAUSSIAN ELIMINATION ON REPEAT
C       CALLS IF FILTER LENGTH REMAINS UNCHANGED.
        IF (N.EQ.NSAVE) ISKIP=.TRUE.
        IF (N.NE.NSAVE) ISKIP=.FALSE.
        NSAVE = N
C
        NM1 = N-1
        NM2 = N-2
        NM3 = N-3
        NM4 = N-4
C        
        RHS(1) = 0
        RHS(N) = 0
        RHS(  2) = EPS*(XY(  1)-2.D0*XY(  2)+XY(  3))
        RHS(NM1) = EPS*(XY(NM2)-2.D0*XY(NM1)+XY(  N)) 
        RHS(  3) = EPS*(-1.D0*(XY(  1)+XY(  5))
     +                  +4.D0*(XY(  2)+XY(  4))
     +                  -6.D0* XY(  3)         )
        RHS(NM2) = EPS*(-1.D0*(XY(  N)+XY(NM4))
     +                  +4.D0*(XY(NM1)+XY(NM3))
     +                  -6.D0* XY(NM2)         )
        DO 1000 J=4,NM3
        RHS(J) = EPS*(       (XY(J-3)+XY(J+3))
     +                - 6.D0*(XY(J-2)+XY(J+2))
     +                +15.D0*(XY(J-1)+XY(J+1))
     +                -20.D0* XY(  J)         )
 1000   CONTINUE
C
C       SOLVE THE LINEAR SYSTEM FOR XDASH
        CALL INVLOW(RHS,N,XDASH,EPS,ISKIP)
C
C       ADD CORRECTION TO GET FILTERED VALUES.
        DO 2000 J=1,N
        XY(J) = XY(J) + XDASH(J)
 2000   CONTINUE
C
        RETURN
        END
C
C---------------------------------------------------------------
C---------------------------------------------------------------
C
        SUBROUTINE INVLOW(BB,N,XANS,EP,ISKIP)
C
C       GAUSSIAN ELIMINATION FOR LOW-PASS FILTER.
C
C       SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER.
C       (REF: WILLIAM H RAYMOND, MWR, 116, 2132-2124)
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL ISKIP
C
        PARAMETER (NNMAX=481)
C
        DIMENSION A(NNMAX),B(NNMAX),C(NNMAX),D(NNMAX),E(NNMAX),
     +            DELTA(NNMAX),BETA(NNMAX),W(NNMAX),GAM(NNMAX),
     +            H(NNMAX),XANS(NNMAX),BB(NNMAX),PI(NNMAX),
     +            AP(NNMAX),F(NNMAX),Z(NNMAX)
C
        SAVE
C
C---------------------------------------------------------------
C
C       SKIP INITIALIZATION OF MATRIX ON REPEAT CALLS.
        IF ( ISKIP ) GO TO 100
C
C       INITIALIZE THE MATRIX
C
        DO 10 I=4,N-3
        Z(I) = 1.D0-EP
        A(I) = 6.D0*(1.D0+EP)
        B(I) = 15.D0*(1.D0-EP)
        C(I) = 20.D0*(1.D0+EP)
        D(I) = B(I)
        E(I) = A(I)
        F(I) = Z(I)
   10   CONTINUE
C
        Z(1) = 0
        Z(2) = 0
        Z(3) = 0
C
        A(1) = 0
        A(2) = 0
        A(3) = 1.D0+EP
C
        B(1) = 0
        B(2) = 1.D0-EP
        B(3) = 4.D0*(1.D0-EP)
C
        C(1) = 1.D0
        C(2) = 2.D0*(1.D0+EP)
        C(3) = 6.D0*(1.D0+EP)
C
        D(1) = 0
        D(2) = 1.D0-EP
        D(3) = 4.D0*(1.D0-EP)
C
        E(1) = 0
        E(2) = 0
        E(3) = 1.D0+EP
C
        F(1) = 0
        F(2) = 0
        F(3) = 0
C
C
        Z(N-2) = 0 
        Z(N-1) = 0
        Z(N) = 0
C
        A(N-2) = 1.D0+EP
        A(N-1) = 0
        A(N) = 0
C
        B(N-2) = 4.D0*(1.D0-EP)
        B(N-1) = 1.D0-EP
        B(N) = 0
C
        C(N-2) = 6.D0*(1.D0+EP)
        C(N-1) = 2.D0*(1.D0+EP)
        C(N) = 1.D0
C
        D(N-2) = 4.D0*(1.D0-EP)
        D(N-1) = 1.D0-EP
        D(N) = 0
C
        E(N-2) = 1.D0+EP
        E(N-1) = 0
        E(N) = 0
C
        F(N-2) = 0 
        F(N-1) = 0
        F(N) = 0
C
C       Step One.
        BETA(1) = D(1)/C(1)
        DELTA(2) = B(2)
        W(1) = C(1)
        PI(1) = F(1)/W(1)
        AP(1) = 0
        AP(2) = 0
        AP(3) = A(3)
        W(2) = C(2)-DELTA(2)*BETA(1)
        GAM(1) = E(1)/C(1)
        BETA(2) = (D(2)-DELTA(2)*GAM(1))/W(2)
        GAM(2) = (E(2)-PI(1)*DELTA(2))/W(2)
        PI(2) = F(2)/W(2)
        DELTA(3) = (B(3)-AP(3)*BETA(1))
        W(3) = C(3)-DELTA(3)*BETA(2)-AP(3)*GAM(1)
        BETA(3) = (D(3)-AP(3)*PI(1)-DELTA(3)*GAM(2))/W(3)
        GAM(3) = (E(3)-DELTA(3)*PI(2))/W(3)
        PI(3) = F(3)/W(3)
C
C       Step Two
        DO 20 I=4,N
        AP(I) = A(I)-Z(I)*BETA(I-3)
        DELTA(I) = B(I)-AP(I)*BETA(I-2)-Z(I)*GAM(I-3)
        W(I) = C(I)-AP(I)*GAM(I-2)-DELTA(I)*BETA(I-1)
     +             -Z(I)*PI(I-3)
        BETA(I) = (D(I)-AP(I)*PI(I-2)-DELTA(I)*GAM(I-1))/W(I)
        GAM(I) = (E(I)-DELTA(I)*PI(I-1))/W(I)
        PI(I) = F(I)/W(I)
   20   CONTINUE
C
  100   CONTINUE
C
C       Step Three
        H(1) = BB(1)/W(1)
        H(2) = (BB(2)-DELTA(2)*H(1))/W(2)
        H(3) = (BB(3)-DELTA(3)*H(2)-AP(3)*H(1))/W(3)
        DO 30 I=4,N
        H(I) = (BB(I)-DELTA(I)*H(I-1)-AP(I)*H(I-2)
     +               -Z(I)*H(I-3))/W(I)
   30   CONTINUE
C
C       Step Four
        XANS(N) = H(N)
        XANS(N-1) = H(N-1)-BETA(N-1)*XANS(N)
        XANS(N-2) = H(N-2)-BETA(N-2)*XANS(N-1)-GAM(N-2)*XANS(N)
        DO 40 I=N-3,1,-1
        XANS(I) = H(I)-BETA(I)*XANS(I+1)-GAM(I)*XANS(I+2)
     +                                  -PI(I)*XANS(I+3)
   40   CONTINUE
        RETURN
        END
c---------------------------------------------------------------
