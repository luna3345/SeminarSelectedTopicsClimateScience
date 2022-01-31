        PROGRAM CALLDOLPH
C
C       Simple test program to call Subroutine DOLPH
C
        REAL HH (0:500 ), HH2(0:1000)
C
C-----------------------------------------------------
C       Set up the parameters for the CALL
C
C       Time Step
          DT = 300.
C       Specify Stop-Band Edge, let Ripple emerge
          NCHOICE = 1
          TAUS = 3.
C       Filter Length
          TSPAN = 3
          NSTEPS = TSPAN*3600./DT +.5
          NSHALF = NSTEPS/2
          IF ( 2*NSHALF.ne.NSTEPS ) STOP ' Use EVEN NSTEPS '
          N = NSTEPS + 1
          M = NSHALF
C
          TYPE *,' DOLPH. Input parameters '
          TYPE *, '        DT           Nchoice      Taus       TSPAN '  
          TYPE *,  DT, Nchoice, Taus, Tspan
          TYPE *, '           M           N '
          TYPE *,  M     , N       

C
C-----------------------------------------------------
C
          CALL  DOLPH(DT, NCHOICE, TauS, ripple, M, HH2)
C
C         DOLPH fills full h-array; we only need one side.
          DO NNN=0,M
            HH(NNN) = HH2(NNN+M)
          ENDDO
C
          TYPE *,' DOLPH. Outputs '
          TYPE *,' DOLPH. Ripple Ratio: ', ripple
          TYPE *,' DOLPH. Coefficients HH(n) : '
          DO NNN=0,M
            TYPE *, NNN, HH(NNN)
          ENDDO
C
C-----------------------------------------------------
        STOP
        END
C-----------------------------------------------------
C************************************************************
C------------------------------------------------------------
C
        SUBROUTINE  DOLPH(DeltaT,NCHOICE,TauS,r,M, WINDOW)
C
C       Calculation of Dolph-Chebyshev window or, for short,
C       Dolph Window, using the expression in the reference:
C
C       Antoniou, Andreas, 1993: Digital Filters: Analysis,
C       Design and Applications. McGraw-Hill, Inc., 689pp.
C
C       The Dolph window is optimal in the following sense:
C       For a given main-lobe width, the stop-band attenuation
C       is minimal; for a given stop-band level, the main-lobe
C       width is minimal.
C
C       It is possible to specify either the ripple-ratio r
C       or the stop-band edge THETA1.
C
C------------------------------------------------------------
C
        REAL WINDOW(0:2*M)
C
        PARAMETER (NMAX = 1000)
        REAL T(0:NMAX)
        REAL w(0:NMAX),TIME(0:NMAX)
        REAL w2(0:2*NMAX),TIME2(0:2*NMAX)
C
        COMPLEX HH(0:1000)
        CHARACTER*10 STRING
        integer system
C
C------------------------------------------------------------
C
C       Set up the parameters to specify the window.
C          There are two options: 
C            (1) specify the cutoff period, TauS
C            (2) specify the Ripple ratio, r.
C
        PI = 4*ATAN(1.D0)
C
        N = 2*M+1
        NM1 = N-1 
C
        IF      ( NCHOICE.eq.1 ) THEN
          TauS = TauS * 3600
          ThetaS = 2*PI*DeltaT/TauS
          X0 = 1/COS(ThetaS/2)  
          TERM1 = (X0 + SQRT(X0**2-1))**(FLOAT(N-1))  
          TERM2 = (X0 - SQRT(X0**2-1))**(FLOAT(N-1))  
          RR = 0.5*(TERM1+TERM2)
          r = 1/RR
          dB = 20*LOG10(r)
C
        ELSE IF ( NCHOICE.eq.2 ) THEN
          R = 0.1
            dB = 20*LOG10(R)
          RR = 1/R
          COSHINV = LOG(RR+SQRT(RR**2-1))
          X0 = COSH(COSHINV/(N-1))
          THETA1 = (FLOAT(N)/FLOAT(M))*ACOS(1/X0)
          THETA1 =                   2*ACOS(1/X0)
C
        ELSE 
          TYPE *, ' Invalid value of NCHOICE ', NCHOICE
          STOP
        ENDIF
C
C------------------------------------------------------------
C
        DO nT=0,M
          SUM = RR
          DO i=1,M
            arg = X0*cos(i*PI/N)
            CALL CHEBY(T,NM1,ARG)
            TERM1 = T(NM1)
            TERM2 = cos(2*nT*PI*i/N)
            SUM = SUM + 2*TERM1*TERM2
          ENDDO
          w(nT) = SUM/N
          TIME(nT) = nT 
C           TYPE *, ' TIME, w ',  TIME(nT), w(nT)
        ENDDO
C
C       Fill in the negative-time values by symmetry.
        DO nT=0,M
          w2(M+nT) = w(nT)
          w2(M-nT) = w(nT)
          TIME2(M+nT) = TIME(nT)
          TIME2(M-nT) = -TIME(nT)
        ENDDO
C
C----   Fill up the array for return. NORMILIZE WINDOW.
        WINSOM = 0
        DO nT=0,2*M
          WINDOW(nT) = w2(nT)
          WINSOM = WINSOM + WINDOW(nT)
        ENDDO
        DO nT=0,2*M
          WINDOW(nT) = WINDOW(nT) / WINSOM
        ENDDO
C
C----------------------------------------------------------
C
        RETURN
        END
C***********************************************************
        SUBROUTINE CHEBY(T,N,X)
C
C       Calculate all Chebyshev polynomials up to order N
C       for the argument value X.
C
C       Reference: Numerical Recipes, Page 184, recurrence
C           T_n(x) = 2xT_{n-1}(x) - T_{n-2}(x) ,  n>=2.
C
        REAL T(0:N)
C
        T(0) = 1
        T(1) = X
        IF(N.lt.2) RETURN
        DO NN=2,N
          T(NN) = 2*X*T(NN-1) - T(NN-2)
        ENDDO
C
        RETURN
        END
C***********************************************************
