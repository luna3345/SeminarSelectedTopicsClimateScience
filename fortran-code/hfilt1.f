      SUBROUTINE HFILT2 (NMAX,NSTEPS,DT,TAUC,IWINDOW,WPARAM,HH)
c
c     Routine to call HFILT1 and fill up the array HH.
c     HFILT1 gets the values HH(0) to HH(NSTEPS), and it
c     is completed under the assumption that it is symmetric.

      REAL HH(-NMAX:NMAX)
c
      CALL HFILT1 (NSTEPS,DT,TAUC,IWINDOW,WPARAM,HH(0))
c
      DO 100 N=1,NSTEPS
        HH(-N) = HH(N)
  100 CONTINUE
      return
      end
c
c--------------------------------------------------------
c
      SUBROUTINE HFILT1 (NSTEPS,DT,TAUC,IWINDOW,WPARAM,HH)
c    
c     Calculate filter weights with selected window.
c
c       Ref: see Hamming, R.W., 1989: Digital Filters,
c                Prentice-Hall International. 3rd Edition.
c
c       INPUT:      NSTEPS  -  Number of timesteps
c                              forward or backward.
c
c                   DT      -  Time step in seconds.
c
c                   TAUC    -  Cut-off period in hours.
c
c                   IWINDOW -  Indicator for selected window.
c
c                   WPARAM  -  Parameter (possibly) needed
c                              to specify window.
c
c
c       OUTPUT:     HH      -  Array(0:NSTEPS) with the
c                              required filter weights
c
c------------------------------------------------------------
c
        real HH(0:NSTEPS)
c
        data pi / 3.14159265358979 /

c------------------------------------------------------------
c
c       Windows are defined by a call and held in HH. 
        IF ( IWINDOW .eq. 0 ) CALL UNIFORM  (NSTEPS,dummy ,HH)
        IF ( IWINDOW .eq. 1 ) CALL LANCZOS  (NSTEPS,WPARAM,HH)
        IF ( IWINDOW .eq. 2 ) CALL HAMMING  (NSTEPS,WPARAM,HH)
        IF ( IWINDOW .eq. 3 ) CALL BLACKMAN (NSTEPS,dummy ,HH)
        IF ( IWINDOW .eq. 4 ) CALL KAISER   (NSTEPS,WPARAM,HH)
        IF ( IWINDOW .eq. 5 ) CALL POTTER2  (NSTEPS,WPARAM,HH)
C
        IF ( IWINDOW .eq. 50 ) CALL CRAZY1   (NSTEPS,WPARAM,HH)
        IF ( IWINDOW .eq. 60 ) CALL CRAZY2   (NSTEPS,WPARAM,HH)
        IF ( IWINDOW .eq. 70 ) CALL CRAZY3   (NSTEPS,WPARAM,HH)
        IF ( IWINDOW .eq. 80 ) CALL CRAZY4   (NSTEPS,WPARAM,HH)
c
c------------------------------------------------------------
c
c       calculate the cutoff frequency
        omegac = 2.*pi/(TAUC*3600.)
c
        NN = IABS(NSTEPS)
        DELTAT = ABS(DT)
c
        do 100 n=0,NN
          window = HH(n)
          if ( n .EQ. 0 ) then
            h = (omegac*DELTAT/pi)
          else
            h = sin(n*omegac*DELTAT)/(n*pi)
          endif
          HH(n) = h*window
  100   continue
c
c       normalize the sums to be unity
        call NORMLZ(HH,NN)
c
        return
        end
c
c===========================================================
c
        SUBROUTINE NORMLZ(HH,NMAX)
c
c       normalize the sum of HH to be unity
        REAL HH(0:NMAX)
c
        sumhh = HH(0) 
        do 100 n=1,NMAX
          sumhh = sumhh + 2*HH(n)
  100   continue
        do 200 n=0,NMAX
          HH(n)  = HH(n)/sumhh
  200   continue
c
        return
        end
c
c===========================================================
c
      SUBROUTINE UNIFORM(NSTEPS,dummy,WW)
c
c     define UNIFORM or RECTANGULAR window function.
c
      real WW(0:NSTEPS)
c
      DO 100 n=0,nsteps
        WW(n) = 1. 
  100 continue
      return
      end
c-----------------------------------------------------------
c
      SUBROUTINE LANCZOS(NSTEPS,power,WW)
c
c     define (Genaralised) LANCZOS window function.
c     (For the usual Lanczos window, power = 1 )
c
      real WW(0:NSTEPS)
c
        pi=4*atan(1.)
        do 100 n=0,nsteps
          if ( n .EQ. 0 ) then
            w = 1.0
          else
            w = sin(n*pi/(NSTEPS+1)) / ( n*pi/(NSTEPS+1))
          endif
          WW(n) = w**power
  100   continue
      return
      end
c-----------------------------------------------------------
c
      SUBROUTINE HAMMING(NSTEPS,alpha,WW)
c
c     define (Genaralised) HAMMING window function.
c     (For the usual Hamming window, alpha=0.54,
c          for the Hann window, alpha=0.50). 
c
      real WW(0:NSTEPS)
c
        pi=4*atan(1.)
        do 100 n=0,nsteps
          if ( n .EQ. 0 ) then
            w = 1.0
          else
            w = alpha + (1-alpha)*cos(n*pi/(NSTEPS))
          endif
          WW(n) = w
  100   continue
      return
      end
c-----------------------------------------------------------
c
      SUBROUTINE BLACKMAN(NSTEPS,dummy,WW)
c
c     define BLACKMAN window function.
c
      real WW(0:NSTEPS)
c
        pi=4*atan(1.)
        do 100 n=0,nsteps
          if ( n .EQ. 0 ) then
            w = 1.0
          else
            w = 0.42 + 0.50*cos(  n*pi/(NSTEPS))
     X                + 0.08*cos(2*n*pi/(NSTEPS))
          endif
          WW(n) = w
  100   continue
      return
      end
c-----------------------------------------------------------
c
      SUBROUTINE KAISER(NSTEPS,alpha,WW)
c
c     define KAISER window function.
c
      real WW(0:NSTEPS)
c
        XI0A =  bessi0(alpha)
        do 100 n=0,NSTEPS
          xn = n
          as = alpha*SQRT(1.-(xn/NSTEPS)**2)
          WW(n) = bessi0(as) / XI0A 
  100   continue
      return
      end
c
      FUNCTION BESSI0(X)
ccc   From NUMERICAL RECIPES (Press, et al.)
c
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,
     *    3.0899424D0,1.2067492D0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
c
c-----------------------------------------------------------
c
      SUBROUTINE POTTER2(NSTEPS,dummy,WW)
c
c     define POTTER window function.
c     Modified (by me) to fall off over twice the range.
c
      real WW(0:NSTEPS), D(0:3)
      DATA D / 0.35577019, 0.2436983, 0.07211497, 0.00630165 /
      DATA PI / 3.14159265 /
c
        CK = 1.0
        do 100 n=0,NSTEPS
          IF(n.EQ.NSTEPS) CK = 0.5
          arg = PI*FLOAT(N)/FLOAT(NSTEPS)
C              MODIFICATION IN NEXT STATEMENT 
               arg = arg/2.
C              END OF MODIFICATION
          sum = D(0)
          DO 50 IP=1,3
             sum = sum + 2.*D(IP)*COS(ARG*FLOAT(IP))
   50     CONTINUE
          WW(n) = CK*sum 
  100   continue
      return
      end
c
c-----------------------------------------------------------
c
      SUBROUTINE CRAZY1 (NSTEPS,alpha,WW)
c
c     define (Genaralised) HAMMING window function.
c     (For the usual Hamming window, alpha=0.54,
c          for the Hann window, alpha=0.50). 
c     CRAZY1 window is inverse of this.
c
      real WW(0:NSTEPS)
c
        pi=4*atan(1.)
        do 100 n=0,nsteps
          if ( n .EQ. 0 ) then
            w = 1.0
          else
            w = alpha + (1-alpha)*cos(n*pi/(NSTEPS))
          endif
          WW(n) = 1./w
  100   continue
      return
      end
c
      SUBROUTINE CRAZY2 (NSTEPS,gamma,WW)
c
c     define (Genaralised) HAMMING window function.
c     (For the usual Hamming window, alpha=0.54,
c          for the Hann window, alpha=0.50). 
c     CRAZY2 window is derived from this. 
c
      real WW(0:NSTEPS)
c
        pi=4*atan(1.)
        do 100 n=0,nsteps
          if ( n .EQ. 0 ) then
            w = 1.0
          else
            w =(1+gamma)-gamma*cos(2*n*pi/(NSTEPS))
          endif
          WW(n) = w
  100   continue
      return
      end
c
      SUBROUTINE CRAZY3 (NSTEPS,gamma,WW)
c
c     define (Genaralised) HAMMING window function.
c     (For the usual Hamming window, alpha=0.54,
c          for the Hann window, alpha=0.50). 
c     CRAZY3 window is derived from this. 
c
      real WW(0:NSTEPS)
c
        pi=4*atan(1.)
        do 100 n=0,nsteps
          if ( n .EQ. 0 ) then
            w = 1.0
          else
            w =(1+gamma)-gamma*cos(2*n*pi/(NSTEPS))
          endif
          WW(n) = 1./w
  100   continue
      return
      end
c
      SUBROUTINE CRAZY4 (NSTEPS,gamma,WW)
c
c     define (Genaralised) HAMMING window function.
c     (For the usual Hamming window, alpha=0.54,
c          for the Hann window, alpha=0.50). 
c     CRAZY4 window is derived from this. 
c
      real WW(0:NSTEPS)
c
        pi=4*atan(1.)
        do 100 n=0,nsteps
          if ( n .EQ. 0 ) then
            w = 1.0
          else
            w =(1+gamma)-gamma*cos(2*n*pi/(NSTEPS))
          endif
          WW(n) = w + 1./w
  100   continue
      return
      end
c-----------------------------------------------------------
