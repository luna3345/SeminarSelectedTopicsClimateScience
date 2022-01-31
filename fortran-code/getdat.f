        SUBROUTINE GETDAT
C
c       Prepare Initial Data for RICHIE.
C 
C   (1) Read the basic data files:
c       PS.DAT --- sea level pressure on the E-grid
c       Z.DAT  --- heights of the pressure surfaces
c                  200, 400, 600 and 800 mb,
c                  on a coarse A_Grid
c       UV.DAT --- Wind Direction/Speed at levels
c                  100, 300, 500, 700 and 900mb,
c                  on a coarse A_Grid
c       T.DAT  --- 200-100mb thickness for calculation of
c                  temperature in the stratosphere
c       OROG.DAT - height of orography on fine A-grid.
c
c   (2) Fill up the E-grid from the coarse A-grid
c       where necessary (for Z, u, v and TTOP ),
c       in S/Rs ATOEZ and ATOEV.
c
c   (3) If required, fill up a fime A-grid to make
c       plotting easier, with the same grid for all.
c
c   (4) Convert sea-level pressure and heights to pressures
c       at ground level and conventional interface levels.
c
c   (5) Convert wind data to momenta in layers.
c
c=================================================================
c
c       size of grid
        PARAMETER (IM=9, JM=5, LM=5)
        PARAMETER (LMM1=LM-1)
c
c=================================================================
c
c      2-D Work fields common to all modules.
       COMMON /WORK/ WORK1(IM,JM), WORK2(IM,JM), WORK3(IM,JM)
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
c       size of fine A-grid
        PARAMETER (JM2=10,JM2M1=JM2-1,JM2M2=JM2-2,JM2M3=JM2-3)
C
c       space for various DIAGNOSTIC variables on FINE grid
        real P2 (IM,JM2,LM),U2(IM,JM2,LM),V2(IM,JM2,LM)
        COMMON / WDIAG2 / P2, U2, V2
c
c       space for the basic PROGNOSTIC variables on FINE grid
        real PP2(IM,JM2,LM),MU2(IM,JM2,LM),
     +       MV2(IM,JM2,LM),TTOP2(IM,JM2)
        COMMON / WPROG2 / PP2, MU2, MV2, TTOP2
C
c=================================================================
c
        LUNPS= 21
        LUNZ = 22
        LUNV = 23
        LUNT = 24
        LUNO = 25
        open(unit=LUNPS,file='PS.DAT')
        open(unit=LUNZ ,file='Z.DAT')
        open(unit=LUNV ,file='UV.DAT')
        open(unit=LUNT ,file='T.DAT')
        open(unit=LUNO ,file='OROG.DAT')
c
c-----------------------------------------------------------------
c
c       read the sea-level pressure.
        CALL READPS(LUNPS,PS,IM,JM)
c
C       read the geopotential heights.
        CALL READZ(LUNZ,Z,IM,JM,LMM1)
c
c       read the wind data.
        CALL READV(LUNV,U,V,DDDFF,IM,JM,LM)
c
c       read the stratospheric temperature.
        CALL READT(LUNT,TTOP,IM,JM)
c
c       read the height of orography
        CALL READO(LUNO,OROG,IM,JM)
c
c       Get the initial pressures.
        CALL ZTOP(PS,Z,TTOP,OROG,P,IM,JM,LM)
c
c       Get the initial momenta.
        CALL UVTOMM(U,V,P,MU,MV,IM,JM,LM)
c
C + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
c
C       WRITE OUT THE INITIAL FIELDS.
c       NDISC=88
c       OPEN(UNIT=NDISC,FILE='RICHIE.DAT',FORM='FORMATTED')
c       WRITE(NDISC,777) P,MU,MV,TTOP,OROG
c777    format( 5E20.8 )
C
C + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
C
c       fill the fine A-grid for plotting, and write out.
c       CALL FILLZ(PS,ZPLOT,IM,JM)
c                LUNPLT = 30
c                open(unit=LUNPLT,file='PSPLOT.OUT')
c                WRITE(LUNPLT,9) ((ZPLOT(I,J),I=1,IM),J=2*JM,1,-1)
c   9            FORMAT(13f6.0)
c                close(unit=LUNPLT)
c
        RETURN
        end
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE READPS(LUN,ARRAY,IM,JM)
c
c       read a scalar field on a staggered E-Grid, 
c       as it appears on the printed page.
C       THIS VERSION CONVERTS SURFACE PRESSURE FROM 
C                10 X mm HG to Pascals.
c
        REAL ARRAY(IM,JM)
        CHARACTER*80 HEADER
c
c       FLUSH OUT the HEADER RECORD
        READ(LUN,9901) HEADER 
        print9901, HEADER
 9901   FORMAT(A80)
c
        IMM1 = IM-1
c
        DO J=JM,1,-1
c 		*there is a sneaky "implied DO" here, which loops over I
          READ(LUN,*) ( ARRAY(I,J) , I=2,IMM1,2 )
          READ(LUN,*) ( ARRAY(I,J) , I=1,IM  ,2 )
              print98, ( ARRAY(I,J) , I=2,IMM1,2 )
              print99, ( ARRAY(I,J) , I=1,IM  ,2 )
   98            FORMAT(13(6x,f6.0))
   99            FORMAT(13(f6.0,6x))
        ENDDO
c
        DO J=1,JM
        DO I=1,IM
           ARRAY(I,J) = (1.33322 * ARRAY(I,J)/10.) * 100.
        ENDDO
        ENDDO
c
        DO J=JM,1,-1
          print998, ( ARRAY(I,J) , I=2,IMM1,2 )
          print999, ( ARRAY(I,J) , I=1,IM  ,2 )
  998            FORMAT(13(6x,-2pf6.0))
  999            FORMAT(13(-2pf6.0,6x))
        ENDDO
C
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE READO(LUN,OROG,IM,JM)
c
c       read a scalar field on a fine A-Grid, 
c
        REAL OROG(IM,2*JM)
        CHARACTER*80 HEADER
c
        IMM1 = IM-1
        JM2 = 2*JM
c
c       FLUSH OUT the HEADER RECORD
        READ(LUN,9901) HEADER 
        print9901, HEADER
 9901   FORMAT(A80)
        DO J=JM2,1,-1
          READ(LUN,*) ( OROG(I,J) , I=1,IM )
                 print99, ( OROG(I,J) , I=1,IM )
   99            FORMAT(13(f6.0))
        ENDDO
c
C       Smooth the Orography from the Fine A-Grid to the E-Grid.
        print*,' OROG(7,6), OROG(7,7) Before Smoothing ', 
     X            OROG(7,6),OROG(7,7)
        call SMEARZ(OROG,OROG,IM,JM)
        print*,' OROG(7,6), OROG(7,7) After  Smoothing ', 
     X            OROG(7,6),OROG(7,7)
C
C       printout the Smoothed Orography.
        print9902
 9902   FORMAT(' Orography After Smoothing')
        DO J=JM2,1,-1
           print99, ( OROG(I,J) , I=1,IM )
        ENDDO
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE READZ(LUN,ARRAY,IM,JM,LM)
c
c       read a scalar field on a Sparse A-Grid,
c       as it appears on the printed page,
c       and fill it out to an E-grid.
c
        REAL ARRAY(IM,JM,LM)
        CHARACTER*80 HEADER
c
        IMM1 = IM-1
c
c       Loop for each level
        DO 1000 L=1,LM
c
c       FLUSH OUT the HEADER RECORD
        READ(LUN,9901) HEADER 
        print9901, HEADER
 9901   FORMAT(A80)
        DO J=JM,1,-1
          READ(LUN,*) ( ARRAY(I,J,L) , I=1,IM  ,2 )
                 print99, ( ARRAY(I,J,L) , I=1,IM  ,2 )
   99            FORMAT(13(f6.0,6x))
        ENDDO
c
c       Fill the E-grid by interpolation.
c       print9991
c       DO J=JM,1,-1
c       print9993, (array(I,J,L),I=2,IMM1,2) 
c       print9994, (array(I,J,L),I=1,IM,2)
c       ENDDO
c
        CALL ATOEZ(ARRAY(1,1,L),IM,JM)
c
c       print9992
c       DO J=JM,1,-1
c       print9993, (array(I,J,L),I=2,IMM1,2)
c       print9994, (array(I,J,L),I=1,IM,2)
c       ENDDO
c
c9991   FORMAT(' Z BEFOR ATOEZ ')
c9992   FORMAT(' Z AFTER ATOEZ ')
c9993   FORMAT(3X,6(1X,F8.0))
c9994   FORMAT(7(F8.0,1X))
c
 1000   CONTINUE
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE READT(LUN,TTOP,IM,JM)
c
c       read the thickness 200-100mb height field on a 
c       Sparse A-Grid,and convert thickness to 
c       stratospheric temperature.
c
        REAL TTOP(IM,JM)
        CHARACTER*80 HEADER
c
        GRAV = 9.80
        RGAS = 287.
        DELTAP = (200.-100.) * 100.
        PBAR = 150. * 100.
        FACTOR = ( GRAV*PBAR ) / ( RGAS*DELTAP )
C
        IMM1 = IM-1
c
c       FLUSH OUT the HEADER RECORD
        READ(LUN,9901) HEADER 
        print9901, HEADER
 9901   FORMAT(A80)
        DO J=JM,1,-1
          READ(LUN,*) ( TTOP(I,J) , I=1,IM  ,2 )
                 print99, ( TTOP(I,J) , I=1,IM  ,2 )
   99            FORMAT(13(f6.0,6x))
        ENDDO
c
c       Fill the E-grid by interpolation.
        CALL ATOEZ(TTOP,IM,JM)
c
C       Convert thickness to temperature
        DO J=1,JM
        DO I=1,IM
          TTOP(I,J) = FACTOR*TTOP(I,J)
        ENDDO
        ENDDO
c 
               DO J=JM,1,-1
                 print9993, (TTOP(I,J),I=2,IMM1,2)
                 print9994, (TTOP(I,J),I=1,IM,2)
 9993            FORMAT(3X,6(1X,F8.0))
 9994            FORMAT(7(F8.0,1X))
               ENDDO
c
        RETURN
        END
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE READV(LUN,U,V,DDDFF,IM,JM,LM)
c
c       Read a vector field on a Sparse A-Grid,
c       as it appears on the printed page,
c       and fill it out to an E-grid.
C       Values input as ddd.ff and split into u,v.
C 		*ddd .. directions in degrees (0=N, 90=E, 180=S, 270=W)
C 		*ff .. wind speed in knots (1kn = 0.514m/s)
c
        REAL U(IM,JM,LM),V(IM,JM,LM), DDDFF(IM,JM)
        CHARACTER*80 HEADER
c
        IMM1 = IM-1
c
c       Loop over each level
        DO 1000 L=1,LM
c
c       FLUSH OUT the HEADER RECORD
        READ(LUN,9901) HEADER 
        print9901, HEADER
 9901   FORMAT(A80)
        DO J=JM,1,-1
          READ(LUN,*) ( DDDFF(I,J) , I=1,IM,2 )
                 print98, ( DDDFF(I,J), I=1,IM,2 )
   98            FORMAT(7f7.2)
        ENDDO
c
        drad = 3.14159265/180.
        DO J=JM,1,-1
        DO I=1,IM,2
          df = DDDFF(I,J)
C		  *INT(df) kappt die Kommastellen .ff weg und lässt nur die direction übrig
          dir= INT(df)
C		  * 100*(ddd.ff-ddd)=100*0.ff = ff lässt nur speed übrig (in knoten???)
          speed = 100*(df-dir)
          angle = 270 - dir
          if (angle.lt.0) angle = angle+360
          rangle = angle*drad
          uu = speed*cos(rangle)
          vv = speed*sin(rangle)
          U(I,J,L) = uu
          V(I,J,L) = vv
c          print9, i,j, dir,angle,speed,uu,vv
    9     format('ij dir angle speed uv',2I3,3F6.0,4X,2F4.0)
        ENDDO
        ENDDO
c
c       Fill the E-grid by interpolation.
c       print9991
c       DO J=JM,1,-1
c       print9995, (U(I,J,L),I=1,IM,2) , (U(I,J,L),I=2,IMM1,2)
c       ENDDO
c       print9992
c       DO J=JM,1,-1
c       print9995, (V(I,J,L),I=1,IM,2) , (V(I,J,L),I=2,IMM1,2)
c       ENDDO
c
        CALL ATOEV(U(1,1,L),V(1,1,L),IM,JM)
c
c       print9993
c       DO J=JM,1,-1
c       print9995, (U(I,J,L),I=1,IM,2) , (U(I,J,L),I=2,IMM1,2)
c       ENDDO
c       print9994
c       DO J=JM,1,-1
c       print9995, (V(I,J,L),I=1,IM,2) , (V(I,J,L),I=2,IMM1,2)
c       ENDDO
c
c9991   FORMAT(' U BEFOR ATOEV ')
c9992   FORMAT(' V BEFOR ATOEV ')
c9993   FORMAT(' U AFTER ATOEV ')
c9994   FORMAT(' V AFTER ATOEV ')
c9995   FORMAT( (6(F5.1,5X),F5.1/6(5X,F5.1)) )
c
 1000   CONTINUE
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE ATOEZ(Z,IM,JM)
c
c       From a scalar field on a sparse A-Grid, 
c       construct it on a denser A-grid,
c       using bilinear interpolation.
c
        REAL Z(IM,JM)
c
        IMM1 = IM-1
        JMM1 = JM-1
c
c       Internal Points.
        DO J=1,JMM1
        DO I=2,IMM1,2
          Z(I,J) = ( Z(I-1,J  )+Z(I+1,J  ) 
     X             + Z(I-1,J+1)+Z(I+1,J+1) ) / 4.
        ENDDO
        ENDDO
c
c       Top Row
        DO I=2,IMM1,2
          Z(I,JM) = ( Z(I-1,JM)+Z(I+1,JM)-Z(I,JMM1) )
        ENDDO
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE ATOEV(U,V,IM,JM)
c
c       From a vector field on a sparse A-Grid, 
c       construct it on a denser A-grid,
c       using bilinear interpolation.
c
        REAL U(IM,JM), V(IM,JM)
c
        IMM1 = IM-1
        JMM1 = JM-1
c
c       Internal Points.
        DO J=2,JM
        DO I=2,IMM1,2
          U(I,J) = ( U(I-1,J  )+U(I+1,J  ) 
     X             + U(I-1,J-1)+U(I+1,J-1) ) / 4.
          V(I,J) = ( V(I-1,J  )+V(I+1,J  ) 
     X             + V(I-1,J-1)+V(I+1,J-1) ) / 4.
        ENDDO
        ENDDO
c
c       Bottom Row
        DO I=2,IMM1,2
          U(I,1) = ( U(I-1,1)+U(I+1,1)-U(I,2) )
          V(I,1) = ( V(I-1,1)+V(I+1,1)-V(I,2) )
        ENDDO
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE FILLZ(Z,Z2,IM,JM)
c
c       From a scalar field on a staggered E-Grid, 
c       construct another field on a denser A-grid.
c       Interpolation is linear at edges
c       and bilinear for internal points.
c
        REAL Z(IM,JM), Z2(IM,2*JM)
c
        IMM1 = IM-1
        IMM2 = IM-2
        JMM1 = JM-1
        JMX2 = 2*JM
c
c       First, fill in the known points.
        DO J=1,JM
          DO I=1,IM,2
            Z2(I,J*2-1) =  Z(I  ,J)
          ENDDO
          DO I=2,IMM1,2
            Z2(I,J*2  ) =  Z(I  ,J)
          ENDDO
        ENDDO
c
C       Next, interpolate the unknown points.
c
c       Bottom Row.
        DO I=2,IMM1,2
          Z2(I,1) = ( Z(I-1,1)+Z(I+1,1) ) /2.
        ENDDO
c
c       Top Row.
        Z2(1 ,JMX2) = ( Z(1 ,JM)+Z(2   ,JM) ) /2.
        DO I=3,IMM2,2
          Z2(I,JMX2) = ( Z(I-1,JM)+Z(I+1,JM) ) /2.
        ENDDO
        Z2(IM,JMX2) = ( Z(IMM1,JM)+Z(IM,JM) ) /2.
c
c       Left and Right Sides
        DO J=1,JMM1
          Z2(1 ,J*2) = ( Z(1 ,J)+Z(1 ,J+1) ) /2.
          Z2(IM,J*2) = ( Z(IM,J)+Z(IM,J+1) ) /2.
        ENDDO
c
c       Internal Points.
        DO J=2,JM
        DO I=2,IMM1,2
          Z2(I,J*2-1) = ( Z(I  ,J)+Z(I,J-1) 
     X                     + Z(I+1,J)+Z(I-1,J) ) / 4.
        ENDDO
        ENDDO
        DO J=1,JMM1
        DO I=3,IMM2,2
          Z2(I,J*2  ) = ( Z(I  ,J)+Z(I,J+1) 
     X                     + Z(I+1,J)+Z(I-1,J) ) / 4.
        ENDDO
        ENDDO
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE FILLV(V,V2,IM,JM)
c
c       From a vector field on a staggered E-Grid, 
c       construct another field on a denser A-grid.
c       Interpolation is linear at edges
c       and bilinear for internal points.
c
        REAL V(IM,JM), V2(IM,2*JM)
c
        IMM1 = IM-1
        IMM2 = IM-2
        JMM1 = JM-1
        JMX2 = 2*JM
c
c       First, fill in the known points.
        DO J=1,JM
          DO I=1,IM,2
            V2(I,J*2  ) =  V(I  ,J)
          ENDDO
          DO I=2,IMM1,2
            V2(I,J*2-1) =  V(I  ,J)
          ENDDO
        ENDDO
c
C       Next, interpolate the unknown points.
c
c       Top Row.
        DO I=2,IMM1,2
          V2(I,JMX2) = ( V(I-1,JM)+V(I+1,JM) ) /2.
        ENDDO
c
c       Bottom Row.
        V2(1 ,1) = ( V(1 ,1)+V(2   ,1) ) /2.
        DO I=3,IMM2,2
         V2(I,1) = ( V(I-1,1)+V(I+1,1) ) /2.
        ENDDO
        V2(IM,1) = ( V(IMM1,1)+V(IM,1) ) /2.
c
c       Left and Right Sides
        DO J=2,JM
          V2(1 ,J*2-1) = ( V(1 ,J)+V(1 ,J-1) ) /2.
          V2(IM,J*2-1) = ( V(IM,J)+V(IM,J-1) ) /2.
        ENDDO
c
c       Internal Points.
        DO J=1,JMM1
        DO I=2,IMM1,2
          V2(I,J*2  ) = ( V(I  ,J)+V(I,J+1) 
     X                  + V(I+1,J)+V(I-1,J) ) / 4.
        ENDDO
        ENDDO
        DO J=2,JM
        DO I=3,IMM2,2
          V2(I,J*2-1) = ( V(I  ,J)+V(I,J-1) 
     X                  + V(I+1,J)+V(I-1,J) ) / 4.
        ENDDO
        ENDDO
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE ZTOP(PS,Z,TTOP,OROG,P,IM,JM,LM)
C
c       Get the initial pressures. 
C
        REAL PS(IM,JM),TTOP(IM,JM)
        REAL OROG(IM,2*JM)
        REAL Z(IM,JM,LM-1)
        REAL P(IM,JM,LM)
c
        integer VMETHOD
        DATA VMETHOD /3/
C
c       Conventional interface levels (km).
        real PBAR(5),TBAR(5), ZBAR(5), ZDYNAM(5)
        DATA PBAR /   200.,   400.,   600.,   800.,  1013./
        DATA TBAR /   -55.,   -32.,   -12.,    +2.,   +15./
        DATA ZBAR /   11.8,    7.2,   4.2 ,   2.0 ,   0.0 /
c       change to 'dynamic kilometers" (LFR, p 181).
        DATA ZDYNAM / 11.543,  7.048,  4.113,  1.959,   0.0 /
c
        IMM1 = IM-1
        LMM1 = LM-1
c
        GRAV=9.80
        RGAS=287.
c
c       convert upper heights to pressures.
        DO 100 L=1,LMM1
        ZL = ZDYNAM(L)*1000.
        TL = TBAR(L)+273.
        PL = PBAR(L)*100.
        HL = RGAS*TL/GRAV         
          print*,' ZDYNAM, PBAR, TBAR : ', ZL,PL,TL
        DO 100 J=1,JM
        DO 100 I=1,IM 
          P(I,J,L) = PL * ( 1. + (Z(I,J,L)-ZL)/HL )
  100   CONTINUE
c
c       interpolate the surface pressure.
           ZL = ZDYNAM(LM)*1000.
           TL = TBAR(LM)+273.
           PL = PBAR(LM)*100.
           print*,' ZDYNAM, PBAR, TBAR : ', ZL,PL,TL
      IF ( VMETHOD.EQ.1 ) THEN
           DO 150 J=1,JM
           DO 150 I=1,IM 
             H9 = Z(I,J,LMM1) / LOG(PS(I,J)/P(I,J,LMM1))
             if(I/2*2.EQ.I) JH=2*J
             if(I/2*2.NE.I) JH=2*J-1
             ZG = OROG(I,JH)
             P(I,J,LM) = PS(I,J) * EXP(-ZG/H9) 
  150      CONTINUE
      ELSE IF ( VMETHOD.EQ.2 ) THEN
           DO 151 J=1,JM
           DO 151 I=1,IM 
             PMID = (PS(I,J)+P(I,J,LMM1))/2.
             PDIF = (PS(I,J)-P(I,J,LMM1))
             ZDIF = (Z(I,J,LMM1)-0.)
             H9 = (PMID/PDIF)*ZDIF
             if(I/2*2.EQ.I) JH=2*J
             if(I/2*2.NE.I) JH=2*J-1
             ZG = OROG(I,JH)
             DELZ = Z(I,J,LMM1)-ZG
             QQ=DELZ/(2*H9)
             P(I,J,LM) = P(I,J,LMM1) * (1+QQ)/(1-QQ)
  151      CONTINUE
      ELSE IF ( VMETHOD.EQ.3 ) THEN
           TZERO = 273. + 15.
           PZERO = 1013.25 * 100.
           GG = 9.80665
           RR = 287.
           GAM = 0.0065
           PWR = GG/(GAM*RR)
           DO 1511 J=1,JM
           DO 1511 I=1,IM 
             if(I/2*2.EQ.I) JH=2*J
             if(I/2*2.NE.I) JH=2*J-1
             ZG = OROG(I,JH)
             DELZ = ZG
             P(I,J,LM) = PS(I,J) * ( 1-(GAM/TZERO)*DELZ )**(+PWR)
CCC          print*,'P CHECK ', P(I,J,LM) , PS(I,J) 
 1511      CONTINUE
      ENDIF
c
        DO 200 L=1,LM
        print97,L
   97   FORMAT(' INTERFACE PRESSURES, Level: ',I4)
        DO J=JM,1,-1
           print98, ( P(I,J,L) , I=2,IMM1,2 )
           print99, ( P(I,J,L) , I=1,IM  ,2 )
   98      FORMAT(13(6x,-2pf6.0))
   99      FORMAT(13(-2pf6.0,6x))
        ENDDO
  200   CONTINUE
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE PTOZ(P2,TTOP2,OROG, Z2,PSEA2,IM,JM2,LM)
C
c       Get the sea-level pressure and upper heights. 
C
        REAL PSEA2(IM,JM2),TTOP2(IM,JM2)
        REAL OROG(IM,JM2)
        REAL Z2(IM,JM2,LM-1)
        REAL P2(IM,JM2,LM)
c
        integer VMETHOD
        DATA VMETHOD /3/
C
c       Conventional interface levels (km).
        real PBAR(5),TBAR(5), ZBAR(5), ZDYNAM(5)
        DATA PBAR /   200.,   400.,   600.,   800.,  1013./
        DATA TBAR /   -55.,   -32.,   -12.,    +2.,   +15./
        DATA ZBAR /   11.8,    7.2,   4.2 ,   2.0 ,   0.0 /
c       change to 'dynamic kilometers" (LFR, p 181).
        DATA ZDYNAM / 11.543,  7.048,  4.113,  1.959,   0.0 /
c
        IMM1 = IM-1
        LMM1 = LM-1
        JM = JM2/2
c
        GRAV=9.80
        RGAS=287.
c
c       convert upper pressures to heights 
        DO 100 L=1,LMM1
        ZL = ZDYNAM(L)*1000.
        TL = TBAR(L)+273.
        PL = PBAR(L)*100.
        HL = RGAS*TL/GRAV         
CCC       print*,' PTOZ:: ZDYNAM, PBAR, TBAR : ', ZL,PL,TL
        DO 100 J=1,JM2
        DO 100 I=1,IM 
          Z2(I,J,L) = ZL + (P2(I,J,L)/PL - 1.)*HL
  100   CONTINUE
c
c       interpolate the surface pressure.
        ZL = ZDYNAM(LM)*1000.
        TL = TBAR(LM)+273.
        PL = PBAR(LM)*100.
CCC     print*,' PTOZ:: ZDYNAM, PBAR, TBAR : ', ZL,PL,TL
      IF ( VMETHOD.EQ.1 ) THEN
          DO 150 J=1,JM2
          DO 150 I=1,IM 
            ZG = OROG(I,J)
            H9 = (Z2(I,J,LMM1)-ZG) / LOG(P2(I,J,LM)/P2(I,J,LMM1))
            PSEA2(I,J) = P2(I,J,LM) * EXP(+ZG/H9) 
  150     CONTINUE
      ELSE IF ( VMETHOD.EQ.2 ) THEN
          DO 151 J=1,JM2
          DO 151 I=1,IM 
             PMID = (P2(I,J,LM)+P2(I,J,LMM1))/2.
             PDIF = (P2(I,J,LM)-P2(I,J,LMM1))
             ZG = OROG(I,J)
             ZDIF = (Z2(I,J,LMM1)-ZG)
             H9 = (PMID/PDIF)*ZDIF
             DELZ = Z2(I,J,LMM1)-0.
             QQ=DELZ/(2*H9)
             PSEA2(I,J) = P2(I,J,LMM1) * (1+QQ)/(1-QQ)
  151     CONTINUE
      ELSE IF ( VMETHOD.EQ.3 ) THEN
          TZERO = 273. + 15.
          PZERO = 1013.25 * 100.
          GG = 9.80665
          RR = 287.
          GAM = 0.0065
          PWR = GG/(GAM*RR)
          DO 1511 J=1,JM2
          DO 1511 I=1,IM 
             ZG = OROG(I,J)
             DELZ = ZG
             PSEA2(I,J) = P2(I,J,LM) * (1-(GAM/TZERO)*DELZ)**(-PWR)
CCC          print*,'P CHECK ', P2(I,J,LM) , PSEA2(I,J) 
 1511     CONTINUE
      ENDIF
c
C       INTERPOLATE FROM E- TO A-GRID.
        call SMEARZ(PSEA2,PSEA2,IM,JM)
C
c       Calculate the 1000 mb geopotential. 
        PL = 1000. * 100.
        TL = TBAR(LM)+273.
        HS = RGAS*TL/GRAV
        DZ = HS/PL 
CCC     print*,' Z1000:: PL,TL,HS,DZ:', PL,TL,HS,DZ
          DO 152 J=1,JM2
          DO 152 I=1,IM 
            Z2(I,J,LM) = ( PSEA2(I,J)-PL ) * DZ 
  152     CONTINUE
c
        print98
   98   FORMAT(' SEA LEVEL PRESSURE ')
        DO J=JM2,1,-1
           print91, ( PSEA2(I,J) , I=1,IM)
   91      FORMAT(13(-2pf6.0))
        ENDDO
c
        RETURN
c
        DO 200 L=1,LM
        print97,L
   97   FORMAT(' INTERFACE HEIGHTS, Level: ',I4)
        DO J=JM2,1,-1
           print99, ( Z2(I,J,L) , I=1,IM)
   99      FORMAT(13(f6.0))
        ENDDO
  200   CONTINUE
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE UVTOMM(U,V,P,MU,MV,IM,JM,LM)
C
c       Get the initial momenta. 
c
c=================================================================
c
c      2-D Work fields common to all modules.
       PARAMETER (IX=13, JX=7, LX=5)
       COMMON /WORK/ WORK1(IX,JX), WORK2(IX,JX), WORK3(IX,JX)
c
c=================================================================
C
        REAL U(IM,JM,LM), V(IM,JM,LM)
        REAL P(IM,JM,LM), MU(IM,JM,LM), MV(IM,JM,LM)
c
        GRAV=9.80
C HACK of initalization of indices
        IMM1 = IM-1
c       convert velocities to momenta.
c
        DO 300 L=1,LM
c
c       Get dp/g and interpolate to wind points.
        DO 100 J=1,JM
        DO 100 I=1,IM 
          IF ( L.EQ.1 ) THEN
             RR = P(I,J,1)/GRAV
          ELSE
             RR = ( P(I,J,L)-P(I,J,L-1) )/GRAV
          ENDIF
          WORK1(I,J) = RR
  100   CONTINUE
C
        CALL ZTOV(WORK1,WORK2,IM,JM)
c
        DO 200 J=1,JM
        DO 200 I=1,IM 
          RR = WORK2(I,J)
          MU(I,J,L) = RR * U(I,J,L) 
          MV(I,J,L) = RR * V(I,J,L) 
  200   CONTINUE
c
  300   CONTINUE
c
        DO 400 L=1,LM
        print96,L
 96     FORMAT(' LAYER MOMENTA (MU) , Level: ',I4)
        DO J=JM,1,-1
           print99, (MU(I,J,L) , I=1,IM  ,2 )
           print98, (MU(I,J,L) , I=2,IMM1,2 )
   99      FORMAT(13(-2pf6.0,6x))
   98      FORMAT(13(6x,-2pf6.0))
CCCCCCCC   these scaling factors make printed values
CCCCCCCC   "look the same size" as in LFR's table.
        ENDDO
        print97,L
   97   FORMAT(' LAYER MOMENTA (MV) , Level: ',I4)
        DO J=JM,1,-1
           print99, (MV(I,J,L) , I=1,IM  ,2 )
           print98, (MV(I,J,L) , I=2,IMM1,2 )
        ENDDO
  400   CONTINUE
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE MMTOUV(P2,MU2,MV2, U2,V2,IM,JM2,LM)
C
c       Get the velocities from momenta. 
c
c=================================================================
c
c      2-D Work fields common to all modules.
       PARAMETER (IX=13, JX=7, LX=5)
       COMMON /WORK/ WORK1(IX,JX), WORK2(IX,JX), WORK3(IX,JX)
c
c=================================================================
C
        REAL P2(IM,JM2,LM), MU2(IM,JM2,LM), MV2(IM,JM2,LM)
        REAL U2(IM,JM2,LM), V2(IM,JM2,LM)
c
        GRAV=9.80
c
c       convert momenta to velocities.
c
        DO 300 L=1,LM
c
        DO 100 J=1,JM2
        DO 100 I=1,IM 
          IF ( L.EQ.1 ) THEN
             RR = P2(I,J,1)/GRAV
          ELSE
             RR = ( P2(I,J,L)-P2(I,J,L-1) )/GRAV
          ENDIF
          U2(I,J,L) = MU2(I,J,L) / RR
          V2(I,J,L) = MV2(I,J,L) / RR
  100   CONTINUE
c
  300   CONTINUE
c
        RETURN
c
        DO 400 L=1,LM
        print96,L
   96   FORMAT(' LAYER VELOCITY (U) , Level: ',I4)
        DO J=JM2,1,-1
           print99, (U2(I,J,L) , I=1,IM)
   99      FORMAT(13(f6.0))
        ENDDO
        print97,L
   97   FORMAT(' LAYER VELOCITY (V) , Level: ',I4)
        DO J=JM2,1,-1
           print99, (V2(I,J,L) , I=1,IM)
        ENDDO
  400   CONTINUE
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE ZTOV(Z,ZV,IM,JM)
c
c       From a scalar field on a staggered E-Grid, 
c       construct values at the vector points.
c       Interpolation is linear at edges
c       and bilinear for internal points.
c
        REAL Z(IM,JM), ZV(IM,JM)
c
        IMM1 = IM-1
        IMM2 = IM-2
        JMM1 = JM-1
c
c       Bottom Row.
        DO I=2,IMM1,2
          ZV(I,1) = ( Z(I-1,1)+Z(I+1,1) ) /2.
        ENDDO
c
c       Top Row.
        ZV(1 ,JM) = ( Z(1 ,JM)+Z(2   ,JM) ) /2.
        DO I=3,IMM2,2
          ZV(I,JM) = ( Z(I-1,JM)+Z(I+1,JM) ) /2.
        ENDDO
        ZV(IM,JM) = ( Z(IMM1,JM)+Z(IM,JM) ) /2.
c
c       Left and Right Sides
        DO J=1,JMM1
          ZV(1 ,J) = ( Z(1 ,J)+Z(1 ,J+1) ) /2.
          ZV(IM,J) = ( Z(IM,J)+Z(IM,J+1) ) /2.
        ENDDO
c
c       Internal Points.
        DO J=2,JM
        DO I=2,IMM1,2
          ZV(I,J) = ( Z(I  ,J)+Z(I,J-1) 
     X                + Z(I+1,J)+Z(I-1,J) ) / 4.
        ENDDO
        ENDDO
        DO J=1,JMM1
        DO I=3,IMM2,2
          ZV(I,J) = ( Z(I  ,J)+Z(I,J+1) 
     X                + Z(I+1,J)+Z(I-1,J) ) / 4.
        ENDDO
        ENDDO
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE SMEAR
C
C       CHECK routine to test SMEARZ and SMEARV.
c
c      2-D Work fields on A-grid.
       PARAMETER ( IM=9, JM=5, JM2=14 )
       COMMON /WORKK/ WORKK1(IM,JM2), WORKK2(IM,JM2), 
     +                WORKK3(IM,JM2), WORKK4(IM,JM2)
C
C       CHECK ACTION OF SMEARZ.
        DO I=1,IM
        DO J=1,JM*2
        WORKK1(I,J) =0.0 
        WORKK2(I,J) =0.0 
        ENDDO
        ENDDO
c
        DO I=1,IM,2
        DO J=1,JM*2-1,2
        WORKK1(I,J) = FLOAT(I)+FLOAT(J)/100.
        ENDDO
        ENDDO
        DO I=2,IM-1,2
        DO J=2,JM*2,2
        WORKK1(I,J) = FLOAT(I)+FLOAT(J)/100.
        ENDDO
        ENDDO
c
        print* , ' WORKK1 '
        print9929, (( WORKK1(I,J), I=1,IM ), J=JM2,1,-1 )
        CALL SMEARZ(WORKK1,WORKK1,IM,JM)
        print* , ' SMEARZ '
        print9929, (( WORKK1(I,J), I=1,IM ), J=JM2,1,-1 )
 9929   FORMAT( 13f6.2 )
C
c
C       CHECK ACTION OF SMEARV.
        DO I=1,IM
        DO J=1,JM*2
        WORKK3(I,J) =0.0 
        WORKK4(I,J) =0.0 
        ENDDO
        ENDDO
c
        DO J=1,JM*2-1,2
        DO I=2,IM-1,2
        WORKK3(I,J) = FLOAT(I)+FLOAT(J)/100.
        ENDDO
        ENDDO
        DO J=2,JM*2,2
        DO I=1,IM,2
        WORKK3(I,J) = FLOAT(I)+FLOAT(J)/100.
        ENDDO
        ENDDO
c
        print* , ' WORKK3 '
        print9929, (( WORKK3(I,J), I=1,IM ), J=JM2,1,-1 )
        CALL SMEARV(WORKK3,WORKK3,IM,JM)
        print* , ' SMEARV '
        print9929, (( WORKK3(I,J), I=1,IM ), J=JM2,1,-1 )
C
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE SMEARZ(Z1,Z2,IM,JM)
c
c       From a scalar field on a dense A-grid,
c       construct another field on the same grid 
c       where values at vector points are averaged
c       Interpolation is linear at edges and bilinear 
c       for internal points.
c
        REAL Z1(IM,2*JM), Z2(IM,2*JM)
c
        IMM1 = IM-1
        IMM2 = IM-2
        JMM1 = JM-1
        JM2 = 2*JM
        JM2M1 = 2*JM-1
        JM2M2 = 2*JM-2
c
c       First, fill in the E-grid scalar points.
        DO J=1,JM2M1,2
          DO I=1,IM,2
            Z2(I,J) =  Z1(I,J)
          ENDDO
        ENDDO
        DO J=2,JM2,2
          DO I=2,IMM1,2
            Z2(I,J) =  Z1(I,J)
          ENDDO
        ENDDO
c
C       Next, interpolate the intermediate points.
c
c       Bottom Row.
        DO I=2,IMM1,2
          Z2(I,1) = ( Z1(I-1,1)+Z1(I+1,1) ) /2.
        ENDDO
c
c       Top Row.
        Z2(1 ,JM2) = ( Z1(1 ,JM2M1)+Z1(2   ,JM2) ) /2.
        DO I=3,IMM2,2
          Z2(I,JM2) = ( Z1(I-1,JM2)+Z1(I+1,JM2) ) /2.
        ENDDO
        Z2(IM,JM2) = ( Z1(IMM1,JM2)+Z1(IM,JM2M1) ) /2.
c
c       Left and Right Sides
        DO J=2,JM2M2,2
          Z2(1 ,J) = ( Z1(1 ,J-1)+Z1(1 ,J+1) ) /2.
          Z2(IM,J) = ( Z1(IM,J-1)+Z1(IM,J+1) ) /2.
        ENDDO
c
c       Internal Points.
        DO J=2,JM2M2,2
        DO I=3,IMM2,2
          Z2(I,J) = ( Z1(I-1,J)+Z1(I+1,J) 
     X              + Z1(I,J-1)+Z1(I,J+1) ) / 4.
        ENDDO
        ENDDO
        DO J=3,JM2M1,2
        DO I=2,IMM1,2
          Z2(I,J) = ( Z1(I-1,J)+Z1(I+1,J) 
     X              + Z1(I,J-1)+Z1(I,J+1) ) / 4.
        ENDDO
        ENDDO
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE SMEARV(Z1,Z2,IM,JM)
c
c       From a vector field on a dense A-grid,
c       construct another field on the same grid 
c       where values at scaler points are averaged
c       Interpolation is linear at edges and bilinear 
c       for internal points.
c
        REAL Z1(IM,2*JM), Z2(IM,2*JM)
c
        IMM1 = IM-1
        IMM2 = IM-2
        JMM1 = JM-1
        JM2 = 2*JM
        JM2M1 = 2*JM-1
        JM2M2 = 2*JM-2
c
c       First, fill in the E-grid vector points.
        DO J=1,JM2M1,2
          DO I=2,IMM1,2
            Z2(I,J) =  Z1(I,J)
          ENDDO
        ENDDO
        DO J=2,JM2,2
          DO I=1,IM,2
            Z2(I,J) =  Z1(I,J)
          ENDDO
        ENDDO
c
C       Next, interpolate the intermediate points.
c
c       Top Row.
        DO I=2,IMM1,2
          Z2(I,JM2) = ( Z1(I-1,JM2)+Z1(I+1,JM2) ) /2.
        ENDDO
c
c       Bottom Row.
        Z2(1 ,1) = ( Z1(1 ,2)+Z1(2   ,1) ) /2.
        DO I=3,IMM2,2
          Z2(I,1) = ( Z1(I-1,1)+Z1(I+1,1) ) /2.
        ENDDO
        Z2(IM,1) = ( Z1(IMM1,1)+Z1(IM,2) ) /2.
c
c       Left and Right Sides
        DO J=3,JM2M1,2
          Z2(1 ,J) = ( Z1(1 ,J-1)+Z1(1 ,J+1) ) /2.
          Z2(IM,J) = ( Z1(IM,J-1)+Z1(IM,J+1) ) /2.
        ENDDO
c
c       Internal Points.
        DO J=2,JM2M2,2
        DO I=2,IMM1,2
          Z2(I,J) = ( Z1(I-1,J)+Z1(I+1,J) 
     X              + Z1(I,J-1)+Z1(I,J+1) ) / 4.
        ENDDO
        ENDDO
        DO J=3,JM2M1,2
        DO I=3,IMM2,2
          Z2(I,J) = ( Z1(I-1,J)+Z1(I+1,J) 
     X              + Z1(I,J-1)+Z1(I,J+1) ) / 4.
        ENDDO
        ENDDO
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE PSFILT(PS2,OROG, PSEA2,IM,JM2,
     X            ISWEST,ISEAST,JSSOUTH,JSNORTH, EpSpace)
C
c       Filter the pressure at the lowest level: 
C          (1) Convert Surface pressure to MSL Pressure
C          (2) Filter that field
C          (3) Convert MSL pressure to Surface Pressure
C      (Note: PS2 is P2 at the level LM: P2(1,1,LM) in CALL.
C
        REAL PS2(IM,JM2)
        REAL PSEA2(IM,JM2)
        REAL OROG(IM,JM2)
C
c       print95
   95   FORMAT(' UP-N-DOWN SURF PRESSURES, BEFOR  ')
        DO J=JM2,1,-1
c          print98, ( PS2(I,J) , I=1,IM )
   98      FORMAT(13(-2pf6.0))
        ENDDO
c
c       interpolate the surface pressure.
          TZERO = 273. + 15.
          PZERO = 1013.25 * 100.
          GG = 9.80665
          RR = 287.
          GAM = 0.0065
          PWR = GG/(GAM*RR)
c
          DO 150 J=1,JM2
          DO 150 I=1,IM 
             ZG = OROG(I,J)
             DELZ = ZG
             PSEA2(I,J) = PS2(I,J) * (1-(GAM/TZERO)*DELZ)**(-PWR)
c            print*,'P CHECK ', PS2(I,J) , PSEA2(I,J) 
 150      CONTINUE
c
c       print995
  995   FORMAT(' UP-N-DOWN MSL. PRESSURES, BEFOR  ')
        DO J=JM2,1,-1
c          print998, ( PSEA2(I,J) , I=1,IM )
  998      FORMAT(13(-2pf6.0))
        ENDDO
C
C - - -  Call the Filtering Routine.
             IMM1 = IM-1
             JM2M2 = JM2-2 
             CALL IMPFIL(PSEA2,PSEA2,IM,JM2,IM,JM2,
     X            ISWEST,ISEAST,JSSOUTH,JSNORTH,EpSpace)
c
c       interpolate the surface pressure.
           DO 151 J=1,JM2
           DO 151 I=1,IM 
             ZG = OROG(I,J)
             DELZ = ZG
             PS2(I,J) = PSEA2(I,J) * ( 1-(GAM/TZERO)*DELZ )**(+PWR)
CCC          print*,'P CHECK ', PS2(I,J) , PSEA2(I,J) 
 151       CONTINUE
c
c       print97
   97   FORMAT(' UP-N-DOWN SURF PRESSURES, AFTER  ')
        DO J=JM2,1,-1
c          print99, ( PS2(I,J) , I=1,IM )
   99      FORMAT(13(-2pf6.0))
        ENDDO
c       print997
  997   FORMAT(' UP-N-DOWN MSL. PRESSURES, AFTER  ')
        DO J=JM2,1,-1
c          print999, ( PSEA2(I,J) , I=1,IM )
  999      FORMAT(13(-2pf6.0))
        ENDDO
c
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
