C---------------------------------------------------------------------
        SUBROUTINE ENERGY(NSTEP, FI,U,V, FIBAR,IM,JM, DX,DY,
     X                    EK,EA,EP,AK,PK,EM, CGRID)
C
C      CALCULATE THE KINETIC AND POTENTIAL ENERGIES.
C      ( ALL ENERGIES NORMALIZED TO BE PER UNIT AREA)
C
C      (DATA ON A C- OR D- GRID SPECIFIED BY CGRID)
C
C       EK = KINETIC ENERGY
C       EA = AVAILABLE POTENTIAL ENERGY
C       EP = TOTAL POTENTIAL ENERGY
C       AK = EA + EK
C       PK = EP + EK
C       EM = DEVIATION FROM TOTAL MASS (PER UNIT AREA)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
        DIMENSION FI(IM,JM), U(IM,JM), V(IM,JM)
        DIMENSION DX(JM,2)
        REAL MASS, VOLUME, AREA
        LOGICAL CGRID
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      return
C
C     RHO = 1.E+03
      RHO = 1.0
C
      GRAV = 9.81
      AREA = 0.
      VOLUME = 0.
C
      EK1 = 0.
      EK2 = 0.
      EA = 0.
      EP = 0.
      EM = 0.
C
      DO 100 II=2,IM-1
      DO 100 JJ=2,JM-1
        DAREA = DX(JJ,2)*DY
        AREA = AREA + DAREA
        FIIJ = FI(II,JJ)
        FIDASH = FIIJ - FIBAR
        HIJ = FIIJ/GRAV
        DVOLUM = DAREA*HIJ
        VOLUME = VOLUME + DVOLUM
        IF ( CGRID ) THEN
           UIJ = ( U(II,JJ)+ U(II-1,JJ) ) / 2.
           VIJ = ( V(II,JJ)+ V(II,JJ-1) ) / 2.
           FIU = (FI(II,JJ)+FI(II+1,JJ) ) / 2.
           FIV = (FI(II,JJ)+FI(II,JJ+1) ) / 2.
        ELSE 
           UIJ = ( U(II,JJ)+ U(II,JJ-1) ) / 2.
           VIJ = ( V(II,JJ)+ V(II-1,JJ) ) / 2.
           FIU = (FI(II,JJ)+FI(II,JJ+1) ) / 2.
           FIV = (FI(II,JJ)+FI(II+1,JJ) ) / 2.
        END IF 
        EK1 = EK1 + DAREA * ( UIJ**2+VIJ**2 ) * FIIJ
        EK2 = EK2 + DAREA * ( U(II,JJ)**2*FIU + V(II,JJ)**2*FIV )
        EA = EA + DAREA * FIDASH**2
        EP = EP + DAREA * FIIJ**2
C       EM = EM + DAREA * FIIJ
        EM = EM + DAREA * FIDASH
100   CONTINUE
C
C     EK = EK1 * ( RHO/(2.*GRAV) )
      EK = EK2 * ( RHO/(2.*GRAV) )
      EA = EA * ( RHO/(2.*GRAV) )
      EP = EP * ( RHO/(2.*GRAV) )
      EM = EM * ( RHO/GRAV )
C
      AK = EK + EA
      PK = EK + EP
C
C     GET THE MEAN ENERGIES PER UNIT AREA
      EK = EK / AREA
      EA = EA / AREA
      EP = EP / AREA
      AK = AK / AREA
      PK = PK / AREA
      EM = EM / AREA
C
      TYPE 9876, NSTEP, EK,EA,EP, AK,PK, EM
 9876 FORMAT(/' *** ENERGY ***  STEP:',I4,
     X        ' KE APE PE APLUSK PPLUSK MASS '/
     X                     1P6E12.3/)
C
      MASS = RHO*VOLUME
      TYPE 9877, MASS, VOLUME, RHO, AREA
 9877 FORMAT(' MASS(KG), VOLUME(M**3), DENSITY(KG/M**3), AREA(M**2)'/
     X         1P4E12.3/)
C
      RETURN
      END
C-------------------------------------------------------------------
