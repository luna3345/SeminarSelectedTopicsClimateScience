      SUBROUTINE PLTMRK(X,Y,SIZE,ICODE,IPEN)
c --- [9-May-1985]
C
C *** PLOT A MARK AT POSITION (X,Y) **********************************
c
c-----------------------------------------------------------------------
c  This routine plots a mark (a symbol) at position (X,Y) of size
c  SIZE. ICODE determines the type of mark that is drawn. For ICODE
c  negative a circle is drawn as well.
c-----------------------------------------------------------------------
C
C#IF 'ICODE' IS ZERO RETURN
      IF(ICODE.EQ.0) GOTO 50
C
C#IF ICODE IS NEGATIVE DRAW CIRCLE OF DIAMETER 'SIZE'
      IF(ICODE.GT.0)  GOTO 15
      RADIUS=0.5*SIZE
      NBITS =32
      DANGLE=6.28318530718/FLOAT(NBITS)
      ANGLE =DANGLE
      XCIRCL=RADIUS
      YCIRCL=0.0
      Call Move(X+XCIRCL,Y+YCIRCL,0)
      DO 10 I=1,NBITS
         XCIRCL=RADIUS*COS(ANGLE)
         YCIRCL=RADIUS*SIN(ANGLE)
         Call Move(X+XCIRCL,Y+YCIRCL,IPEN)
         ANGLE=ANGLE+DANGLE
   10    CONTINUE
C
C#DRAW THE APPROPIATE SYMBOL
   15 CONTINUE
      IGOTO=ICODE
      IF(IGOTO.LT.0) IGOTO=-IGOTO
      IF(IGOTO.GT.6) GOTO 50
      GOTO (20,25,30,35,40,45),IGOTO
C
C#A DIAMOND SHAPE
   20 CONTINUE
      RADIUS=0.5*SIZE
      Call Move(X+RADIUS,Y,0)
      Call Move(X,Y+RADIUS,IPEN)
      Call Move(X-RADIUS,Y,IPEN)
      Call Move(X,Y-RADIUS,IPEN)
      Call Move(X+RADIUS,Y,IPEN)
      GOTO 50
C
C#A SQUARE BOX
   25 CONTINUE
      SIDE = SIZE / (2.0*SQRT(2.0))
      Call Move(X+SIDE,Y+SIDE,0)
      Call Move(X-SIDE,Y+SIDE,IPEN)
      Call Move(X-SIDE,Y-SIDE,IPEN)
      Call Move(X+SIDE,Y-SIDE,IPEN)
      Call Move(X+SIDE,Y+SIDE,IPEN)
      GOTO 50
C
C#A CROSS
   30 CONTINUE
      RADIUS=0.5*SIZE
      Call Move(X,Y+RADIUS,0)
      Call Move(X,Y-RADIUS,IPEN)
      Call Move(X+RADIUS,Y,0)
      Call Move(X-RADIUS,Y,IPEN)
      GOTO 50
C
C#AN X
   35 CONTINUE
      SIDE = SIZE / (2.0*SQRT(2.0))
      Call Move(X+SIDE,Y+SIDE,0)
      Call Move(X-SIDE,Y-SIDE,IPEN)
      Call Move(X-SIDE,Y+SIDE,0)
      Call Move(X+SIDE,Y-SIDE,IPEN)
      GOTO 50
C
C#A DOT (ACTUALLY A SMALL DIAMOND)
   40 CONTINUE
      RADIUS=0.1*SIZE
      Call Move(X+RADIUS,Y,0)
      Call Move(X,Y+RADIUS,IPEN)
      Call Move(X-RADIUS,Y,IPEN)
      Call Move(X,Y-RADIUS,IPEN)
      Call Move(X+RADIUS,Y,IPEN)
      GOTO 50
C
C#AN EQUILATERAL TRIANGLE
   45 CONTINUE
      RADIUS=0.5*SIZE
      ANGL30=3.14159265359/6.0
      RCOS30=RADIUS*COS(ANGL30)
      RSIN30=RADIUS*SIN(ANGL30)
      Call Move(X,Y+RADIUS,0)
      Call Move(X+RCOS30,Y-RSIN30,IPEN)
      Call Move(X-RCOS30,Y-RSIN30,IPEN)
      Call Move(X,Y+RADIUS,IPEN)
      GOTO 50
C
C#END OF PLOTTING - RETURN
   50 CONTINUE
C
      RETURN
      END
