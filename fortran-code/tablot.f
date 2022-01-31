C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE TABLOT(P,MU,MV,TTOP,OROG,IM,JM,LM,MESAGE)
c
c       PRINT A TABLE OF CENTRAL VALUES.
c
        CHARACTER*10 MESAGE
        REAL P(IM,JM,LM), MU(IM,JM,LM), MV(IM,JM,LM)
        REAL TTOP(IM,JM)
        REAL OROG(IM,2*JM)
c
c       Print a HEADER RECORD
        print9800, MESAGE
        print9901, ' Table of CENTRAL Values '
c       print9900
 9800   FORMAT(A10)
 9900   FORMAT(1X)
 9901   FORMAT(A50)
 9802   FORMAT('_______________________________',
     +          '_____________________________________________')
 9902   FORMAT('|______________|______________|',
     +          '______________|______________|______________|')
c
        print9802
        DO J=5,3,-1
c
          print901, ( TTOP(I,J) , I=6,8,2 )
          DO L=1,5
             print902, (MU(I,J,L)/100,MV(I,J,L)/100.0 , I=5,9,2 )
             print903, ( P(I,J,L) , I=6,8,2 )
          ENDDO
          JU = 2*J
          print916, (OROG(I,JU),I=5,9,1)
  901 FORMAT('| ',12x,' | ',6x,f6.0,' | ',12x,' | ',
     +            6x,f6.0,' | ',12x,' |')
  902 FORMAT('| ',   2f6.0,' | ',12x,' | ',2f6.0,' | ',
     +            12x,' | ',2f6.0,' |')
  903 FORMAT('| ',12x,' | ',   f10.1,2x,' | ',12x,' | ',
     +            f10.1,2x,' | ',12x,' |')
  916 FORMAT('| ',f12.0,' | ',f12.0,' | ',f12.0,' | ',
     +            f12.0,' | ',f12.0,' |')
c
        print9902
          print904, ( TTOP(I,J) , I=5,9,2 )
          DO L=1,5
             print905, (MU(I,J,L)/100,MV(I,J,L)/100 , I=6,8,2 )
             print906, ( P(I,J,L) , I=5,9,2 )
          ENDDO
          JL = 2*J-1
          print916, (OROG(I,JL),I=5,9,1)
  904 FORMAT('| ',6x,f6.0,' | ',12x,' | ',6x,
     +            f6.0,' | ',12x,' | ',6x,f6.0,' |')
  905 FORMAT('| ',12x,' | ',   2f6.0,' | ',12x,' | ',
     +            2f6.0,' | ',12x,' |')
  906 FORMAT('| ',   f10.1,2x,' | ',12x,' | ',
     +             f10.1,2x,' | ',12x,' | ',f10.1,2x,' |')
c
          print9902
        ENDDO
C
c        CALL TEXTAB(P,MU,MV,TTOP,OROG,IM,JM,LM)
c
        RETURN
        END
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
        SUBROUTINE TEXTAB(P,MU,MV,TTOP,OROG,IM,JM,LM)
c
c==============  TABLE Output for TeX  ====================
c
        REAL P(IM,JM,LM), MU(IM,JM,LM), MV(IM,JM,LM)
        REAL TTOP(IM,JM)
        REAL OROG(IM,2*JM)

        PARAMETER( IMM= 13, JMM=07, LMM=5 )
        INTEGER IP(IMM,JMM,LMM),IMU(IMM,JMM,LMM),IMV(IMM,JMM,LMM)
        INTEGER ITTOP(IMM,JMM)
        INTEGER IOROG(IMM,2*JMM)
c
        DO I=1,IM
        DO J=1,JM
          DO L=1,LM
            IP (I,J,L) = NINT( P (I,J,L)/100. )
            IMU(I,J,L) = NINT( MU(I,J,L)/100. )
            IMV(I,J,L) = NINT( MV(I,J,L)/100. )
          ENDDO
          ITTOP(I,J) = NINT( TTOP(I,J) )
        ENDDO
        ENDDO

        DO I=1,IM
        DO J=1,2*JM
          IOROG(I,J) = NINT( OROG(I,J) )
        ENDDO
        ENDDO
c
C       "Blank" Line.
c       print99802
99802   FORMAT('|&|&|&|&|& \\nr ')
99842   FORMAT(' ')
c
c------------------------- Header Block
c
        print99802
        print90052
        print99803
90052   FORMAT('|$ 5\\Deg E$&|$ 8\\Deg E$&|$11\\Deg E$',
     X        '&|$14\\Deg E$&|$17\\Deg E$     & \\nr ')
99803   FORMAT('|&|&|&|&|& \\cr ')
c
c------------------------- First Two Blocks
        DO J=5,5,-1
c
          JU = 2*J
          JL = 2*J-1
c
c         print8401, ( TTOP(I,J) , I=6,8,2 )
          print99802
          DO L=1,5
             IF(L.EQ.3) THEN
             print8683,  100*IMU(7,J,L),100*IMV(7,J,L) 
             ELSE
             print8613,  100*IMU(7,J,L),100*IMV(7,J,L) 
             ENDIF
          ENDDO
C         JU = 2*J
C         print8416, IOROG(7,JU)
          print99803
 8401 FORMAT('|&|& ',I5,'|&|& ',I5,'|& \\nr ')
 8416 FORMAT('|&|&|&',I5,'|&|& \\cr ')
 8613 FORMAT('|&|&|\\bf ',   I5,'&\\bf ',I5,'|&|& \\nr ')
 8683 FORMAT('$54.0\\Deg N$ |&|&|\\bf ',   I5,
     X                         '&\\bf ',   I5,'|&|& \\nr ')
c
        print99842
c
          print8504, ( ITTOP(I,J) , I=7,7,1 )
          DO L=1,5
             IF(L.EQ.3) THEN
             print8685, 
     +          100*IMU(6,J,L), 100*IP(7,J,L), 100*IMU(8,J,L)
             ELSE
             print8615, 
     +          100*IMU(6,J,L), 100*IP(7,J,L), 100*IMU(8,J,L)
             ENDIF
          ENDDO
C         JL = 2*J-1
          print8516, (IOROG(I,JL),I=7,7,1)
 8504 FORMAT('|&|&|&',I5,'|&|&\\nr ')
 8516 FORMAT('|&|&|& ',I5,'|&|& \\cr ')
 8615 FORMAT('|&|\\bf ',   I5,'&|',
     X          I5,'&|\\bf ',I5,'&|& \\nr ')
 8685 FORMAT('$52.2\\Deg N$ |&|\\bf ',   I5,'&|',
     X          I5,'&|\\bf ',I5,'&|& \\nr ')
c
        print99842
c
        ENDDO
C
c------------------------- Next Two Blocks
        DO J=4,4,-1
c
          JU = 2*J
          JL = 2*J-1
c
          print9401, ( ITTOP(I,J) , I=6,8,2 )
          DO L=1,5
c            print9402, (MU(I,J,L),MV(I,J,L) , I=5,9,2 )
c            print9403, ( P(I,J,L) , I=6,8,2 )
           IF(L.EQ.10) THEN
             print9603, 
     X        MU(5,J,L),MV(5,J,L) , P(6,J,L) , 100*TTOP(6,J),
     X        MU(7,J,L),MV(7,J,L) , P(8,J,L) , 100*TTOP(8,J),
     X        MU(9,J,L),MV(9,J,L)  
           ELSE IF(L.EQ.50) THEN
             print9693, 
     X        MU(5,J,L),MV(5,J,L) , P(6,J,L) , 100*OROG(6,JU),
     X        MU(7,J,L),MV(7,J,L) , P(8,J,L) , 100*OROG(8,JU),
     X        MU(9,J,L),MV(9,J,L)  
           ELSE IF(L.EQ.3) THEN
             print9683, 
     X        IMU(5,J,L),IMV(5,J,L) , IP(6,J,L) ,
     X        IMU(7,J,L),IMV(7,J,L) , IP(8,J,L) ,
     X        IMU(9,J,L),IMV(9,J,L)  
           ELSE 
             print9613, 
     X        IMU(5,J,L),IMV(5,J,L) , IP(6,J,L) ,
     X        IMU(7,J,L),IMV(7,J,L) , IP(8,J,L) ,
     X        IMU(9,J,L),IMV(9,J,L)  
           ENDIF
          ENDDO
C         JU = 2*J
          print9416, (IOROG(I,JU),I=6,8,2)
 9401 FORMAT('|&|& ',I5,'\\hfill|&|& ',I5,'\\hfill|& \\nr ')
C9402 FORMAT('|',   I5,'&',I5,'|&|',I5,'&',I5,
C    X                        '|&|',2x,I5,'&',I5,' \\nr ')
C9403 FORMAT('|&|',   f5.1,'&|&|',f6.1,'&|&\\nr ')
 9416 FORMAT('|&|&',I5,'|&|&',I5,'|& \\cr ')
 9603 FORMAT('|',   I5,'&',I5,'|',I5,'& \\it ',I5,'|',
     X                   I5,'&',I5,
     X                        '|',I5,'& \\it ',I5,'|',
     X                   I5,'&',I5, ' \\nr ')
 9693 FORMAT('|',   I5,'&',I5,'|',I5,'& \\bf ',I5,'|',
     X                   I5,'&',I5,
     X                        '|',I5,'& \\bf ',I5,'|',
     X                   I5,'&',I5, ' \\cr ')
 9613 FORMAT('|\\bf ',   I5,'&\\bf ',I5,'|',I5,
     X            '&|\\bf ',I5,'&\\bf ',I5,
     X   '|',I5,'&|\\bf ',I5,'&\\bf ',I5,' \\nr ')
 9683 FORMAT('$50.4\\Deg N$ |\\bf ',   I5,'&\\bf ',I5,'|',
     X             I5,'&|\\bf ',I5,'&\\bf ',I5,
     X            '|',I5,'&|\\bf ',I5,'&\\bf ',I5,' \\nr ')
c
C       print99802
        print99842
c
          print9504, ( ITTOP(I,J) , I=5,9,2 )
          DO L=1,5
C            print9505, (MU(I,J,L),MV(I,J,L) , I=6,8,2 )
C            print9506, ( P(I,J,L) , I=5,9,2 )
           IF(L.EQ.10) THEN
             print9605, 
     X        P(5,J,L), 100*TTOP(5,J), MU(6,J,L),MV(6,J,L) ,
     X        P(7,J,L), 100*TTOP(7,J), MU(8,J,L),MV(8,J,L) ,
     X        P(9,J,L), 100*TTOP(9,J)
           ELSE IF(L.EQ.50) THEN
             print9695, 
     X        P(5,J,L), 100*OROG(5,JL), MU(6,J,L),MV(6,J,L) ,
     X        P(7,J,L), 100*OROG(7,JL), MU(8,J,L),MV(8,J,L) ,
     X        P(9,J,L), 100*OROG(9,JL)
           ELSE IF(L.EQ.3) THEN
             print9685, 
     X        IP(5,J,L),  IMU(6,J,L),
     X        IP(7,J,L),  IMU(8,J,L),
     X        IP(9,J,L)
           ELSE
             print9615, 
     X        IP(5,J,L),  IMU(6,J,L),
     X        IP(7,J,L),  IMU(8,J,L),
     X        IP(9,J,L)
           ENDIF
          ENDDO
C         JL = 2*J-1
          print9516, (IOROG(I,JL),I=5,9,2)
 9504 FORMAT('|&',I5,'|&|&',I5,
     X                '|&|&',I5,'\\nr ')
C9505 FORMAT('|&|',   I5,'&',I5,'|&|',I5,'&',I5,'|& \\nr ')
C9506 FORMAT('|&',   f5.1,'&|&|',f6.1,' \\nr ',
C    X               '&|&|',f5.1,' \\nr ')
C9516 FORMAT('|&\\bf ',I5,'|&\\bf ',I5,'|&\\bf ',I5,
C    X                '|&\\bf ',I5,'|&\\bf ',I5,' \\cr ')
 9516 FORMAT('|& ',I5,'|&|& ',I5,'|&|& ',I5,' \\cr ')
 9605 FORMAT('|',   I5,'& \\it ',I5,'|',
     X          I5,'&',I5,'|',
     X          I5,'& \\it ',I5,'|',
     X          I5,'&',I5,'|',
     X          I5,'& \\it ',I5,'\\nr ')
 9695 FORMAT('|',   I5,'& \\bf ',I5,'|',
     X          I5,'&',I5,'|',
     X          I5,'& \\bf ',I5,'|',
     X          I5,'&',I5,'|',
     X          I5,'& \\bf ',I5,'\\cr ')
 9615 FORMAT('|',   I5,'&|\\bf ',I5,'&|',
     X          I5,'&|\\bf ',I5,'&|',I5,'& \\nr ')
 9685 FORMAT('$48.6\\Deg N$|',   I5,'&|\\bf ',I5,'&|',
     X          I5,'&|\\bf ',I5,'&|',I5,'& \\nr ')
c
C       print99802
        print99842
c
        ENDDO
C
c------------------------- LAST Two Blocks
        DO J=3,3,-1
c
          JU = 2*J
          JL = 2*J-1
c
          print7401, ( ITTOP(I,J) , I=6,8,2 )
          DO L=1,5
             IF(L.EQ.3) THEN
             print7683, 
     X         IP(6,J,L), IMU(7,J,L),IMV(7,J,L) , IP(8,J,L)
             ELSE
             print7613, 
     X         IP(6,J,L), IMU(7,J,L),IMV(7,J,L) , IP(8,J,L)
             ENDIF
          ENDDO
C         JU = 2*J
          print7416, (IOROG(I,JU),I=6,8,2)
 7401 FORMAT('|&|& ',I5,'\\hfill|&|& ',I5,'\\hfill|& \\nr ')
 7416 FORMAT('|&|&',I5,'|&|&',I5,'|& \\cr ')
 7613 FORMAT('|&|',   I5,'&|\\bf ',I5,'&\\bf ',
     X                I5,'|',I5,'&|& \\nr ')
 7683 FORMAT('$46.8\\Deg N$|&|',   I5,'&|\\bf ',I5,'&\\bf ',
     X                           I5,'|',I5,'&|& \\nr ')
c
        print99842
c
          print7504, ( ITTOP(I,J) , I=7,7,1 )
          DO L=1,5
             IF(L.EQ.3) THEN
             print7685, IP(7,J,L)
             ELSE
             print7615, IP(7,J,L)
             ENDIF
          ENDDO
C         JL = 2*J-1
          print7516, (IOROG(I,JL),I=7,7,1)
 7504 FORMAT('|&|&|&',I5,'|&|& \\nr ')
 7516 FORMAT('|&|&|&',I5,'|&|& %%% \\cr  ')
 7615 FORMAT('|&|&|',    I5,'&|&|& \\nr ')
 7685 FORMAT('$45.0\\Deg N$|&|&|',    I5,'&|&|& \\nr ')
c
        print99842
c
        ENDDO
C
        RETURN
        END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
