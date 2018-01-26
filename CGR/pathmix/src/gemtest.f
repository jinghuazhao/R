C              G E M I N I    T E S T    R O U T I N E
C              -----------    -------    -------------
C
C DRIVING PROCEDURE FOR  G E M I N I, VERSION 1, FEB. 79, J.M.LALOUEL
C
C REVISED FOR USE WITH ALMINI LIBRARY BY S.R. (3/85)
C
C *********************************************************************
C
      SUBROUTINE GEMTEST
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MAXPAR=5, MAXITP=MAXPAR, MAXCC=1)
      DIMENSION   XALL(MAXPAR),ITP(MAXPAR),SE(MAXPAR),BND(MAXPAR,2)
      DIMENSION   TL(MAXITP,MAXITP),WORK(MAXITP,16),WRKNLC(MAXCC,4)
C
C ******* TST BLOCK IS TEMPORARY FOR MULTIPLE RUNS *******
      COMMON /TST/   IRUN
C
C DEFINE LOGICAL FILES (IQ=TERSE FILE, IP=PROLIX FILE)
C
      IQ = 21
      OPEN (IQ, FILE='GEMTEST.OUT')
C
      IP = 22
      OPEN (IP, FILE='GEMTEST.PLX')
C
C LOOP ON SEVERAL TESTS
C
      DO 10 IRUN=1,10
         CALL INPUTgem(NALL,XALL,ITP,BND,MAXPAR,
     +      IBND,IHESS,IVAR,IHX,IQOB,IVERB,MAXIT,H,TOL,TRUPB,INLC,IQ)
         CALL ALMINI(FX,NALL,XALL,ITP,SE,BND,
     +      IBND,IHESS,IVAR,IHX,IQOB,IVERB,MAXIT,IP,H,TOL,TRUPB,
     +      NIT,NFE,PTG,IDG,TL,MAXPAR,MAXITP,WORK,
     +      INLC,MAXCC,ISCAL,RATE,SCL,EPSCC,FC,FXC,NK,CC,IAC,WRKNLC,2)
         CALL OUTPUTgem(FX,NALL,XALL,SE,IVAR,NIT,NFE,PTG,IDG,IQ)
10    CONTINUE
C
      RETURN
      END
      SUBROUTINE INPUTgem(NALL,XALL,ITP,BND,MAXPAR,
     +   IBND,IHESS,IVAR,IHX,IQOB,IVERB,MAXIT,H,TOL,TRUPB,INLC,IQ)
C---
C--- READ INPUT DATA AND PARAMETERS
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION   XALL(*),ITP(*),BND(MAXPAR,2)
C
      CHARACTER TITLE*80
C
      COMMON /TST/   IRUN
C
C DEFINE MACHINE PRECISION AS 2**-P WHERE P IS SIZE OF MANTISSA
C
CCC      PRECIS = 2.**(-38)   ! 38 BIT MANTISSA ON HARRIS COMPUTERS
C
      IPOWER = 0
10    PRECIS = 2.0D0**IPOWER
      IF (1.0D0+PRECIS .NE. 1.0D0) THEN
         IPOWER = IPOWER - 1
         GO TO 10
      END IF
C
C DEFINE INPUT PARAMETERS FOR TEST PROBLEMS
C ONE MAY DECIDE TO READ THE FOLLOWING FROM AN INPUT JOB FILE
C
      H = SQRT(PRECIS)
      TRUPB = SQRT(H)
      TOL = 0.0
C
      IBND  = 0
      IHESS = 0
      IVAR  = 0
      IHX   = 0
      IQOB  = 0
      IVERB = 3
      INLC  = 0
C
      IF (IRUN .EQ. 1) THEN
         TITLE = 'ROSENBROCK PARABOLIC VALLEY P1'
         NALL = 2
         XALL(1) = -1.2
         XALL(2) = 1.0
C
      ELSE IF (IRUN .EQ. 2) THEN
         TITLE = 'FLETCHER HELICAL VALLEY P2'
         NALL = 3
         XALL(1) = -1.0
         XALL(2) = 0.0
         XALL(3) = 0.0
C
      ELSE IF (IRUN .EQ. 3) THEN
         TITLE = 'POWELL FUNCTION OF FOUR VARIABLES'
         NALL = 4
         XALL(1) = 3.0
         XALL(2) = -1.0
         XALL(3) = 0.0
         XALL(4) = 1.0
C
      ELSE IF (IRUN .EQ. 4) THEN
         TITLE = 'A SYSTEM NEARLY LINEAR NEAR THE ROOTS P4'
         NALL = 2
         XALL(1) = 0.1
         XALL(2) = 2.0
C
      ELSE IF (IRUN .EQ. 5) THEN
         TITLE = 'A GAUSSIAN FUNCTION P5'
         NALL = 3
         XALL(1) = 0.4
         XALL(2) = 1.0
         XALL(3) = 0.0
C
      ELSE IF (IRUN .EQ. 6) THEN
         TITLE = 'BOX FUNCTION P7'
         NALL = 3
         XALL(1) = 0.0
         XALL(2) = 20.0
         XALL(3) = 20.0
C
      ELSE IF (IRUN .EQ. 7) THEN
         TITLE = 'WOOD FUNCTION P8'
         NALL = 4
         XALL(1) = -3.0
         XALL(2) = -1.0
         XALL(3) = -3.0
         XALL(4) = -1.0
C
      ELSE IF (IRUN .EQ. 8) THEN
         TITLE = 'A VARIATION OF POWELL SECOND FUNCTION P9'
         NALL = 3
         XALL(1) = 0.0
         XALL(2) = 1.0
         XALL(3) = 2.0
C
      ELSE IF (IRUN .EQ. 9) THEN
         TITLE = 'A VARIABLE DIMENSIONED PROBLEM (HERE WITH N=5) P10'
         NALL = 5
         XALL(1) = 0.1
         XALL(2) = 0.1
         XALL(3) = 0.1
         XALL(4) = 0.1
         XALL(5) = 0.1
C
      ELSE IF (IRUN .EQ. 10) THEN
         TITLE = 'CHEBYQUAD ( USED HERE WITH N=4 )'
         NALL = 4
         XALL(1) = 0.2
         XALL(2) = 0.4
         XALL(3) = 0.6
         XALL(4) = 0.8
      END IF
C
C ALLOW FOR SPECIFICATION OF ITERATED PARAMETERS, SUCH THAT
C ITP(I)=1 IF I-TH PARAMETER ITERATED
C
      N = 0
      DO 100 I=1,NALL
         ITP(I) = 1
         N = N + 1
100   CONTINUE
C
C DEFINE MAXIT BASED ON NUMBER OF ITERATED PARAMETERS
C
      MAXIT = 20*N
C
C WRITE INPUT DATA TO TERSE FILE
C
      call dblepr (TITLE,-1,0,0)
      WRITE (IQ,1010) TITLE
      WRITE (IQ,1020) NALL,IBND,IHESS,H,TRUPB,TOL
      WRITE (IQ,1030) (XALL(I),I=1,NALL)
      WRITE (IQ,1040) (ITP(I),I=1,NALL)
C
 1010 FORMAT(/1X,100('*'),
     +      //10X,A)
 1020 FORMAT(/10X,'NALL=',I3,'  IBND=',I3,'  IHESS=',I3,
     +          '  H=',E10.3,'  TRUPB=',E10.3,'  TOLERANCE=',E10.3)
 1030 FORMAT(/10X,'INITIAL VALUES OF PARAMETERS:',
     +      /(10X,10F8.5))
 1040 FORMAT(/10X,'IT. INDICATORS  ',50I2)
C
      RETURN
      END
      SUBROUTINE OUTPUTgem(FX,NALL,XALL,SE,IVAR,NIT,NFE,PTG,IDG,IQ)
C---
C--- WRITE FINAL OUTPUT ON LOGICAL UNIT IQ
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION   XALL(NALL),SE(NALL)
C
      WRITE (IQ,1000)
      WRITE (IQ,1010) NIT,NFE,FX,PTG
      WRITE (IQ,1020) (XALL(I),I=1,NALL)
      IF (IVAR.GT.0)  WRITE (IQ,1030) (SE(I),I=1,NALL)
C
 1000 FORMAT(/10X,33('-'),' FINAL OUTPUT ',33('-'))
 1010 FORMAT(/10X,'NIT=',I4,'  NFE=',I5,'  FINAL F=',E18.11,
     +       /35X,'PTG=',E18.11)
 1020 FORMAT(/10X,'FINAL XALL=', /(10X,6E18.11))
 1030 FORMAT(/10X,'STANDARD-ERRORS=', /(10X,6E18.11))
C
C --------------------------------------------------------------------
C
C INDICATE REASON FOR TERMINATION
C
      WRITE (IQ,1500) IDG
 1500 FORMAT(/10X,'CONVERGENCE STATUS:',
     +      //10X,'ITERATIVE PROCESS TERMINATED BECAUSE:   ',
     +        '(FINAL IDG =',I2,')')
C
      GO TO (900,901,902,903,904,905,906,907,908,909), IDG+1
      call intpr('INVALID IDG',-1,0,0)
      STOP
C
900   WRITE (IQ,1900)
      GO TO 999
901   WRITE (IQ,1901)
      GO TO 999
902   WRITE (IQ,1902)
      GO TO 999
903   WRITE (IQ,1903)
      GO TO 999
904   WRITE (IQ,1904)
      GO TO 999
905   WRITE (IQ,1905)
      GO TO 999
906   WRITE (IQ,1906)
      GO TO 999
907   WRITE (IQ,1907)
      GO TO 999
908   WRITE (IQ,1908)
      GO TO 999
909   WRITE (IQ,1909)
      GO TO 999
C
 1900 FORMAT(10X,'*** MAXIMUM POSSIBLE ACCURACY REACHED ***')
 1901 FORMAT(10X,'*** SEARCH DIRECTION NOT DOWNHILL ANYMORE ***')
 1902 FORMAT(10X,'*** ACCUMULATION OF ROUNDING ERRORS PREVENTS ',
     +               'FURTHER PROGRESS ***')
 1903 FORMAT(10X,'*** ALL SIGNIFICANT DIGITS LOST THROUGH ',
     +               'CANCELLATION IN OPTIMAL CONDITIONING ***')
 1904 FORMAT(10X,'*** SPECIFIED TOLERANCE ON NORMALIZED GRADIENT ',
     +               'WAS MET ***')
 1905 FORMAT(10X,'*** EXCESSIVE CANCELLATION IN GRADIENT ***')
 1906 FORMAT(10X,'*** MAXIMUM NUMBER OF ITERATIONS REACHED ***')
 1907 FORMAT(10X,'*** EXCESSIVE CANCELLATION IN GRADIENT ',
     +               'CALCULATION ***')
 1908 FORMAT(10X,'*** NO ITERATED PARAMETERS ***')
 1909 FORMAT(10X,'*** PARAMETER SET TO A BOUND ***')
C
999   CONTINUE
C
      RETURN
      END
      SUBROUTINE FUNgem(FX,V,NALL,CC,NC)
C---
C--- TEST FUNCTIONS TO BE MINIMIZED
C---
C--- WHAT IS XALL ELSEWHERE IS RENAMED V IN THIS SUBROUTINE
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION   V(NALL)
C
      PARAMETER (MAXPAR=5)
      DIMENSION   Y(15),Z(15),YY(MAXPAR),TI(MAXPAR),TIM(MAXPAR)
      LOGICAL     IEVEN
C
      COMMON /TST/ IRUN
C
      DATA PI /3.1415926535897932385/
C
      GO TO (1,2,3,4,5,6,7,8,9,10), IRUN
C
C ROSENBROCK PARABOLIC VALLEY P1
C
1     FX = (1.-V(1))*(1.-V(1)) + 100.*(V(2)-V(1)*V(1))*(V(2)-V(1)*V(1))
      GO TO 999
C
C FLETCHER HELICAL VALLEY P2
C
2     R = SQRT(V(1)*V(1)+V(2)*V(2))
      IF (V(1)) 21,22,23
   21 U = ATAN(V(2)/V(1))/(2.*PI) + 0.5
      GOTO 24
   22 U = 0.25
      GOTO 24
   23 U = ATAN(V(2)/V(1))/(2.*PI)
   24 CONTINUE
      F1 = 10.*(V(3)-10.*U)
      F2 = 10.*(R-1.)
      F3 = V(3)
      FX = F1*F1 + F2*F2 + F3*F3
      GO TO 999
C
C POWELL FUNCTION OF 4 VARIABLES P3
C
3     FX = (V(1)+10.*V(2)) * (V(1)+10.*V(2)) +
     +    5.*(V(3)-V(4)) * (V(3)-V(4)) +
     +    ((V(2)-2.*V(3))*(V(2)-2.*V(3))) *
     +    ((V(2)-2.*V(3))*(V(2)-2.*V(3))) +
     +    10.*((V(1)-V(4))*(V(1)-V(4))) * ((V(1)-V(4))*(V(1)-V(4)))
      GO TO 999
C
C A SYSTEM NEARLY LINEAR NEAR THE ROOTS P4
C
4     F1 = V(1)*V(1) - V(2) - 1.0
      F2 = (V(1)-2.)*(V(1)-2.) + (V(2)-0.5)*(V(2)-0.5) - 1.0
      FX = F1*F1 + F2*F2
      GO TO 999
C
C A GAUSSIAN FUNCTION P5
C
5     Z(1) = 3.5
      DO 51 I=2,15
         Z(I) = Z(I-1)-0.5
   51 CONTINUE
      Y(1) =  0.0009
      Y(15) = 0.0009
      Y(2)  = 0.0044
      Y(14) = 0.0044
      Y(3)  = 0.0175
      Y(13) = 0.0175
      Y(4)  = 0.0540
      Y(12) = 0.0540
      Y(5)  = 0.1295
      Y(11) = 0.1295
      Y(6)  = 0.2420
      Y(10) = 0.2420
      Y(7)  = 0.3521
      Y(9)  = 0.3521
      Y(8)  = 0.3989
      FX = 0.0
      DO 52 I=1,15
         U = V(1) * EXP(-(Z(I)-V(3))*(Z(I)-V(3))*V(2)/2.0) - Y(I)
         FX = FX + U*U
   52 CONTINUE
      FX = FX - 1.12793277E-08
      GO TO 999
C
C BOX FUNCTION P7
C
6     FX = 0.0
      DO 61 I=1,10
         DI = 0.1*I
         FI = EXP(-V(1)*DI)-EXP(-V(2)*DI)-V(3)*(EXP(-DI)-EXP(-10.*DI))
         FX = FX + FI*FI
   61 CONTINUE
      GO TO 999
C
C WOOD FUNCTION P8
C
7     FX = 100.*(V(2)-V(1)*V(1))**2 + (1.-V(1))**2 +
     +     90.*(V(4)-V(3)*V(3))**2 + (1.-V(3))**2 +
     +     10.1*((V(2)-1.)**2  + (V(4)-1.)**2) +
     +     19.8*(V(2)-1.)*(V(4)-1.)
      GO TO 999
C
C A VARIATION OF POWELL SECOND FUNCTION P9
C
8     F1 = -1.0 + 1.0/(1.0+(V(1)-V(2))*(V(1)-V(2)))
      F2 = SIN(.5*PI*V(2)*V(3)) - 1.0
      F3 = EXP(-((V(1)+V(3))/V(2)-2.)*((V(1)+V(3))/V(2)-2.)) - 1.0
      FX = F1*F1 + F2*F2 + F3*F3
      GO TO 999
C
C A VARIABLE DIMENSIONED PROBLEM (USED HERE WITH N=5) P10
C
9     F1 = 0.0
      F2 = 0.0
      DO 91 I = 1,NALL
         F1 = F1+V(I)*V(I)
         RI = I
         F2 = F2+SQRT(RI)*V(I)
   91 CONTINUE
      F3 = F2*F2
      FX = F1+F3+F3*F3
      GO TO 999
C
C CHEBYQUAD
C
10    FN = NALL
      DELTA = 0.0
      DO 110 J = 1,NALL
         YY(J) = 2.*V(J)-1.
         DELTA = DELTA+YY(J)
         TI(J) = YY(J)
         TIM(J) = 1.0
  110 CONTINUE
      FX = DELTA*DELTA
      IEVEN = .FALSE.
      DO 130 I = 2,NALL
         IEVEN = .NOT.IEVEN
         DELTA = 0.0
         DO 120 J = 1,NALL
            TIP = 2.*YY(J)*TI(J)-TIM(J)
            DELTA = DELTA+TIP
            TIM(J) = TI(J)
            TI(J) = TIP
  120    CONTINUE
         FI = 0.0
         IF (IEVEN) FI = -1.0/(I*I-1.0)
         DELTA = DELTA/FN-FI
         FX = FX + DELTA*DELTA
  130 CONTINUE
      GO TO 999
C
999   RETURN
      END
