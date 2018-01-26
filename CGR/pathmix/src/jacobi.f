      SUBROUTINE JACOBI(AA,N,IDIM,MODE,EIGVCT,EIGVAL,COND)
C
C          DIAGONALIZES SYMMETRIC MATRIX AA
C          --------------------------------
C
C     FINDS SPECTRAL DECOMPOSITION OF AA BY JACOBI METHOD, AS
C     PRESENTED BY GREENSTADT IN RALSTON AND WILF FIRST BOOK.
C
C INPUT:   AA      -- INPUT MATRIX
C          N       -- ORDER OF MATRIX AA
C          IDIM    -- DIMENSION OF SQUARE MATRICES AA AND EIGVCT
C          MODE    -- FLAG REQUESTING ORDERED EIGENVALUES:
C                      0=NO ORDERING
C                      1=EIGVAL AND EIGVCT ORDERED
C
C OUTPUT:  EIGVCT  -- COLUMN EIGENVECTOR MATRIX
C          EIGVAL  -- EIGENVALUES
C          COND    -- CONDITION NUMBER (MODE=1 ONLY)
C
C NOTE:    ONLY THE PORTION OF AA(I,J), WHERE I.LE.J, IS USED AS INPUT
C
C
C WRITTEN BY JML
C REVISED  3/85 BY SR -- RESTRUCTURED FOR PORTABLE LIBRARY
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION   AA(IDIM,IDIM), EIGVCT(IDIM,IDIM), EIGVAL(*)
      LOGICAL     AGAIN
C
C HANDLE THE TRIVIAL CASE OF N=1
C
      IF (N .LE. 1) THEN
         EIGVAL(1) = AA(1,1)
         EIGVCT(1,1) = 1.0
         COND = 1.0
         GO TO 999
      END IF
C
C INITIALIZE THR, EIGVAL, AND EIGVCT
C
      ANORM = 0.0
      DO 20 I=1,N
         EIGVAL(I) = AA(I,I)
         DO 10 J=1,I
            ANORM = ANORM + 2.*AA(J,I)*AA(J,I)
            AA(I,J) = AA(J,I)
10       CONTINUE
20    CONTINUE
      ANORM = SQRT(ANORM)
      THR = ANORM
      ALMMAX = THR
C
      DO 40 I=1,N
         DO 30 J=1,N
            EIGVCT(J,I) = 0.0
30       CONTINUE
         EIGVCT(I,I) = 1.0
40    CONTINUE
C
C FIND NUMBER OF ITERATIONS
C
      PRECIS = 8.0
      NITER = PRECIS / LOG10( FLOAT(N) ) + 1.5
C
      DO 150 ITER=1,NITER
C
C REPEAT THE FOLLOWING UNTIL ALL OFF-DIAGONALS ARE LESS THAN THR
C
         THR = THR / N
         IF (ALMMAX .LE. THR) GO TO 150
C
100      AGAIN = .FALSE.
         ALMMAX = 0.0
C
         DO 130 L=1,N
            DO 120 M=L+1,N
C
               ALM = AA(L,M)
               TEST = ABS(ALM)
               IF (TEST .GT. ALMMAX) ALMMAX = TEST
               IF (TEST .LE. THR) GO TO 120
               AGAIN = .TRUE.
C
               ALL = EIGVAL(L)
               AMM = EIGVAL(M)
               TEMP = .5 * (ALL-AMM)
               TSQ = ALM*ALM + TEMP*TEMP
               IF (ALL .LT. AMM) THEN
                  Y   = ALM / SQRT(TSQ)
               ELSE
                  Y   = -ALM / SQRT(TSQ)
               END IF
               YSQ = ALM*ALM / TSQ
               RSQ = 2. + 2.*SQRT(1.-YSQ)
               R   = SQRT(RSQ)
               XSQ = RSQ - YSQ
               X   = SQRT(XSQ)
C
               SIN2 = YSQ / RSQ
               SIN  = Y / R
               COS2 = XSQ / RSQ
               COS  = X / R
               SINCOS = X*Y / RSQ
C
               DO 110 I=1,N
                  IF (I.LT.L) THEN
                     TEMP      = AA(I,L)*COS - AA(I,M)*SIN
                     AA(I,M)   = AA(I,L)*SIN + AA(I,M)*COS
                     AA(I,L)   = TEMP
                  ELSE IF (I.GT.L .AND. I.LT.M) THEN
                     TEMP      = AA(L,I)*COS - AA(I,M)*SIN
                     AA(I,M)   = AA(L,I)*SIN + AA(I,M)*COS
                     AA(L,I)   = TEMP
                  ELSE IF (I.GT.M) THEN
                     TEMP      = AA(L,I)*COS - AA(M,I)*SIN
                     AA(M,I)   = AA(L,I)*SIN + AA(M,I)*COS
                     AA(L,I)   = TEMP
                  END IF
C
                  TEMP         = EIGVCT(I,L)*COS - EIGVCT(I,M)*SIN
                  EIGVCT(I,M)  = EIGVCT(I,L)*SIN + EIGVCT(I,M)*COS
                  EIGVCT(I,L)  = TEMP
110            CONTINUE
C
               EIGVAL(L) = ALL*COS2 + AMM*SIN2 - 2.*ALM*SINCOS
               EIGVAL(M) = ALL*SIN2 + AMM*COS2 + 2.*ALM*SINCOS
               AA(L,M)   = (ALL-AMM)*SINCOS + ALM*(COS2-SIN2)
C
120         CONTINUE
130      CONTINUE
C
         IF (AGAIN) GO TO 100
C
150   CONTINUE
C
C RETURN INPUT MATRIX AA INTACT
C
      DO 220 I=1,N
         DO 210 J=1,I-1
            AA(J,I) = AA(I,J)
210      CONTINUE
220   CONTINUE
C
C ---------------------------------------------------------------------
C
C ORDERING IS PERFORMED ONLY IF MODE=1 WAS SPECIFIED
C
      IF (MODE .NE. 1) GO TO 999
C
C ORDER EIGENVALUES IN DECREASING MAGNITUDE
C
      DO 530 I=1,N
C
C FIND I-TH GREATEST EIGENVALUE:
C
         JMAX = I
         DO 510 J=I+1,N
            IF (ABS(EIGVAL(J)) .GT. ABS(EIGVAL(JMAX))) JMAX = J
510      CONTINUE
C
C EXCHANGE EIGENVALUES:
C
         IF (I .EQ. JMAX) GO TO 530
C
         TEMP = EIGVAL(I)
         EIGVAL(I) = EIGVAL(JMAX)
         EIGVAL(JMAX) = TEMP
C
C EXCHANGE CORRESPONDING COLUMN EIGENVECTORS:
C
         DO 520 J=1,N
            TEMP = EIGVCT(J,I)
            EIGVCT(J,I) = EIGVCT(J,JMAX)
            EIGVCT(J,JMAX) = TEMP
520      CONTINUE
530   CONTINUE
C
C GET CONDITION NUMBER
C
      COND = EIGVAL(N)
      IF (COND .NE. 0.0) COND = ABS( EIGVAL(1) / COND )
CCCC
CCCC FOLLOWING CODE IS FOR TESTING ONLY
CCCC
CCC      DO 620 I=1,N
CCC         DEBUG = 0.0
CCC         DO 610 J=1,N
CCC            DEBUG = DEBUG + EIGVCT(J,I)*EIGVCT(J,I)
CCC610      CONTINUE
CCC         WRITE (*,*) I, ') DEBUG:', DEBUG, 'SHOULD BE 1.0'
CCC620   CONTINUE
C
999   RETURN
      END
