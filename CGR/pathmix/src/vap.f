      SUBROUTINE VAP(NALL,XALL,ITP,COV,NCOV,NC,CC,LFNTER,H)
C
C VARIANCE COMPONENTS SUBROUTINE FOR PATHMIX, PARTS 1B, 2B, AND 3B
C ----------------------------------------------------------------
C
C PRINTS VARIANCE COMPONENTS WITH STANDARD ERRORS TO TERSE FILE
C
C ESSENTIAL VARIABLES:
C
C VV()     = VARIANCE COMPONENTS
C UVEC()   = DERIVATIVES OF VARIANCE COMPONENTS FOR ALL PARAMETERS
C SV()     = STANDARD ERRORS
C COV()    = COVARIANCE MATRIX (FROM GEMINI)
C H        = DIFFERENTIATION INTERVAL
C
      PARAMETER (MXV=8, MXP=14)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION   XALL(NALL),ITP(NALL),COV(NCOV,NCOV),CC(NC)
      DIMENSION   VV(MXV),SV(MXV),UVEC(MXP)

      IF (NALL .NE. MXP) THEN
         call intpr('VAP: NALL INCONSISTENT',-1,0,0)
         STOP
      END IF
C
C WRITE HEADER
C
      CALL EJECT(LFNTER)
      WRITE (LFNTER,1000)
 1000 FORMAT (/33X,'VARIANCE COMPONENTS WITH STANDARD ERRORS:')
C
C EVALUATE VAPFUN FUNCTIONS FOR VARIANCE COMPONENTS, INCLUDING
C         COVARIANCES AND RESIDUALS
C
      DO 100 I=1,MXV
         VV(I) = VAPFUN(I,XALL,NALL)
100   CONTINUE
C
C ------- GET STANDARD ERRORS FOR VARIANCE COMPONENTS -------
C
      DO 500 I=1,MXV
C
C GET UVEC() BY FORWARD DIFFERENCES
C
         N = 0
         DO 200 J=1,NALL
            UVEC(J) = 0.0
            IF (ITP(J) .EQ. 1) THEN
               N = N + 1

               X       = XALL(J)
               XALL(J) = X + H
               U = VAPFUN(I,XALL,NALL)
               XALL(J) = X

               UVEC(N) = (U - VV(I)) / H
            END IF
200      CONTINUE
C
C GET THE VARIANCE: UVEC' SIGMA-INV UVEC
C
         V = 0.0
         DO 320 J=1,N
            V = V + UVEC(J)*COV(J,J)*UVEC(J)
            DO 310 K=1,J-1
               V = V + 2.*UVEC(K)*COV(K,J)*UVEC(J)
310         CONTINUE
320      CONTINUE
C
C NOW WE HAVE THE STANDARD ERROR
C
         IF (V .GT. 0.0) THEN
            SV(I) = SQRT(V)
         ELSE
            SV(I) = 0.0
         END IF
500   CONTINUE
C
C ----------- END OF STANDARD ERROR CALCULATIONS ------------
C
C GET VALUE OF PARAMETER A
C
      A = VAPFUN(MXV+1,XALL,NALL)
C
C CONVERT CONSTRAINTS BACK TO VARIANCES, ETC.
C
      DO 600 I=1,9
         CC(I) = 1.0 - CC(I)
600   CONTINUE
C
C WRITE THE RESULTS
C
900   WRITE (LFNTER,1910)
      WRITE (LFNTER,1920) (VV(I),SV(I),VV(I+4),SV(I+4),I=1,4)
      WRITE (LFNTER,1930)
      WRITE (LFNTER,1940) A, (CC(I),I=1,NC)

 1910 FORMAT (/33X,46('*'),
     +        /33X,'SOURCE      ',15X,'PHENOTYPE',10X,
     +        /33X,'            ',6X,'CHILDREN',11X,'ADULTS   ',
     +        /33X,46('*'))
 1920 FORMAT (/33X,'GENETIC     ',F8.3,'(',F5.3,')',F12.3,'(',F5.3,')',
     +        /33X,'CULTURAL    ',F8.3,'(',F5.3,')',F12.3,'(',F5.3,')',
     +        /33X,'COVARIANCE  ',F8.3,'(',F5.3,')',F12.3,'(',F5.3,')',
     +        /33X,'RESIDUAL    ',F8.3,'(',F5.3,')',F12.3,'(',F5.3,')')
 1930 FORMAT (/33X,46('*'))
 1940 FORMAT(//33X,'GENOTYPE ENVIRONMENT CORRELATION =',F10.4,
     +       //33X,'VALUES OF NON-LINEAR CONSTRAINTS:',
     +       //33X,'VARIANCE OF CHILD''S PHENOTYPE    =',F10.4,
     +        /33X,'VARIANCE OF ADULT''S PHENOTYPE    =',F10.4,
     +        /33X,'VARIANCE OF CHILD''S ENVIRONMENT  =',F10.4,
     +        /33X,'VARIANCE OF ADULT''S ENVIRONMENT  =',F10.4,
     +        /33X,'VARIANCE OF CHILD''S INDEX        =',F10.4,
     +        /33X,'VARIANCE OF ADULT''S INDEX        =',F10.4,
     +        /33X,'CORRELATION BETWEEN MARITAL C''S  =',F10.4,
     +        /33X,'CORRELATION BETWEEN MARITAL G''S  =',F10.4,
     +        /33X,'MARITAL C  MARITAL G CORRELATION =',F10.4,
     +        /33X,'U * M                            =',F10.4)
      END

      FUNCTION VAPFUN(NUM,XALL,NALL)
C---
C--- FUNCTION TO RETURN VALUE OF FUNCTION 'NUM'
C---
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER     NUM,NALL
      DIMENSION   XALL(NALL)
C
C GET THE PARAMETERS:
C
      H = XALL(1)
      C = XALL(2)
      Y = XALL(3)
      Z = XALL(4)
      M = XALL(5)
      U = XALL(6)
      P = XALL(7)
      FF= XALL(8)
      FM= XALL(9)
      B = XALL(10)
      X = 1.0
      I = XALL(11)
      V = XALL(12)
      J = XALL(13)
      W = XALL(14)
      T = 0.0
      RG= 1.0
      RC= 1.0
      R = 0.0
      IF (NALL .NE. 14) THEN
         call intpr('VAPFUN: NALL INCONSISTENT',-1,0,0)
         STOP
      END IF
C
C DEFINE SPECIAL TERMS
C
      HZ  = H*Z
      CY  = C*Y
      S = SQRT(M*U)
C
C DEFINE PARAMETER A
C
      A = 0.0
      IF (2.0-(FF+FM) .NE. 0.0) THEN
         ALA = (FF+FM) / (2.0-(FF+FM))
      ELSE
         ALA = 0.0
      END IF
      IF (P .EQ. 0.0) THEN
         A = S * ALA
      ELSE
         PHZCY = P*HZ*CY
         IF (ALA*PHZCY .NE. 0.0) THEN
            AAA = 1.0
            BBB = (ALA*P*(HZ*HZ + CY*CY) - 1.0) / (ALA*PHZCY)
            CCC = 1.0 + S/PHZCY
            DISC = BBB*BBB - 4.0*AAA*CCC
            IF (DISC.LE.0.0) THEN
               A = (-BBB) / (2.0*AAA)
            ELSE
               A = (-BBB - SQRT(DISC)) / (2.0*AAA)
            END IF
         END IF
      END IF
C
C CHILD'S PHENOTYPE VARIANCE COMPONENTS:
C
      IF (NUM .EQ. 1) THEN
C       GENETIC
         VAPFUN = H ** 2

      ELSE IF (NUM .EQ. 2) THEN
C       CULTURAL
         VAPFUN = C ** 2

      ELSE IF (NUM .EQ. 3) THEN
C       COVARIANCE
         VAPFUN = 2.0*H*C*A*RG*RC

      ELSE IF (NUM .EQ. 4) THEN
C       RESIDUAL
         VAPFUN = 1.0 - (H**2 + C**2 + 2.0*H*C*A*RG*RC)
C
C ADULT'S PHENOTYPE:
C
      ELSE IF (NUM .EQ. 5) THEN
C       GENETIC
         VAPFUN = HZ**2

      ELSE IF (NUM .EQ. 6) THEN
C       CULTURAL
         VAPFUN = CY**2

      ELSE IF (NUM .EQ. 7) THEN
C       COVARIANCE
         VAPFUN = 2.0*HZ*CY*A

      ELSE IF (NUM .EQ. 8) THEN
C       RESIDUAL
         VAPFUN = 1.0 - (HZ**2 + CY**2 + 2.0*HZ*CY*A)
C
C GENOTYPE ENVIRONMENT CORRELATION
C
      ELSE IF (NUM .EQ. 9) THEN
         VAPFUN = A
C
C OTHERWISE
C
      ELSE
         call intpr('VAPFUN: NUM INVALID',-1,0,0)
         STOP
      END IF
      END
