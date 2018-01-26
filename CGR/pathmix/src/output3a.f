      SUBROUTINE OUTPUT3a(NALL,PAR,XALL,ITP,SE,U,SS,COV,NCOV,NC,CC,
     +   IVAR,LFNTER,NIT,NFE,PTG,IDG,FX,FC,NK,MSG)
C---
C--- WRITES FINAL OUTPUT TO TERSE
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER   PAR(*)*(*),MSG*(*)
      DIMENSION   XALL(NALL),ITP(NALL),SE(NALL),U(NALL),SS(NALL)
      DIMENSION   COV(NCOV,NCOV),CC(NC)

#include "path3.h"

      CHARACTER   IPAR(9)*3

      PARAMETER   (MXCORR=18)
      INTEGER     IGROUP(MXCORR,2)

      DATA IPAR   /'PF','IF','PM','IM','PC','IC','PC1','IC1','IC2'/

      DATA IGROUP /1,1,1,1,1,2,2,2,2,3,3,3,4,4,7,5,6,7,
     +             2,3,4,5,6,3,4,5,6,4,5,6,5,6,9,5,6,8/
C
C INDICATE IF THE LIKELIHOOD FUNCTION IS UNDEFINED AT THE FINAL VALUE
C
      IF (IPD .NE. 0) WRITE (LFNTER,1110) IPD
 1110 FORMAT (/'  ** WARNING:',I6,' SIGMA-MATRICES WERE NOT POSITIVE-',
     + 'DEFINITE IN THE LAST FUNCTION EVALUATION **')
C
C WRITE COVARIANCE MATRIX TO TERSE
C
      IF (IVAR .NE. 0) THEN
         IF (IVAR .EQ. 1) THEN
            WRITE (LFNTER,1200)
            WRITE (LFNTER,1210) ((PAR(I), J=1,ITP(I)), I=1,NALL)
            K = 0
            DO 210 I=1,NALL
               IF (ITP(I) .GT. 0) THEN
                  K = K + 1
                  WRITE (LFNTER,1220) PAR(I), (COV(J,K),J=1,K)
               END IF
210         CONTINUE
         ELSE
            WRITE (LFNTER,1250)
         END IF
      END IF

 1200 FORMAT (/' ','FACTORIZATION OF B-MATRIX HAS SUCCEEDED')
 1210 FORMAT (/' ','COVARIANCE MATRIX',
     +        /3X,(18(4X,A)))
 1220 FORMAT (' ',A3,18F7.3, /(4X,18F7.3))
 1250 FORMAT (/' ','FACTORIZATION OF B-MATRIX HAS FAILED')
C
C WRITE CPU TIME
C
      CALL TIMER(CPUMIN)

      WRITE (LFNTER,1300) CPUMIN
 1300 FORMAT (//' ','OPTIMIZATION TOOK:',F12.2,' MIN. CPU TIME')
C
C WRITE TERMINATION MESSAGE
C
      NMSG = LENSTR(MSG)
      WRITE (LFNTER,1400) IDG, MSG(1:NMSG)
 1400 FORMAT (/' ITERATIVE PROCESS TERMINATED BECAUSE: (IDG =',I2,')',
     +        /10X,'*** ',A,' ***')
C
C ESTIMATE SAMPLE SIZE
C
      K = 0
      DO 410 I=1,NALL
         IF (ITP(I) .GT. 0) THEN
            K = K + 1
            VAR = COV(K,K)
            IF (IVAR .EQ. 0 .OR. VAR .EQ. 0.0) THEN
               SS(I) = 0.0
            ELSE
               SS(I) = (1.-XALL(I)**2)**2 / VAR
            END IF
         END IF
410   CONTINUE
C
C WRITE SUMMARY
C
      CALL EJECT(LFNTER)
      WRITE (LFNTER,1500) -PTG, 2.*FX, FC, NK, NIT, NFE

 1500 FORMAT (/' ',29('- '),' SUMMARY OUTPUT ',28(' -'),
     +       //' ','UK-1U =',F16.10,'    -2 LN L =',1PG19.10,
     +             '    FC =',G8.2,
     +             '    NK =',I3,'    NIT =',I4,'    NFE =',I5)
C
C WRITE FINAL PARAMETER ESTIMATES TO TERSE
C
      WRITE (LFNTER,1600)
      WRITE (LFNTER,1605)

C    MEANS:
      K = 0
      DO 510 I=1,NGROUP
         K = K + 1
         IF (ITP(K) .GT. 0) THEN
            WRITE (LFNTER,1610) PAR(K),IPAR(I),XALL(K),U(K),SE(K)
         END IF
510   CONTINUE

C    VARIANCES:
      DO 520 I=1,NGROUP
         K = K + 1
         IF (ITP(K) .GT. 0) THEN
            WRITE (LFNTER,1620) PAR(K),IPAR(I),XALL(K),U(K),SE(K)
         END IF
520   CONTINUE

C    CORRELATIONS:
      K = K + 1
      DO 530 I=K,NALL
         IG1 = IGROUP(I-K+1,1)
         IG2 = IGROUP(I-K+1,2)
         IF (ITP(I) .GT. 0) THEN
            WRITE (LFNTER,1630) PAR(I),IPAR(IG1),IPAR(IG2),
     +         XALL(I),U(I),SE(I),SS(I)
         ELSE IF (ITP(I) .EQ. 0) THEN
            WRITE (LFNTER,1631) PAR(I),IPAR(IG1),IPAR(IG2),XALL(I)
         ELSE
            J = -ITP(I)
            WRITE (LFNTER,1632) PAR(I),IPAR(IG1),IPAR(IG2),XALL(I),
     +         PAR(J)
         END IF
530   CONTINUE
      WRITE (LFNTER,1605)

 1600 FORMAT (/' ',101X,'EFFECTIVE',
     +        /' PARAMETER',11X,'GROUP(S) ',
     +          10X,'  ESTIMATE',10X,'   U-SCORE',
     +          10X,'STANDARD ERROR',6X,' SAMPLE SIZE')
 1605 FORMAT (' ',111('-'))
 1610 FORMAT (1X,A,' - MEAN           ',3X,'(',A3,')',3F20.5)
 1620 FORMAT (1X,A,' - VARIANCE       ',3X,'(',A3,')',3F20.5)
 1630 FORMAT (1X,A,' - CORRELATION   (',A3,',',A3,')',3F20.5,F20.1)
 1631 FORMAT (1X,A,' - CORRELATION   (',A3,',',A3,')',F20.5,' (FIXED)')
 1632 FORMAT (1X,A,' - CORRELATION   (',A3,',',A3,')',F20.5,
     +   ' (=',A,')')
C
C WRITE INDIVIDUAL GOODNESS OF FIT TO PROLIX
C
      GOFOUT = .TRUE.
      CALL FUN3a(FX,XALL,NALL,CC,NC)
      GOFOUT = .FALSE.
      END
