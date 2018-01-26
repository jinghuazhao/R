      SUBROUTINE FUN3a(FX,XALL,NALL,CC,NC)
C---
C--- DEFINES FUNCTION TO BE MINIMIZED BY GEMINI -- FUNCTION VALUE IS FX
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION   XALL(NALL), CC(NC)

#include "path3.h"

      PARAMETER   (MXFMSZ=30)
      DIMENSION   NN(NGROUP),UUU(NGROUP),VAR(NGROUP)
      DIMENSION   DVEC(MXFMSZ),SIGMA(MXFMSZ,MXFMSZ)
      DIMENSION   EIGVCT(MXFMSZ,MXFMSZ),EIGVAL(MXFMSZ)
      INTEGER     NOLD(NGROUP)
      LOGICAL     DIFSIG

C READ XALL INTO UUU(), VAR(), AND RHO()
C
C XALL IS MAPPED AS FOLLOWS:
C
C     XALL(1-6):    MEANS FOR PF,IF,PM,IM,PC,IC
C     XALL(7-12):   VARIANCES FOR PF,IF,PM,IM,PC,IC
C     XALL(13-27):  CORRELATIONS FOR PF-IF,PF-PM...PC-IC
C     XALL(28,29):  CORRELATIONS FOR PC-PC, AND IC-IC
C     XALL(30):     SIBLING PHENOTYPE - SELF INDEX CORRELATION
C
      L = 0
      DO 110 I=1,NGROUP
         L = L + 1
         UUU(I) = XALL(L)
110   CONTINUE
      DO 120 I=1,NGROUP
         L = L + 1
         VAR(I) = XALL(L)
120   CONTINUE
      DO 140 I=1,NGROUP
         RHO(I,I) = 1.0
         DO 130 J=I+1,NGROUP
            L = L + 1
            RHO(I,J) = XALL(L)
130      CONTINUE
140   CONTINUE
      L = L + 1
      RHO(5,5) = XALL(L)
      L = L + 1
      RHO(6,6) = XALL(L)
      L = L + 1
      RHO(6,5) = XALL(L)
      IF (L .NE. NALL) THEN 
         call intpr('FUN: INTERNAL ERROR',-1,0,0)
         STOP
      END IF

C      WRITE (LFNTER,1110)
C      WRITE (LFNTER,1150) (UUU(I),I=1,6)
C      WRITE (LFNTER,1120)
C      WRITE (LFNTER,1150) (VAR(I),I=1,6)
C      WRITE (LFNTER,1130)
C      DO 150 I=1,6
C         WRITE (LFNTER,1150) (RHO(J,I),J=1,6)
C150   CONTINUE
C 1110 FORMAT (/' MEANS:')
C 1120 FORMAT (/' VARIANCES:')
C 1130 FORMAT (/' RHO MATRIX:')
C 1150 FORMAT (1X,6F14.6)

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C           L N    L I K E L I H O O D    F U N C T I O N
C
C  DESCRIPTION OF ESSENTIAL VARIABLES IN THE LIKELIHOOD CALCULATION
C
C   NN(6)  = THE NUMBER OF OBS. IN THIS FAMILY FOR PF,IF,PM,IM,PC,IC
C   UUU(6) = THE MEANS FOR PF...IC
C   VAR(6) = THE VARIANCES FOR PF...IC
C   XOBS() = THE DATA FOR PF...IC FOR ALL FAMILIES
C   DIFSIG = FLAG INDICATING A NEW STRUCTURE FOR THE SIGMA MATRIX
C   NOLD() = USED TO CHECK STRUCTURE FOR THE SIGMA MATRIX
C
C INITIALIZE THE LN L FUNCTION
C
      FX = 0.0
      IPD = 0
      IOBS = 0

      DO 210 IROW=1,NGROUP
         NOLD(IROW) = 0
210   CONTINUE
C
C GET INFO FOR NEXT FAMILY  [ NN(), XOBS() ]
C
      DO 490 IFAM=1,NFAM
         IX = 0
         DO 310 I=1,NGROUP
            NN(I) = NLINK(I,IFAM)
310      CONTINUE
C
C FIND OUT IF OUR LAST DEFINITION OF SIGMA STILL HOLDS
C
         DIFSIG = .FALSE.
         DO 320 IROW=1,NGROUP
            IF (NOLD(IROW) .NE. NN(IROW)) THEN
               DIFSIG = .TRUE.
               NOLD(IROW) = NN(IROW)
            END IF
320      CONTINUE
         IF (DIFSIG) THEN
C
C DEFINE SIGMA MATRIX
C
         IX = 0
         DO 390 IROW=1,NGROUP
            IF (NN(IROW) .GT. 0) THEN
             SX = SQRT(VAR(IROW))
             IY0 = IX
             DO 350 I=1,NN(IROW)
                IY = IY0
                IX = IX + 1
                DO 340 ICOL=IROW,NGROUP
                   IF (NN(ICOL) .GT. 0) THEN
                    SY = SQRT(VAR(ICOL))
                       DO 330 J=1,NN(ICOL)
                       IY = IY + 1
                       IF (IX .EQ. IY) THEN
C                         DIAGONAL ELEMENTS:
                          SIGMA(IX,IX) = VAR(IROW)
                       ELSE IF ((IROW.EQ.5 .AND. ICOL.EQ.6)
     +                 .AND. (I.EQ.J)) THEN
C                         SELF PHE.-VS-INDICES:
                          SIGMA(IX,IY) = RHO(ICOL,IROW)*SX*SY
                          SIGMA(IY,IX) = SIGMA(IX,IY)
                       ELSE IF (IX .LT. IY) THEN
C                         OTHER COVARIANCES:
                          SIGMA(IX,IY) = RHO(IROW,ICOL)*SX*SY
                          SIGMA(IY,IX) = SIGMA(IX,IY)
                       END IF
330                  CONTINUE
                   END IF
340             CONTINUE
350          CONTINUE
            END IF
390      CONTINUE
C      HOW ABOUT
C     THAT, EH?
C
         CALL OPTSUM(IX,SIGMA,MXFMSZ,SIGSUM)
C
C GET SIGMA-INVERSE AND LOG DET. SIGMA
C
         CALL PDINV(SIGMA,IX,MXFMSZ,SIGMA,EIGVCT,EIGVAL,DETLOG,ISTAT)

      END IF

      IF (ISTAT .NE. 0) IPD = IPD + 1
C
C DEFINE DVEC()
C
         IX = 0
         DO 420 IROW=1,NGROUP
            IF (NN(IROW).GT.0) THEN
             DO 410 I=1,NN(IROW)
                IX = IX + 1
                IOBS = IOBS + 1
                DVEC(IX) = XOBS(IOBS)
                DVEC(IX) = DVEC(IX) - UUU(IROW)
410          CONTINUE
            END IF
420      CONTINUE
C
C DVEC' SIGMA-INVERSE DVEC:
C
         DSD = 0.0
         DO 440 I=1,IX
            DO 430 J=1,IX
               DSD = DSD + DVEC(I)*SIGMA(I,J)*DVEC(J)
430         CONTINUE
440      CONTINUE
C
C UPDATE LN L
C
         FX = FX + .5*(DETLOG+DSD)

         CALL OPTOUT(LFNGOF,IFAM,IX,DSD,DVEC,SIGSUM)

490   CONTINUE
      END
