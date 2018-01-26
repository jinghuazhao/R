      SUBROUTINE FUN3b(FX,XALL,NALL,CC,NC)
C---
C--- DEFINES FUNCTION TO BE MINIMIZED BY GEMINI -- FUNCTION VALUE IS FX
C---
      PARAMETER   (MXFMSZ=30)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION   XALL(NALL), CC(NC)

#include "path3.h"

C   DECLARATIONS FOR PART A:
      DIMENSION   MAPEQN(NGROUP,NGROUP)
C   DECLARATIONS FOR PART B:
      DIMENSION   NN(NGROUP),UUU(NGROUP),VAR(NGROUP)
      DIMENSION   DVEC(MXFMSZ),SIGMA(MXFMSZ,MXFMSZ)
      DIMENSION   EIGVCT(MXFMSZ,MXFMSZ),EIGVAL(MXFMSZ)
      INTEGER     NOLD(NGROUP)
      LOGICAL     DIFSIG
      DOUBLE PRECISION  M,II,JJ
      DOUBLE PRECISION  IV,JW
      common /parms/ h,c,y,z,m,u,p,ff,fm,b,x,ii,v,jj,w,t,rg,rc,r,
     +               hz,cy,iv,jw,s,a,ph,ga

      DATA MAPEQN
     +     /    0,  501,  100,  600,  210,  710,
     +          0,    0,  600, 1000,  810, 1110,
     +          0,    0,    0,  501,  220,  720,
     +          0,    0,    0,    0,  820, 1120,
     +          0,    0,    0,    0,  300,  900,
     +          0,    0,    0,    0,  500, 1200/
      data dvec /mxfmsz*0./,
     +     sigma /(mxfmsz*mxfmsz)*0./, nold /ngroup*-99/

C PART A:        D E F I N I T I O N   O F   R H O ()
C
C INITIALIZE EXPECTED CORRELATION CALCULATIONS
C
      CALL RHINIT(XALL,NPARM)
C
C GET EXPECTED CORRELATIONS
C NOTE: RHO(IC,PC) REFERS TO SELF-INDEX, RHO(PC,IC) TO SIB'S
C
      ESP = 0.0

      DO 120 I=1,NGROUP
         rho(i,i) = 1.0
         DO 110 J=1,NGROUP
            rho(i,j) = 0.
            IEQNO = MAPEQN(J,I)
            IF (IEQNO .GT. 0) THEN
               CALL RHOEXP(RHO(I,J),IEQNO,ESP)
            END IF
110      CONTINUE
120   CONTINUE
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C PART B:        L N    L I K E L I H O O D    F U N C T I O N
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
C GET THE MEANS AND VARIANCES:
C
      DO 200 I=1,NGROUP
         J = I + NGROUP
         UUU(I) = XALL(NPARM+I)
         VAR(I) = XALL(NPARM+J)
200   CONTINUE
      IF (NPARM+J .NE. NALL) THEN
         call intpr('FUN: NALL INCONSISTENT',-1,0,0)
         STOP
      END IF
C
C INITIALIZE THE LN L FUNCTION
C
      FX = 0.0
      IPD = 0
      IOBS = 0

      DO 210 IROW=1,6
         NOLD(IROW) = 0
210   CONTINUE
C
C GET INFO FOR NEXT FAMILY  [ NN(), XOBS() ]
C
      DO 490 IFAM=1,NFAM
         IX = 0
         DO 310 I=1,6
            NN(I) = NLINK(I,IFAM)
310      CONTINUE
C
C FIND OUT IF OUR LAST DEFINITION OF SIGMA STILL HOLDS
C
         DIFSIG = .FALSE.
         DO 320 IROW=1,6
            IF (NOLD(IROW).NE.NN(IROW)) THEN
               DIFSIG = .TRUE.
               NOLD(IROW) = NN(IROW)
            END IF
320      CONTINUE
         IF (DIFSIG) THEN
C
C DEFINE SIGMA MATRIX
C
         IX = 0
         DO 390 IROW=1,6
            IF (NN(IROW).GT.0) THEN
             SX = SQRT(VAR(IROW))
             IY0 = IX
             DO 350 I=1,NN(IROW)
                IY = IY0
                IX = IX + 1
                DO 340 ICOL=IROW,6
                   IF (NN(ICOL).GT.0) THEN
                    SY = SQRT(VAR(ICOL))
                       DO 330 J=1,NN(ICOL)
                       IY = IY + 1
                       IF (IX.EQ.IY) THEN
C                         DIAGONAL ELEMENTS:
                          SIGMA(IX,IX) = VAR(IROW)
                       ELSE IF ((IROW.EQ.5 .AND. ICOL.EQ.6)
     +                 .AND. (I.EQ.J)) THEN
C                         SELF PHE.-VS-INDICES:
                          SIGMA(IX,IY) = RHO(ICOL,IROW)*SX*SY
                          SIGMA(IY,IX) = SIGMA(IX,IY)
                       ELSE IF (IX.LT.IY) THEN
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

      IF (ISTAT .NE. 0) IPD = IPD + ISTAT
C
C DEFINE DVEC()
C
         IX = 0
         DO 420 IROW=1,6
            IF (NN(IROW).GT.0) THEN
             DO 410 I=1,NN(IROW)
                IX = IX + 1
                IOBS = IOBS + 1
                DVEC(IX) = XOBS(IOBS) - UUU(IROW)
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

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C DEFINE THE CONSTRAINTS
C
      CALL CONSTR(CC)
      END
