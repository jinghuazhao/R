      SUBROUTINE MERLIN(NVAR,UU,VV,NN,MX,RR,NRR,IDIM,ISTAN,
     +                  NFAM,NLINK,XOBS,MAXVAR,MAXFAM)
C
C GENERALIZED ROUTINE TO RETURN MEANS, VARIANCES, AND CORRELATION
C ESTIMATES FOR FAMILY DATA, FOR USE AS STARTING VALUES IN
C MAXIMUM LIKELIHOOD CALCULATIONS.
C
C INPUT ARGUMENTS:
C     NVAR   -- NUMBER OF VARIABLES WITHIN EACH FAMILY GROUP
C     IDIM   -- DIMENSION OF SQUARE MATRICES RR(,) AND NRR(,)
C     ISTAN  -- IF 1, INDICATES DATA IN XOBS() IS TO BE STANDARDIZED
C     NFAM   -- NUMBER OF FAMILIES
C     NLINK()-- NUMBER OF OBS. IN EACH FAMILY BY VARIABLE
C     XOBS() -- THE DATA
C     MAXVAR -- MAXIMUM NUMBER OF VARIABLES AS DIMENSIONED IN NLINK(,)
C     MAXFAM -- MAXIMUM NUMBER OF FAMILIES AS DIMENSIONED IN NLINK(,)
C
C RETURNED:
C     UU()   -- SAMPLE MEANS (OVER ALL PAIRS)
C     VV()   -- SAMPLE VARIANCES (OVER ALL PAIRS)
C     NN()   -- SAMPLE SIZE FOR GIVEN VARIABLE
C     MX()   -- MAXIMUM SIZE FOR GIVEN VARIABLE
C     RR(,)  -- MATRIX OF CORRELATION ESTIMATES;
C     NRR(,) -- SAMPLE SIZES FOR GIVEN CORRELATION IN RR(,)
C
C     RR AND NRR STORED AS (I,J) WHERE:
C     -- DIAGONAL ELEMENTS (I=J) ARE INTRACLASS CORRELATIONS,
C     -- LOWER TRIANGULAR ELEMENTS (I<J) ARE INTER-CLASS CORRELATIONS,
C     -- UPPER TRIANGULAR ELEMENTS (I>J) ARE 'SELF' CORRELATIONS.
C
C  ESSENTIAL VARIABLES USED:
C
C     S()    -- SUM X   (FOR MEAN/VARIANCES)
C     S2()   -- SUM X^2 (FOR MEAN/VARIANCES)
C     TEMP() -- INTERMEDIATE FOR ALL CORRELATIONS
C     DIAG   -- INTRACLASS CORRELATIONS
C     SSW()  -- INTERMEDIATE FOR INTRACLASS CORRELATIONS
C     NF()   -- INTERMEDIATE FOR INTRACLASS CORRELATIONS
C
C
C REVISION HISTORY:
C
C --- REVISION IMPLEMENTED JAN, 1985 BY S.R. ---
C
C THE ESTIMATION OF INTER-CLASS CORRELATIONS INCORPORATES X AND Y
C VARIANCE TERMS, WHICH ARE DERIVED FROM ALL FAMILIES CONTAINING
C EITHER X OR Y DATA VALUES; NOT JUST THOSE FAMILIES WHICH CONTAIN
C DATA FROM BOTH.  THIS OCCASIONALLY RESULTS IN CORRELATION ESTIMATES
C WHICH FALL OUTSIDE OF THE RANGE -1 < R < 1, (I.E. WHEN X, Y AND XY
C INFORMATION ARE DERIVED FROM SUBSTANTIALLY DIFFERENT SUBSETS OF THE
C POPULATION).
C
C THE FIX IS TO IGNORE FAMILIES WHICH, FOR A GIVEN CORRELATION, DO
C NOT HAVE APPROPRIATELY PAIRED DATA VALUES, WHEN WE CALCULATE THE
C CORRESPONDING VARIANCE TERMS.
C
C --- REVISION IMPLEMENTED AUG, 1985 BY S.R. ---
C
C INTER-CLASS (LOWER TRIANGULAR) CORRELATIONS ARE NOW REPLACED BY
C PAIRWISE 'SELF' CORRELATIONS ANY TIME THAT NO MORE THAN ONE PAIR
C OF OBSERVATIONS APPEARS IN ANY GIVEN GROUP.
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CCC      REAL        N
      DIMENSION   RR(IDIM,IDIM),NRR(IDIM,IDIM)
      DIMENSION   UU(NVAR),VV(NVAR),NN(NVAR),MX(NVAR)
      REAL        XOBS
      DIMENSION   NLINK(MAXVAR,*),XOBS(*)
      DIMENSION   DIAG(12),TEMP(12),S(12),S2(12),SSW(12),NF(12)
C
C INITIALIZE CORRELATIONS, SAMPLE SIZES, INTERMEDIATES:
C
      DO 20 I=1,NVAR
         DO 10 J=1,NVAR
           RR(J,I) = 0.0
           NRR(J,I) = 0
10       CONTINUE
         S(I)   = 0.0
         S2(I)  = 0.0
         NN(I)  = 0
         MX(I)  = 0
         SSW(I) = 0.0
         NF(I)  = 0
20    CONTINUE
      SXY = 0.0

CCC---
CCC---  M E A N S   A N D   V A R I A N C E S
CCC---

100   IOBS = 0
      DO 130 IFAM=1,NFAM
         DO 120 I=1,NVAR
            N = NLINK(I,IFAM)
            IF (N.GT.0) THEN
               Z = 0.0
               DO 110 J=1,N
                  IOBS = IOBS + 1
                  X = XOBS(IOBS)
                  Z     = Z     + X
                  S(I)  = S(I)  + X
                  S2(I) = S2(I) + X*X
                  NN(I) = NN(I) + 1
110            CONTINUE
C             UPDATE MAXIMUM COUNT
               IF ( N .GT. MX(I) ) MX(I) = N
C             AND SUM MORE TERMS FOR LATER ON...
               SSW(I) = SSW(I) + Z*Z/N
               NF(I)  = NF(I) + 1
            END IF
120      CONTINUE
130   CONTINUE
C
      DO 150 I=1,NVAR
         IF (NN(I).GT.0) THEN
            UU(I) = S(I) / NN(I)
            VV(I) = (S2(I)-S(I)*S(I)/NN(I)) / NN(I)
         ELSE
            UU(I) = 0.0
            VV(I) = 0.0
         END IF
150   CONTINUE

CCC---
CCC---  I N T R A C L A S S   C O R R E L A T I O N S
CCC---

C
C IMPLEMENTS "C.A.B SMITH'S ITERATIVE METHOD" TO ALLOW FOR VARIABLE
C GROUP SIZES.  THESE WILL BE SAVED IN THE DIAGONAL OF RR.
C
200   DO 260 I=1,NVAR
         IF (NN(I).GT.NF(I)) THEN
            SSW(I) = (S2(I)-SSW(I)) / (NN(I)-NF(I))
         ELSE
            SSW(I) = 0.0
         END IF
         R = .5
         DO 250 IMAXIT=1,100
            PAS = 0.0
            PWK11 = 0.0
            PWK12 = 0.0
            PWK22 = 0.0
            PWSQ1 = 0.0
            PWSQ2 = 0.0
            BIGN  = 0.0
            PASRHO = R
C
            IOBS = 0
            DO 240 IFAM=1,NFAM
               DO 210 J=1,I-1
                  N = NLINK(J,IFAM)
                  IOBS = IOBS + N
210            CONTINUE
               N = NLINK(I,IFAM)
               IF (N .GT. 0.0) THEN
                  IF (N.GT.BIGN) BIGN = N
                  Z = 0.0
                  DO 220 J=1,N
                     IOBS = IOBS + 1
                     X = XOBS(IOBS)
                     Z = Z + X
220               CONTINUE
                  WXK   = Z/N
                  ET1   = 1. / (1.+(N-1.)*R)
                  PAS   = PAS   + ET1*N
                  PWK11 = PWK11 + ET1
                  PWK12 = PWK12 + ET1*ET1*N
                  PWK22 = PWK22 + ET1*ET1*N*N
                  PWSQ1 = PWSQ1 + ET1*WXK*N
                  PWSQ2 = PWSQ2 + ET1*WXK*WXK*N
               END IF
               I1 = I+1
               DO 230 J=I1,NVAR
                  IOBS = IOBS + NLINK(J,IFAM)
230            CONTINUE
240         CONTINUE
C
            R = 0.0
            IF (PAS.NE.0.0) THEN
               TEMP1 = PWSQ2 - PWSQ1*PWSQ1/PAS
               TEMP2 = PWK11 - PWK12/PAS
               TEMP3 = PAS   - PWK22/PAS
               IF (TEMP3.NE.0.0) THEN
                  TEMP4 = SSW(I) + (TEMP1-TEMP2*SSW(I)) / TEMP3
                  IF (TEMP4.NE.0.0) R = (TEMP4 - SSW(I)) / TEMP4
               END IF
            END IF
            IF (ABS(PASRHO-R) .LT. .0001) GO TO 251
250      CONTINUE
251      CONTINUE
C
C       TRAP LARGE NEGATIVE CORRELATION:
         RHOBND = -1.0
         IF (BIGN.GT.1.) RHOBND = -1. / (BIGN - 1.)
         IF (R.LE.RHOBND) R = RHOBND + .05
C
C       NOW WE HAVE ESTIMATED INTRACLASS CORRELATION:
C       (UPDATE RR(,) LATER, SO WE CAN USE IT FOR WORKING STORAGE)
         DIAG(I) = R
260   CONTINUE

CCC---
CCC---  I N T E R C L A S S   C O R R E L A T I O N S
CCC---

C
C IMPLEMENTS UNPUBLISHED METHOD OF D.C. RAO  (FOR MOTIVATION, SEE
C P. 47 OF GERRARD ET AL, AMER. JOURNAL OF HUMAN GENET., 1978).
C THESE WILL BE SAVED IN RR(I,J), I < J.
C
300   IOBS = 0
      DO 350 IFAM=1,NFAM
C       USE TEMP() AS INTERMEDIATE TERM FOR CURRENT FAMILY:
         DO 320 I=1,NVAR
            N = NLINK(I,IFAM)
            IF (N.GT.0) THEN
               Z = 0.0
               DO 310 J=1,N
                  IOBS = IOBS + 1
                  X = XOBS(IOBS)
                  Z = Z + X
310            CONTINUE
               TEMP(I) = Z/N - UU(I)
            ELSE
               TEMP(I) = 0.0
            END IF
320      CONTINUE
C
C       UPDATE INTERMEDIATE SUM RR():
         DO 340 I=1,NVAR
            N = NLINK(I,IFAM)
            IF (N.GT.0) THEN
               DO 330 J=I,NVAR
                  N = NLINK(J,IFAM)
                  IF (N.GT.0) THEN
                     IF (I.EQ.J) THEN
                      IF (N .GT. 1) THEN
                        RR(I,I) = RR(I,I) +
     +                    TEMP(I)*TEMP(I)*N / (1.+(N-1.)*DIAG(I))
                        NRR(I,I) = NRR(I,I) + 1
                      END IF
                     ELSE
                        RR(I,J) = RR(I,J) + TEMP(I)*TEMP(J)
                        NRR(I,J) = NRR(I,J) + 1
                     END IF
                  END IF
330            CONTINUE
            END IF
340      CONTINUE
350   CONTINUE
C
C GET FINAL RR() OFF-DIAGONALS:
C
      DO 370 I=1,NVAR
         IF (NRR(I,I).GT.0) THEN
            T1 = SQRT( RR(I,I) / NRR(I,I) )
         ELSE
            IF (NN(I) .GT. 0) RR(I,I) = 1.0
            T1 = 1.0
         END IF
         I1 = I+1
         DO 360 J=I1,NVAR
            IF (NRR(J,J).GT.0) THEN
               T2 = T1 * SQRT( RR(J,J) / NRR(J,J) )
            ELSE
               T2 = 1.0
            END IF
            IF (T2.NE.0.0 .AND. NRR(I,J).GT.0) THEN
               RR(I,J) = RR(I,J) / NRR(I,J) / T2
            ELSE
               RR(I,J) = 0.0
            END IF
360      CONTINUE
370   CONTINUE

CCC---
CCC---  S E L F   I N T E R C L A S S   C O R R E L A T I O N S
CCC---

C
C PRODUCT-MOMENT CORRELATIONS OF ALL 'MATCHING' PAIRS, (PAIRS FOR
C WHICH ORDINALS ARE THE SAME, IN GIVEN GROUP).  ALL UNMATCHED
C PAIRS ARE IGNORED.  ALLOWS FOR 'PHENOTYPE - SELF INDEX' TYPE CORRS.
C THESE WILL BE SAVED IN RR(I,J), I > J.
C
400   IOBS = 0
      DO 440 IFAM=1,NFAM
         DO 430 I=1,NVAR
            NX = NLINK(I,IFAM)
C          GET SUM(XY) IGNORING GROUPING:
            JOBS = IOBS + NX
            I1 = I+1
            DO 420 J=I1,NVAR
               NY = NLINK(J,IFAM)
               DO 410 K=1,MIN0(NX,NY)
                  X = XOBS(IOBS+K) - UU(I)
                  Y = XOBS(JOBS+K) - UU(J)
                  RR(J,I) = RR(J,I) + X*Y
                  NRR(J,I) = NRR(J,I) + 1
410            CONTINUE
               JOBS = JOBS + NY
420         CONTINUE
            IOBS = IOBS + NX
430      CONTINUE
440   CONTINUE
C
C FINAL RR(,)
C
      DO 460 I=1,NVAR
         I1 = I+1
         DO 450 J=I1,NVAR
            T1 = SQRT( VV(I) * VV(J) )
            IF (T1.NE.0.0 .AND. NRR(J,I).NE.0) THEN
               RR(J,I) = RR(J,I) / NRR(J,I) / T1
C           -- FOLLOWING LINE ADDED 8/85 BY SR
               IF (MX(J) .LE. 1 .AND. MX(I) .LE. 1) RR(I,J) = RR(J,I)
            ELSE
               RR(J,I) = 0.0
            END IF
450      CONTINUE
460   CONTINUE

CCC---
CCC---  S T A N D A R D I Z E   D A T A
CCC---

700   IF (ISTAN .EQ. 1) THEN
         IOBS = 0
         DO 730 IFAM=1,NFAM
            DO 720 I=1,NVAR
               SD = 1.
               IF (VV(I).GT.0.0) SD = SQRT(VV(I))
               N = NLINK(I,IFAM)
               DO 710 J=1,N
                  IOBS = IOBS + 1
                  XOBS(IOBS) = (XOBS(IOBS) - UU(I)) / SD
710            CONTINUE
720         CONTINUE
730      CONTINUE
      END IF
C
C SAVE INTRACLASS CORRELATIONS IN THE DIAGONAL OF RR.
C
      DO 610 I=1,NVAR
         RR(I,I) = DIAG(I)
610   CONTINUE
C
      RETURN
      END
