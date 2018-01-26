C     PATH 3 - PART B                                         FORTRAN 77
C     WRITTEN BY SKIP RUSSELL                                APRIL, 1982
C
C     REVISION HISTORY:
C     11/82 -- BYPASS DUPLICATE SIGMA INVERSIONS WHEN ADJACENT FAMILES
C              HAVE THE SAME STRUCTURE (TO REDUCE EXECUTION TIME)
C     11/82 -- ADDED TO PATHMIX PACKAGE (AS ENTRY POINT 5)
C     11/82 -- ADDED CONSTRAINT 10 (M*U FOR PARAMETER S)
C      8/83 -- CHANGED PROCESSED DATA FILE FROM BINARY TO SYMBOLIC FORM
C      4/84 -- ADDED ALGEBRAIC "TX" TRANSFORMATIONS TO PATHIN PROGRAM
C     10/84 -- ADDED GOODNESS OF FIT CS, BY FAMILY, TO PROLIX OUTPUT
C      6/85 -- CONVERTED FROM HARRIS FORTRAN TO ANSI FORTRAN 77;
C              ALSO CREATED FIRST VERSION TO RUN ON IBM PC
C     12/85 -- REMOVED FORMS CONTROL FROM COL 1 OF OUTPUT FOR
C              PORTABLILTY; ADDED HEADER LINE TO TOP OF EACH PAGE
C     11/87 -- CODE TO GET AROUND 256 CHARACTER LIMIT IN PA/IT CARDS
C      6/88 -- MERGED CODE FROM PATHIN TO CREATE A STAND-ALONE PROGRAM
C      1/89 -- FIXED BUG IN SUBROUTINE PXFAMR
C      2/89 -- CHANGED IHX to IHX=0 IN SUBROUTINE INPUT (R.W.)
C
C     (IF FUTURE REVISIONS ARE MADE, PLEASE CHANGE VERSION DATE BELOW IN
C     CALL TO SUBROUTINE HEADER)
C
C    *******************************************************************
C    *                                                                 *
C    *                           PATHMIX                               *
C    *                           =======                               *
C    *            THIS IS METHOD 3, PART B, OF PATHMIX                 *
C    *                                                                 *
C    *   PATH3B:  FITTING OF PATH MODELS DIRECTLY TO THE FAMILY        *
C    *   ------   DATA IS PROGRAMMED IN THIS PART.  NO CORRELATIONS    *
C    *            ARE CALCULATED EXPLICITLY.  THIS PROGRAM INVOLVES    *
C    *            ONLY NUCLEAR FAMILIES:  THE PARAMETERS 'RC' AND      *
C    *            'R' ARE NOT INCORPORATED                             *
C    *                                                                 *
C    *   INPUT:   TWO FILES ARE REQUIRED:                              *
C    *   ------   1) A "DATA FILE" WITH QUANTITATIVE MEASURES OF       *
C    *            PHENOTYPES AND ENVIRONMENTAL INDICES IN NUCLEAR      *
C    *            FAMILIES.  EACH DATA RECORD (ONE RECORD PER          *
C    *            INDIVIDUAL) CONSISTS OF A FAMILY ID CODE, A POSITION *
C    *            CODE, THE TRAIT TO WHICH PATH MODELS ARE TO BE       *
C    *            FITTED, AND THE ASSOCIATED INDEX VALUE.  RECORDS     *
C    *            WHICH ARE BLANK IN ANY OF THE ABOVE FIELDS EXCEPT    *
C    *            INDEX ARE REJECTED.                                  *
C    *            2) A "JOB FILE" WITH INFORMATION ON DATA FORMATS,    *
C    *            OPTIONAL DATA TRANSFORMATIONS, AND THE PARAMETERS    *
C    *            OF THE MODEL TO ESTIMATE.  SEE THE ENCLOSED PATHMIX  *
C    *            DOCUMENTATION FOR DETAILS.                           *
C    *                                                                 *
C    *                                                                 *
C    *            DIVISION OF BIOSTATISTICS                            *
C    *            WASHINGTON UNIVERSITY SCHOOL OF MEDICINE             *
C    *            BOX 8067, 660 S. EUCLID AVE                          *
C    *            ST. LOUIS, MISSOURI 63110                            *
C    *                                                                 *
C    *                                                                 *
C    *******************************************************************

      SUBROUTINE PATH3B(FILNAM1, FILNAM2, FILNAM3, FILNAM4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (NALL=26, NC=10)
      CHARACTER   PAR(NALL)*3, MSG*62
      character   filnam1*60, filnam2*60, filnam3*60, filnam4*60
      DIMENSION   XINIT(NALL), XALL(NALL), ITP(NALL), ITPSAV(NALL)
      DIMENSION   SE(NALL), USCORE(NALL)
      DIMENSION   BND(NALL,2), TL(NALL,NALL), WORK(NALL,16)
      DIMENSION   CC(NC),IAC(NC),WRKNLC(NC,4)
      EQUIVALENCE (USCORE,WORK(1,2))

C INITIALIZE PARAMETERS, DEFAULT VALUES AND BOUNDS

      DATA PAR /   'H', 'C', 'Y', 'Z', 'M', 'U', 'P', 'FF', 'FM',
     +             'B', 'I', 'V', 'J', 'W',
     +           'UPF', 'UIF', 'UPM', 'UIM', 'UPC', 'UIC',
     +           'VPF', 'VIF', 'VPM', 'VIM', 'VPC', 'VIC'/

      DATA XINIT /  0.,  0.,  1.,  1.,  0.,  0.,  0.,   0.,   0.,
     +              0.,  0.,  1.,  0.,  1.,
     +              0.,    0.,    0.,    0.,    0.,    0.,
     +              1.,    1.,    1.,    1.,    1.,    1./

      DATA BND /    0.,  0.,  0.,  0., -1., -1., -1.,   0.,   0.,
     +              0.,  0.,  0.,  0.,  0.,
     +          -1000.,-1000.,-1000.,-1000.,-1000.,-1000.,
     +              0.,    0.,    0.,    0.,    0.,    0.,

     +              1.,  1.,  5.,  5.,  1.,  1.,  1.,   1.,   1.,
     +              1.,  1.,  5.,  1.,  5.,
     +           1000., 1000., 1000., 1000., 1000., 1000.,
     +           1000., 1000., 1000., 1000., 1000., 1000./

C DEFINE PROGRAM NAME AND VERSION NUMBER

      CALL HEADER('PATH-3B','VERSION 2/89')

C INITIALIZE FILES:

C     LFNJOB -- JOB FILE, CONTAINS CONTROL CARDS (INPUT)
C     LFNDAT -- CONTAINS RAW FAMILY DATA (INPUT)
C     LFNTER -- SUMMARY OUTPUT ("TERSE") FILE
C     LFNPLX -- DIAGNOSTIC OUTPUT ("PROLIX") FILE

      LFNJOB = 50
      LFNDAT = 51
      LFNTER = 52
      LFNPLX = 53

      CALL FILES( LFNJOB, LFNDAT, LFNTER, LFNPLX,
     +            FILNAM1,FILNAM2,FILNAM3,FILNAM4)

C READ AND TRANSFORM INPUT DATA

      CALL PATHIN( LFNJOB, LFNDAT, LFNTER, LFNPLX, ISTAT )
      IF ( ISTAT .NE. 0 ) then
         call intpr('PATH-3B: ERROR IN JOB FILE',-1,istat,1)
         STOP
      end if

C GENERAL INITIALIZATION AND STARTING ESTIMATES FOR MEANS AND VARIANCES

      CALL DEFALT3b(NALL,XINIT,LFNTER,LFNPLX)

C POSITION JOB FILE IN ORDER TO READ 'PA' AND 'IT' CARDS (FILE 3)

      REWIND LFNJOB
      CALL ADVFIL(LFNJOB,2)

C---- LOOP THROUGH ALL PA(),IT() CARDS --------------------------------

100   CONTINUE

C        READ PA & IT CARDS FROM JOB FILE TO OVERRIDE DEFAULTS
         CALL INPUT3b(NALL,PAR,XINIT,XALL,ITP,LFNJOB,LFNTER,
     +      H,TOL,TRUPB,INLC,IVAR,IHESS,IHX,IQOB,IVERB,ISTAT)
         IF (ISTAT .NE. 0) GO TO 900

C        GET THE MAXIMUM LIKELIHOOD SOLUTION
         call intpr('ALMINIZING...',-1,0,0)
         CALL ALMIZE(FX,NALL,PAR,XALL,ITP,ITPSAV,SE,USCORE,BND,
     +      IBND,IHESS,IVAR,IHX,IQOB,IVERB,LFNTER,LFNPLX,H,TOL,TRUPB,
     +      NIT,NFE,PTG,IDG,TL,NALL,NALL,WORK,
     +      INLC,NC,FC,NK,CC,IAC,WRKNLC,N,MSG,4)

C        WRITE FINAL OUTPUT TO TERSE
         call intpr('WRITING OUTPUT...',-1,0,0)
         CALL OUTPUT3b(NALL,PAR,XALL,ITP,SE,USCORE,TL,NALL,NC,CC,
     +      IVAR,IVERB,LFNTER,H,NIT,NFE,PTG,IDG,FX,FC,NK,MSG)

      GO TO 100

C----------------------------------------------------------------------

C END OF JOB FILE, COPY JOB FILE TO END OF TERSE FILE

900   CALL EJECT(LFNTER)
      WRITE (LFNTER,1900)
 1900 FORMAT (' ', 20('> '), 'JOBFILE', 20(' >') /)
      REWIND LFNJOB
      CALL CPYFIL(3,LFNJOB,LFNTER)

c     STOP 'PATH-3B: NORMAL EOJ'
      RETURN
      END
