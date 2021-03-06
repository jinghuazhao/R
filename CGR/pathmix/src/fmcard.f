      SUBROUTINE FMCARD( STRING, LFNTER, LFNPLX, ISTAT )
C---
C--- READS 'FM' CARD, SAVING RELEVANT INFORMATION IN COMMON
C---
C--- ARGUMENTS:
C---     STRING() -- 'FM' CARD (CHARACTER FORMAT)
C---     LFNTER   -- LOGICAL UNIT FOR LIST OUTPUT
C---     LFNPLX   -- LOGICAL UNIT FOR PROLIX OUTPUT
C---     ISTAT    -- SET TO -1 IF ERROR, 0 OTHERWISE
C---
C--- INPUT:
C---     A STRING OF THE FORM: (FORTRAN FORMAT) VAR1, VAR2,..., VARn
C---
C--- OUTPUT:
C---     THE FOLLOWING VARIABLES IN COMMON BLOCK /FIELD/:
C---        NFLDIN, FLDTYP, FLDNAM, FLDNCH, FMTALP, FMTX
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C NOTE: PARAMETER MAXDP IS THE NUMBER OF BYTES IN A DOUBLE PRECISION
C DATA ELEMENT FOR THE CURRENT MACHINE.  IT IS USED TO DETERMINE
C HOW MUCH CHARACTER DATA CAN BE STORED IN FLDOBS(), AND CAN SAFELY
C BE INCREASED FROM 6 TO 8 ON MANY MACHINES.
      PARAMETER (MAXDP=6)

#include "field.h"

      CHARACTER   STRING*(*)

C ======================================================================

C GET LIST OF FIELD NAMES

C     BYPASS FORMAT SECTION TEMPORARILY
      ICOL = 1
      CALL PAR( STRING, ICOL, ILEFT, IRIGHT, ISTAT )
      IF ( ISTAT .NE. 0 ) GO TO 800

      ICOL = NSPACE( STRING, IRIGHT+1, ISTAT )
      IF ( ISTAT .NE. 0 ) GO TO 190

      NNN = 0

100   CONTINUE

         MODE = 2
         CALL MAKVAR( STRING, ICOL, MODE,
     @      MAXFLD, FLDOBS, FLDBLK, FLDTYP, FLDNAM, FLDNCH,
     @      IVAR, ISTAT )
         IF ( ISTAT .NE. 0 ) GO TO 810

         NNN = NNN + 1
         IF ( NNN .NE. IVAR ) THEN
            WRITE (LFNTER,*) 'ERROR : VAR', NNN,
     @         ' .NE. TABLE ENTRY', IVAR
            call intpr
     +      ('INTERNAL ERROR IN FMCARD: VAR. TABLE NOT INITIALIZED',
     +      -1,0,0)
            STOP
         END IF

C        CHECK FOR END OF LIST

         ICOL = NSPACE( STRING, ICOL, ISTAT )
         IF ( ISTAT .NE. 0 ) GO TO 190

C        CHECK THE DELIMITER

         IF ( STRING(ICOL:ICOL) .LT. 'A' ) GO TO 820

      GO TO 100

190   CONTINUE

C     SAVE A COUNT OF THE NUMBER OF INPUT FIELDS

      NFLDIN = NNN

C ======================================================================

C PARSE FORMAT PORTION OF STRING

      FMTALP = ' '

C     COPY INTO ALPHA FORMAT:
      NCOUNT = IRIGHT - ILEFT + 1
      IOFFS = ILEFT - 1
      FMTALP = STRING(ILEFT:IRIGHT)

C LOOP THROUGH FORMAT FIELDS

      ICOL = 1
      NNN = 0

200   CONTINUE

         ICOL = NSPACE( FMTALP(:NCOUNT), ICOL+1, ISTAT )
         IF ( ISTAT .NE. 0 ) GO TO 830

         NREP = IXNO( FMTALP(:NCOUNT), ICOL, ISTAT )
         IF ( NREP .EQ. 0 ) NREP = 1

         ICOL = NSPACE( FMTALP(:NCOUNT), ICOL, ISTAT )
         IF ( ISTAT .NE. 0 ) GO TO 830

C     HANDLE CURRENT FORMAT FIELD:
C     FIELDS ARE: 'F', 'A', 'I', 'X', 'T', '/'

C----
C---- "F"
C----
         IF ( FMTALP(ICOL:ICOL) .EQ. 'F' ) THEN

C        MAKE SURE A NUMBER FOLLOWS:
         ITEMP = NSPACE( FMTALP(:NCOUNT), ICOL+1, ISTAT )
         IF ( ISTAT .NE. 0 ) GO TO 830
         J = INT( XNO(FMTALP(:NCOUNT), ITEMP, ISTAT ) )
         IF ( J.LE.0 .OR. J.GT.MAXFLD ) GO TO 830

C        INITIALIZE FIELD SPECIFICATIONS:
         DO 310 I=1,NREP
            FLDTYP( NNN+I ) = 2
            FMTX( NNN+I ) = '(F'
310      CONTINUE
         JCOL = 2

C        COPY FORMAT FIELD, REPLACING "F" WITH "A" IN ALPHA FORMAT:
         FMTALP(ICOL:ICOL) = 'A'
320      CONTINUE
            ICOL = ICOL + 1
            JCOL = JCOL + 1
            DO 325 I=1,NREP
               FMTX(NNN+I)(JCOL:JCOL) = FMTALP(ICOL:ICOL)
325         CONTINUE
            IF (ICOL .GE. ITEMP) GO TO 830
         IF ( FMTALP(ICOL:ICOL) .NE. '.' ) GO TO 320

C        AND COPY PAST DECIMAL, BLANKING ALPHA FORMAT:
         DO 340 ICOL=ICOL,ITEMP-1
            DO 330 I=1,NREP
               FMTX(NNN+I)(JCOL:JCOL) = FMTALP(ICOL:ICOL)
330         CONTINUE
            JCOL = JCOL + 1
            FMTALP(ICOL:ICOL) = ' '
340      CONTINUE

C        FINISH UP INPUT FORMAT
         DO 350 I=1,NREP
            FMTX(NNN+I)(JCOL:) = ')'
350      CONTINUE

         NNN = NNN + NREP
         ICOL = ITEMP
C----
C---- "A"
C----
      ELSE IF (FMTALP(ICOL:ICOL) .EQ. 'A') THEN

C        MAKE SURE A NUMBER FOLLOWS:
         ITEMP = NSPACE( FMTALP(:NCOUNT), ICOL+1, ISTAT )
         IF ( ISTAT  .NE. 0) GO TO 830
         J = IXNO( FMTALP(:NCOUNT), ITEMP, ISTAT )
         IF ( J.LE.0 .OR. J.GT.MAXDP ) GO TO 830

C        DEFINE FIELD SPECIFICATIONS:
         DO 410 I=1,NREP
            FLDTYP(NNN+I) = 1
            FMTX(NNN+I) = '(A' // FMTALP(ICOL+1:ITEMP-1) // ')'
410      CONTINUE

         NNN = NNN + NREP
         ICOL = ITEMP
C----
C---- "I"
C----
      ELSE IF ( FMTALP(ICOL:ICOL) .EQ. 'I' ) THEN

C        MAKE SURE A NUMBER FOLLOWS:
         ITEMP = NSPACE( FMTALP(:NCOUNT), ICOL+1, ISTAT )
         IF ( ISTAT  .NE. 0) GO TO 830
         J = IXNO( FMTALP(:NCOUNT), ITEMP, ISTAT )
         IF ( J.LE.0 .OR. J.GT.MAXFLD ) GO TO 830

C        REPLACE "I" WITH "A" IN ALPHA FORMAT:
         FMTALP(ICOL:ICOL) = 'A'

C        DEFINE FIELD SPECIFICATIONS:
         DO 510 I=1,NREP
            FLDTYP(NNN+I) = 2
            FMTX(NNN+I) = '(F' // FMTALP(ICOL+1:ITEMP-1) // '.0)'
510      CONTINUE

         NNN = NNN + NREP
         ICOL = ITEMP
C----
C---- IGNORE "X"
C----
         ELSE IF ( FMTALP(ICOL:ICOL) .EQ. 'X' ) THEN
            ICOL = ICOL + 1
C----
C---- "T"   IGNORE TABS
C----
         ELSE IF ( FMTALP(ICOL:ICOL) .EQ. 'T' ) THEN

C           MAKE SURE A NUMBER FOLLOWS:
            ICOL = NSPACE( FMTALP(:NCOUNT), ICOL+1, ISTAT )
            IF ( ISTAT .NE. 0 ) GO TO 830
            J = IXNO( FMTALP(:NCOUNT), ICOL, ISTAT )
            IF ( J .EQ. 0 ) GO TO 830
C----
C---- IGNORE "/"
C----
         ELSE IF ( FMTALP(ICOL:ICOL) .EQ. '/' ) THEN
            ICOL = ICOL + 1
C----
C---- ANYTHING ELSE IS AN ERROR
C----
         ELSE
            GO TO 830
         END IF

C SKIP TO NEXT DELIMITER

         ICOL = NSPACE( FMTALP(:NCOUNT), ICOL, ISTAT )
         IF ( ISTAT .NE. 0) GO TO 830

C CONTINUE LOOP IF DELIMITER IS A COMMA

      IF ( FMTALP(ICOL:ICOL) .EQ. ',' ) GO TO 200
      IF ( FMTALP(ICOL:ICOL) .NE. ')' ) GO TO 830

C MAKE SURE THE FORMAT IS CONSISTENT WITH THE NUMBER OF VARIABLE NAMES

      IF ( NNN .NE. NFLDIN ) GO TO 840

C ======================================================================

C WRITE FORMAT SUMMARY TO PROLIX

      WRITE (LFNPLX,1700) STRING(ILEFT:IRIGHT), FMTALP(1:NCOUNT)
      DO 700 I=1,NFLDIN
         WRITE (LFNPLX,1710) FLDNAM(I), FMTX(I), FLDTYP(I)
700   CONTINUE

 1700 FORMAT (' ',20X,'INPUT FORMAT:       ',A,
     +       /' ',20X,'BLANK-CHECK FORMAT: ',A,
     +      //' ',20X,'FIELD FORMATS AND TYPES:'
     +      //' ',20X,'VARIABLE NAME     FORMAT',
     +                     '    TYPE (1=CHAR 2=REAL 3=INT)')
 1710 FORMAT (' ',20X,4X,A12,2X,A,I4)

      GO TO 900

C BAD 'FM' CARD

800   WRITE (LFNTER,*) '** ERROR IN "FM" CARD **'
      GO TO 890

810   WRITE (LFNTER,*) '** ILLEGAL VARIABLE NAME **'
      GO TO 890

820   WRITE (LFNTER,*) '** ILLEGAL DELIMITER BETWEEN VARIABLE NAMES **'
      GO TO 890

830   ICOL = IOFFS + ICOL
      WRITE (LFNTER,*) '** ILLEGAL FORMAT SPECIFICATION **'
      GO TO 890

840   WRITE (LFNTER,*) '** VARIABLE LIST DOES NOT AGREE WITH FORMAT **'
      GO TO 890

890   ISTAT = -1
      CALL ERRLOC( LFNTER, STRING, ICOL )
      WRITE (LFNTER,*) 'FORMAT IS:'
      WRITE (LFNTER,*) '  FM (FORTRAN FORMAT) LIST OF VARIABLE NAMES'
      RETURN

900   ISTAT = 0
      END

      BLOCK DATA INITFM
C---
C--- INITIALIZE THE COMMON BLOCK WHICH CONTAINS DATA STRUCTURES
C--- FOR ACCESSING INPUT DATA
C---

#include "field.h"

      DATA  NFLDIN   / 0 /
      DATA  FLDTYP   / 20 * 0 /
      DATA  FLDNAM   / 20 * '?' /
      DATA  FLDNCH   / 20 * 0 /
      DATA  FLDBUF   / 20 * ' ' /
      DATA  FLDOBS   / 20 * 0.0 /
      DATA  FLDBLK   / 20 * 999.0 /
      DATA  FMTALP   / '()' /
      DATA  FMTX     / 20 * '()' /

      END
