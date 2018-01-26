      SUBROUTINE EVAL(MODE,STRING,IQ,IEXNO,
     @   MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @   IVAR,ANS,ISTAT)

C---------------------------------------------------------------------+
C                                                                     !
C                             E  V  A  L                              !
C                             ----------                              !
C                                                                     !
C for evaluating arithmetic expressions                               !
C written 3/84 in Fortran 77 by Skip Russell                          !
C                                                                     !
C revision history:                                                   !
C  6/88  -- moved variable table from common to argument list for     !
C           modularity; implemented random number function ranu();    !
C           implemented "integer" variable type (SR)                  !
C                                                                     !
C                                                                     !
C WHAT IT DOES:                                                       !
C                                                                     !
C     This subroutine parses and evaluates arithmetic expressions.    !
C     Several forms of expressions are currently recognized:          !
C                                                                     !
C     form                          result of evaluation              !
C     ----                          --------------------              !
C                                                                     !
C     <expr>                        the arithmetic expression <expr>  !
C                                   is evaluated and the result       !
C                                   returned.                         !
C                                                                     !
C     <dest> = <expr>               the destination variable <dest>   !
C                                   is assigned the value of the      !
C                                   evaluated expression.             !
C                                                                     !
C     IF (<cond>) <dest> = <expr>   the destination variable is       !
C                                   assigned a value only if the      !
C                                   conditional expression <cond>     !
C                                   evaluates to true.                !
C                                                                     !
C     IF (<cond>) REJECT            the rejection status code is      !
C                                   returned if the conditional       !
C                                   expression evaluates to true.     !
C                                                                     !
C                                                                     !
C HOW TO USE IT:                                                      !
C                                                                     !
C INITIALIZATION:                                                     !
C     define appropriatly the values of MAXTOK, MAXCON, and MAXOPT    !
C     in common block /TX/, and call this subroutine with MODE=0 to   !
C     initialize.  Alternatively, one may use a BLOCK DATA            !
C     subprogram to initialize the tokenized expression space.        !
C                                                                     !
C INPUT ARGUMENTS:                                                    !
C     MODE     -- specifies the function to perform, as follows:      !
C        0 = clear: erases variables & deallocates saved expressions  !
C        1 = read and evaluate expression (tokenized form not saved)  !
C        2 = read expression only, saving tokenized form              !
C        3 = evaluate expression previously saved (in tokenized form) !
C     STRING   -- the expression to read (if MODE=1 or 2)             !
C     IQ       -- logical unit to write error messages to             !
C     IEXNO    -- indicates which expression to read or evaluate, as  !
C        follows:                                                     !
C                                                                     !
C             if MODE=0 or 1, the value of IEXNO is ignored.          !
C                                                                     !
C             if MODE=2, a value of 0 or 1 indicates that token space !
C        is to be cleared, and the given expression saved as the      !
C        first entry.  A larger value for IEXNO, say 3 for example,   !
C        indicates that the first 2 expression are to be bypassed,    !
C        and the given expression saved as the 3rd entry.  In either  !
C        case, the new entry becomes the last in the table.  If a     !
C        value is given which is more than the current number of      !
C        entries, the given entry will be placed at the end of the    !
C        list (see also below).                                       !
C                                                                     !
C             if MODE=3, the IEXNO-th expression is evaluated, and    !
C        the result returned.  If the value of IEXNO is more than     !
C        the number of saved expressions, an 'end of expression'      !
C        code is returned in ISTAT (see below).                       !
C                                                                     !
C                                                                     !
C        The remainder of the input arguments comprise the variable   !
C        table, containing variable names and their values.           !
C                                                                     !
C     MAXVAR   -- the dimension of the arrays X, BLANK, ITYPE,        !
C                 NAMVAR and LENVAR (defines the limit on the         !
C                 number of variables in the variable table)          !
C     X()      -- the current value of each variable                  !
C     BLANK()  -- the missing value code to use for each variable     !
C     ITYPE()  -- the type codes for each variable as follows:        !
C                 0 = variable not yet defined                        !
C                 1 = alphanumeric, represented as hollerith within   !
C                     double precision scalers (cannot be used in     !
C                     numeric operations)                             !
C                 2 = floating point numeric (double precision)       !
C                 3 = integer, represented in double precision        !
C     NAMVAR   -- name assigned to each variable                      !
C     LENVAR   -- lengths of names above                              !
C                                                                     !
C                                                                     !
C RETURNED ARGUMENTS:                                                 !
C     IEXNO    -- entry number of saved expression (if MODE=2)        !
C     IVAR     -- index number of destination variable (if > 0)       !
C     ANS      -- result of expression evaluation                     !
C     ISTAT    -- error status as follows:                            !
C        definition time errors (MODE 1 or 2):                        !
C         0 = operation performed as requested                        !
C        -1 = syntax error in input string                            !
C                                                                     !
C        evaluation time errors (MODE 1 or 3)                         !
C         0 = operation performed as requested                        !
C        -2 = evaluation time error (e.g. division by zero)           !
C        -3 = evaluation time error in condition part of expression   !
C         1 = condition false                                         !
C         2 = blank variable encountered (yields missing value code)  !
C         3 = rejection                                               !
C         4 = end of expression list (IEXNO > # saved expressions)    !
C                                                                     !
C                                                                     !
C HOW IT WORKS:                                                       !
C                                                                     !
C     EVAL uses several "sub-"subroutines to accomplish its goal.     !
C The process involves two steps:                                     !
C                                                                     !
C     1) a call to the parser 'MAKENT'                                !
C     2) a call to the evaluator 'GETENT'                             !
C                                                                     !
C The two processes are separated in order to allow one to define     !
C and store expressions which are functions of one or more variables. !
C Once defined, the function can be evaluated as often as needed      !
C through repeated calls to GETENT.                                   !
C                                                                     !
C                                                                     !
C PARSING:                                                            !
C                                                                     !
C     Input to MAKENT is a character string containing an arithmetic  !
C "sentence" of one of the forms described above.  The sentence is    !
C parsed into its component parts, containing a conditional           !
C expression and/or predicate expression.  These arithmetic           !
C expressions are in turn passed to the subprogram MAKEXP, one at a   !
C time.                                                               !
C                                                                     !
C     MAKEXP parses the expression into 'tokens' consisting of        !
C numeric constants (operands) and operators.  When all tokens are    !
C processed, they are converted from infix to postfix notation,       !
C reflecting the heirarchy of the operators (e.g. * is evaluated      !
C before +).  Associativity to the left is assumed for most binary    !
C operators, and to the right for unary operators.                    !
C                                                                     !
C     MAKENT encodes the sentence structure with some tokens of       !
C its own, between the tokens for the arithmetic expressions.  The    !
C completed token stream is then returned to to EVAL.  This           !
C represents the saved, executable form of an arithmetic sentence.    !
C                                                                     !
C     Three arrays are required to store the tokens:  one containing  !
C constants (operands), another for operators, and a third with       !
C pointers into the two lists in the order needed by GETENT.  It is   !
C up to the applications programmer to define the arrays and specify  !
C the corresponding size limits.                                      !
C                                                                     !
C                                                                     !
C EVALUATION:                                                         !
C                                                                     !
C     GETENT reads the tokens created during the parsing phase and    !
C decodes the sentence structure.  GETENT then calls GETEXP which     !
C performs the actual evaluation of an arithmetic expression.  The    !
C tokens are processed in sequence by GETEXP, stacking operands and   !
C perfoming operations as operators are encountered.  The result,     !
C if all went well, is a single number returned to GETENT.  The       !
C final operation is performed as dictated by the sentence structure. !
C                                                                     !
C     GETEXP and MAKEXP are constructed in such a way that they       !
C can be used independently of GETENT and MAKENT.  That is, the       !
C token stream output by MAKEXP is used as input to GETEXP, so one    !
C could, in fact, streamline the process of evaluating arithmetic     !
C expressions, so long as the added "sentence" constructs are not     !
C required.                                                           !
C                                                                     !
C                                                                     !
C AUXILIARY SUBPROGRAMS:                                              !
C                                                                     !
C     All operators are implemented in the subprograms GETOP and      !
C DOOP.  GETOP contains the names of the operators, and returns       !
C to MAKEXP the information saved in the token list.  DOOP defines    !
C the functions of those operators, and returns the result of a       !
C given operation to GETEXP.  To add, delete, or get more             !
C information on operators, see details there.                        !
C                                                                     !
C     Variables are implemented in the subprograms MAKVAR, SETVAR     !
C and GETVAR.  MAKVAR adds a variable name to the variable list if    !
C it hasn't already been defined, SETVAR defines the value of a       !
C given variable, and GETVAR returns the value of a variable.         !
C (These three subroutines may have to be modified to accomodate the  !
C specific storage scheme used in a particular application.)          !
C                                                                     !
C                                                                     !
C                      DIVISION OF BIOSTATISTICS                      !
C                        WASHINGTON UNIVERSITY                        !
C                              ST. LOUIS                              !
C                                                                     !
C                                                                     !
C---------------------------------------------------------------------+

      IMPLICIT INTEGER (A-W)

      INTEGER           MODE
      CHARACTER         STRING*(*)
      INTEGER           IEXNO
      INTEGER           MAXVAR
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)
      INTEGER           LENVAR(*)
      CHARACTER         NAMVAR(*)*(*)
      INTEGER           IVAR
      DOUBLE PRECISION  ANS
      INTEGER           ISTAT

#include "tx.h"

C
C DEFINE END-OF-EXPRESSION AND REJECTION CODES:
C
      DATA  EOX /999/
      DATA  REJECT /-999/
C
C MODE 0 : CLEAR VARIABLES AND TABLE SPACE
C
      IF (MODE .EQ. 0) THEN
C       CLEAR VARIABLES:
         DO 100 I=1,MAXVAR
            X(I) = BLANK(I)
            ITYPE(I) = 0
            NAMVAR(I) = ' '
            LENVAR(I) = 0
100      CONTINUE

C       CLEAR TOKEN SPACE:
         NTOKEN = 0
         NOPTOR = 0
         NCONST = 0
         ITOKEN(1) = EOX

C       EXIT:
         ISTAT = 0
         RETURN
C
C OTHERWISE, GET NUMBER OF ENTRIES TO BYPASS
C
      ELSE IF (MODE .EQ. 1) THEN
         BYPASS = MAXTOK
      ELSE IF (MODE .EQ. 2) THEN
         BYPASS = IEXNO - 1
      ELSE IF (MODE .EQ. 3) THEN
         BYPASS = IEXNO - 1
      ELSE
         call intpr('EVAL: INVALID CALL TO EVAL (MODE)',-1,0,0)
         STOP
      END IF
C
C FIND BEGINNING OF CURRENT ENTRY
C
      SAVE = 0
      IEXNO = 1
      DO 200 I=1,BYPASS
         IF (SAVE .GE. NTOKEN) GO TO 201
         TYPE = ITOKEN(SAVE+1)
         IF (TYPE .EQ. EOX) GO TO 201
         SAVE = SAVE + ITOKEN(SAVE+2) + 2
         IF (TYPE.LT.0 .AND. TYPE.NE.REJECT)
     @      SAVE = SAVE + ITOKEN(SAVE+1) + 1
         IEXNO = IEXNO + 1
200   CONTINUE
201   CONTINUE
C
C MODE 1 : TOKENIZE AND EVALUATE EXPRESSION, WITHOUT SAVING IT
C
      IF (MODE .EQ. 1) THEN

C       MAKE TEMPORARY COPIES OF COUNTERS:
         NTOK0 = NTOKEN
         NCON0 = NCONST
         NOPT0 = NOPTOR

C       TOKENIZE:
         CALL MAKENT (STRING,IQ, ITOKEN,NTOK0,MAXTOK,
     @      RCONST,NCON0,MAXCON, IOPTOR,NOPT0,MAXOPT,
     @      MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @      IVAR,ISTAT)
         IF (ISTAT .NE. 0) GO TO 800

C       EVALUATE:
         CALL GETENT (IQ,ITOKEN,NTOKEN,RCONST,IOPTOR,
     @      X,BLANK,ITYPE,
     @      ANS,ISTAT)
         IF (ISTAT .NE. 0) GO TO 800

C       RE-MARK END OF PREEXISTING EXPRESSION:
         ITOKEN(NTOKEN+1) = EOX
C
C MODE 2 : TOKENIZE AND SAVE EXPRESSION ONLY:
C
      ELSE IF (MODE .EQ. 2) THEN

C       RESET TOKEN COUNT:
         NTOKEN = SAVE

C       TOKENIZE:
         CALL MAKENT (STRING,IQ, ITOKEN,NTOKEN,MAXTOK,
     @      RCONST,NCONST,MAXCON, IOPTOR,NOPTOR,MAXOPT,
     @      MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @      IVAR,ISTAT)
         IF (ISTAT .NE. 0) GO TO 800

C       MARK END OF NEW EXPRESSION:
         ITOKEN(NTOKEN+1) = EOX
C
C MODE 3 : EVALUATE SAVED EXPRESSION:
C
      ELSE IF (MODE .EQ. 3) THEN

C       CHECK FOR END OF EXPRESSION LIST
         IF (ITOKEN(SAVE+1) .EQ. EOX) GO TO 810

C       EVALUATE SAVED EXPRESSION:
         CALL GETENT (IQ,ITOKEN,SAVE,RCONST,IOPTOR,
     @      X,BLANK,ITYPE,
     @      ANS,ISTAT)
         IF (ISTAT .NE. 0) GO TO 800

      END IF
      GO TO 900
C
C UNSUCCESSFUL RETURN (ISTAT ALREADY DEFINED)
C
800   GO TO 999
C
C END OF SAVED EXPRESSION LIST
C
810   ISTAT = 4
CD      WRITE (IQ,*) 'DEBUG: <END OF EXPRESSION LIST>'
      GO TO 999
C
C SUCCESSFUL RETURN
C
900   ISTAT = 0
999   RETURN
      END

      SUBROUTINE MAKVAR(STRING,ICOL,MODE,
     @   MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @   IVAR,ISTAT)
C
C LOOKS UP VARIABLE NAME IN 'STRING' AND RETURNS THE INDEX INTO THE
C VARIABLE TABLE.  VARIABLE 'MODE' CONTROLS HOW UNMATCHED NAMES ARE
C TO BE HANDLED (I.E. WHETHER NEW VARIABLES SHOULD BE ADDED TO THE
C TABLE.)
C
C INPUT ARGUMENTS:
C     STRING   -- INPUT STRING CONTAINING VARIABLE NAME
C     ICOL     -- COLUMN IN 'STRING' TO BEGIN PATTERN MATCH
C     MODE     -- IF THE NAME DOES NOT MATCH ANY IN THE VARIABLE TABLE:
C                 MODE=1 MEANS RETURN ERROR FLAG;
C                 MODE=2 MEANS CREATE A NEW ENTRY IN THE TABLE
C     MAXVAR   -- DIMENSION OF ARRAYS IN THE VARIABLE LIST
C     X()      -- LIST OF VARIABLE VALUES
C     BLANK()  -- LIST OF MISSING VALUE CODES BY VARIABLE
C     ITYPE()  -- LIST OF TYPE CODES BY VARIABLE AS FOLLOWS:
C                 0 = UNDEFINED VARIABLE
C                 1 = CHARACTER (ALPHANUMERIC)
C                 2 = FLOATING POINT NUMERIC
C                 3 = INTEGER
C     NAMVAR   -- LIST OF VARIABLE NAMES
C     LENVAR   -- LENGTHS OF ABOVE NAMES
C
C OUTPUT ARGUMENTS:
C     ICOL     -- PUSHED TO NEXT COLUMN FOLLOWING VARIABLE NAME
C     IVAR     -- INDEX INTO THE VARIABLE TABLE
C     ISTAT    -- 0=OK ; -1=ERROR
C
      CHARACTER         STRING*(*)
      INTEGER           ICOL
      INTEGER           MODE
      INTEGER           MAXVAR
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)
      INTEGER           LENVAR(*)
      CHARACTER         NAMVAR(*)*(*)
      INTEGER           IVAR
      INTEGER           ISTAT
      CHARACTER         C*1
      data jcol /0/
C
C SEE IF VARIABLE NAME EXISTS
C
      IUNDEF = 0
      IVAR = 0
      DO 100 I=1,MAXVAR
         IF (ITYPE(I) .GT. 1) THEN
            J = INDEX( STRING(ICOL:), NAMVAR(I)(1:LENVAR(I)) )
            IF (J .EQ. 1) THEN
               IVAR = I
               ICOL = ICOL + LENVAR(I)
               GO TO 900
            END IF
         ELSE IF (IUNDEF.EQ.0 .AND. ITYPE(I).EQ.0) THEN
            IUNDEF = I
         END IF
100   CONTINUE
      IF (MODE .NE. 2) GO TO 800
C
C CREATE IT IF NOT
C
C    FIND LENGTH OF VARIABLE NAME:
      N = 0
      DO 200 I=ICOL,LEN(STRING)
         C = STRING(I:I)
         IF ( (C.GE.'A' .AND. C.LE.'Z') .OR.
     @        (C.GE.'0' .AND. C.LE.'9' .AND. I.NE.ICOL) ) THEN
            N = N + 1
            JCOL = I
         ELSE
            GO TO 201
         END IF
200   CONTINUE
201   IF (N .LE. 0) GO TO 800

      IF (IUNDEF .LE. 0) THEN
         call intpr('EVAL: *ERROR* TOO MANY VARIABLES',-1,0,0)
         STOP
      END IF

C    ADD IT TO THE TABLE:
      IVAR = IUNDEF
      NAMVAR(IVAR) = STRING(ICOL:JCOL)
      LENVAR(IVAR) = N
      ITYPE(IVAR) = 2
      X(IVAR) = BLANK(IVAR)
CD      WRITE (*,1600) IVAR
CD 1600 FORMAT (' ---> CREATING NEW X(',I2,')')

      ICOL = ICOL + N
      GO TO 900
C
C UNSUCCESSFUL RETURN
C
800   ISTAT = -1
      GO TO 999
C
C SUCCESSFUL RETURN
C
900   ISTAT = 0
999   RETURN
      END

      SUBROUTINE MAKVA2(IVAR,
     @   MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @   ISTAT)
C
C THIS ALTERNATE FORM OF MAKVAR IS REQUIRED FOR COMPATABILITY WITH
C THE HAWAII-STYLE PROGRAMS IN WHICH VARIABLES ARE REFERENCED BY
C INDEX, RATHER THAN BY NAME.  THIS SUBROUTINE LOOKS UP VARIABLE
C 'IVAR' IN THE VARIABLE TABLE, AND CREATES THE ENTRY IF IT WAS
C UNDEFINED.
C
C INPUT ARGUMENTS:
C     IVAR     -- INDEX INTO THE VARIABLE TABLE
C     MAXVAR   -- DIMENSION OF ARRAYS IN THE VARIABLE LIST
C     X()      -- LIST OF VARIABLE VALUES
C     BLANK()  -- LIST OF MISSING VALUE CODES BY VARIABLE
C     ITYPE()  -- LIST OF TYPE CODES BY VARIABLE AS FOLLOWS:
C                 0 = UNDEFINED VARIABLE
C                 1 = CHARACTER (ALPHANUMERIC)
C                 2 = FLOATING POINT NUMERIC
C                 3 = INTEGER
C     NAMVAR   -- LIST OF VARIABLE NAMES
C     LENVAR   -- LENGTHS OF ABOVE NAMES
C
C OUTPUT ARGUMENTS:
C     ISTAT    -- 0=OK ; -1=ERROR
C
      INTEGER           IVAR
      INTEGER           MAXVAR
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)
      INTEGER           LENVAR(*)
      CHARACTER         NAMVAR(*)*(*)
      INTEGER           ISTAT

      IF (IVAR.LT.1 .OR. IVAR.GT.MAXVAR) GO TO 800

      IF (LENVAR(IVAR) .LE. 0) THEN
         NAMVAR(IVAR) = 'NONAME'
         LENVAR(IVAR) = 6
         ITYPE(IVAR) = 2
         X(IVAR) = BLANK(IVAR)
CD         WRITE (*,1600) IVAR
CD 1600    FORMAT (' ---> CREATING NEW X(',I2,')')
      END IF

      GO TO 900
C
C UNSUCCESSFUL RETURN
C
800   ISTAT = -1
      GO TO 999
C
C SUCCESSFUL RETURN
C
900   ISTAT = 0
999   RETURN
      END

      SUBROUTINE GETVAR(IVAR, X,BLANK,ITYPE,ANS, ISTAT)
C
C RETURNS VALUE OF VARIABLE IVAR
C
C INPUT ARGUMENTS:
C     IVAR     -- ORDINAL OF VARIABLE
C     X()      -- LIST OF VARIABLE VALUES
C     BLANK()  -- LIST OF MISSING VALUE CODES BY VARIABLE
C     ITYPE()  -- LIST OF TYPE CODES BY VARIABLE AS FOLLOWS:
C                 0 = UNDEFINED VARIABLE
C                 1 = CHARACTER (ALPHANUMERIC)
C                 2 = FLOATING POINT NUMERIC
C                 3 = INTEGER
C
C OUTPUT ARGUMENTS:
C     ANS   - VALUE RETURNED
C     ISTAT - SET TO -1 IF MISSING, 0 OTHERWISE
C
      DOUBLE PRECISION  ANS
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)

      ANS = 0.0
C
C WRONG VARIABLE TYPE (OR UNDEFINED) IS A FATAL ERROR
C
      IF ( ITYPE(IVAR) .NE. 2 .AND. ITYPE(IVAR) .NE. 3 ) THEN
         WRITE (*,1800) IVAR
 1800    FORMAT (' VARIABLE X(',I2,') IS UNINITIALIZED')
         GO TO 800
      END IF
C
C RETURN VALUE IN ANS
C
      ANS = X(IVAR)
C
C CHECK FOR MISSING VALUE
C
      IF (X(IVAR) .EQ. BLANK(IVAR)) GO TO 800
C
C ROUND TO NEAREST INTEGER IF VARIABLE TYPE CODE IS INTEGER
C
      IF ( ITYPE(IVAR) .EQ. 3 ) ANS = AINT( ANS + SIGN( 0.5D0, ANS ) )
      GO TO 900
C
C MISSING VALUE
C
800   ISTAT = -1
      GO TO 999
C
C SUCCESSFUL RETURN
C
900   ISTAT = 0
999   RETURN
      END

      SUBROUTINE SETVAR(IVAR,VAL,MODE, X,BLANK,ITYPE, ISTAT)
C
C DEFINES THE VALUE OF VARIABLE 'IVAR' IN THE VARIABLE TABLE
C
C INPUT ARGUMENTS:
C     IVAR     -- INDEX INTO THE VARIABLE TABLE
C     VAL      -- NEW VALUE OF VARIABLE
C     MODE     -- SPECIFIES WHETHER TO ASSIGN A MISSING VALUE CODE
C                 MODE=1 MEANS ASSIGN THE VALUE OF 'VAL'
C                 MODE=2 MEANS ASSIGN A MISSING VALUE CODE
C     X()      -- LIST OF VARIABLE VALUES
C     BLANK()  -- LIST OF MISSING VALUE CODES BY VARIABLE
C     ITYPE()  -- LIST OF TYPE CODES BY VARIABLE AS FOLLOWS:
C                 0 = UNDEFINED VARIABLE
C                 1 = CHARACTER (ALPHANUMERIC)
C                 2 = FLOATING POINT NUMERIC
C                 3 = INTEGER
C
C OUTPUT ARGUMENTS:
C     ISTAT    -- 0=OK ; -1=UNDEFINED VARIABLE
C
      DOUBLE PRECISION  VAL
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)
C
C CHECK MODE FOR VALIDITY
C
      IF (MODE .NE. 1 .AND. MODE .NE. 2) THEN
         call intpr('SETVAR: INVALID MODE',-1,0,0)
         STOP
      END IF
C
C WRONG VARIABLE TYPE (OR UNDEFINED) IS A FATAL ERROR
C
      IF ( ITYPE(IVAR) .NE. 2 .AND. ITYPE(IVAR) .NE. 3 ) THEN
         WRITE (*,1800) IVAR
 1800    FORMAT (' VARIABLE X(',I2,') IS OF WRONG TYPE FOR OPERATION')
         call intpr('SETVAR: ERROR IN VARIABLE TYPE',-1,0,0)
         STOP
      END IF
C
C SET NEW VALUE
C
      IF (MODE .EQ. 1) THEN
         IF ( ITYPE(IVAR) .EQ. 3 ) THEN
C           IF TYPE IS INTEGER, ROUND TO NEAREST WHOLE NUMBER
            X(IVAR) = AINT( VAL + SIGN( 0.5D0, VAL ) )
         ELSE
C           OTHERWISE SAVE FLOATING POINT NUMBER
            X(IVAR) = VAL
         END IF
      ELSE
         X(IVAR) = BLANK(IVAR)
      END IF
CD      WRITE (*,1600) IVAR, X(IVAR)
CD 1600 FORMAT (' X(',I2,')=',G20.12)
      GO TO 900
C
C UNSUCCESSFUL RETURN
C
800   ISTAT = -1
      GO TO 999
C
C SUCCESSFUL RETURN
C
900   ISTAT = 0
999   RETURN
      END

      SUBROUTINE MAKENT(STRING,IQ, TOKEN,NTOKEN,MTOKEN,
     @   CONST,NCONST,MCONST, OPTOR,NOPTOR,MOPTOR,
     @   MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @   IVAR,ISTAT)
C
C read an expression and save it as a tokenized entry
C
C INPUT ARGUMENTS:
C     IQ       -- logical unit to write error messages to
C     STRING   -- input character string containing expression
C     TOKEN    -- list of pointers to operators and operands
C     NTOKEN   -- current number of elements in TOKEN
C     MTOKEN   -- maximum number of elements in TOKEN
C     CONST    -- list of operands (numeric constants)
C     NCONST   -- current number of elements in CONST
C     MCONST   -- macimum number of elements in CONST
C     OPTOR    -- list of operators
C     NOPTOR   -- current number of elements in OPTOR
C     MOPTOR   -- maximum number of elements in OPTOR
C     MAXVAR   -- DIMENSION OF ARRAYS IN THE VARIABLE LIST
C     X()      -- LIST OF VARIABLE VALUES
C     BLANK()  -- LIST OF MISSING VALUE CODES BY VARIABLE
C     ITYPE()  -- LIST OF TYPE CODES BY VARIABLE AS FOLLOWS:
C                 0 = UNDEFINED VARIABLE
C                 1 = CHARACTER (ALPHANUMERIC)
C                 2 = FLOATING POINT NUMERIC
C                 3 = INTEGER
C     NAMVAR   -- LIST OF VARIABLE NAMES
C     LENVAR   -- LENGTHS OF ABOVE NAMES
C
C RETURNED:
C     TOKEN    -- previous TOKEN list with new information appended
C     NTOKEN   -- updated number of elements in TOKEN
C     CONST    -- previous CONST list with new operands appended
C     NCONST   -- updated number of elements in CONST
C     OPTOR    -- previous OPTOR list with new operators appended
C     NOPTOR   -- updated number of elements in OPTOR
C     IVAR     -- index number of destination variable (0 if none)
C     ISTAT    -- 0=ok; others: see below
C
C OVERVIEW:
C     This subroutine saves an expression in a tokenized form
C     in the arrays TOKEN, CONST and OPTOR, each with its respective
C     counter.  TOKEN is the master array, with pointers into the
C     other two.  An entry in 'TOKEN' beginning at 'SAVE' is coded
C     as follows:
C
C     TOKEN(SAVE+0)   = type of expression entry (see below)
C     TOKEN(SAVE+1)   = start of tokenized conditional expression
C     TOKEN(SAVE+c)   = end of conditional with 'c' tokens
C     TOKEN(SAVE+c+1) = start of tokenized predicate expression
C     TOKEN(SAVE+c+p) = end of predicate with 'p' tokens
C
C where the first token indicates the expression type as follows:
C     > 0   -- index of destination variable, unconditional
C     < 0   -- negative index of destination variable, conditional
C     = 0   -- no destination variable
C     =-999 -- conditional rejection
C     = 999 -- end of token list
C
      IMPLICIT INTEGER (A-W)

      CHARACTER         STRING*(*)
      INTEGER           TOKEN(*)
      DOUBLE PRECISION  CONST(*)
      INTEGER           OPTOR(3,*)
      INTEGER           MAXVAR
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)
      INTEGER           LENVAR(*)
      CHARACTER         NAMVAR(*)*(*)
C
C DEFINE REJECTION CODE:
C
      DATA  REJECT /-999/, jcol1 /0/, jcol2 /0/
C
C INITIALIZE
C
      NTOK0 = NTOKEN
      NCON0 = NCONST
      NOPT0 = NOPTOR

      IVAR = 0

      TYPE = 1
C
C RESERVE FIRST TOKEN FOR TYPE OF ENTRY
C
      NTOK0 = NTOK0 + 1
      HOLD = NTOK0
      IF (NTOK0 .GE. MTOKEN) THEN
         call intpr('EVAL: TOO MANY SAVED EXPRESSIONS',-1,0,0)
         STOP
      END IF
C
C HANDLE CONDITIONAL PORTION OF EXPRESSION
C
      ICOL = 1
      ICOL = NSPACE(STRING,ICOL,ISTAT)
      IF (ISTAT .NE. 0) GO TO 810

      IF (STRING(ICOL:ICOL+1) .EQ. 'IF') THEN

C       FIND BEGINNING AND END OF CONDITIONAL EXPRESSION:
         ICOL = ICOL + 2
         CALL PAR(STRING,ICOL,ILEFT,IRIGHT,ISTAT)
         IF (ISTAT .NE. 0) GO TO 820

C       TOKENIZE CONDITIONAL, ADDING LIST OF TOKENS TO 'TOKEN':
         CALL MAKEXP( STRING(ILEFT+1:IRIGHT-1),IQ, TOKEN,NTOK0,MTOKEN,
     @      CONST,NCON0,MCONST, OPTOR,NOPT0,MOPTOR,
     @      MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @      ISTAT)
         IF (ISTAT .NE. 0) GO TO 810

C       PUSH ICOL PAST CONDITIONAL:
         ICOL = IRIGHT + 1
         ICOL = NSPACE(STRING,ICOL,ISTAT)
         IF (ISTAT .NE. 0) GO TO 820

C       INDICATE THAT THIS IS A CONDITIONAL:
         TYPE = -1
      END IF
C
C PROCESS PREDICATE EXPRESSION
C
      IF (STRING(ICOL:) .EQ. 'REJECT') THEN
C       INDICATE CONDITIONAL REJECTION:
         IF (TYPE .NE. -1) GO TO 820
         IVAR = -1
         TYPE = REJECT
      ELSE
C
C IF AN EQUAL SIGN APPEARS IN THE PREDICATE, ASSUME THAT
C A DESTINATION VARIABLE IS SPECIFIED AND PROCESS THAT LAST
C
         JCOL = INDEX(STRING(ICOL:),'=')
         IF (JCOL .GT. 0) THEN
            JCOL1 = ICOL
            JCOL2 = ICOL + JCOL - 1
            ICOL  = JCOL2 + 1
         ELSE
            IF (TYPE .NE. 1) GO TO 830
            TYPE = 0
         END IF
C
C PROCESS ARITHMETIC PORTION OF PREDICATE
C
         CALL MAKEXP ( STRING(ICOL:),IQ, TOKEN,NTOK0,MTOKEN,
     @      CONST,NCON0,MCONST, OPTOR,NOPT0,MOPTOR,
     @      MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @      ISTAT)

         IF (ISTAT .NE. 0) GO TO 810
C
C PROCESS DESTINATION VARIABLE
C
         IF (JCOL .GT. 0) THEN
            ICOL = JCOL1

C          HANDLE TWO FORMS OF VARIABLE REFERNECES:
            IF (STRING(ICOL:ICOL) .EQ. 'X') THEN

C             VARIABLE IS IN INDEX FORM, 'Xn' :
               ICOL = ICOL + 1
               IVAR = IXNO(STRING,ICOL,ISTAT)

               ICOL = NSPACE(STRING,ICOL,ISTAT)
               IF (ISTAT .NE. 0) GO TO 830
               IF (ICOL .NE. JCOL2) GO TO 830

               CALL MAKVA2(IVAR,
     @            MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @            ISTAT)
               IF (ISTAT .NE. 0) GO TO 830

            ELSE
C             VARIABLE IS IN NAME FORM:
               MODE = 2
               CALL MAKVAR(STRING,ICOL,MODE,
     @            MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @            IVAR,ISTAT)
               IF (ISTAT .NE. 0) GO TO 830

               ICOL = NSPACE(STRING,ICOL,ISTAT)
               IF (ISTAT .NE. 0) GO TO 830
               IF (ICOL .NE. JCOL2) GO TO 830
            END IF

C          ENCODE DESTINATION VARIABLE IN TYPE:
            TYPE = TYPE * IVAR

         END IF

      END IF
C
C FINALIZE SAVED EXPRESSION
C
      NTOKEN = NTOK0
      NCONST = NCON0
      NOPTOR = NOPT0

      TOKEN(HOLD) = TYPE
CD      WRITE (IQ,*) 'DEBUG: ENTRY TYPE',TYPE,' ->',HOLD,' THROUGH',NTOKEN

      GO TO 900
C
C UNSUCCESSFUL RETURN
C
810   WRITE (IQ,*) '** SYNTAX ERROR IN EXPRESSION **'
      GO TO 999

820   WRITE (IQ,*) '** INVALID FORM FOR CONDITIONAL EXPRESSION **'
      GO TO 890

830   WRITE (IQ,*) '** INVALID DESTINATION VARIABLE SPECIFICATION **'
      GO TO 890

840   WRITE (IQ,*) '** EMPTY EXPRESSION **'
      GO TO 999

890   ISTAT = -1
      CALL ERRLOC(IQ,STRING,ICOL)
      GO TO 999
C
C SUCCESSFUL RETURN
C
900   ISTAT = 0
999   RETURN
      END

      SUBROUTINE GETENT(IQ,TOKEN,NTOKEN,CONST,OPTOR,
     @   X,BLANK,ITYPE,
     @   ANS,ISTAT)
C
C evaluate a saved entry in the tokenized expression list
C
C INPUT ARGUMENTS:
C     IQ       -- logical unit to write error messages to
C     TOKEN    -- tokenized expression
C     NTOKEN   -- number of tokens already processed (0 if none)
C     CONST    -- list of operands (numeric constants)
C     OPTOR    -- list of operators
C     X()      -- LIST OF VARIABLE VALUES
C     BLANK()  -- LIST OF MISSING VALUE CODES BY VARIABLE
C     ITYPE()  -- LIST OF TYPE CODES BY VARIABLE AS FOLLOWS:
C                 0 = UNDEFINED VARIABLE
C                 1 = CHARACTER (ALPHANUMERIC)
C                 2 = FLOATING POINT NUMERIC
C                 3 = INTEGER
C
C OUTPUT ARGUMENTS:
C     ANS      -- evaluated result
C     ISTAT    -- 0=ok; others: see below
C
      IMPLICIT INTEGER (A-W)

      INTEGER           TOKEN(*)
      DOUBLE PRECISION  CONST(*)
      INTEGER           OPTOR(3,*)
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)
      DOUBLE PRECISION  ANS

      DOUBLE PRECISION  XC
C
C DEFINE REJECTION CODE
C
      DATA  REJECT /-999/
C
C INITIALIZE RESULT
C
      ANS = 0.0
C
C FIND OUT WHICH TYPE OF ENTRY THIS IS
C (SEE SUBPROGRAM MAKENT FOR DETAILS ON TYPE CODES)
C
      IF (NTOKEN.LT.0) THEN
         call intpr('EVAL: INVALID CALL TO GETENT (NTOKEN)',-1,0,0)
         STOP
      END IF

      SAVE = NTOKEN + 1
      TYPE = TOKEN(SAVE)
CD      WRITE (IQ,*) 'DEBUG: ENTRY TYPE',TYPE,' EVALUATION BEGINS AT',SAVE
C
C EVALUATE CONDITIONAL, VAMOSE IF FALSE
C
      IF (TYPE .LT. 0) THEN
         CALL GETEXP (TOKEN,SAVE,CONST,OPTOR,
     @      X,BLANK,ITYPE,
     @      XC,ISTAT)
         IF (ISTAT .NE. 0) GO TO 820
         IF (ABS(XC) .LT. 0.5) GO TO 830
      END IF
C
C HANDLE REJECTION
C
      IF (TYPE .EQ. REJECT) GO TO 850
C
C EVALUATE PREDICATE
C
      CALL GETEXP (TOKEN,SAVE,CONST,OPTOR,
     @   X,BLANK,ITYPE,
     @   ANS,ISTAT)
      IF (ISTAT .NE. 0 .AND. ISTAT .NE. -2) GO TO 810
C
C SAVE RESULT IN DESTINATION VARIABLE IF INDICATED
C
      IVAR = IABS(TYPE)
      IF (IVAR .GT. 0) THEN
         IF (ISTAT .NE. -2) THEN
CD            WRITE (IQ,*) 'DEBUG: SAVING',ANS,' IN VARIABLE',IVAR
            MODE = 1
            CALL SETVAR(IVAR,ANS,MODE, X,BLANK,ITYPE, ISTAT)
            IF (ISTAT .NE. 0) GO TO 810
         ELSE
CD            WRITE (IQ,*) 'DEBUG: SAVING MISSING CODE IN VARIABLE',IVAR
            MODE = 2
            CALL SETVAR(IVAR,ANS,MODE, X,BLANK,ITYPE, ISTAT)
            IF (ISTAT .NE. 0) GO TO 810
            GO TO 840
         END IF
      END IF

      GO TO 900
C
C UNSUCCESSFUL RETURN
C
810   ISTAT = -2
      WRITE (IQ,*) '** ERROR IN EXPRESSION EVALUATION **'
      GO TO 890

820   ISTAT = -3
      WRITE (IQ,*) '** ERROR IN EVALUATION OF CONDITION **'
      GO TO 890

830   ISTAT = 1
CD      WRITE (IQ,*) 'DEBUG: <CONDITION FALSE>'
      GO TO 890

840   ISTAT = 2
CD      WRITE (IQ,*) 'DEBUG: <MISSING VALUE>'
      GO TO 890

850   ISTAT = 3
CD      WRITE (IQ,*) 'DEBUG: <REJECTION>'
      GO TO 890

890   GO TO 999
C
C SUCCESSFUL RETURN
C
900   ISTAT = 0
999   RETURN
      END

      SUBROUTINE MAKEXP(STRING,IQ, TOKEN,NTOKEN,MTOKEN,
     @   CONST,NCONST,MCONST, OPTOR,NOPTOR,MOPTOR,
     @   MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @   ISTAT)
C
C parse input string into operator/operand codes, handling parentheses
C
C This subroutine updates three arrays and their associated
C counters.  The arrays are updated by appending the required
C information to them, and by incrementing their counters.
C The current expression is coded as follows:
C
C     TOKEN(NTOKEN+1)   = number of tokens (t) in expression
C     TOKEN(NTOKEN+2)   = start of tokenized expression
C     TOKEN(NTOKEN+t+1) = end of tokenized expression
C
C     CONST(NCONST+1)   = start of (c) new constants
C     CONST(NCONST+c)   = end of new constants
C
C     OPTOR(NOPTOR+1)   = start of (o) new operators
C     OPTOR(NOPTOR+o)   = end of new operators
C
C
C INPUT ARGUMENTS:
C     IQ       -- logical unit to write error messages to
C     STRING   -- input character string containing expression
C     TOKEN    -- list of pointers to operators and operands
C     NTOKEN   -- current number of elements in TOKEN
C     MTOKEN   -- maximum number of elements in TOKEN
C     CONST    -- list of operands (numeric constants)
C     NCONST   -- current number of elements in CONST
C     MCONST   -- maximum number of elemetts in CONST
C     OPTOR    -- list of operators
C     NOPTOR   -- current number of elements in OPTOR
C     MOPTOR   -- maximum number of elements in OPTOR
C     MAXVAR   -- DIMENSION OF ARRAYS IN THE VARIABLE LIST
C     X()      -- LIST OF VARIABLE VALUES
C     BLANK()  -- LIST OF MISSING VALUE CODES BY VARIABLE
C     ITYPE()  -- LIST OF TYPE CODES BY VARIABLE AS FOLLOWS:
C                 0 = UNDEFINED VARIABLE
C                 1 = CHARACTER (ALPHANUMERIC)
C                 2 = FLOATING POINT NUMERIC
C                 3 = INTEGER
C     NAMVAR   -- LIST OF VARIABLE NAMES
C     LENVAR   -- LENGTHS OF ABOVE NAMES
C
C RETURNED:
C     TOKEN    -- previous TOKEN list with new information appended
C     NTOKEN   -- updated number of elements in TOKEN
C     CONST    -- previous CONST list with new operands appended
C     NCONST   -- updated number of elements in CONST
C     OPTOR    -- previous OPTOR list with new operators appended
C     NOPTOR   -- updated number of elements in OPTOR
C     ISTAT    -- 0=OK; -1=ERROR
C
      IMPLICIT INTEGER (A-W)

      PARAMETER   (MAXSTK=20)

      CHARACTER         STRING*(*)
      INTEGER           TOKEN(*)
      DOUBLE PRECISION  CONST(*)
      INTEGER           OPTOR(3,*)
      INTEGER           MAXVAR
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)
      INTEGER           LENVAR(*)
      CHARACTER         NAMVAR(*)*(*)
      INTEGER           STACK(MAXSTK)
      DOUBLE PRECISION  XTEMP
      DOUBLE PRECISION  XNO
C
C PARENS IS CURRENT NESTING LEVEL OF PARENTHESES
C FLAG = 1 EXPECTS NEXT TOKEN TO BE AN OPERAND (OR UNARY OPERATOR)
C FLAG =-1 EXPECTS NEXT TOKEN TO BE A BINARY OPERATOR
C
      PARENS = 0
      FLAG = 1
C
C RESERVE FIRST TOKEN ENTRY FOR NUMBER OF TOKENS IN THIS EXPRESSION
C
      IF (NTOKEN .GE. MTOKEN) THEN
         call intpr('EVAL: TOKEN LIST OVERFLOW',-1,0,0)
         STOP
      END IF
      NTOKEN = NTOKEN + 1
      HOLD   = NTOKEN
C
C PARSE STRING INTO TOKENS
C
      ICOL = 1
100   CONTINUE
C
C     ADVANCE TO NEXT NON-BLANK
C
CCC         I = VERIFY(STRING(ICOL:),' ')
CCC         IF (I.LE.0) GO TO 199
CCC         ICOL = ICOL + I - 1
         ICOL = NSPACE(STRING,ICOL,ISTAT)
         IF (ISTAT .NE. 0) GO TO 199
C
C     HANDLE PARENTHESES:
C
         IF      ( STRING(ICOL:ICOL) .EQ. '(') THEN
            IF (FLAG .NE. 1) THEN
               WRITE (IQ,*) 'ILLEGAL PLACEMENT OF ('
               GO TO 800
            END IF
            PARENS = PARENS + 1
            ICOL = ICOL + 1

         ELSE IF ( STRING(ICOL:ICOL) .EQ. ')') THEN
            IF (FLAG .NE. -1) THEN
               WRITE (IQ,*) 'ILLEGAL PLACEMENT OF )'
               GO TO 800
            END IF
            PARENS = PARENS - 1
            IF (PARENS .LT. 0) THEN
               WRITE (IQ,*) 'UNMATCHED )'
               GO TO 800
            END IF
            ICOL = ICOL + 1

         ELSE
C
C     IF WE'RE EXPECTING AN OPERAND, CHECK FOR THAT FIRST.  IF WE
C     FAIL, CHECK SECONDLY FOR A UNARY OPERATOR.
C
            IF (FLAG .EQ. 1) THEN
C             CHECK FOR A NUMBER:
               XTEMP = XNO(STRING,ICOL,ISTAT)
               IF (ISTAT .EQ. 0) THEN
C                ADD IT TO CONSTANT LIST:
                  IF (NCONST .GE. MCONST) THEN
                     call intpr('EVAL: CONSTANT LIST OVERFLOW',
     @                    -1,0,0)
                     STOP
                  END IF
                  NCONST = NCONST + 1
                  CONST(NCONST) = XTEMP
                  NEWTOK = NCONST
                  FLAG = -1

C             CHECK FOR A UNARY OPERATOR:
               ELSE
                  NANDS = 1
                  CALL GETOP(STRING,ICOL,NANDS,
     @               MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @               CODE,PRI,ISTAT)
                  IF (ISTAT.EQ.0) THEN
C                   GET PRIORITY (RIGHT ASSOCIATIVE)
                     PRI = PRI + PARENS*100 + NOPTOR
C                   ADD IT TO OPERATOR LIST:
                     IF (NOPTOR .GE. MOPTOR) THEN
                        call intpr('EVAL: OPERATOR LIST OVERFLOW',
     @                       -1,0,0)
                        STOP
                     END IF
                     NOPTOR = NOPTOR + 1
                     OPTOR(1,NOPTOR) = CODE
                     OPTOR(2,NOPTOR) = PRI
                     OPTOR(3,NOPTOR) = NANDS
                     NEWTOK = -NOPTOR
                     IF (NANDS .EQ. 0) FLAG = -1
                  ELSE
                     WRITE (IQ,*) 'ILLEGAL OPERAND'
                     GO TO 800
                  END IF
               END IF
C
C     CHECK FOR A BINARY OPERATOR
C
            ELSE
               NANDS = 2
               CALL GETOP(STRING,ICOL,NANDS,
     @            MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @            CODE,PRI,ISTAT)
               IF (ISTAT.EQ.0) THEN
C                GET PRIORITY (LEFT ASSOCIATIVE UNLESS PRI=10)
                  IF (PRI .EQ. 10) PRI = PRI + NOPTOR
                  PRI = PRI + PARENS*100
C                ADD IT TO OPERATOR LIST:
                  IF (NOPTOR .GE. MOPTOR) THEN
                     call intpr('EVAL: OPERATOR LIST OVERFLOW',
     @                    -1,0,0)
                     STOP
                  END IF
                  NOPTOR = NOPTOR + 1
                  OPTOR(1,NOPTOR) = CODE
                  OPTOR(2,NOPTOR) = PRI
                  OPTOR(3,NOPTOR) = NANDS
                  NEWTOK = -NOPTOR
               ELSE
                  WRITE (IQ,*) 'ILLEGAL OPERATOR'
                  GO TO 800
               END IF
               FLAG = 1
            END IF

            IF (NTOKEN .GE. MTOKEN) THEN
               call intpr('EVAL: TOKEN LIST OVERFLOW',-1,0,0)
               STOP
            END IF
            NTOKEN = NTOKEN + 1
            TOKEN(NTOKEN) = NEWTOK

         END IF
      GO TO 100
199   CONTINUE

      IF (PARENS .NE. 0) THEN
         WRITE (IQ,*) 'UNMATCHED ('
         GO TO 800
      END IF

      IF (FLAG .EQ. 1) THEN
         WRITE (IQ,*) 'MISSING OPERAND'
         GO TO 800
      END IF
C
C SAVE NUMBER OF TOKENS USED IN THIS EXPRESSION
C
      TOKEN(HOLD) = NTOKEN - HOLD

CD      WRITE (IQ,*) 'DEBUG: TOKENS IN INFIX FORM:'
CD      WRITE (IQ,*) (TOKEN(I), I=HOLD,NTOKEN)

C ---------------------------------------------------------------------
C     CONVERT INFIX EXPRESSION TO POSTFIX ('REVERSE POLISH')
C ---------------------------------------------------------------------
C
C INITIALIZE THE STACK AND OUTPUT TOKEN COUNTER
C
      STKPTR = 0
      NOUT = HOLD
C
C LOOP THROUGH ALL INPUT TOKENS
C
      DO 400 I=HOLD+1,NTOKEN
C
C     COPY NUMERIC OPERANDS DIRECTLY TO OUTPUT LIST
C
         IF (TOKEN(I) .GT. 0) THEN
            NOUT = NOUT + 1
            TOKEN(NOUT) = TOKEN(I)
C
C     UNSTACK OPERATORS OF EQUAL OR GREATER PRECEDENCE
C
         ELSE
            OURPRI = OPTOR(2, -TOKEN(I) )
200         CONTINUE
               IF (STKPTR .LE. 0) GO TO 299
               STKPRI = OPTOR(2, -STACK(STKPTR) )
               IF (STKPRI .LT. OURPRI) GO TO 299
               NOUT = NOUT + 1
               TOKEN(NOUT) = STACK(STKPTR)
               STKPTR = STKPTR - 1
            GO TO 200
299         CONTINUE
C
C ADD CURRENT OPERATOR TO THE STACK
C
            IF (STKPTR .GE. MAXSTK) THEN
               call intpr('EVAL: STACK OVERFLOW',-1,0,0)
               STOP
            END IF
            STKPTR = STKPTR + 1
            STACK(STKPTR) = TOKEN(I)
         END IF

400   CONTINUE
C
C EMPTY THE STACK
C
500   CONTINUE
         IF (STKPTR .LE. 0) GO TO 599
         NOUT = NOUT + 1
         TOKEN(NOUT) = STACK(STKPTR)
         STKPTR = STKPTR - 1
      GO TO 500
599   CONTINUE
      IF (NOUT .NE. NTOKEN) THEN
         call intpr('CONVRT: INTERNAL ERROR',-1,0,0)
         STOP
      END IF

CD      WRITE (IQ,*) 'DEBUG: TOKENS IN POSTFIX ORDER:'
CD      WRITE (IQ,*) (TOKEN(I), I=HOLD,NTOKEN)
      GO TO 900
C
C UNSUCCESSFUL RETURN
C
800   ISTAT = -1
      CALL ERRLOC (IQ,STRING,ICOL)
      GO TO 999
C
C SUCCESSFUL RETURN
C
900   ISTAT = 0
999   RETURN
      END

      SUBROUTINE GETEXP(TOKEN,NTOKEN,CONST,OPTOR,
     @   X,BLANK,ITYPE,
     @   ANS,ISTAT)
C
C evaluate the postfix expression
C
C INPUT ARGUMENTS:
C     TOKEN    -- tokenized expression
C     NTOKEN   -- number of tokens already processed (0 if none)
C     CONST    -- list of operands (numeric constants)
C     OPTOR    -- list of operators
C     X()      -- LIST OF VARIABLE VALUES
C     BLANK()  -- LIST OF MISSING VALUE CODES BY VARIABLE
C     ITYPE()  -- LIST OF TYPE CODES BY VARIABLE AS FOLLOWS:
C                 0 = UNDEFINED VARIABLE
C                 1 = CHARACTER (ALPHANUMERIC)
C                 2 = FLOATING POINT NUMERIC
C                 3 = INTEGER
C
C OUTPUT ARGUMENTS:
C     ANS      -- evaluated result
C     NTOKEN   -- updated number of tokens processed
C     ISTAT    -- 0=ok; -1=error, -2=missing value
C
      IMPLICIT INTEGER (A-W)

      PARAMETER   (MAXSTK=20)

      INTEGER           TOKEN(*)
      DOUBLE PRECISION  CONST(*)
      INTEGER           OPTOR(3,*)
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)
      DOUBLE PRECISION  ANS
      DOUBLE PRECISION  STACK(MAXSTK)
      DOUBLE PRECISION  XBLANK
      DOUBLE PRECISION  DOOP, XARG, YARG

      DATA  XBLANK  /-99999./
C
C INITIALIZE THE STACK
C
      STKPTR = 0
C
C FIND FIRST AND LAST TOKENS TO PROCESS
C
      FIRST = NTOKEN + 2
      LAST  = NTOKEN + TOKEN(NTOKEN+1) + 1

      IF (NTOKEN .LT. 0) THEN
         call intpr('INVALID CALL TO GETEXP (NTOKEN)',-1,0,0)
         STOP
      END IF
      IF (TOKEN(NTOKEN+1) .LE. 0) THEN
         call intpr('INVALID CALL TO GETEXP (TOKEN)',-1,0,0)
         STOP
      END IF
C
C LOOP THROUGH ALL TOKENS IN CURRENT EXPRESSION
C
      DO 100 I=FIRST,LAST
C
C     STACK PREVIOUS RESULT
C
         IF (I .GT. FIRST) THEN
            STKPTR = STKPTR + 1
            IF (STKPTR .GE. MAXSTK) THEN
               call intpr('EVAL: STACK OVERFLOW',-1,0,0)
               STOP
            END IF
            STACK(STKPTR) = XARG
         END IF
C
C     GET OPERAND
C
         IF (TOKEN(I) .GT. 0) THEN
            XARG = CONST( TOKEN(I) )
C
C     GET OPERATOR: PERFORM INDICATED OPERATION
C
         ELSE
            CODE  = OPTOR(1, -TOKEN(I) )
            NANDS = OPTOR(3, -TOKEN(I) )

C        GET REQUIRED NUMBER OF OPERANDS FROM STACK:
            YARG = 0.0
            XARG = 0.0
            IF (NANDS .GT. 1) THEN
               YARG = STACK(STKPTR)
               STKPTR = STKPTR - 1
            END IF
            IF (NANDS .GT. 0) THEN
               XARG = STACK(STKPTR)
               STKPTR = STKPTR - 1
            END IF

C        PERFORM OPERATION:
            XARG = DOOP(CODE,XARG,YARG,
     @         X,BLANK,ITYPE,
     @         ISTAT)
            IF (ISTAT .NE. 0) GO TO 810
            IF (I.NE.LAST) THEN
CD               WRITE (*,*) 'DEBUG:', STKPTR+1,
CD     @            ') INTERMEDIATE X =', XARG
            END IF
         END IF

100   CONTINUE
      IF (STKPTR .NE. 0) THEN
         call intpr('POSTFX: STACK ERROR',-1,0,0)
         STOP
      END IF
C
C CLEAN UP AND RETURN
C
      ANS = XARG
      NTOKEN = LAST
      IF (ANS .EQ. XBLANK) GO TO 820
      GO TO 900
C
C UNSUCCESSFUL RETURN
C
C    ARITHMETIC ERROR:
810   ISTAT = -1
      GO TO 999

C    MISSING VALUE:
820   ISTAT = -2
      GO TO 999
C
C SUCCESSFUL RETURN
C
900   ISTAT = 0
999   RETURN
      END

      SUBROUTINE GETOP(OPER,ICOL,NANDS,
     @   MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @   CODE,PRI,ISTAT)
C
C looks up operator 'OPER' and returns the operation code,
C priority, and number of operands required
C
C
C INPUT ARGUMENTS:
C     OPER     -- input string
C     ICOL     -- column in 'OPER' to begin pattern match
C     NANDS    -- number of operands required for operator
C     MAXVAR   -- DIMENSION OF ARRAYS IN THE VARIABLE LIST
C     X()      -- LIST OF VARIABLE VALUES
C     BLANK()  -- LIST OF MISSING VALUE CODES BY VARIABLE
C     ITYPE()  -- LIST OF TYPE CODES BY VARIABLE AS FOLLOWS:
C                 0 = UNDEFINED VARIABLE
C                 1 = CHARACTER (ALPHANUMERIC)
C                 2 = FLOATING POINT NUMERIC
C                 3 = INTEGER
C     NAMVAR   -- LIST OF VARIABLE NAMES
C     LENVAR   -- LENGTHS OF ABOVE NAMES
C
C OUTPUT ARGUMENTS:
C     ICOL     -- pushed to next column following operator
C     CODE     -- internal code
C     PRI      -- evaluation precedence
C     ISTAT    -- 0=ok ; -1=error
C
C
C ENTERING NEW OPERATORS:
C
C the table of operator attributes is defined as follows:
C
C     OPSTR    -- string representation
C     OPLEN    -- length of string
C     OPCODE   -- internal codes
C     OPPRI    -- evaluation precedence
C     OPNUM    -- number of operands
C     -- (all operators of priority 10 are right associative)
C
C     Enter new operators into table in order of length (OPSTR)
C     if there is a possibly of ambiguity in the names (e.g. LOG
C     must preceed LOG10 in order to be properly interpreted).
C     please revise the version date below if changes are made
C
C
C NOTE: Exponentiation as defined, associates to the left, rather
C       than to the right because of confusion in unparenthesized
C       unary expressions; so "LOG(X)**2" becomes "(LOG X)**2"
C       instead of LOG(X**2).  Be careful to explicitly parenthesize
C       expressions of the form "X**Y**Z".

      PARAMETER (NOPS=30)

      CHARACTER         OPER*(*)
      INTEGER           ICOL
      INTEGER           NANDS
      INTEGER           MAXVAR
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)
      INTEGER           LENVAR(*)
      CHARACTER         NAMVAR(*)*(*)
      INTEGER           CODE
      INTEGER           PRI
      INTEGER           ISTAT

      CHARACTER         OPSTR(NOPS)*5
      INTEGER           OPLEN(NOPS)
      INTEGER           OPCODE(NOPS)
      INTEGER           OPPRI(NOPS)
      INTEGER           OPNUM(NOPS)

      DATA OPSTR /
     @ '+','-','*','/','-', 'X',
     @ 'EQ','LT','GT','LE','GE','NE','OR','AND','BLANK',
     @ '**','PI','LOG','EXP','MOD','INT','RND','SQRT','LOG10',
     @ 'SIN','COS','TAN','ASIN','ACOS','ATAN'/

      DATA OPLEN /
     @   1,  1,  1,  1,  1,   1,
     @    2,   2,   2,   2,   2,   2,   2,    3,      5,
     @    2,   2,    3,    3,    3,    3,    3,     4,      5,
     @     3,    3,    3,     4,     4,     4/

      DATA OPCODE /
     @   1,  2,  3,  4,  5,  99,
     @   51,  52,  53,  54,  55,  56,  58,   57,     59,
     @    6,   7,    8,    9,   10,   11,   12,    13,     14,
     @    15,   16,   17,    18,    19,    20/

      DATA OPPRI /
     @   5,  5,  6,  6, 10,  10,
     @    3,   3,   3,   3,   3,   3,   1,    2,     10,
     @    7,  10,   10,   10,    7,   10,   10,    10,     10,
     @    10,   10,   10,    10,    10,    10/

      DATA OPNUM /
     @   2,  2,  2,  2,  1,   1,
     @    2,   2,   2,   2,   2,   2,   2,    2,      1,
     @    2,   0,    1,    1,    2,    1,    0,     1,     1,
     @     1,    1,    1,     1,     1,     1/
C
C INITIALIZE
C
      ISTAT = 0
      NLEFT = LEN(OPER) - ICOL + 1
C
C SEARCH OPSTR TABLE FOR A MATCH
C
      DO 100 I=NOPS,1,-1
         NMATCH = OPLEN(I)
         ICOL2 = ICOL + NMATCH - 1
         IF (NLEFT .GE. NMATCH .AND.
     @       OPER(ICOL:ICOL2) .EQ. OPSTR(I)) THEN
            IF ((NANDS .EQ. OPNUM(I)) .OR.
     @          (NANDS .EQ. 1 .AND. OPNUM(I) .LE. 1)) THEN
               CODE  = OPCODE(I)
               PRI   = OPPRI(I)
               NANDS = OPNUM(I)
               ICOL  = ICOL + NMATCH
               RETURN
            END IF
         END IF
100   CONTINUE
C
C SEE IF IT'S A VARIABLE NAME
C
      IF (NANDS .EQ. 1) THEN
         MODE = 1
         CALL MAKVAR(OPER,ICOL,MODE,
     @      MAXVAR,X,BLANK,ITYPE,NAMVAR,LENVAR,
     @      IVAR,ISTAT)
         IF (ISTAT .EQ. 0) THEN
            CODE  = 100 + IVAR
            PRI   = 10
            NANDS = 0
            RETURN
         END IF
      END IF
C
C SPECIAL OPERATOR "?" DUMPS TABLE OF VALID OPERATORS
C
      IF (OPER(ICOL:ICOL) .EQ. '?') THEN
         WRITE (*,*) 'LIST OF OPERATORS -- VERSION DATE: 6/88'
         WRITE (*,*) (OPSTR(I),I=1,NOPS)
         STOP
      END IF
C
C OPERATOR NOT FOUND
C
      ISTAT = -1

      RETURN
      END

      DOUBLE PRECISION FUNCTION DOOP(CODE,XARG,YARG,
     @   X,BLANK,ITYPE,
     @   ISTAT)
C
C performs operation 'CODE' on XARG and YARG
C returns ISTAT=0 if ok; =-1 if error
C
C AUX INPUT VARIABLES:
C     X()      -- LIST OF VARIABLE VALUES
C     BLANK()  -- LIST OF MISSING VALUE CODES BY VARIABLE
C     ITYPE()  -- LIST OF TYPE CODES BY VARIABLE AS FOLLOWS:
C                 0 = UNDEFINED VARIABLE
C                 1 = CHARACTER (ALPHANUMERIC)
C                 2 = FLOATING POINT NUMERIC
C                 3 = INTEGER
C
C missing values are handled as follows:
C -- in arithmetic expressions, the numeric blank code is returned;
C -- in logical expressions, the numeric code is substituted
C    and the true/false result appropriately returned.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER           CODE
      DOUBLE PRECISION  X(*)
      DOUBLE PRECISION  BLANK(*)
      INTEGER           ITYPE(*)

      LOGICAL           NOX, NOY

      DATA TRUE, FALSE /1.0,0.0/
      DATA XBLANK /-99999./

      ANS = XBLANK

      IF (XARG .EQ. XBLANK) THEN
         NOX = .TRUE.
      ELSE
         NOX = .FALSE.
      END IF

      IF (YARG .EQ. XBLANK) THEN
         NOY = .TRUE.
      ELSE
         NOY = .FALSE.
      END IF
C
C ARITHMETIC OPERATORS
C
C ADD
      IF      (CODE .EQ. 1) THEN
         IF (NOX .OR. NOY) GO TO 820
         ANS = XARG + YARG
C SUBTRACT
      ELSE IF (CODE .EQ. 2) THEN
         IF (NOX .OR. NOY) GO TO 820
         ANS = XARG - YARG
C MULTIPLY
      ELSE IF (CODE .EQ. 3) THEN
         IF (NOX .OR. NOY) GO TO 820
         ANS = XARG * YARG
C DIVIDE
      ELSE IF (CODE .EQ. 4) THEN
         IF (NOX .OR. NOY) GO TO 820
         ANS = XARG / YARG
C UNARY MINUS
      ELSE IF (CODE .EQ. 5) THEN
         IF (NOX) GO TO 820
         ANS = -XARG
C POWER
      ELSE IF (CODE .EQ. 6) THEN
         IF (NOX .OR. NOY) GO TO 820
         ANS = ABS(XARG) ** YARG
C
C STANDARD FUNCTIONS
C
C VALUE OF PI
      ELSE IF (CODE .EQ. 7) THEN
         ANS = 3.14159265358D0
C NATURAL LOG
      ELSE IF (CODE .EQ. 8) THEN
         IF (NOX) GO TO 820
         ANS = LOG(XARG)
         IF (XARG .EQ. 1.0) ANS = 0.0
C EXPONENTIATE
      ELSE IF (CODE .EQ. 9) THEN
         IF (NOX) GO TO 820
         ANS = EXP(XARG)
C MODULUS
      ELSE IF (CODE .EQ. 10) THEN
         IF (NOX .OR. NOY) GO TO 820
         ANS = MOD(XARG,YARG)
C TRUNCATE TO INTEGER
      ELSE IF (CODE .EQ. 11) THEN
         IF (NOX) GO TO 820
         ANS = AINT(XARG)
C UNIFORM RANDOM NUMBER IN THE CLOSED RANGE 0,1
      ELSE IF (CODE .EQ. 12) THEN
         IF (NOX) GO TO 820
         ANS = RANU()
C SQUARE ROOT
      ELSE IF (CODE .EQ. 13) THEN
         IF (NOX) GO TO 820
         ANS = SQRT(XARG)
C COMMON LOG
      ELSE IF (CODE .EQ. 14) THEN
         IF (NOX) GO TO 820
         ANS = LOG10(XARG)
C
C TRIGONOMETRIC FUNCTIONS
C
C SINE
      ELSE IF (CODE .EQ. 15) THEN
         IF (NOX) GO TO 820
         ANS = SIN(XARG)
C COSINE
      ELSE IF (CODE .EQ. 16) THEN
         IF (NOX) GO TO 820
         ANS = COS(XARG)
C TANGENT
      ELSE IF (CODE .EQ. 17) THEN
         IF (NOX) GO TO 820
         ANS = TAN(XARG)
C ARC-SINE
      ELSE IF (CODE .EQ. 18) THEN
         IF (NOX) GO TO 820
         ANS = ASIN(XARG)
C ARC-COSINE
      ELSE IF (CODE .EQ. 19) THEN
         IF (NOX) GO TO 820
         ANS = ACOS(XARG)
C ARC-TANGENT
      ELSE IF (CODE .EQ. 20) THEN
         IF (NOX) GO TO 820
         ANS = ATAN(XARG)
C
C BOOLEAN OPERATORS
C
C EQUALITY
      ELSE IF (CODE .EQ. 51) THEN
         ANS = FALSE
         IF (XARG .EQ. YARG) ANS = TRUE
C LESS THAN
      ELSE IF (CODE .EQ. 52) THEN
         ANS = FALSE
         IF (XARG .LT. YARG) ANS = TRUE
C GREATER THAN
      ELSE IF (CODE .EQ. 53) THEN
         ANS = FALSE
         IF (XARG .GT. YARG) ANS = TRUE
C LESS THAN OR EQUAL TO
      ELSE IF (CODE .EQ. 54) THEN
         ANS = FALSE
         IF (XARG .LE. YARG) ANS = TRUE
C GREATER THAN OR EQUAL TO
      ELSE IF (CODE .EQ. 55) THEN
         ANS = FALSE
         IF (XARG .GE. YARG) ANS = TRUE
C NOT EQUAL TO
      ELSE IF (CODE .EQ. 56) THEN
         ANS = FALSE
         IF (XARG .NE. YARG) ANS = TRUE
C LOGICAL AND
      ELSE IF (CODE .EQ. 57) THEN
         ANS = FALSE
         IF (XARG.NE.FALSE .AND. YARG.NE.FALSE) ANS = TRUE
         IF (NOX .OR. NOY) ANS = FALSE
C LOGICAL OR
      ELSE IF (CODE .EQ. 58) THEN
         ANS = FALSE
         IF (XARG.NE.FALSE .AND. (.NOT. NOX)) ANS = TRUE
         IF (YARG.NE.FALSE .AND. (.NOT. NOY)) ANS = TRUE
C BLANK CHECK FUNCTION
      ELSE IF (CODE .EQ. 59) THEN
         ANS = FALSE
         IF (NOX) ANS = TRUE
C
C VARIABLES
C
C VALUE OF VARIABLE INDEXED BY ENTRY NUMBER (DYNAMIC REFERENCE)
      ELSE IF (CODE .EQ. 99) THEN
         IF (NOX) GO TO 820
         IVAR = XARG + .5
         CALL GETVAR(IVAR, X,BLANK,ITYPE,ANS, ISTAT)
         IF (ISTAT .NE. 0) GO TO 820
C VALUE OF VARIABLE INDEXED BY NAME
      ELSE IF (CODE .GE. 100) THEN
         IVAR = CODE - 100
         CALL GETVAR(IVAR, X,BLANK,ITYPE,ANS, ISTAT)
         IF (ISTAT .NE. 0) GO TO 820

      ELSE
         WRITE (*,*) 'CODE =',CODE
         call intpr('DOOP: UNRECOGNIZED INTERNAL CODE',-1,0,0)
         STOP
      END IF
      GO TO 900
C
C UNSUCCESSFUL RETURN
C
C    ARITHMETIC ERROR:
810   ISTAT = -1
      DOOP = 0.0
      GO TO 999

C    MISSING VALUE:
820   ISTAT = 0
      DOOP = XBLANK
      GO TO 999
C
C SUCCESSFUL RETURN
C
900   ISTAT = 0
      DOOP = ANS

999   RETURN
      END

      BLOCK DATA INITTX
C---
C--- INITIALIZE THE COMMON BLOCK WHICH CONTAINS DATA STRUCTURES
C--- FOR TRANSFORMATIONS
C---

#include "tx.h"

      DATA  ITOKEN   / 200*999 /
      DATA  NTOKEN   / 0 /
      DATA  NCONST   / 0 /
      DATA  NOPTOR   / 0 /

      END
