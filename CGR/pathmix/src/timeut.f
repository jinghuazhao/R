C     +---------------------------------------------------------------+
C     !         ULTRIX VERSION OF TIME AND DATE UTILITIES             !
C     +---------------------------------------------------------------+
C
      SUBROUTINE TIMER(ELAPSE)
C---
C--- RETURNS CPU-TIME (IN MINUTES) ELAPSED SINCE THE LAST CALL
C---
C--- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C--- !! NOTE: IMPLEMENTATION OF THIS SUBROUTINE IS MACHINE DEPENDENT !!
C--- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C---
      DOUBLE PRECISION  ELAPSE
      REAL              TARRAY(2), result

      CALL DTIME( TARRAY, result )
      ELAPSE = TARRAY(1) / 60.0
      END
      SUBROUTINE GETDT(STRING)
C---
C--- RETURNS THE DATE AND TIME IN CHARACTER STRING FORM
C---
C--- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C--- !! NOTE: IMPLEMENTATION OF THIS SUBROUTINE IS MACHINE DEPENDENT !!
C--- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C---
      CHARACTER   STRING*(*)
      INTEGER     SCLOCK, T(9), TIME
      CHARACTER   STRTMP*18, MNTH(0:11)*3
      DATA MNTH /'JAN','FEB','MAR','APR','MAY','JUN',
     +           'JUL','AUG','SEP','OCT','NOV','DEC'/

      SCLOCK = TIME()
      CALL LTIME( SCLOCK, T )

      WRITE (STRTMP,1000) T(4),MNTH( T(5) ),T(6), T(3),T(2),T(1)
 1000 FORMAT (I2.2,'-',A3,'-',I2.2, ' ',I2,':',I2.2,':',I2.2)

      STRING = STRTMP
      END
