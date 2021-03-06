      SUBROUTINE FILES( LFNJOB, LFNDAT, LFNTER, LFNPLX,
     +                  FILNAM1,FILNAM2,FILNAM3,FILNAM4 )
C---
C--- INITIALIZE ALL FILES USED IN THIS PROGRAM
C---
      CHARACTER   FILNAM1*60,FILNAM2*60,FILNAM3*60,FILNAM4*60
      LOGICAL     NAMED

c     WRITE (*,*)
c     WRITE (*,*) 'ENTER NAME OF INPUT JOB FILE '
c     READ (*,'(A)') FILNAM1
      OPEN (UNIT=LFNJOB, FILE=FILNAM1, STATUS='OLD')

c     WRITE (*,*)
c     WRITE (*,*) 'ENTER NAME OF INPUT FAMILY DATA FILE '
c     READ (*,'(A)') FILNAM2
      OPEN (UNIT=LFNDAT, FILE=FILNAM2, STATUS='OLD')

c     WRITE (*,*)
c     WRITE (*,*) 'ENTER NAME OF TERSE (SUMMARY) OUTPUT FILE '
c     READ (*,'(A)') FILNAM3
      OPEN (UNIT=LFNTER, FILE=FILNAM3, STATUS='UNKNOWN')

c     WRITE (*,*)
c     WRITE (*,*) 'ENTER NAME OF PROLIX OUTPUT FILE (BLANK=STD OUT) '
c     READ (*,'(A)') FILNAM4
C     IF FILE NAME IS BLANK, USE STANDARD OUTPUT
      IF (FILNAM4 .EQ. ' ') THEN
C ======================================================================
C NOTE:  THE FOLLOWING LINE SETS THE LOGICAL FILE UNIT OF PROLIX TO THE
C        "STANDARD OUTPUT" UNIT.  IN MOST IMPLEMENTATIONS OF FORTRAN 77,
C        UNIT NUMBER SIX IS CUSTOMARILY USED.  HOWEVER, THIS NUMBER IS
C        MACHINE SPECIFIC AND COULD POSSIBLY BE DIFFERENT ON OTHER
C        MACHINES.  (STANDARD OUTPUT REPRESENTS THE FILE WRITTEN TO IN
C        THE "WRITE (*,*)" STATEMENT.)
C ======================================================================
         LFNPLX = 6
      ELSE
         OPEN (UNIT=LFNPLX, FILE=FILNAM4, STATUS='UNKNOWN')
      END IF

C WRITE TERSE AND PROLIX FILE HEADINGS

      CALL EJECT( LFNTER )
      CALL EJECT( LFNPLX )

C     JOB FILE
      INQUIRE (UNIT=LFNJOB,NAMED=NAMED,NAME=FILNAM1)
      IF (NAMED) WRITE (LFNTER,1010) FILNAM1(1:LENSTR(FILNAM1) )

C     DATA FILE
      INQUIRE (UNIT=LFNDAT,NAMED=NAMED,NAME=FILNAM2)
      IF (NAMED) WRITE (LFNTER,1020) FILNAM2(1:LENSTR(FILNAM2) )

C     TERSE FILE
      INQUIRE (UNIT=LFNTER,NAMED=NAMED,NAME=FILNAM3)
      IF (NAMED) WRITE (LFNTER,1030) FILNAM3(1:LENSTR(FILNAM3) )

C     PROLIX FILE
      INQUIRE (UNIT=LFNPLX,NAMED=NAMED,NAME=FILNAM4)
      IF (NAMED) WRITE (LFNTER,1040) FILNAM4(1:LENSTR(FILNAM4) )

 1010 FORMAT(/' ','JOB FILE:           ',A)
 1020 FORMAT (' ','DATA FILE:          ',A)
 1030 FORMAT (' ','TERSE OUTPUT FILE:  ',A)
 1040 FORMAT (' ','PROLIX OUTPUT FILE: ',A)
      call flush(lfnter)
      call flush(lfnplx)
      END
