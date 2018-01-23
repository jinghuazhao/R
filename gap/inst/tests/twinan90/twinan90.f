C------------------------------------------------------------------C
C     NEWTWIN4.FOR                                                 C
C     BY CHRIS WILLIAMS,   January 1993                            C
C------------------------------------------------------------------C
C     NEWTWIN.FOR PERFORMS A TWIN ANALYSIS USING THE METHODS       C
C     OF CHRISTIAN, KANG AND NORTON (1974) AM J HUM GENETICS       C
C     AND ALSO PERFORMS A SIMPLE PATH MODEL ANALYSIS.              C
C------------------------------------------------------------------C
      subroutine newtw5(d1,d2,nmz,ndz,mp,xlamb,const,vmiss,path,ped,
     &names,nchar,nvar,nf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION D1(nmz,2),D2(ndz,2),DAT1(nmz,2),DAT2(nmz,2),N(2)
      DIMENSION DIFF1(nmz),DIFF2(ndz),SUM1(nmz),SUM2(ndz)
      DIMENSION ADIFF1(nmz,2),ADIFF2(ndz,2),SDAT1(nmz*2),SDAT2(ndz*2)
      DIMENSION VN(2),TCOV(2,2,2),BETA1(3,1),SOL1(3,1)
      DIMENSION A1(3,3),STEP(3),AIN1(3,3)
      DIMENSION BETA2(3,1),BETA3(3,1),SOL2(3,1),SOL3(3,1)
      DIMENSION A2(3,3),A3(3,3)
      DIMENSION AIN2(3,3),AIN3(3,3),SOL4(3,1),AIN4(3,3)
      CHARACTER*2 PKSMZ,PKSDZ,FF
      CHARACTER*16 NAMA,NAMC,NAMD,NAMP
cc      CHARACTER*80 FORM
      character*255 form2
      CHARACTER*55 NAMB
      CHARACTER*65 RES1,RES2,RES3,RESULT
      character*255 names,form1
      DATA FF/'12'/
      RES1 = 'THE EQUAL VARIANCE HYPOTHESIS IS NOT REJECTED.'
      RES2 =
     &  'THE EQUAL VAR HYPOTHESIS IS REJECTED, WITH DZ VAR > MZ VAR.'
      RES3 =
     &  'THE EQUAL VAR HYPOTHESIS IS REJECTED, WITH MZ VAR > DZ VAR.'
C-----------------------------------------------------------C
C     INPUT THE DATA FILES                                  C
C-----------------------------------------------------------C
10    CONTINUE
*     WRITE(*,*) 'TWINAN90 TWIN ANALYSIS PROGRAM, VERSION 1.0'
c     WRITE(*,*) 'ENTER THE NAME OF THE DATA FILE.'
c     READ(*,50) NAMA
      nama=names((4*nchar+13):(5*nchar+16))
cc50    FORMAT(A16)
c     WRITE(*,*) 'ENTER THE NAME OF THE VARIABLE.'
c     READ(*,60) NAMB
      namb=names(1:nchar)
cc60    FORMAT(A55)
c     WRITE(*,*) 'ENTER THE NAME OF THE OUTPUT FILE.'
c     READ(*,50) NAMC
      namc=names((nchar*2+5):(3*nchar+8))
c     WRITE(*,*) 'ENTER THE NAME OF THE DIAGNOSTIC OUTPUT FILE.'
c     READ(*,50) NAMD
      namd=names((nchar+1):(2*nchar+4))
      namp=names((3*nchar+9):(4*nchar+12))
      OPEN(UNIT=2,FILE=NAMD,STATUS='NEW')
*     OPEN(UNIT=3,FILE=NAMA,STATUS='OLD')
c     WRITE(*,*) 'ENTER THE INPUT FORMAT (A FORTRAN FORMAT).'
c     READ(*,70) FORM
cc70    FORMAT(A80)
      form1=names((5*nchar+17):(5*nchar+16+nf))
      form2=form1
c     WRITE(*,*) '   ENTER THE POWER TRANSFORMATION PARAMETER.'
c     WRITE(*,*) 'ENTER 1. IF YOU DO NOT WISH TO TRANSFORM YOUR DATA.'
c     READ(*,*)  XLAMB
c     WRITE(*,*) 'ENTER A VALUE TO ADD TO EACH DATA POINT.'
c     READ(*,*) CONST
c     WRITE(*,*) 'ENTER THE VALUE FOR MISSING DATA.'
c     READ(*,*) VMISS
c     WRITE(*,*) '   ENTER 1. IF YOU WISH TO OBTAIN'
c     WRITE(*,*) 'PATH MODEL RESULTS, OTHERWISE ENTER 0.'
c     READ(*,*) PATH
c     WRITE(*,*) '    ENTER 1. IF YOU WISH TO PRODUCE'
c     WRITE(*,*) 'AN OUTPUT PEDIGREE FILE, OTHERWISE ENTER 0.'
c     READ(*,*) PED
      INEG = 0
      KT1 = 0
      KT2 = 0
c     READ(3,*) X1,X2
      x1=nmz
      x2=ndz
      N(1) = DINT(X1)
      N(2) = DINT(X2)
      WRITE(2,*) 'READING MZ DATA '
      WRITE(2,*) ' TWIN PAIR,    ORIGINAL DATA,    TRANSFORMED DATA'
      DO 100 I = 1,N(1)
        DAT1(I,1) = 0.
        DAT1(I,2) = 0.
c       READ(3,FORM,END=105) D1(I,1),D1(I,2)
c       write(*,'(2hMZ,2F4.0)') d1(i,1), d1(i,2)
        IF (D1(I,1) .EQ. VMISS .OR. D1(I,2) .EQ. VMISS) GO TO 100
        KT1 = KT1 + 1
        IF (XLAMB .EQ. 1.) THEN
          DAT1(KT1,1) = D1(I,1) +CONST
          DAT1(KT1,2) = D1(I,2) +CONST
        ELSE IF (XLAMB .EQ. 0.) THEN
          DAT1(KT1,1) = DLOG(D1(I,1)+CONST)
          DAT1(KT1,2) = DLOG(D1(I,2)+CONST)
        ELSE
          DAT1(KT1,1) = ((D1(I,1)+CONST)**XLAMB -1.)/XLAMB
          DAT1(KT1,2) = ((D1(I,2)+CONST)**XLAMB -1.)/XLAMB
        END IF
        IF (DAT1(KT1,1) .LT. 0.0 .OR. DAT1(KT1,2) .LT. 0.0)
     &   INEG = 1
        WRITE(2,110) KT1,D1(I,1),D1(I,2),DAT1(KT1,1),DAT1(KT1,2)
        DIFF1(KT1) = DAT1(KT1,1) -DAT1(KT1,2)
        SUM1(KT1) = DAT1(KT1,1) +DAT1(KT1,2)
        ADIFF1(KT1,1) = ABS(DIFF1(KT1))
        ADIFF1(KT1,2) = ADIFF1(KT1,1)
100   CONTINUE
cc105   CONTINUE
110   FORMAT(3X,I5,3X,4(g9.3,1X))
      WRITE(2,*) 'READING DZ DATA '
      WRITE(2,*) ' TWIN PAIR,    ORIGINAL DATA,    TRANSFORMED DATA'
      DO 120 I = 1,N(2)
        DAT2(I,1) = 0.
        DAT2(I,2) = 0.
c       READ(3,FORM,END=125) D2(I,1),D2(I,2)
c       write(*,'(2HDZ,2F4.0)') d2(i,1), d2(i,2)
        IF (D2(I,1) .EQ. VMISS .OR. D2(I,2) .EQ. VMISS) GO TO 120
        KT2 = KT2 + 1
        IF (XLAMB .EQ. 1.) THEN
          DAT2(KT2,1) = D2(I,1) +CONST
          DAT2(KT2,2) = D2(I,2) +CONST
        ELSE IF (XLAMB .EQ. 0.) THEN
          DAT2(KT2,1) = DLOG(D2(I,1)+CONST)
          DAT2(KT2,2) = DLOG(D2(I,2)+CONST)
        ELSE
          DAT2(KT2,1) = ((D2(I,1)+CONST)**XLAMB -1.)/XLAMB
          DAT2(KT2,2) = ((D2(I,2)+CONST)**XLAMB -1.)/XLAMB
        END IF
        IF (DAT2(KT2,1) .LT. 0.0 .OR. DAT2(KT2,2) .LT. 0.0)
     &   INEG = 1
        WRITE(2,110) KT2,D2(I,1),D2(I,2),DAT2(KT2,1),DAT2(KT2,2)
        DIFF2(KT2) = DAT2(KT2,1) -DAT2(KT2,2)
        SUM2(KT2) = DAT2(KT2,1) +DAT2(KT2,2)
        ADIFF2(KT2,1) = ABS(DIFF2(KT2))
        ADIFF2(KT2,2) = ADIFF2(KT2,1)
120    CONTINUE
cc125    CONTINUE
*     WRITE(*,*) 'THE SAMPLE SIZES ARE ',KT1,' AND ',KT2
      N(1) = KT1
      N(2) = KT2
      X1 = DBLE(KT1)
      X2 = DBLE(KT2)
*     CLOSE(UNIT=3)
      IF (PED .EQ. 1.) CALL FOUT(NAMA,VMISS,NAMP,nvar,form2)
C-----------------------------------------------------------C
C     COMPUTE THE VALUES OF DESCRIPTIVE STATISTICS          C
C     AND SUMS OF SQUARES BY CALLING SUBROUTINE DESCRIP,    C
C     ADAPTED FROM A SUBROUTINE FROM NUMERICAL RECIPES      C
C-----------------------------------------------------------C
      CALL DESCRIP(DAT1,N(1),AVE1,VAR1,VMZ1,VMZ2,SK1,
     &  CRT1,SSAMZ,SSWMZ,mp)
C     WRITE(2,*) 'AVE1,VAR1,VMZ1,VMZ2,SK1,CRT1,SSAMZ,SSWMZ'
C     WRITE(2,130) AVE1,VAR1,VMZ1,VMZ2,SK1,CRT1,SSAMZ,SSWMZ
cc130   FORMAT(2X,8(g8.3,1X))
      CALL DESCRIP(DAT2,N(2),AVE2,VAR2,VDZ1,VDZ2,SK2,
     &  CRT2,SSADZ,SSWDZ,mp)
C     WRITE(2,130) AVE2,VAR2,VDZ1,VDZ2,SK2,CRT2,SSADZ,SSWDZ
      CALL DESCRIP(ADIFF1,N(1),AVE3,VAR3,VJ1,VJ2,SK3,
     &  CRT3,SADF1,SWDF1,mp)
C     WRITE(2,130) AVE3,VAR3,VJ1,VJ2,SK3,CRT3,SADF1,SWDF1
      CALL DESCRIP(ADIFF2,N(2),AVE4,VAR4,VJ3,VJ4,SK4,
     &  CRT4,SADF2,SWDF2,mp)
C     WRITE(2,130) AVE4,VAR4,VJ3,VJ4,SK4,CRT4,SADF2,SWDF2
      C1 = (2.*X1 -1.)/(2.*X1 -2.)
      VAR3 = VAR3/C1
      SK3 = C1**1.5*SK3
      CRT3 = C1**2*CRT3
      C2 = (2.*X2 -1.)/(2.*X2 -2.)
      VAR4 = VAR4/C2
      SK4 = C2**1.5*SK4
      CRT4 = C2**2*CRT4
      SESKMZ = DSQRT(6./(2*X1))
      SESKDZ = DSQRT(6./(2*X2))
      SEKUMZ = DSQRT(24./(2*X1))
      SEKUDZ = DSQRT(24./(2*X2))
      AMSMZ = SSAMZ/(X1 -1)
      WMSMZ = SSWMZ/X1
      AMSDZ = SSADZ/(X2 -1)
      WMSDZ = SSWDZ/X2
      TSSMZ = SSAMZ + SSWMZ
      TSSDZ = SSADZ + SSWDZ
      DIFMN = AVE1 -AVE2
      DIFTST = DIFMN/DSQRT( AMSMZ/(2*X1) +AMSDZ/(2*X2) )
        DFTDF = ( AMSMZ/(2*X1) +AMSDZ/(2*X2) )**2/
     &  ( (AMSMZ/(2*X1))**2/(X1-1) +
     &    (AMSDZ/(2*X2))**2/(X2-1) )
        DFTARG = DFTDF/(DFTDF +DIFTST**2)
        PDFT = BETAI(DFTDF/2.D0,.5D0,DFTARG)
      CALL CONVKS(DAT1,SDAT1,N(1),AVE1,VAR1,DMZ,PRDMZ,mp)
       PKSMZ = '= '
       IF (PRDMZ .EQ. .15) PKSMZ = '> '
       IF (PRDMZ .EQ. .01) PKSMZ = '< '
       WRITE(2,*) 'BEGINNING MZ BCOX ANALYSIS'
        VLAMB1 = 1.
        IF (PRDMZ .LT. .05 .AND. XLAMB .EQ. 1. .AND. INEG .NE. 1)
     &     VLAMB1 = BCOX(DAT1,N(1),mp)
       WRITE(2,150) VLAMB1
      CALL CONVKS(DAT2,SDAT2,N(2),AVE2,VAR2,DDZ,PRDDZ,mp)
       PKSDZ = '= '
       IF (PRDDZ .EQ. .15) PKSDZ = '> '
       IF (PRDDZ .EQ. .01) PKSDZ = '< '
       WRITE(2,*) 'BEGINNING DZ BCOX ANALYSIS'
        VLAMB2 = 1.
        IF (PRDDZ .LT. .05 .AND. XLAMB .EQ. 1. .AND. INEG .NE. 1)
     &     VLAMB2 = BCOX(DAT2,N(2),mp)
        VLAMB = (VLAMB1 + VLAMB2)/2.
       WRITE(2,152) VLAMB2
       WRITE(2,154) VLAMB
150    FORMAT(3X,'MZ LAMBDA = ',F8.4)
152    FORMAT(3X,'DZ LAMBDA = ',F8.4)
154    FORMAT(3X,'AVERAGE LAMBDA = ',F8.4)
C-----------------------------------------------------------C
C     CONDUCTING TWIN ANALYSIS                              C
C-----------------------------------------------------------C
      RATIO = (AMSDZ +WMSDZ)/(AMSMZ +WMSMZ)
        DFMZ = (AMSMZ +WMSMZ)**2/(AMSMZ**2/(X1-1) +WMSMZ**2/X1)
        DFDZ = (AMSDZ +WMSDZ)**2/(AMSDZ**2/(X2-1) +WMSDZ**2/X2)
        RATARG = DFMZ/(DFMZ +DFDZ*RATIO)
        PRAT = BETAI(DFMZ/2.,DFDZ/2.,RATARG)
        IF (PRAT .LT. .5) THEN
          PRAT = 2.*PRAT
        ELSE IF (PRAT .GE. .5) THEN
          PRAT = 2.*(1. - PRAT)
        END IF
C      WRITE(2,*) 'PRAT = ',PRAT
        IF (PRAT .GT. .2) THEN
          RESULT = RES1
        ELSE IF (RATIO .GT. 1.) THEN
          RESULT = RES2
        ELSE
          RESULT = RES3
        ENDIF
      GWP = 2*(WMSDZ - WMSMZ)
        VARGWP = 8*( WMSDZ**2/(X2+2) +WMSMZ**2/(X1+2) )
        FGWP = WMSDZ/WMSMZ
        GWPARG = X1/(X1 +X2*FGWP)
        PGWP = BETAI(X1/2.,X2/2.,GWPARG)
      GAP = 2*(AMSMZ - AMSDZ)
        VARGAP = 8*( AMSDZ**2/(X2+1) +AMSMZ**2/(X1+1) )
      GAC = (AMSMZ -AMSDZ +WMSDZ -WMSMZ)
        VARGAC = (VARGWP +VARGAP)
        FGAC = (AMSMZ +WMSDZ)/(AMSDZ +WMSMZ)
        DFNAC = (AMSMZ +WMSDZ)**2/(AMSMZ**2/(X1-1) +WMSDZ**2/X2)
        DFDAC = (AMSDZ +WMSMZ)**2/(AMSDZ**2/(X2-1) +WMSMZ**2/X1)
        GACARG = DFDAC/(DFDAC +DFNAC*FGAC)
        PGAC = BETAI(DFDAC/2.,DFNAC/2.,GACARG)
      GMIN = (VARGWP*GAP +VARGAP*GWP)/(VARGWP +VARGAP)
        VARGM = (VARGAP*VARGWP**2 +VARGWP*VARGAP**2)/
     &          (VARGAP+VARGWP)**2
      VICMZ = (AMSMZ -WMSMZ)/(AMSMZ +WMSMZ)
        FICMZ = AMSMZ/WMSMZ
        ARGMIC = X1/(X1 +(X1-1.)*FICMZ)
        DFMIC = (X1 -1.)/2.
        PZMZ = BETAI(X1/2.,DFMIC,ARGMIC)
      VICDZ = (AMSDZ -WMSDZ)/(AMSDZ +WMSDZ)
        FICDZ = AMSDZ/WMSDZ
        ARGDIC = X2/(X2 +(X2-1.)*FICDZ)
        DFDIC = (X2 -1.)/2.
        PZDZ = BETAI(X2/2.,DFDIC,ARGDIC)
      VHAC = 4*GAC/(AMSMZ +AMSDZ +WMSMZ +WMSDZ)
      VHWP = 4*GWP/(AMSMZ +AMSDZ +WMSMZ +WMSDZ)
      VHIC = 2*(VICMZ -VICDZ)
        ZICMZ = .5*DLOG( (1.+VICMZ)/(1.-VICMZ) ) +.5*DLOG(X1/(X1-1))
        ZICDZ = .5*DLOG( (1.+VICDZ)/(1.-VICDZ) ) +.5*DLOG(X2/(X2-1))
        VZMZ = 1./(X1 - 1.5)
        VZDZ = 1./(X2 - 1.5)
        ZHIC = (ZICMZ -ZICDZ)/DSQRT(VZMZ +VZDZ)
        PZHIC = .5*CHI123(ZHIC**2,1)
C     WRITE(2,*) 'ZICMZ,ZICDZ,VZDZ,VZMZ,ZHIC '
C     WRITE(2,160) ZICMZ,ZICDZ,VZDZ,VZMZ,ZHIC
cc160   FORMAT(2X,5(g10.6,2X))
      VDTMZ = .72676046*WMSMZ/X1
      VDTDZ = .72676046*WMSDZ/X2
      ADT2 = (AVE4 -AVE3)/DSQRT(VDTMZ +VDTDZ)
        DFADT2 = (VDTMZ +VDTDZ)**2/(VDTMZ**2/X1 +VDTDZ**2/X2)
        AD2ARG = DFADT2/(DFADT2 +ADT2**2)
        PADT2 = .5*BETAI(DFADT2/2.,.5D0,AD2ARG)
C-----------------------------------------------------------C
C      PATH MODEL ANALYSIS                                  C
C-----------------------------------------------------------C
      IF (PATH .NE. 1.) GO TO 190
      VN(1) = X1
      VN(2) = X2
      TCOV(1,1,1) = VMZ1
      TCOV(1,1,2) = VAR1*VICMZ
      TCOV(1,2,1) = TCOV(1,1,2)
      TCOV(1,2,2) = VMZ2
      TCOV(2,1,1) = VDZ1
      TCOV(2,1,2) = VAR2*VICDZ
      TCOV(2,2,1) = TCOV(2,1,2)
      TCOV(2,2,2) = VDZ2
      CALL CE(TCOV,VN,SOL4,AIN4,VLCE)
      VAVG = (X1*VAR1 +X2*VAR2)/(X1 +X2)
      INO = 3
      MONIT = 2
      STEP(1) = 1.
      STEP(2) = 1.
      STEP(3) = 1.
      DO 185 I = 1,3
      DO 180 J = 1,3
      A1(I,J) = 0.
      A2(I,J) = 0.
      A3(I,J) = 0.
      AIN3(I,J) = 0.
180   CONTINUE
      SOL1(I,1) = 0.
      SOL2(I,1) = 0.
      SOL3(I,1) = 0.
      BETA1(I,1) = VAVG/3.
      BETA2(I,1) = VAVG/3.
      BETA3(I,1) = VAVG/2.
185   CONTINUE
      DETMZ = TCOV(1,1,1)*TCOV(1,2,2) -TCOV(1,1,2)**2.
      DETDZ = TCOV(2,1,1)*TCOV(2,2,2) -TCOV(2,1,2)**2.
      VARE = ((X1-1.)/(X1+X2-2.))*.5*(TCOV(1,1,1)+TCOV(1,2,2))
     &      +((X2-1.)/(X1+X2-2.))*.5*(TCOV(2,1,1)+TCOV(2,2,2))
      VLE = -.5*(X1-1.)*(2.*DLOG(VARE) +(TCOV(1,1,1)+TCOV(1,2,2))/VARE)
     &  -.5*(X2-1.)*(2.*DLOG(VARE) +(TCOV(2,1,1)+TCOV(2,2,2))/VARE)
      VL0 = -.5*(X1-1.)*(DLOG(DETMZ) +2.)
     &  -.5*(X2-1.)*(DLOG(DETDZ) +2.)
      VL1 = 0.
      VL2 = 0.
      VL3 = 0.
      IMOD = 1
      CALL QHOUSE(VN,IMOD,TCOV,BETA1,SOL1,A1,VL1,INO,MONIT,STEP)
      GOF1 = 2*(VL0 -VL1)
      GENVAR1 = SOL1(1,1) +SOL1(2,1)
      PHERIT1 = GENVAR1/(GENVAR1 +SOL1(3,1))
      ERR = 0.
      CALL SINV3(A1,AIN1,ERR)
      AVGENVAR1 = AIN1(1,1) +AIN1(2,2) +2*AIN1(1,2)
      ACOVGE1 = AIN1(1,3) +AIN1(2,3)
      SEGEN1 = DSQRT(AVGENVAR1)
      SEERR1 = DSQRT(AIN1(3,3))
      IMOD = 2
      CALL QHOUSE(VN,IMOD,TCOV,BETA2,SOL2,A2,VL2,INO,MONIT,STEP)
      GOF2 = 2*(VL0 -VL2)
      ERR = 0.
      CALL SINV3(A2,AIN2,ERR)
      IMOD = 3
      CALL QHOUSE(VN,IMOD,TCOV,BETA3,SOL3,A3,VL3,INO,MONIT,STEP)
      GOF3 = 2*(VL0 -VL3)
      GOF4 = 2*(VL3 -VLE)
C      GOFE = 2*(VL0 -VLE)
      GOFCE = 2*(VL0 -VLCE)
      GOFC = 2*(VLCE -VLE)
*     WRITE(*,*) 'VARE, VLE, VL3, VL1, VL2, AND VL0 ARE'
*     WRITE(*,*) VARE,VLE,VL3,VL1,VL2,VL0
      DET3 = ( A3(1,1)*A3(3,3) -A3(1,3)**2 )
      AIN3(1,1) = A3(3,3)/DET3
      AIN3(1,3) = -A3(1,3)/DET3
      AIN3(3,3) = A3(1,1)/DET3
      AIN3(3,1) = AIN3(1,3)
C-----------------------------------------------------------C
C     OUTPUT RESULTS                                        C
C-----------------------------------------------------------C
190   CONTINUE
      OPEN(UNIT=3,FILE=NAMC,STATUS='NEW')
C***********************************************************C
      WRITE(3,200)
200   FORMAT(78('_')/)
      WRITE(3,*) 'TWINAN90 TWIN ANALYSIS PROGRAM, VERSION 1.0'
      WRITE(3,250) NAMB
250   FORMAT(21H   FOR THE VARIABLE: ,A55)
      IF (XLAMB .NE. 1.) WRITE(3,260) XLAMB
260   FORMAT('   WITH THE POWER TRANSFORMATION LAMBDA = ',F4.2)
      WRITE(3,270) N(1),N(2)
270   FORMAT('    NUMBER OF TWIN PAIRS: MZ = ',I5,'    DZ = ',I5)
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'ASSESSMENT OF MODEL ASSUMPTIONS, ',
     & 'DESCRIPTIVE STATISTICS'
      WRITE(3,320)
320   FORMAT(/23X,'Individuals',16X,'Twin absolute differences')
      WRITE(3,340)
340   FORMAT(20X,'MZ',14X,'DZ',15X,'MZ',14X,'DZ')
      WRITE(3,360) AVE1,AVE2,AVE3,AVE4
360   FORMAT(1X,'Mean',10X,g11.4,5X,g11.4,4X,2X,2(g11.4,5X))
      WRITE(3,380) VAR1,VAR2,VAR3,VAR4
380   FORMAT(1X,'Variance',6X,g11.4,5X,g11.4,4X,2X,2(g11.4,5X))
      WRITE(3,400) SK1,SESKMZ,SK2,SESKDZ
400   FORMAT(1X,'Skew;SE',6X,2(f6.2,';',f5.2,4X))
C     &       ,1X,2(f6.2,';',f5.2,4X))
      WRITE(3,420) CRT1,SEKUMZ,CRT2,SEKUDZ
420   FORMAT(1X,'Ex Kurt;SE',3X,2(f6.2,';',f5.2,4X))
C     &       ,1X,2(f6.2,';',f5.2,4X))
      WRITE(3,*)
      WRITE(3,*)'Kolmogorov-Smirnov test of normality (valid for N > ',
     &'25 twin pairs)'
      WRITE(3,440) DMZ,PKSMZ,PRDMZ,DDZ,PKSDZ,PRDDZ
440   FORMAT(2X,'MZ: D" = ',F5.2,'; P ',A2,F5.3,5X,
     &       'DZ: D" = ',F5.2,'; P ',A2,F5.3)
      IF (VLAMB .NE. 1.0) WRITE(3,450) VLAMB
450   FORMAT(4X,'A recommended power transformation is lambda = ',F5.2)
      WRITE(3,*)
      WRITE(3,*) 'Test for difference in MZ and DZ means :'
      WRITE(3,460) DIFMN,DIFTST,DFTDF,PDFT
460   FORMAT('  MZ - DZ = ',g12.4,2X,' T" = ',f6.2,2X,
     &  ' DF = ',F6.1,2X,' P = ',F5.3)
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'ANALYSIS OF VARIANCE'
      WRITE(3,480)
480   FORMAT(37X,'MZ',10X,'DF',13X,'DZ',10X,'DF')
      WRITE(3,500) AMSMZ,N(1)-1,AMSDZ,N(2)-1
500   FORMAT(' Among mean squares   :',1X,2(3X,g17.4,3X,I4))
      WRITE(3,550) WMSMZ,N(1),WMSDZ,N(2)
550   FORMAT(' Within mean squares  :',1X,2(3X,g17.4,3X,I4))
      WRITE(3,580)
580   FORMAT(34X,10('-'),17X,10('-'))
      WRITE(3,600) AMSMZ+WMSMZ,DFMZ,AMSDZ+WMSDZ,DFDZ
600   FORMAT(' Sum of mean squares  :',3X,2(1X,g17.4,2X,F7.1))
C      WRITE(3,650) TSSMZ,2*N(1)-1,TSSDZ,2*N(2)-1
cc650   FORMAT(' Total sum of squares :',1X,2(3X,g17.2,3X,I4))
      WRITE(3,*)
      WRITE(3,680)
680   FORMAT(14('- '),'TEST OF TWIN MODEL:',1X,15('- '))
      WRITE(3,700) RATIO,PRAT
700   FORMAT(' Sum of mean squares ratio :'
     &       ,7X,'F" = ',F6.2,5X,'Approx P = ',F5.3)
      WRITE(3,710) RESULT
710   FORMAT(5X,A65)
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'ESTIMATES OF GENETIC VARIANCE'
      WRITE(3,750)
750   FORMAT(/27X,'Estimate',8X,'Variance',9X,'F-ratio',4X,'Prob')
      WRITE(3,800) GWP,VARGWP,FGWP,PGWP
800   FORMAT(1X,'Within pair     (WP):',
     &       1X,g14.6,2X,g14.6,6X,F7.2,4X,F5.3)
      WRITE(3,850) GAC,VARGAC,FGAC,PGAC
850   FORMAT(1X,'Among component (AC):',
     &       1X,g14.6,2X,g14.6,6X,F7.2,4X,F5.3)
C      WRITE(3,900) GMIN,VARGM
C900   FORMAT(/1X,'MINIMUM VARIANCE    :'1X,F14.6,2X,F19.6)
      WRITE(3,950) VICMZ,PZMZ,VICDZ,PZDZ
950   FORMAT(/1X,'Intraclass correlations:',2X,
     & 'MZ = ',F6.3,'; P = ',F5.3,4X,'DZ = ',F6.3,'; P = ',F5.3)
      WRITE(3,970) ADT2,DFADT2,PADT2
970   FORMAT(/1X,'Average abs difference test :  T" = ',F6.2,
     &  2X,' DF = ',F6.1,2X,' P = ',F5.3)
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'HERITABILITY ESTIMATES'
      WRITE(3,1000)
1000  FORMAT(/24X,'Formula',17X,'Estimate',4X,'Probability')
      WRITE(3,1050) VHWP
1050  FORMAT(1X,'Within pair    :  2*(WDZ-WMZ)/((SMZ+SDZ)/4)'
     &       ,4X,F8.3,4X,'SEE PROB(WP)')
      WRITE(3,1100) VHAC
1100  FORMAT(1X,'Among component:  (2*AC)/((SMZ+SDZ)/4)'
     &       ,9X,F8.3,4X,'SEE PROB(AC)')
      WRITE(3,1250) VHIC,PZHIC
1250  FORMAT(/1X,'Intraclass correlation: 2*(RMZ-RDZ)',12X,
     &        F8.3,8X,F5.3)
C***********************************************************C
      WRITE(3,200)
      IF (PATH .NE. 1.) GO TO 1650
      WRITE(3,1290) FF
1290  FORMAT(A2)
C***********************************************************C
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'TWINAN90 TWIN ANALYSIS PROGRAM, VERSION 1.0'
      WRITE(3,*) 'TWINAN90 PATH ANALYSIS OUTPUT, page 1'
      WRITE(3,250) NAMB
      IF (XLAMB .NE. 1.) WRITE(3,260) XLAMB
      WRITE(3,270) N(1),N(2)
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'FOR THE COVARIANCE MATRICES'
      WRITE(3,*)
      WRITE(3,1350) VMZ1,VAR1*VICMZ,VDZ1,VAR2*VICDZ,
     &              VAR1*VICMZ,VMZ2,VAR2*VICDZ,VDZ2
1350  FORMAT(3X,'MZ = ',2X,g11.4,5X,g11.4,3X,'DZ = ',2X,g11.4,
     &        5X,g11.4/10X,2(g11.4,5X,g11.4,10X))
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'ADE MODEL MAXIMUM LIKELIHOOD (ML) ESTIMATES'

      call sehml(sol1(1,1),sol1(2,1),sol1(3,1),DSQRT(AIN1(1,1)),
     &   DSQRT(AIN1(2,2)),DSQRT(AIN1(3,3)),
     & AIN1(1,2)/DSQRT(AIN1(1,1)*AIN1(2,2)),
     & AIN1(1,3)/DSQRT(AIN1(1,1)*AIN1(3,3)),
     & AIN1(2,3)/DSQRT(AIN1(3,3)*AIN1(2,2)),hadehat,sehade,
     &  dadehat,sedade,eadehat,seeade)
c     hadehat = sol1(1,1)/( SOL1(1,1)+SOL1(2,1)+SOL1(3,1) )
c     dadehat = sol1(2,1)/( SOL1(1,1)+SOL1(2,1)+SOL1(3,1) )
c     eadehat = sol1(3,1)/( SOL1(1,1)+SOL1(2,1)+SOL1(3,1) )
c     setest = .5
      write(3,1370)
1370  format(/32x,'Estimate',6x,'S. E.',6x,'Proportion',4x,'S. E.'/)
      write(3,1372) sol1(1,1),DSQRT(AIN1(1,1)),hadehat,sehade
1372  format(' Additive genetic variance',2x,g15.4,2x,g11.4,5x,f5.3,
     &         4x,f5.3/)
      write(3,1374) sol1(2,1),DSQRT(AIN1(2,2)),dadehat,sedade
1374  format(' Dominance genetic variance',1x,g15.4,2x,g11.4,5x,f5.3,
     &         4x,f5.3/)
      write(3,1376) sol1(3,1),DSQRT(AIN1(3,3)),eadehat,seeade
1376  format(' Error variance',13x,g15.4,2x,g11.4,5x,f5.3,
     &         4x,f5.3//)

c      WRITE(3,1370) SOL1(1,1),SOL1(2,1),SOL1(3,1)
c1370  FORMAT(/2X,' (ADD GEN VAR, DOM GEN VAR, ERR VAR) = ',
c     &        /24X,'(',g15.4,', ',g15.4,', ',g15.4,')')
c      WRITE(3,1375) DSQRT(AIN1(1,1)),DSQRT(AIN1(2,2)),DSQRT(AIN1(3,3))
c1375  FORMAT(3X,' WITH STANDARD ERRORS OF ',5X,
c     &  g11.4,4X,g11.4,4X,g11.4)

      WRITE(3,1380) 1.,AIN1(1,2)/DSQRT(AIN1(1,1)*AIN1(2,2)),
     & AIN1(1,3)/DSQRT(AIN1(1,1)*AIN1(3,3)),
     & AIN1(1,2)/DSQRT(AIN1(1,1)*AIN1(2,2)),1.,
     & AIN1(2,3)/DSQRT(AIN1(2,2)*AIN1(3,3)),
     & AIN1(3,1)/DSQRT(AIN1(1,1)*AIN1(3,3)),
     & AIN1(3,2)/DSQRT(AIN1(2,2)*AIN1(3,3)),1.
1380  FORMAT(/,3X,' The asymptotic correlation matrix  = ',
     &  3X,3(F7.4,5X)/44X,3(F7.4,5X)/44X,3(F7.4,5X))
      WRITE(3,1385) GOF1,CHI123(GOF1,3)
1385  FORMAT(/,3X,' Goodness-of-fit likelihood ratio statistic = ',
     &  F7.3,' ON 3 DF, P = ',F5.3)
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'ACE MODEL ML ESTIMATES'

      call sehml(sol2(1,1),sol2(2,1),sol2(3,1),DSQRT(AIN2(1,1)),
     &   DSQRT(AIN2(2,2)),DSQRT(AIN2(3,3)),
     & AIN2(1,2)/DSQRT(AIN2(1,1)*AIN2(2,2)),
     & AIN2(1,3)/DSQRT(AIN2(1,1)*AIN2(3,3)),
     & AIN2(2,3)/DSQRT(AIN2(3,3)*AIN2(2,2)),hacehat,sehace,
     &  cacehat,secace,eacehat,seeace)
      write(3,1470)
1470  format(/32x,'Estimate',6x,'S. E.',6x,'Proportion',4x,'S. E.'/)
      write(3,1472) sol2(1,1),DSQRT(AIN2(1,1)),hacehat,sehace
1472  format(' Additive genetic variance',2x,g15.4,2x,g11.4,5x,f5.3,
     &         4x,f5.3/)
      write(3,1474) sol2(2,1),DSQRT(AIN2(2,2)),cacehat,secace
1474  format(' Common environmenal variance',1x,g13.4,2x,g11.4,5x,f5.3,
     &         4x,f5.3/)
      write(3,1476) sol2(3,1),DSQRT(AIN2(3,3)),eacehat,seeace
1476  format(' Error variance',13x,g15.4,2x,g11.4,5x,f5.3,
     &         4x,f5.3//)

c      WRITE(3,1470) SOL2(1,1),SOL2(2,1),SOL2(3,1)
c1470  FORMAT(/2X,' (ADD GEN VAR, COMMON ENV VAR, ERR VAR) = ',
c     &       /24X,'(',g15.4,', ',g15.4,', ',g15.4,')')
c      WRITE(3,1480) DSQRT(AIN2(1,1)),DSQRT(AIN2(2,2)),DSQRT(AIN2(3,3))
c1480  FORMAT(3X,' WITH STANDARD ERRORS OF ',5X,
c     &  g11.4,4X,g11.4,4X,g11.4)


      WRITE(3,1495) 1.,AIN2(1,2)/DSQRT(AIN2(1,1)*AIN2(2,2)),
     & AIN2(1,3)/DSQRT(AIN2(1,1)*AIN2(3,3)),
     & AIN2(1,2)/DSQRT(AIN2(1,1)*AIN2(2,2)),1.,
     & AIN2(2,3)/DSQRT(AIN2(2,2)*AIN2(3,3)),
     & AIN2(3,1)/DSQRT(AIN2(1,1)*AIN2(3,3)),
     & AIN2(3,2)/DSQRT(AIN2(2,2)*AIN2(3,3)),1.
1495  FORMAT(/,3X,' The asymptotic correlation matrix  = ',
     &  3X,3(F7.4,5X)/44X,3(F7.4,5X)/44X,3(F7.4,5X))
      WRITE(3,1500) GOF2,CHI123(GOF2,3)
1500  FORMAT(/,3X,' Goodness-of-fit likelihood ratio statistic = ',
     &  F7.3,' ON 3 DF, P = ',F5.3)
1650  CONTINUE
      WRITE(3,200)
      WRITE(3,1290) FF
C***********************************************************C
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'TWINAN90 TWIN ANALYSIS PROGRAM, VERSION 1.0'
      WRITE(3,*) 'TWINAN90 PATH ANALYSIS OUTPUT, page 2'
      WRITE(3,250) NAMB
      IF (XLAMB .NE. 1.) WRITE(3,260) XLAMB
      WRITE(3,270) N(1),N(2)
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'FOR THE COVARIANCE MATRICES'
      WRITE(3,*)
      WRITE(3,1350) VMZ1,VAR1*VICMZ,VDZ1,VAR2*VICDZ,
     &              VAR1*VICMZ,VMZ2,VAR2*VICDZ,VDZ2
C***********************************************************C
      WRITE(3,200)
      WRITE(3,*) 'AE MODEL ML ESTIMATES'

      call sehml(sol3(1,1),sol3(3,1),0.d0,DSQRT(AIN3(1,1)),
     &   DSQRT(AIN3(3,3)),0.d0,
     & AIN3(1,3)/DSQRT(AIN3(1,1)*AIN3(3,3)),0.d0,0.d0,
     & haehat,sehae,eaehat,seeae,xaehat,sexae)
      write(3,1510)
1510  format(/32x,'Estimate',6x,'S. E.',6x,'Proportion',4x,'S. E.'/)
      write(3,1512) sol3(1,1),DSQRT(AIN3(1,1)),haehat,sehae
1512  format(' Additive genetic variance',2x,g15.4,2x,g11.4,5x,f5.3,
     &         4x,f5.3/)
      write(3,1516) sol3(3,1),DSQRT(AIN3(3,3)),eaehat,seeae
1516  format(' Error variance',13x,g15.4,2x,g11.4,5x,f5.3,
     &         4x,f5.3/)

c      WRITE(3,1510) SOL3(1,1),SOL3(3,1)
c1510  FORMAT(/2X,' (ADD GEN VAR, ERR VAR) = ',10X,
c     &        '(',g11.4,', ',g11.4,')')
c      WRITE(3,1530) DSQRT(AIN3(1,1)),DSQRT(AIN3(3,3))
c1530  FORMAT(3X,' WITH STANDARD ERRORS OF ',8X,g11.4,5X,g11.4)



      WRITE(3,1550) 1.,AIN3(1,3)/DSQRT(AIN3(1,1)*AIN3(3,3)),
     & AIN3(3,1)/DSQRT(AIN3(1,1)*AIN3(3,3)),1.
1550  FORMAT(/,3X,' The asymptotic correlation matrix  = ',
     &  3X,2(F7.4,5X)/44X,2(F7.4,5X))
      WRITE(3,1570) GOF3,CHI123(GOF3,4)
1570  FORMAT(/,3X,' Goodness-of-fit likelihood ratio statistic = ',
     &  F7.3,' ON 4 DF, P = ',F5.3)
      WRITE(3,*)
      WRITE(3,200)
C***********************************************************C
      WRITE(3,*) 'CE MODEL ML ESTIMATES'

      call sehml(sol4(1,1),sol4(3,1),0.d0,DSQRT(AIN4(1,1)),
     &   DSQRT(AIN4(3,3)),0.d0,
     & AIN4(1,3)/DSQRT(AIN4(1,1)*AIN4(3,3)),0.d0,0.d0,
     & hcehat,sehce,ecehat,seece,xaehat,sexae)
      write(3,1590)
1590  format(/32x,'Estimate',6x,'S. E.',6x,'Proportion',4x,'S. E.'/)
      write(3,1592) sol4(1,1),DSQRT(AIN4(1,1)),hcehat,sehce
1592  format(' Common environmental variance',1x,g12.4,2x,g11.4,5x,f5.3,
     &         4x,f5.3/)
      write(3,1596) sol4(3,1),DSQRT(AIN4(3,3)),ecehat,seece
1596  format(' Error variance',13x,g15.4,2x,g11.4,5x,f5.3,
     &         4x,f5.3/)

      WRITE(3,1600) 1.,AIN4(1,3)/DSQRT(AIN4(1,1)*AIN4(3,3)),
     & AIN4(3,1)/DSQRT(AIN4(1,1)*AIN4(3,3)),1.
1600  FORMAT(/,3X,' The asymptotic correlation matrix  = ',
     &  3X,2(F7.4,5X)/44X,2(F7.4,5X))
      WRITE(3,1610) GOFCE,CHI123(GOFCE,4)
1610  FORMAT(/,3X,' Goodness-of-fit likelihood ratio statistic = ',
     &  F7.3,' ON 4 DF, P = ',F5.3)
      WRITE(3,*)
      WRITE(3,200)
C***********************************************************C
      WRITE(3,1655) GOF4,CHI123(GOF4,1)
1655  FORMAT(3X,'Model E vs model AE likelihood ratio statistic = ',
     &  F7.3,' ON 1 DF, P = ',F5.3)
      WRITE(3,1660) GOFC,CHI123(GOFC,1)
1660  FORMAT(3X,'Model E vs model CE likelihood ratio statistic = ',
     &  F7.3,' ON 1 DF, P = ',F5.3)
      WRITE(3,1670) GOF3 -GOF2,CHI123(GOF3-GOF2,1)
1670  FORMAT(3X,'Model AE vs model ACE likelihood ratio statistic = ',
     &  F7.3,' ON 1 DF, P = ',F5.3)
      WRITE(3,1680) GOFCE -GOF2,CHI123(GOFCE-GOF2,1)
1680  FORMAT(3X,'Model CE vs model ACE likelihood ratio statistic = ',
     &  F7.3,' ON 1 DF, P = ',F5.3)
      WRITE(3,1690) GOF3 -GOF1,CHI123(GOF3-GOF1,1)
1690  FORMAT(3X,'Model AE vs model ADE likelihood ratio statistic = ',
     &  F7.3,' ON 1 DF, P = ',F5.3)

C***********************************************************C
c      WRITE(3,1510) SOL3(1,1),SOL3(3,1)
c1510  FORMAT(/2X,' (ADD GEN VAR, ERR VAR) = ',10X,
c     &        '(',g11.4,', ',g11.4,')')
c      WRITE(3,1530) DSQRT(AIN3(1,1)),DSQRT(AIN3(3,3))
c1530  FORMAT(3X,' WITH STANDARD ERRORS OF ',8X,g11.4,5X,g11.4)
c      WRITE(3,1550) 1.,AIN3(1,3)/DSQRT(AIN3(1,1)*AIN3(3,3)),
c     & AIN3(3,1)/DSQRT(AIN3(1,1)*AIN3(3,3)),1.
c1550  FORMAT(/,3X,' THE ASYMPTOTIC CORRELATION MATRIX  = ',
c     &  3X,2(F7.4,5X)/44X,2(F7.4,5X))
c      WRITE(3,1570) GOF3,CHI123(GOF3,4)
c1570  FORMAT(/,3X,' GOODNESS-OF-FIT LIKELIHOOD RATIO STATISTIC = ',
c     &  F7.3,' ON 4 DF, P = ',F5.3)
c      WRITE(3,*)
c      WRITE(3,1580) GOF4,CHI123(GOF4,1)
c1580  FORMAT(3X,'MODEL E VS MODEL 3 LIKELIHOOD RATIO STATISTIC = ',
c     &  F7.3,' ON 1 DF, P = ',F5.3)
c      WRITE(3,1600) GOF3 -GOF1,CHI123(GOF3-GOF1,1)
c1600  FORMAT(3X,'MODEL 3 VS MODEL 1 LIKELIHOOD RATIO STATISTIC = ',
c     &  F7.3,' ON 1 DF, P = ',F5.3)
c      WRITE(3,1630) GOF3 -GOF2,CHI123(GOF3-GOF2,1)
c1630  FORMAT(3X,'MODEL 3 VS MODEL 2 LIKELIHOOD RATIO STATISTIC = ',
c     &  F7.3,' ON 1 DF, P = ',F5.3)
      WRITE(3,200)
      IF (PED .EQ. 1.) WRITE(3,1700) NAMP
      IF (PED .NE. 1.) WRITE(3,1750)
1700  FORMAT(37H   THE OUTPUT PEDIGREE FILE IS NAMED ,A16)
1750  FORMAT(39H   NO OUTPUT PEDIGREE FILE WAS CREATED.,A16)
      CLOSE(UNIT=3)
c     WRITE(*,*) 'IF YOU WOULD LIKE TO RUN THE PROGRAM AGAIN,'
c     WRITE(*,*) '    ENTER A 1, OTHERWISE ENTER A 0 .  '
c     READ(*,*) DONE
      done=0
      IF (DONE .EQ. 1.) GO TO 10
c     STOP
      return
      END
C-----------------------------------------------------------C
C             SUBROUTINES AND FUNCTIONS                     C
C-----------------------------------------------------------C
      SUBROUTINE DESCRIP(DATA,NUM,AVE,VAR,V2,V3,SKEW,CURT,SSA,
     & SSW,mp)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(mp,2)
      NN = DBLE(2*NUM)
      SUM1 = 0.
      SUM2 = 0.
      SUM3 = 0.
      SUM4 = 0.
      S1 = 0.
      S2 = 0.
      S3 = 0.
      VAR = 0.
      SKEW = 0.
      CURT = 0.
      DO 8 I = 1,NUM
        SUM1 = SUM1 + DATA(I,1) + DATA(I,2)
        SUM3 = SUM3 + DATA(I,1)
        SUM4 = SUM4 + DATA(I,2)
8     CONTINUE
      AVE = SUM1/NN
      AVE2 = SUM3/DBLE(NUM)
      AVE3 = SUM4/DBLE(NUM)
      DO 17 I = 1,NUM
        SUM2 = SUM2 +2*( (DATA(I,1)+DATA(I,2))/2. -AVE)**2
        S2 = S2 +(DATA(I,1) -AVE2)**2
        S3 = S3 +(DATA(I,2) -AVE3)**2
      DO 15 J = 1,2
        S1 = DATA(I,J) - AVE
        P = S1*S1
        VAR = VAR + P
        P = P*S1
        SKEW = SKEW + P
        P = P*S1
        CURT = CURT + P
15    CONTINUE
17    CONTINUE
      SSA = SUM2
      SSW = (VAR - SUM2)
      VAR = VAR/(NN -1.)
      V2 = S2/(DBLE(NUM) -1.)
      V3 = S3/(DBLE(NUM) -1.)
      SKEW = NN*SKEW/((NN-1.)*(NN-2.)*DSQRT(VAR)**3)
      CURT = NN*(NN+1.)*CURT/((NN-1.)*(NN-2.)*(NN-3.)*VAR**2)
     &     - 3.*(NN-1.)**2/((NN-2.)*(NN-3.))
      RETURN
      END
C-----------------------------------------------------------C
C     CHI123 COMPUTES RIGHT TAIL P VALUES FOR CHI-SQUARE    C
C     STATISTICS TO ABOUT 3 DECIMAL ACCURACY                C
C-----------------------------------------------------------C
      DOUBLE PRECISION FUNCTION CHI123(Y,NDF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      X = DSQRT(Y)
      CHI123 = 0
        IF (NDF .EQ. 1) CHI123 = 2*(1. - FO(X))
        IF (NDF .EQ. 2) CHI123 = EXP(-1*X/2)
        IF (NDF .GT. 2) THEN
          DF = DBLE(NDF)
          B1 = (Y/DF)**.333333 -1. + 2/(9*DF)
          B2 = SQRT(2/(9*DF))
          CHI123 = 1. - FO(B1/B2)
        END IF
      RETURN
      END
C-----------------------------------------------------------C
      DOUBLE PRECISION FUNCTION FO(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      A1 = 1. + 0.04417*X**2
      A2 = -1.5976*X*A1
      A3 = 1. + DEXP(A2)
      FO = 1./A3
      END
C-----------------------------------------------------------C
      SUBROUTINE SORT(N,RA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RA(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
C-----------------------------------------------------------C
      SUBROUTINE CONVKS(DATA,SDATA,N,AVE,VAR,D,PROB,mp)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(mp,2),SDATA(2*N)
      DO 10 I = 1,N
      SDATA(I) = DATA(I,1)
10    CONTINUE
      DO 20 I = 1,N
      SDATA(N+I) = DATA(I,2)
20    CONTINUE
      NUMB = 2*N
      CALL KSONE(SDATA,NUMB,AVE,VAR,D,PROB)
      RETURN
      END
C-----------------------------------------------------------C
      SUBROUTINE KSONE(DATA,N,AVE,VAR,DST,PROB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(N)
      CALL SORT(N,DATA)
      EN=N
      D=0.
      FO=0.
      DO 11 J=1,N
        FN=J/EN
        Z = (DATA(J) -AVE)/DSQRT(VAR)
        FF=CDNORM(Z)
        DT=DMAX1(ABS(FO-FF),ABS(FN-FF))
        IF(DT.GT.D)D=DT
        FO=FN
11    CONTINUE
      DST = (DSQRT(EN) -.01 + .85/DSQRT(EN))*D
      IF (DST .LE. .775) THEN
        PROB = .15
      ELSE IF (DST .GT. .775 .AND. DST .LE. .819) THEN
        PROB = .15 - .05*(.775 -DST)/(.775 -.819)
      ELSE IF (DST .GT. .819 .AND. DST .LE. .895) THEN
        PROB = .10 - .05*(.819 -DST)/(.819 -.895)
      ELSE IF (DST .GT. .895 .AND. DST .LE. .955) THEN
        PROB = .05 - .025*(.895 -DST)/(.895 -.955)
      ELSE IF (DST .GT. .955 .AND. DST .LE. 1.035) THEN
        PROB = .025 - .015*(.955 -DST)/(.955 -1.035)
      ELSE IF (DST .GT. 1.035) THEN
        PROB = .01
      END IF
      RETURN
      END
C-----------------------------------------------------------C
      DOUBLE PRECISION FUNCTION CDNORM(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      A1 = 1. + 0.04417*X**2
      A2 = -1.5976*X*A1
      A3 = 1. + DEXP(A2)
      CDNORM = 1./A3
      END
C-----------------------------------------------------------C
      DOUBLE PRECISION FUNCTION BCOX(X,N,mp)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(mp,2),XT(mp,2)
      WRITE(2,*)'V, SIG, SUM LNY, LIKELIHOOD, VMAX, MAX LIKELIHOOD ARE '
      VMAX = -1.
      SUM1 = 0.
      DO 5 I = 1,N
        X(I,1) = X(I,1) + 2.
        X(I,2) = X(I,2) + 2.
        SUM1 = SUM1 + DLOG(X(I,1)) + DLOG(X(I,2))
5     CONTINUE
C      WRITE(*,*) 'SUM1 IS ',SUM1
      DO 30 J = 1,32
        V = -1. + 0.1 * (J - 1)
        SUM2 = 0.
        SUM3 = 0.
        DO 10 I = 1,N
          IF (V .NE. 0.0) XT(I,1) = ( X(I,1)**V -1.)/V
          IF (V .EQ. 0.0) XT(I,1) = DLOG(X(I,1))
          IF (V .NE. 0.0) XT(I,2) = ( X(I,2)**V -1.)/V
          IF (V .EQ. 0.0) XT(I,2) = DLOG(X(I,2))
          SUM2 = SUM2 + XT(I,1) + XT(I,2)
10      CONTINUE
        XTBAR = SUM2/(2.*DBLE(N))
C     WRITE(*,*) 'XTBAR IS ',XTBAR
        DO 20 I = 1,N
          SUM3 = SUM3 +(XT(I,1) -XTBAR)**2 +(XT(I,2) -XTBAR)**2
20      CONTINUE
        SIG = SUM3/(2*DBLE(N))
        ALIKE = -DLOG(SIG)/2. + (V-1)*SUM1/(2*DBLE(N))
        IF (V .EQ. -1.) ALMAX = ALIKE
        IF (ALIKE .GT. ALMAX) VMAX = V
        IF (ALIKE .GT. ALMAX) ALMAX = ALIKE
       WRITE(2,25) V,SIG,SUM1,ALIKE,VMAX,ALMAX
25    FORMAT(2X,6(g11.4,1X))
30    CONTINUE
      BCOX = VMAX
      END
C-----------------------------------------------------------C
      DOUBLE PRECISION FUNCTION BETAI(A,B,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
cc      IF(X.LT.0..OR.X.GT.1.)PAUSE 'bad argument X in BETAI'
      IF(X.EQ.0..OR.X.EQ.1.)THEN
        BT=0.
      ELSE
        BT=DEXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)
     *      +A*DLOG(X)+B*DLOG(1.-X))
      ENDIF
      IF(X.LT.(A+1.)/(A+B+2.))THEN
        BETAI=BT*BETACF(A,B,X)/A
        RETURN
      ELSE
        BETAI=1.-BT*BETACF(B,A,1.-X)/B
        RETURN
      ENDIF
      END
C-----------------------------------------------------------C
      DOUBLE PRECISION FUNCTION BETACF(A,B,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (ITMAX=100,EPS=3.D-7)
      AM=1.
      BM=1.
      AZ=1.
      QAB=A+B
      QAP=A+1.
      QAM=A-1.
      BZ=1.-QAB*X/QAP
      DO 11 M=1,ITMAX
        EM=M
        TEM=EM+EM
        D=EM*(B-M)*X/((QAM+TEM)*(A+TEM))
        AP=AZ+D*AM
        BP=BZ+D*BM
        D=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM))
        APP=AP+D*AZ
        BPP=BP+D*BZ
        AOLD=AZ
        AM=AP/BPP
        BM=BP/BPP
        AZ=APP/BPP
        BZ=1.
        IF(ABS(AZ-AOLD).LT.EPS*ABS(AZ))GO TO 1
11    CONTINUE
cc      PAUSE 'A or B too big, or ITMAX too small'
1     BETACF=AZ
      RETURN
      END
C-----------------------------------------------------------C
      DOUBLE PRECISION FUNCTION GAMMLN(XX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*DLOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+DLOG(STP*SER)
      RETURN
      END
C------------------------------------------------------------------C
C    SUBROUTINES FOR THE PATH ANALYSIS ESTIMATES                   C
C------------------------------------------------------------------C
      SUBROUTINE QHOUSE(VN,IM,T,BETA,SN,A,VL,IN,MON,STEP)
C--------------------------------------------------C
C                ENTER A,D,U                       C
C--------------------------------------------------C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(IN,IN),U(3,1),AIN(3,3),VN(2)
      DIMENSION STEP(3),T(2,2,2),DF(3,1),DGF(2)
      DIMENSION BETA(3,1),SN(3,1),AU(3,1)
      OLDLIK = 0.
      DGF(1) = (VN(1) -1.)/2.
      DGF(2) = (VN(2) -1.)/2.
      IF (IM .EQ. 3) BETA(2,1) = 0.
      DO 7000 III = 1,500
      P1 = BETA(1,1)
      P2 = BETA(2,1)
      P3 = BETA(3,1)
      DO 5 I = 1,3
        U(I,1) = 0.
5     CONTINUE
      DO 10 I = 1,3
      DO 7 J = 1,3
        A(I,J) = 0.
7     CONTINUE
10    CONTINUE
      DO 70 I = 1,2
      CALL INFDER(IM,I,T(I,1,1),T(I,1,2),T(I,2,2),P1,P2,P3,DER1,DER2,
     &            DER3,D2P1,D2P2,D2P3,D2P12,D2P13,D2P23)
      U(1,1) = U(1,1) + DGF(I)*DER1
      U(2,1) = U(2,1) + DGF(I)*DER2
      U(3,1) = U(3,1) + DGF(I)*DER3
      A(1,1) = A(1,1) + DGF(I)*D2P1
      A(2,2) = A(2,2) + DGF(I)*D2P2
      A(3,3) = A(3,3) + DGF(I)*D2P3
      A(1,2) = A(1,2) + DGF(I)*D2P12
      A(1,3) = A(1,3) + DGF(I)*D2P13
      A(2,3) = A(2,3) + DGF(I)*D2P23
*     IF (MON .EQ. 3) WRITE(*,*) 'T1-T3, P1-P3, DER1-DER3'
*     IF (MON .EQ. 3) WRITE(*,*)
*    &    T(I,1,1),T(I,1,2),T(I,2,2),P1,P2,P3,DER1,DER2,DER3
*     IF (MON .EQ. 3) WRITE(*,*)
*    &    D2P1,D2P2,D2P3,D2P12,D2P13,D2P23
      IF (MON .EQ. 3) WRITE(2,*) 'T1-T3, P1-P3, DER1-DER3'
      IF (MON .EQ. 3) WRITE(2,*)
     &    T(I,1,1),T(I,1,2),T(I,2,2),P1,P2,P3,DER1,DER2,DER3
      IF (MON .EQ. 3) WRITE(2,*)
     &    D2P1,D2P2,D2P3,D2P12,D2P13,D2P23
70    CONTINUE
      A(3,2) = A(2,3)
      A(2,1) = A(1,2)
      A(3,1) = A(1,3)
      WRITE(2,*) ' ITERATION ',III,' COMPLETED'
      IF (MON .GE. 2) WRITE(2,*) 'THE FIRST DERIVATIVES ARE'
      IF (MON .GE. 2) WRITE(2,99) U
      IF (MON .EQ. 3) WRITE(2,*) 'A IS'
      IF (MON .EQ. 3) WRITE(2,99) A
*     WRITE(*,*) ' ITERATION ',III,' COMPLETED'
*     IF (MON .GE. 2) WRITE(*,*) 'THE FIRST DERIVATIVES ARE'
*     IF (MON .GE. 2) WRITE(*,99) U
*     IF (MON .EQ. 3) WRITE(*,*) 'A IS'
*     IF (MON .EQ. 3) WRITE(*,99) A
99    FORMAT(3(g14.4,1X))
101   FORMAT(4X,g20.8)
C-------------------------------------------C
C     OBTAIN NEW PARAMETER ESTIMATES        C
C-------------------------------------------C
      ERR = 0.
      IF (IM .EQ. 1 .OR. IM .EQ. 2) THEN
        CALL SINV3(A,AIN,ERR)
        CALL MULT(AIN,U,AU,3,3,1)
      ELSE IF (IM .EQ. 3) THEN
        D = ( A(1,1)*A(3,3) - A(1,3)**2 )
        AU(1,1) = ( A(3,3)*U(1,1) - A(1,3)*U(3,1) )/D
        AU(2,1) = 0.
        AU(3,1) = ( A(1,1)*U(3,1) - A(3,1)*U(1,1) )/D
      END IF
      DO 200 I = 1,3
        IF (ERR .EQ. 0.) SN(I,1) = BETA(I,1) +STEP(I)*AU(I,1)
        IF (ERR .EQ. 1.) SN(I,1) = BETA(I,1) +STEP(I)*U(I,1)
200   CONTINUE
      IF (MON .GE. 2) WRITE(2,*) 'THE NEW ESTIMATES ARE '
      IF (MON .GE. 2) WRITE(2,230) SN
*     IF (MON .GE. 2) WRITE(*,*) 'THE NEW ESTIMATES ARE '
*     IF (MON .GE. 2) WRITE(*,230) SN
230   FORMAT(4X,3(g14.4,1X))
      SUMVL = VLLIK(DGF,IM,T,SN)
      IF (MON .GE. 2) WRITE(2,*) 'THE LOGLIKELIHOOD IS'
      IF (MON .GE. 2) WRITE(2,101) SUMVL
*     IF (MON .GE. 2) WRITE(*,*) 'THE LOGLIKELIHOOD IS'
*     IF (MON .GE. 2) WRITE(*,101) SUMVL
C------------------------------------------------------C
C       CHECK GLOBAL CONVERGENCE CRITERION             C
C------------------------------------------------------C
      SUM1 = 0.
      SUM2 = 0.
      DFLIK = SUMVL - OLDLIK
      DO 6000 II=1,3
      DF(II,1) = SN(II,1) - BETA(II,1)
      SUM1 = SUM1 + ( BETA(II,1)**2  )
      SUM2 = SUM2 + ( DF(II,1)**2 )
6000  CONTINUE
      CHECK = SUM2/SUM1
      IF (CHECK  .LT.  1.0D-12) GO TO 9000
      IF (MON.GE.2) WRITE(2,6050) SUM1,SUM2,CHECK
*     IF (MON.GE.2) WRITE(*,6050) SUM1,SUM2,CHECK
6050  FORMAT(2X,'SUM1 = ',g18.4,' SUM2 = ',g10.4,1X,f18.14)
      DO 6500 II=1,3
      IF (III .EQ. 1 .OR. DFLIK .GE. 0.00)
     &     BETA(II,1) = SN(II,1)
      IF (III .EQ. 1 .OR. DFLIK .GE. 0.00)
     &     OLDLIK = SUMVL
      IF (III .GT. 1 .AND. DFLIK .LT. 0.00)
     &     STEP(II) = .9*STEP(II)
6500  CONTINUE
7000  CONTINUE
9000  CONTINUE
      VL = SUMVL
      RETURN
      END
C----------------------------------------------------------------C
C       SUBROUTINE TO CALCULATE DERIVATIVES                      C
C----------------------------------------------------------------C
      SUBROUTINE INFDER(IM,IZ,T1,T2,T3,P1,P2,P3,DER1,DER2,DER3,
     &                  D2P1,D2P2,D2P3,D2P12,D2P13,D2P23)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION PHI2(2,2),DEL(2,2),O(2,2),OIN(2,2),TW(2,2)
      DIMENSION XM2(2,2),XM3(2,2),XM4(2,2),XM5(2,2)
      DIMENSION XM6(2,2),XM8(2,2),XM9(2,2),XM10(2,2)
      DIMENSION XM11(2,2),XM12(2,2),XM13(2,2)
      DIMENSION XB2(2,2),XB3(2,2),XB5(2,2)
      DIMENSION XB6(2,2),XB8(2,2)
      DIMENSION XB11(2,2),XB12(2,2),XB13(2,2),XB14(2,2)
      DIMENSION XB15(2,2),XB16(2,2),XB17(2,2),XB18(2,2)
      DER1 = 0.
      DER2 = 0.
      DER3 = 0.
      D2P1 = 0.
      D2P2 = 0.
      D2P3 = 0.
      D2P12 = 0.
      D2P13 = 0.
      D2P23 = 0.
      DO 7 I = 1,2
      DO 5 J = 1,2
        PHI2(I,J) = 1.
        DEL(I,J) = 1.
5     CONTINUE
7     CONTINUE
      O(1,1) = P1 +P2 +P3
      O(1,2) = P1 +P2
      IF (IZ .EQ. 2 .AND. IM .EQ. 1) THEN
        PHI2(1,2) = .5
        PHI2(2,1) = .5
        DEL(1,2) = .25
        DEL(2,1) = .25
        O(1,2) = .5*P1 +.25*P2
      ELSE IF (IZ .EQ. 2 .AND. IM .EQ. 2) THEN
        PHI2(1,2) = .5
        PHI2(2,1) = .5
        DEL(1,2) = 1.
        DEL(2,1) = 1.
        O(1,2) = .5*P1 +P2
      ELSE IF (IZ .EQ. 1 .AND. IM .EQ. 3) THEN
        DEL(1,1) = 0.
        DEL(2,2) = 0.
        DEL(1,2) = 0.
        DEL(2,1) = 0.
        O(1,2) = P1
      ELSE IF (IZ .EQ. 2 .AND. IM .EQ. 3) THEN
        PHI2(1,2) = .5
        PHI2(2,1) = .5
        DEL(1,1) = 0.
        DEL(2,2) = 0.
        DEL(1,2) = 0.
        DEL(2,1) = 0.
        O(1,2) = .5*P1
      END IF
      O(2,1) = O(1,2)
      O(2,2) = O(1,1)
      CALL SINV2(O,OIN)
      TW(1,1) = T1
      TW(1,2) = T2
      TW(2,2) = T3
      TW(2,1) = T2
      CALL MULT(OIN,PHI2,XM2,2,2,2)
      CALL MULT(TW,OIN,XB2,2,2,2)
      CALL MULT(XB2,PHI2,XB3,2,2,2)
      CALL MULT(XB3,OIN,XM3,2,2,2)
      DER1 = -TR(XM2,2) +TR(XM3,2)
      CALL MULT(OIN,DEL,XM4,2,2,2)
      CALL MULT(TW,OIN,XB5,2,2,2)
      CALL MULT(XB5,DEL,XB6,2,2,2)
      CALL MULT(XB6,OIN,XM5,2,2,2)
      DER2 = -TR(XM4,2) +TR(XM5,2)
      CALL MULT(TW,OIN,XB8,2,2,2)
      CALL MULT(XB8,OIN,XM6,2,2,2)
      DER3 = -TR(OIN,2) +TR(XM6,2)
      CALL MULT(OIN,PHI2,XB11,2,2,2)
      CALL MULT(XB11,OIN,XB12,2,2,2)
      CALL MULT(XB12,PHI2,XM8,2,2,2)
      D2P1 = TR(XM8,2)
      CALL MULT(OIN,DEL,XB13,2,2,2)
      CALL MULT(XB13,OIN,XB14,2,2,2)
      CALL MULT(XB14,DEL,XM9,2,2,2)
      D2P2 = TR(XM9,2)
      CALL MULT(OIN,OIN,XM10,2,2,2)
      D2P3 = TR(XM10,2)
      CALL MULT(OIN,PHI2,XB15,2,2,2)
      CALL MULT(XB15,OIN,XB16,2,2,2)
      CALL MULT(XB16,DEL,XM11,2,2,2)
      D2P12 = TR(XM11,2)
      CALL MULT(OIN,PHI2,XB17,2,2,2)
      CALL MULT(XB17,OIN,XM12,2,2,2)
      D2P13 = TR(XM12,2)
      CALL MULT(OIN,DEL,XB18,2,2,2)
      CALL MULT(XB18,OIN,XM13,2,2,2)
      D2P23 = TR(XM13,2)
      RETURN
      END
C----------------------------------------------------------------C
      DOUBLE PRECISION FUNCTION VLLIK(VN,IM,TWIN,P)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TWIN(2,2,2),P(3,1),T(2,2),VN(2)
      DIMENSION XA(2,2),O(2,2),OIN(2,2)
C     WRITE(2,*) 'CALCULATING LL VALUE'
      SUM = 0.
      DO 70 I = 1,2
      O(1,1) = P(2,1) +P(3,1) +P(1,1)
      O(1,2) = P(2,1) +P(1,1)
      IF (I .EQ. 2 .AND. IM .EQ. 1) O(1,2) = .5*P(1,1) +.25*P(2,1)
      IF (I .EQ. 2 .AND. IM .EQ. 2) O(1,2) = .5*P(1,1) +P(2,1)
      IF (I .EQ. 1 .AND. IM .EQ. 3) O(1,2) = P(1,1)
      IF (I .EQ. 2 .AND. IM .EQ. 3) O(1,2) = .5*P(1,1)
      O(2,1) = O(1,2)
      O(2,2) = O(1,1)
      CALL SINV2(O,OIN)
      V2 = 0.
      V1 = (O(1,1)*O(2,2) -O(1,2)*O(2,1))
      IF (V1 .GE. 0.00) V2 = DLOG(V1)
      T(1,1) = TWIN(I,1,1)
      T(1,2) = TWIN(I,1,2)
      T(2,2) = TWIN(I,2,2)
      T(2,1) = T(1,2)
      CALL MULT(T,OIN,XA,2,2,2)
C      WRITE(2,*)
C     &    'T, DF, O, OIN, S*OIN, LOG(DET(O)), AND TR(S*OIN) ARE'
C      WRITE(2,66) T
C      WRITE(2,67) VN(I)
C      WRITE(2,66) O
C      WRITE(2,66) OIN
C      WRITE(2,66) XA
C      WRITE(2,65) V2,TR(XA,2)
cc65    FORMAT(2X,2(F10.4,2X))
cc66    FORMAT(2X,4(F15.7,1X))
cc67    FORMAT(2X,F8.1)
      SUM = SUM -VN(I)*(V2 +TR(XA,2))
70    CONTINUE
      VLLIK = SUM
      RETURN
      END
C----------------------------------------------------------------C
      DOUBLE PRECISION FUNCTION TR(A,ID)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(ID,ID)
      SUM = 0.
      DO 5 I = 1,ID
        SUM = SUM + A(I,I)
5     CONTINUE
      TR = SUM
      RETURN
      END
C----------------------------------------------------------------C
      SUBROUTINE MULT(A,B,C,ID,IE,IF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(ID,IE),B(IE,IF),C(ID,IF)
      DO 2200 I=1,ID
       DO 2100 J=1,IF
         C(I,J)=0
         DO 2000 K=1,IE
           C(I,J)=C(I,J) + A(I,K)*B(K,J)
2000     CONTINUE
2100   CONTINUE
2200  CONTINUE
      RETURN
      END
C----------------------------------------------------------------C
      subroutine sehml(SIGA2,SIGO2,SIGE2,SEA2,SEO2,SEE2,CAO,CAE,COE,
     &                 ha,seha,hc,sehc,he,sehe)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SIGTOT = SIGA2 +SIGO2 +SIGE2
      S12 = SEA2**2 +CAO*SEA2*SEO2 +CAE*SEA2*SEE2
      S22 = SEA2**2 +SEO2**2 +SEE2**2
     &      +2*CAO*SEA2*SEO2 +2*CAE*SEA2*SEE2 +2*COE*SEO2*SEE2
      S11 = SEA2**2
      ha = SIGA2/SIGTOT
      PF1 = 1./SIGTOT
      PF2 = -SIGA2/SIGTOT**2
      VMLH2 = S11*PF1**2 +2.*S12*PF1*PF2 +S22*PF2**2
      seha = DSQRT(VMLH2)
C----------------------------------------------------C
      SS12 = SEO2**2 +CAO*SEA2*SEO2 +COE*SEO2*SEE2
      SS22 = SEA2**2 +SEO2**2 +SEE2**2
     &      +2*CAO*SEA2*SEO2 +2*CAE*SEA2*SEE2 +2*COE*SEO2*SEE2
      SS11 = SEO2**2
      hc = SIGO2/SIGTOT
      PFO1 = 1./SIGTOT
      PFO2 = -SIGO2/SIGTOT**2
      VMLO2 = SS11*PFO1**2 +2.*SS12*PFO1*PFO2 +SS22*PFO2**2
      sehc = DSQRT(VMLO2)
C----------------------------------------------------C
      SSS12 = SEE2**2 +CAE*SEA2*SEE2 +COE*SEO2*SEE2
      SSS11 = SEE2**2
      he = SIGE2/SIGTOT
      PFE2 = -SIGE2/SIGTOT**2
      VMLE2 = SSS11*PF1**2 +2.*SSS12*PF1*PFE2 +S22*PFE2**2
      sehe = DSQRT(VMLE2)
C----------------------------------------------------C
      return
      end
C----------------------------------------------------------------C
      SUBROUTINE CE(TCOV,VN,SOL4,AIN4,VLCE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VN(2),TCOV(2,2,2),SOL4(3,1),AIN4(3,3),A4(3,3)
      DO 10 I = 1,3
        SOL4(I,1) = 0.
        DO 8  J = 1,3
          AIN4(I,J) = 0.
          A4(I,J) = 0.
8       CONTINUE
10    CONTINUE
      X1 = VN(1)
      X2 = VN(2)
      degf = (x1 +x2 -2.)/2.
      DETMZ = TCOV(1,1,1)*TCOV(1,2,2) -TCOV(1,1,2)**2.
      DETDZ = TCOV(2,1,1)*TCOV(2,2,2) -TCOV(2,1,2)**2.
      EXPRMZ = TCOV(1,1,1) +TCOV(1,2,2) +2.*TCOV(1,1,2)
      EXPRDZ = TCOV(2,1,1) +TCOV(2,2,2) +2.*TCOV(2,1,2)
      FXPRMZ = TCOV(1,1,1) +TCOV(1,2,2) -2.*TCOV(1,1,2)
      FXPRDZ = TCOV(2,1,1) +TCOV(2,2,2) -2.*TCOV(2,1,2)
      VARE = ((X1-1.)/(X1+X2-2.))*.5*(TCOV(1,1,1)+TCOV(1,2,2))
     &      +((X2-1.)/(X1+X2-2.))*.5*(TCOV(2,1,1)+TCOV(2,2,2))
      VLE = -.5*(X1-1.)*(2.*DLOG(VARE) +(TCOV(1,1,1)+TCOV(1,2,2))/VARE)
     &  -.5*(X2-1.)*(2.*DLOG(VARE) +(TCOV(2,1,1)+TCOV(2,2,2))/VARE)
      VL0 = -.5*(X1-1.)*(DLOG(DETMZ) +2.)
     &  -.5*(X2-1.)*(DLOG(DETDZ) +2.)
      VLAMB1 = ((X1-1.)*(EXPRMZ) +(X2-1.)*(EXPRDZ))/
     &        (2.*(X1 +X2 -2))
      VLAMB2 = ((X1-1.)*(FXPRMZ) +(X2-1.)*(FXPRDZ))/
     &        (2.*(X1 +X2 -2))
      VLCE = -.5*(X1-1.)*(DLOG(VLAMB1*VLAMB2) +EXPRMZ/(2.*VLAMB1)
     &                                     +FXPRMZ/(2.*VLAMB2))
     &       -.5*(X2-1.)*(DLOG(VLAMB1*VLAMB2) +EXPRDZ/(2.*VLAMB1)
     &                                     +FXPRDZ/(2.*VLAMB2))
      SOL4(1,1) = (VLAMB1 -VLAMB2)/2.
      SOL4(3,1) = VLAMB2
      A4(1,1) = 4.D0*degf/(2.D0*SOL4(1,1) +SOL4(3,1))**2.
      A4(3,3) = 2.*degf*(2.D0*SOL4(1,1)**2. +2.D0*SOL4(1,1)*SOL4(3,1)
     &                   +SOL4(3,1)**2. )
     &         /(2.D0*SOL4(1,1)*SOL4(3,1) +SOL4(3,1)**2. )**2.
      A4(1,3) = 2.*degf/(2.D0*SOL4(1,1) +SOL4(3,1))**2.
      A4(3,1) = A4(1,3)
      DET = DABS(A4(1,1)*A4(3,3) -A4(1,3)**2. )
      AIN4(1,1) = A4(3,3)/DET
      AIN4(1,3) = -A4(1,3)/DET
      AIN4(3,3) = A4(1,1)/DET
      AIN4(3,1) = AIN4(1,3)
      RETURN
      END
C----------------------------------------------------------------C
      SUBROUTINE TRANS(A,B,ID,IE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(ID,IE),B(IE,ID)
      DO 2600 I=1,IE
      DO 2500 J=1,ID
        B(I,J)=A(J,I)
2500  CONTINUE
2600  CONTINUE
      RETURN
      END
C----------------------------------------------------------------C
      SUBROUTINE SINV2(A,AIN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(2,2),AIN(2,2)
      DET = A(1,1)*A(2,2) -A(1,2)*A(2,1)
      AIN(1,1) = A(2,2)/DET
      AIN(1,2) = -A(1,2)/DET
      AIN(2,1) = -A(2,1)/DET
      AIN(2,2) = A(1,1)/DET
      RETURN
      END
C----------------------------------------------------------------C
        SUBROUTINE SINV3(SA,AIN,ERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION SA(3,3),SB(3,3),SC(3,3)
        DIMENSION SD(3,3),SE(3,3),SF(3,3)
        DIMENSION SG(3,3),AIN(3,3)
      ERR = 0.
        DO 20 IJK = 1,3
          SB(IJK,IJK)=1.
20      CONTINUE
        SD(1,1) = 1./SA(1,1)
        SE(1,1) = SA(1,1)
        DO 30 I = 1,3
        SC(1,I) = SA(1,I)/SA(1,1)
30      CONTINUE
        DO 200 K = 2,3
        SUMK = 0.
        DO 50 J = 1,K-1
        SUMK = SUMK + SC(J,K)*SC(J,K)*SE(J,J)
50      CONTINUE
        SE(K,K) = SA(K,K) - SUMK
        DO 80 I = K,3
        SUMY = 0.
        DO 60 J = 1,K-1
        SUMY = SUMY + SC(J,K)*SC(J,I)*SE(J,J)
60      CONTINUE
        SC(K,I) = (SA(K,I) - SUMY)/SE(K,K)
80      CONTINUE
        DO 120 I = 1,K
        SUMZ = 0.
        DO 100 J = 1,K-1
        SUMZ = SUMZ + SC(J,K)*SD(J,I)*SE(J,J)
100     CONTINUE
        SD(K,I) = (SB(K,I) - SUMZ)/SE(K,K)
120     CONTINUE
200     CONTINUE
        DO 220 I = 1,3
        DO 210 J = 1,3
      IF (SE(I,J) .LT. 0) ERR = 1.
      IF (SE(I,J) .LT. 0) SE(I,J) = 0.
        SE(I,J) = DSQRT(SE(I,J))
210     CONTINUE
220     CONTINUE
        CALL MULT(SE,SD,SF,3,3,3)
        CALL TRANS(SF,SG,3,3)
        CALL MULT(SG,SF,AIN,3,3,3)
        RETURN
        END
C----------------------------------------------------------C
C    SUBROUTINE TO OUTPUT THE PEDIGREE FILE FOR FISHER     C
C----------------------------------------------------------C
      SUBROUTINE FOUT(NAM1,VMISS,NAM2,var,form2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*1 SEX1,SEX2,IS1,IS2
      CHARACTER*16 NAM1,NAM2
      CHARACTER*255 FORM2
      integer var
cc8     FORMAT(A16)
c     WRITE(*,*) 'ENTER THE NAME OF THE OUTPUT FILE.'
c     READ(*,8) NAM2
      OPEN(UNIT=3,FILE=NAM1,STATUS='OLD')
      OPEN(UNIT=4,FILE=NAM2,STATUS='NEW')
c     WRITE(*,*) 'ENTER THE NUMBER OF VARIABLES TO BE OUTPUT (1-3).'
c     READ(*,*) VAR
c     WRITE(*,*) 'ENTER THE FORTRAN FORMAT STATEMENT FOR EACH TWIN.'
c     WRITE(*,*) 'PAIR.  (INCLUDE AN "A1" CHARACTER IN THE FORMAT'
c     WRITE(*,*) 'FOR EACH PERSON, PRECEDING THE VARIABLES. THIS'
c     WRITE(*,*) 'VARIABLE CAN BE USED TO IDENTIFY TWIN GENDER.'
c     WRITE(*,*) 'AN EXAMPLE IS (1X,A1,5X,F6.2/1X,A1,5X,F6.2) FOR '
c     WRITE(*,*) 'OUTPUTING ONE VARIABLE.  THE "A1" MUST BE INCLUDED'
c     WRITE(*,*) 'EVEN IF IT IS A DUMMY VARIABLE, AND NO INFORMATION'
c     WRITE(*,*) '        IS AVAILABLE ON TWIN GENDER.)'
c     READ(*,10) FORM2
cc10    FORMAT(A80)
      READ(3,*) X1,X2
      WRITE(4,*) '(2I3,A8)'
      WRITE(4,*) '(A6,4(1X,A6),(T36,3(1X,A8),:))'
      KT1 = 0
      DO 15 J = 1,INT(X1)
      SEX1 = 'M'
      SEX2 = 'M'
      IF (VAR .EQ. 1.) READ(3,FORM2,END=20) IS1,V1,IS2,V2
      IF (VAR .EQ. 2.) READ(3,FORM2,END=20) IS1,V1,V2,IS2,V3,V4
      IF (VAR .EQ. 3.) READ(3,FORM2,END=20) IS1,V1,V2,V3,IS2,V4,V5,V6
      IF (V1 .EQ. VMISS .OR. V2 .EQ. VMISS) GO TO 15
      IF (VAR .GE. 2. .AND. V3 .EQ. VMISS .OR. V4 .EQ. VMISS) GO TO 15
      IF (VAR .EQ. 3. .AND. V5 .EQ. VMISS .OR. V6 .EQ. VMISS) GO TO 15
      KT1 = KT1 + 1
      IF (IS1 .EQ. 'F' .OR. IS1 .EQ. '2') SEX1 = 'F'
      IF (IS2 .EQ. 'F' .OR. IS2 .EQ. '2') SEX2 = 'F'
      WRITE(4,40) KT1
      IF (VAR .EQ. 1.) WRITE(4,45) SEX1,KT1,V1
      IF (VAR .EQ. 1.) WRITE(4,46) SEX2,KT1,V2
      IF (VAR .EQ. 2.) WRITE(4,45) SEX1,KT1,V1,V2
      IF (VAR .EQ. 2.) WRITE(4,46) SEX2,KT1,V3,V4
      IF (VAR .EQ. 3.) WRITE(4,45) SEX1,KT1,V1,V2,V3
      IF (VAR .EQ. 3.) WRITE(4,46) SEX2,KT1,V4,V5,V6
15    CONTINUE
20    CONTINUE
      DO 25 J = 1,INT(X2)
      SEX1 = 'M'
      SEX2 = 'M'
      IF (VAR .EQ. 1.) READ(3,FORM2,END=30) IS1,V1,IS2,V2
      IF (VAR .EQ. 2.) READ(3,FORM2,END=30) IS1,V1,V2,IS2,V3,V4
      IF (VAR .EQ. 3.) READ(3,FORM2,END=30) IS1,V1,V2,V3,IS2,V4,V5,V6
      IF (V1 .EQ. VMISS .OR. V2 .EQ. VMISS) GO TO 25
      IF (VAR .GE. 2. .AND. V3 .EQ. VMISS .OR. V4 .EQ. VMISS) GO TO 25
      IF (VAR .EQ. 3. .AND. V5 .EQ. VMISS .OR. V6 .EQ. VMISS) GO TO 25
      KT1 = KT1 + 1
      IF (IS1 .EQ. 'F' .OR. IS1 .EQ. '2') SEX1 = 'F'
      IF (IS2 .EQ. 'F' .OR. IS2 .EQ. '2') SEX2 = 'F'
      WRITE(4,40) KT1
      IF (VAR .EQ. 1.) WRITE(4,47) SEX1,V1
      IF (VAR .EQ. 1.) WRITE(4,48) SEX2,V2
      IF (VAR .EQ. 2.) WRITE(4,47) SEX1,V1,V2
      IF (VAR .EQ. 2.) WRITE(4,48) SEX2,V3,V4
      IF (VAR .EQ. 3.) WRITE(4,47) SEX1,V1,V2,V3
      IF (VAR .EQ. 3.) WRITE(4,48) SEX2,V4,V5,V6
25    CONTINUE
30    CONTINUE
40    FORMAT('  0  4',1X,I6/'DAD',T22,'M'/'MOM',T22,'F')
45    FORMAT('TWIN1',T8,'DAD',T15,'MOM',T22,A1,T29,I3,T36,3(1X,F8.3))
46    FORMAT('TWIN2',T8,'DAD',T15,'MOM',T22,A1,T29,I3,T36,3(1X,F8.3))
47    FORMAT('TWIN1',T8,'DAD',T15,'MOM',T22,A1,T36,3(1X,F8.3))
48    FORMAT('TWIN2',T8,'DAD',T15,'MOM',T22,A1,T36,3(1X,F8.3))
cc50    CONTINUE
      CLOSE(UNIT=3)
      CLOSE(UNIT=4)
c     WRITE(*,*) 'THE OUTPUT PEDIGREE FILE IS CREATED.'
      RETURN
      END
