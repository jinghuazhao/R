      SUBROUTINE RHINIT(XALL, NPARM)
C---------------------------------------------------------------------+
C                                                                     !
C RHINIT --  READS PATH COEFFICIENTS FROM XALL, AND                   !
C            CALCULATES PRELIMINARY TERMS COMMON TO ALL LATER         !
C            CALCULATIONS OF EXPECTED CORRELATIONS                    !
C                                                                     !
C INPUT:  XALL()                                                      !
C OUTPUT: NPARM -- NUMBER OF PARAMETERS READ FROM XALL()              !
C                                                                     !
C---------------------------------------------------------------------+

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C   DECLARATIONS FOR PATH COEFFICIENTS:
      DIMENSION XALL(*)
      DOUBLE PRECISION  M,II,JJ

C   DECLARATIONS FOR PRELIMINARY TERMS:
      DOUBLE PRECISION  IV,JW,IVV,JWW

C   DECLARATIONS FOR CONSTRAINTS:
      DIMENSION CC(*)
      common /parms/ h,c,y,z,m,u,p,ff,fm,b,x,ii,v,jj,w,t,rg,rc,r,
     +               hz,cy,iv,jw,s,a,ph,ga
C
C GET THE PARAMETERS:
C
      H = XALL( 1)
      C = XALL( 2)
      Y = XALL( 3)
      Z = XALL( 4)
      M = XALL( 5)
      U = XALL( 6)
      P = XALL( 7)
      FF= XALL( 8)
      FM= XALL( 9)
      B = XALL(10)
      X = 1.0
      II= XALL(11)
      V = XALL(12)
      JJ= XALL(13)
      W = XALL(14)
      T = 0.0
      RG= 1.0
      RC= 1.0
      R = 0.0

      NPARM = 14
C
C DEFINE PRELIMINARY TERMS WHICH ARE USED THROUGHOUT IN THE
C CALCULATION OF EXPECTED CORRELATIONS.
C
      HZ  = H*Z
      CY  = C*Y
      IV  = II*V
      JW  = JJ*W
      IF (M*U .GT. 0.0) THEN
         S = SQRT(M*U)
      ELSE
         S = 0.0
      END IF
C
C DEFINE PARAMETER A
C
      A = 0.0
      IF (2.0-FF-FM .NE. 0.0) THEN
         ALA = (FF+FM) / (2.0-FF-FM)
      ELSE
         ALA = 0.0
      END IF
      IF (P .EQ. 0.0) THEN
         A = S * ALA
      ELSE
         PHZCY = P*HZ*CY
         IF (ALA*PHZCY .NE. 0.0) THEN
            AAA = 1.0
            BBB = (ALA*P*(HZ*HZ + CY*CY) - 1.0) / (ALA*PHZCY)
            CCC = 1.0 + S/PHZCY
            DISC = BBB*BBB - 4.0*AAA*CCC
            IF (DISC .GT. 0.0) THEN
               A = (-BBB - SQRT(DISC)) / (2.0*AAA)
            ELSE
               A = (-BBB) / (2.0*AAA)
            END IF
         END IF
      END IF
C
C DEFINE TERMS USING A
C
      PH  = CY + HZ*A
      GA  = HZ + CY*A

      RETURN

C---------------------------------------------------------------------+
C                                                                     !
C RHOEXP --  CALCULATES GIVEN EXPECTED CORRELATION.                   !
C            ALL CALCULATIONS HERE USE THE PRELIMINARY TERMS          !
C            OBTAINED IN THE INITIZATION, ABOVE                       !
C                                                                     !
C INPUT:  IEQNO -- EQUATION NUMBER, AS DEFINED IN EQUATION FILE       !
C         SPEF  -- SPECIFIC EFFECT: =1.0 IF PRESENT, =0.0 OTHERWISE   !
C                                                                     !
C OUTPUT: RHO                                                         !
C                                                                     !
C---------------------------------------------------------------------+
      ENTRY RHOEXP(RHO,IEQNO,SPEF)
C
C DEFINE TERMS WHICH DEPEND ON SEX AND GENERATION
C
      IGEN = MOD( IEQNO, 10 )
      IF (IGEN .EQ. 0) THEN
         HZZ = H
         CYY = C
         RGG = RG
         RCC = RC
         IVV = II
         JWW = JJ
         ARG = A*RG*RC
         RHOGG = .5 * (1.0 + RG*RG*(M+P*GA*GA))
         RHOCC = B*B + FF*FF + FM*FM + 2.0*FF*FM*RC*RC*(U+P*PH*PH)
      ELSE IF (IGEN .EQ. 1) THEN
         HZZ = HZ
         CYY = CY
         RGG = 1.0
         RCC = 1.0
         IVV = IV
         JWW = JW
         ARG = A
         RHOGG = .5 * (1.0 + M + P*GA*GA)
         RHOCC = B*B*X*X + FF*FF + FM*FM + 2.0*FF*FM*(U+P*PH*PH)
      ELSE
         call intpr('RHOEXP: INVALID GEN',-1,0,0)
         STOP
      END IF
C
C     (PF,PM) #100
C
      IF (IEQNO .EQ. 100) THEN
         RHO = P + M*HZ*HZ + U*CY*CY + 2.*HZ*CY*S
C
C     (PF,PC) #210
C
      ELSE IF (IEQNO .EQ. 210) THEN
         F1  = FF
         F2  = FM
         RHO = HZZ*RGG*(GA*(1.0+P)+HZ*M+CY*S)*.5 +
     @         CYY*RCC*(F1*PH+F2*(P*PH+HZ*S+CY*U))
C
C     (PM,PC) #220
C
      ELSE IF (IEQNO .EQ. 220) THEN
         F1  = FM
         F2  = FF
         RHO = HZZ*RGG*(GA*(1.0+P)+HZ*M+CY*S)*.5 +
     @         CYY*RCC*(F1*PH+F2*(P*PH+HZ*S+CY*U))
C
C     (PC,PC) #300
C
      ELSE IF (IEQNO .EQ. 300) THEN
         RHO = HZZ*HZZ*RHOGG + CYY*CYY*RHOCC + 2.0*HZZ*CYY*ARG
C
C     (PC,IC) #500
C
      ELSE IF (IEQNO .EQ. 500) THEN
         RHO = IVV*(CYY+HZZ*ARG) + JWW*(HZZ+CYY*ARG)
C
C     (PF,IF)=(PM,IM) #501
C
      ELSE IF (IEQNO .EQ. 501) THEN
         RHO = IV*(CY+HZ*A) + JW*(HZ+CY*A)
C
C     (PF,IM)=(PM,IF) #600
C
      ELSE IF (IEQNO .EQ. 600) THEN
         RHO = IV*(P*PH+HZ*S+CY*U) + JW*(P*GA+HZ*M+CY*S)
C
C     (PF,IC) #710
C
      ELSE IF (IEQNO .EQ. 710) THEN
         F1  = FF
         F2  = FM
         RHO = IVV*RCC*(F1*PH+F2*(P*PH+HZ*S+CY*U)) +
     @         JWW*RGG*(GA*(1.0+P)+HZ*M+CY*S)*.5
C
C     (PM,IC) #720
C
      ELSE IF (IEQNO .EQ. 720) THEN
         F1  = FM
         F2  = FF
         RHO = IVV*RCC*(F1*PH+F2*(P*PH+HZ*S+CY*U)) +
     @         JWW*RGG*(GA*(1.0+P)+HZ*M+CY*S)*.5
C
C     (IF,PC) #810
C
      ELSE IF (IEQNO .EQ. 810) THEN
         F1  = FF
         F2  = FM
         RHO = HZZ*RGG*.5*(IV*(A+S+GA*P*PH)+JW*(1.0+M+P*GA*GA)) +
     @          CYY*RCC*F1*(IV+A*JW) +
     @          CYY*RCC*F2*(IV*(U+P*PH*PH)+JW*(S+GA*P*PH))
C
C     (IM,PC) #820
C
      ELSE IF (IEQNO .EQ. 820) THEN
         F1  = FM
         F2  = FF
         RHO = HZZ*RGG*.5*(IV*(A+S+GA*P*PH)+JW*(1.0+M+P*GA*GA)) +
     @          CYY*RCC*F1*(IV+A*JW) +
     @          CYY*RCC*F2*(IV*(U+P*PH*PH)+JW*(S+GA*P*PH))
C
C     (PC1,IC2) #900
C
      ELSE IF (IEQNO .EQ. 900) THEN
         RHO = IVV*(CYY*RHOCC+HZZ*ARG) + JWW*(HZZ*RHOGG+CYY*ARG)
C
C     (IF,IM) #1000
C
      ELSE IF (IEQNO .EQ. 1000) THEN
         RHO = IV*IV*(U+P*PH*PH) + JW*JW*(M+P*GA*GA) +
     @          2.0*IV*JW*(S+GA*P*PH)
C
C     (IF,IC) #1110
C
      ELSE IF (IEQNO .EQ. 1110) THEN
         F1  = FF
         F2  = FM
         RHO = IVV*RCC*F1*(IV+A*JW) +
     @          IVV*RCC*F2*(IV*(U+P*PH*PH) + JW*(S+GA*P*PH)) +
     @          JWW*RGG*.5*(IV*(A+S+GA*P*PH) + JW*(1.0+M+P*GA*GA))
C
C     (IM,IC) #1120
C
      ELSE IF (IEQNO .EQ. 1120) THEN
         F1  = FM
         F2  = FF
         RHO = IVV*RCC*F1*(IV+A*JW) +
     @          IVV*RCC*F2*(IV*(U+P*PH*PH) + JW*(S+GA*P*PH)) +
     @          JWW*RGG*.5*(IV*(A+S+GA*P*PH) + JW*(1.0+M+P*GA*GA))
C
C     (IC,IC) #1200
C
      ELSE IF (IEQNO .EQ. 1200) THEN
         RHO = IVV*IVV*RHOCC + JWW*JWW*RHOGG + 2.0*IVV*JWW*ARG
C
C TWIN CORRELATIONS
C
C        MZ TWINS TOGETHER #2300
C
      ELSE IF (IEQNO .EQ. 2300) THEN
         RHO = HZZ*HZZ + CYY*CYY*RHOCC + 2.*HZZ*CYY*A
C
C        DZ TWINS TOGETHER #2500
C
      ELSE IF (IEQNO .EQ. 2500) THEN
         RHO = HZZ*HZZ*RHOGG + CYY*CYY*RHOCC + 2.*HZZ*CYY*A

      ELSE
         WRITE (3,1399) IEQNO
 1399    FORMAT (/' ** ILLEGAL EQUATION NUMBER =',I6,' **')
         call intpr('RHOEXP: INTERNAL ERROR EQNO',-1,0,0)
         STOP
      END IF
C
C INCORPORATE SPECIFIC EFFECT
C
      RHO = RHO + SPEF*T

      RETURN

C---------------------------------------------------------------------+
C                                                                     !
C CONSTR --  DEFINE THE CONTRAINTS FOR THIS FUNCTION                  !
C            USING TERMS CALCULATED IN INITIALIZATION                 !
C INPUT: NONE       OUTPUT: CC()                                      !
C                                                                     !
C---------------------------------------------------------------------+
      ENTRY CONSTR(CC)

      CC(1) = 1. - (H*H + C*C + 2.*H*C*A*(RG*RC))
      CC(2) = 1. - (H*H*Z*Z + C*C*Y*Y + 2.*H*Z*C*Y*A)
      CC(3) = 1. - (B*B + FF*FF + FM*FM + 2.*FF*FM*RC*RC*(U+P*PH*PH))
      CC(4) = 1. - (B*B*X*X + FF*FF + FM*FM+2.*FF*FM*(U+P*PH*PH))
      CC(5) = 1. - (II*II + JJ*JJ + 2.*II*JJ*A*(RG*RC))
      CC(6) = 1. - (II*II*V*V + JJ*JJ*W*W + 2.*II*V*JJ*W*A)
      CC(7) = 1. - (U + P*PH*PH)
      CC(8) = 1. - (M + P*GA*GA)
      CC(9) = 1. - (S + GA*P*PH)
      CC(10)= M * U

      RETURN
      END
