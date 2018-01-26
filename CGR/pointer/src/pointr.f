c                                                                       
c         
c *******************************************************************   
c         
c                                                                       
c         
      subroutine auxout(iq)
c -------------------------------------------------------------------   
c         
c                                                                       
c         
c GET ADDITIONAL OUTPUT IF IT()                                         
c         
c                                                                       
c         
c -------------------------------------------------------------------   
c         
      implicit double precision (o-z, a-h)
      integer idebug
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, nev, idebug, np
      common /auxp1/ p(3), g(3), clm2, clm1, cxp2, cxp1(3), rcdg, 
     &cdp, z(10), vv, uu, gg, ca, ck, ea, ek, scals(4), sclq, rsqw
     &, notrs, ng1, nox
      common /odds/ pis(10), fi(10), ititle(27), itypgr, irap, npi, nfi
c GET P(AF/G,Z)                                                         
c         
      dimension pafzg(10, 3), pgaf(2, 3)
      do 10 i = 1, nia
      if (nfi .eq. 0) fi(i) = 1.0
      do 10 k = 1, ng1
      x = (z(i) - g(k)) * rcdg
      ax = abs(x)
      t = 1.0 / (1.0 + (.2316419 * ax))
      d = .3989422804 * exp(- ((.5 * x) * x))
c                                                                       
c         
      t1 = (d * t) * ((((((((1.330274 * t) - 1.821256) * t) + 1.781478)
     & * t) - 0.3565638) * t) + 0.3193815)
      if (x) 1, 2, 2
    1 t1 = 1 - t1
    2 qfn = t1
      pafzg(i,k) = qfn
c-----------------------------------------------------------------------
c-------- 
c PRINT RELEVANT OUTPUT                                                 
c         
   10 continue
      write(unit=iq, fmt=2171) 
 2171 format(44x,21hP(AFFECTION/GENOTYPE)//27x,9hINCIDENCE,3x,
     &9hTHRESHOLD,8x,3h(1),9x,3h(2),9x,3h(3))
      do 30 i = 1, nia
   30 write(unit=iq, fmt=2172) yi(i), z(i), (pafzg(i,k),k = 1, ng1)
 2172 format(24x,5f12.5)
      write(unit=iq, fmt=2173) 
c---------------------------------                                      
c         
c GET P(G/AF)                                                           
c         
c ------------------------------                                        
c         
c FIRST, FOR EACH LIABILITY CLASS                                       
c         
 2173 format(/44x,22hP(GENOTYPE/AFF.STATUS))
      do 100 i = 1, nia
      sua = 0.0
      sun = 0.0
      do 50 k = 1, ng1
      pgaf(1,k) = pafzg(i,k) * p(k)
      pgaf(2,k) = (1. - pafzg(i,k)) * p(k)
      sua = sua + pgaf(1,k)
   50 sun = sun + pgaf(2,k)
      do 55 k = 1, ng1
      pgaf(1,k) = pgaf(1,k) / sua
   55 pgaf(2,k) = pgaf(2,k) / sun
      write(unit=iq, fmt=2174) yi(i), z(i), (pgaf(1,k),k = 1, ng1)
      write(unit=iq, fmt=2175) (pgaf(2,k),k = 1, ng1)
 2174 format(14x,6hP(G/A),4x,5f12.5)
 2175 format(14x,6hP(G/N),28x,3f12.5)
c ------------------------------                                        
c         
c OVER LIABILITY CLASSES                                                
c         
  100 continue
      sua = 0.0
      sun = 0.0
      do 20 k = 1, ng1
      pgaf(1,k) = 0.0
      pgaf(2,k) = 0.0
      do 15 i = 1, nia
      pgaf(1,k) = pgaf(1,k) + (fi(i) * pafzg(i,k))
      pgaf(2,k) = pgaf(2,k) + (fi(i) * (1. - pafzg(i,k)))
   15 continue
      pgaf(1,k) = pgaf(1,k) * p(k)
      pgaf(2,k) = pgaf(2,k) * p(k)
      sua = sua + pgaf(1,k)
   20 sun = sun + pgaf(2,k)
      do 25 k = 1, ng1
      pgaf(1,k) = pgaf(1,k) / sua
   25 pgaf(2,k) = pgaf(2,k) / sun
      write(unit=iq, fmt=2177) (pgaf(1,k),k = 1, ng1)
 2177 format(/14x,6hP(G/A),8x,12hOVER CLASSES,8x,3f12.5)
      write(unit=iq, fmt=2178) (pgaf(2,k),k = 1, ng1)
 2178 format(14x,6hP(G/N),28x,3f12.5)
      write(unit=iq, fmt=2179) (p(k),k = 1, ng1)
 2179 format(14x,6hP(G/U),28x,3f12.5)
      return 
      end
c                                                                       
c         
c                                                                       
c         
c **********************************************************************
c**       
c                                                                       
c         
      subroutine auxpar(xall, nall, itp)
c ----------------------------------------------------------------------
c--       
c                                                                       
c         
c MODIFIED TO CORRECT FOR PREVIOUSLY NEGLECTED POST-SELECTION           
c         
c                                                                       
c         
c GET AUXILIARY PARAMETERS                                              
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c--       
      implicit double precision (o-z, a-h)
      implicit integer (i-n)
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, nev, idebug, np
      common /auxp1/ p(3), g(3), clm2, clm1, cxp2, cxp1(3), rcdg, 
     &cdp, z(10), vv, uu, gg, ca, ck, ea, ek, scals(4), sclq, rsqw
     &, notrs, ng1, nox
      common /intg/ xuu(20), auu(20), x(20), a(20), fzdqz(10), xdev, 
     &tvuu, fcaea, igot, nuu, nl, nm
      common /auxp2/ ntrel, mxrel, clm2s(3), clm1s(3), rrrs(2, 3, 5), 
     &cxc2s(2, 2), cxc1s(2, 2), cdcs(2, 2), cxd2s(2, 2, 5), cxd1s(2, 2, 
     &5), cdds(2, 2, 5), cxa2s(2, 5), cxa1s(2, 5), cdas(2, 5)
      common /trsm/ tr(18), tm(6, 3, 70)
      common /hang/ iwrong
      dimension tpo(3, 3), top(3, 3), s(3, 3)
      dimension tpor(3, 3), sr(3, 3), tmr(6, 3), rf(3), anorm(6)
c FOR NON-MENDELIAN TRANSMISSION                                        
c         
      dimension xall(15), itp(15)
      dimension tau(3), ut(3), utt(3)
c PGC( ) : GENOTYPE PROBABILITIES AT CONCEPTION                         
c         
      dimension trg(3, 2), trm(2, 2), trz(2, 3)
c **** NOTE **** ALL TRANSITION MATRICES ARE ROW-STOKASTIC, THAT IS     
c         
c P(J/I)=P(I,J)                                                         
c         
c --------------------------------------------------------------------  
c         
c INTERPRET PARAMETERS                                                  
c         
      dimension pgc(3)
      vv = xall(1)
      uu = xall(2)
      dd = xall(3)
c DESCALE Q                                                             
c         
      tt = xall(4) * sqrt(sclq)
      qq = xall(5) / sclq
c *** CA ***                                                            
c         
      ck = xall(6)
      ca = xall(7)
      if (itp(7) .eq. 1) goto 4
      ca = ck * zinp
      xall(7) = ca
    4 continue
      sp = xall(8)
      ww = xall(9)
      u = sp * qq
      tau(1) = xall(10)
      tau(2) = xall(11)
      tau(3) = xall(12)
c                                                                       
c         
c ALLOWING FOR NO TRANSMISSION                                          
c         
      u = sp * qq
      if (notrs .eq. 0) goto 668
      if (notrs .eq. 1) tau(1) = 1. - qq
      tau(2) = tau(1)
      tau(3) = tau(1)
      xall(10) = tau(1)
      xall(11) = tau(1)
      xall(12) = tau(1)
  668 continue
c     WRITE(6,667)(XALL(I),I=1,NALL),NIA                                
c         
c 667 FORMAT(1X,'XALL',10F8.5,I2)                                       
c         
c ----------------------------------------------------------------------
c-------- 
c CHARACTERISTICS OF MAJOR LOCUS                                        
c         
      rr = xall(13)
c     -----                                                             
c         
      ng1 = 3
c SP IS X PROPORTION OF SPORADICS                                       
c         
      gg = 0.0
      if (qq) 1, 1, 2
    1 g(1) = uu
      p(1) = 1.0
      pgc(1) = 1.0
      ng1 = 1
      goto 5
    2 continue
      t1 = 1.0 - qq
      t2 = (1. - u) * (1. - u)
      t3 = 1. - sp
      p(1) = (((qq * qq) * t3) * t3) / t2
      p(2) = (((2. * qq) * t1) * t3) / t2
      p(3) = (t1 * t1) / t2
      pgc(1) = qq * qq
      pgc(2) = (2. * qq) * t1
      pgc(3) = t1 * t1
      zu = (uu - ((qq * qq) * tt)) - ((((2. * qq) * t1) * tt) * dd)
      g(1) = zu + tt
      g(2) = zu + (tt * dd)
      g(3) = zu
      do 3 i = 1, 3
    3 gg = gg + ((pgc(i) * (g(i) - uu)) * (g(i) - uu))
c ----------------------------------------------------------------------
c-------- 
c VARIANCE COMPONENTS AND OTHER VARIANCES                               
c         
    5 continue
c REQUIRED FOR TRANSLATION                                              
c         
c VARIABLE LINEAR TRANSFORM                                             
c         
c     IF(IGOT)6,6,7                                                     
c         
      twopi = 6.28318530718
    6 tvuu = uu
      fcaea = ca / (vv - gg)
      if (qq .gt. 1.e-08) fcaea = 0.5 * fcaea
      if ((vv - gg) .eq. 0.) iwrong = 1
    7 continue
      ek = (vv - gg) - ck
c     WRITE(6,6337)GG,VV,CK,CA,EK,EA                                    
c         
c6337 FORMAT(1X,'GG,VV,CK,CA,EK,EA=  ',6F8.5)                           
c         
      ea = (vv - gg) - ca
      cdg = sqrt((vv - gg) + ww)
      if (cdg .eq. 0.) iwrong = 2
c SKIP SOME CONSTANT CALCULATIONS WHEN H=0.0                            
c         
      rcdg = 1.0 / cdg
      rsqz = 1.0
      clm2s(1) = 1.0
      clm2s(2) = 1.0
      clm2s(3) = 1.0
      clm1s(1) = 1.0
      clm1s(2) = 1.0
      clm1s(3) = 1.0
      scals(1) = 1.0
      scals(2) = 1.0
      scals(3) = 1.0
      scals(4) = 1.
      if (ck .eq. 0.0) goto 8
      rsqz = sqrt(ck / ca)
      if (ca .eq. 0.) iwrong = 3
      clm2s(1) = - (1. / (2. * ca))
      clm2s(2) = clm2s(1)
      clm2s(3) = - (1. / ca)
      clm1s(1) = 1. / sqrt(twopi * ca)
      clm1s(2) = clm1s(1)
c REQUIRED FOR SCALING (1,2,3,4: IPH=0,1,2,3)                           
c         
      clm1s(3) = 1. / sqrt((twopi * .5) * ca)
      scals(1) = (2.39 * exp(- (7.451 * ca))) + 1.179
      scals(2) = scals(1)
      scals(3) = 0.70 / sqrt((rr * ca) * (vv - (rr * ca)))
c                                                                       
c         
      scals(4) = 0.70 / sqrt(ca * (vv - ca))
    8 continue
      rsqw = 0.
      if (ww .gt. 1.e-05) rsqw = 1. / sqrt(ww)
      cxp2 = - (1. / (2. * ea))
      if (ea .eq. 0.) iwrong = 4
      temp = sqrt(twopi * ea)
      cxp1(1) = p(1) / temp
      cxp1(2) = p(2) / temp
      cxp1(3) = p(3) / temp
c CONSTANTS FOR CHILDREN ARE FUNCTION OF GENERATION AND TYPE OF INTEGRAL
c         
      cdp = 1. / sqrt(ea + ww)
      omr2 = 1. - (rr * rr)
      om2r2 = 1. - ((2. * rr) * rr)
      cxc2s(2,1) = - (1. / (2. * ((ck * omr2) + ek)))
      cxc1s(2,1) = 1. / sqrt(twopi * ((ck * omr2) + ek))
      cdcs(2,1) = 1. / sqrt(((ck * omr2) + ek) + ww)
      cxc2s(1,1) = - (1. / (2. * ((ca * omr2) + ea)))
      cxc1s(1,1) = 1. / sqrt(twopi * ((ca * omr2) + ea))
      cdcs(1,1) = 1. / sqrt(((ca * omr2) + ea) + ww)
      cxc2s(2,2) = - (1. / (2. * ((ck * om2r2) + ek)))
      cxc1s(2,2) = 1. / sqrt(twopi * ((ck * om2r2) + ek))
      cdcs(2,2) = 1. / sqrt(((ck * om2r2) + ek) + ww)
      cxc2s(1,2) = - (1. / (2. * ((ca * om2r2) + ea)))
      cxc1s(1,2) = 1. / sqrt(twopi * ((ca * om2r2) + ea))
c CONSTANTS FOR DESCENDENT POINTERS ARE FUNSTION OF GENERATION, TYPE OF 
c         
c INTEGRATION AND DEGREE OF RELATIONSHIP                                
c         
c CONSTANTS FOR ANCESTRAL POINTERS ARE FUNCTION OF GENERATION AND DEGREE
c         
c OF RELATIONSHIP                                                       
c         
      cdcs(1,2) = 1. / sqrt(((ca * om2r2) + ea) + ww)
      do 12 irel = 1, 5
c RRRS VARIABLES FUNCTION OF GENERATION, TYPE OF INTEGRATION, AND DEGREE
c         
c OF RELATIONSHIP                                                       
c         
      rel = rr ** irel
      rrrs(1,1,irel) = rel
      rrrs(1,2,irel) = rel
      rrrs(1,3,irel) = 2. * rel
      rrrs(2,1,irel) = rel * rsqz
      rrrs(2,2,irel) = rel * rsqz
      rrrs(2,3,irel) = (2. * rel) * rsqz
      frel = 1. - (rel * rel)
      ftrel = 1. - ((2. * rel) * rel)
      cxd2s(2,1,irel) = - (1. / (2. * ((frel * ck) + ek)))
      cxd1s(2,1,irel) = 1. / sqrt(twopi * ((frel * ck) + ek))
      cdds(2,1,irel) = 1. / sqrt(((frel * ck) + ek) + ww)
      cxd2s(1,1,irel) = - (1. / (2. * ((frel * ca) + ea)))
      cxd1s(1,1,irel) = 1. / sqrt(twopi * ((frel * ca) + ea))
      cdds(1,1,irel) = 1. / sqrt(((frel * ca) + ea) + ww)
      cxd2s(2,2,irel) = - (1. / (2. * ((ftrel * ck) + ek)))
      cxd1s(2,2,irel) = 1. / sqrt(twopi * ((ftrel * ck) + ek))
      cdds(2,2,irel) = 1. / sqrt(((ftrel * ck) + ek) + ww)
      cxd2s(1,2,irel) = - (1. / (2. * ((ftrel * ca) + ea)))
      cxd1s(1,2,irel) = 1. / sqrt(twopi * ((ftrel * ca) + ea))
      cdds(1,2,irel) = 1. / sqrt(((ftrel * ca) + ea) + ww)
      cxa2s(2,irel) = - (1. / (2. * ((frel * ck) + ek)))
      cxa1s(2,irel) = 1. / sqrt(twopi * ((frel * ck) + ek))
      cdas(2,irel) = 1. / sqrt(((frel * ck) + ek) + ww)
      cxa2s(1,irel) = - (1. / (2. * ((frel * ca) + ea)))
      cxa1s(1,irel) = 1. / sqrt(twopi * ((frel * ca) + ea))
      cdas(1,irel) = 1. / sqrt(((frel * ca) + ea) + ww)
c ----------------------------------------------------------------------
c-------- 
c COMPUTE THRESHOLDS                                                    
c         
   12 continue
      if (jli .eq. 0) goto 2220
      if (jli .eq. 3) z(nia + 1) = zi(1)
      do 20 i = 1, nia
      niter = 0
      if (jli .ne. 2) goto 800
      z(i) = zi(i)
      goto 20
  800 aff = 0.0
      zai = yi(i)
      xa = sqrt(vv + ww)
   10 de = 0.0
      af = 0.0
c RECALL THAT WHEN MAJOR GENOTYPE PRESENT, (3) IS NORMAL GENOTYPE       
c         
      do 15 k = 1, ng1
      tempa = xa - g(k)
      tempa = tempa * rcdg
      af = af + (pgc(k) * qf(tempa))
      de = de + ((pgc(k) * .398942280385) * exp(- ((.5 * tempa) * tempa)
     &))
   15 continue
      niter = niter + 1
      af = af - aff
      de = ((af - zai) / de) * cdg
      tt1 = abs(de) - (.00001 * cdg)
      if (tt1) 19, 119, 119
  119 if (niter .gt. 100) goto 1119
      tt1 = abs(de) - cdg
      if (tt1) 18, 18, 17
   17 sgn = abs(de) / de
      de = sgn * cdg
   18 xa = xa + de
      goto 10
 1119 write(unit=61, fmt=1919) i
 1919 format(1x,16h WARNING:  CLASS,i2,19h THRESHOLD LOOPING.)
   19 z(i) = xa
   20 continue
c ---------------------------------------------------------------       
c         
c GET F(Z)/Q(Z) FOR TRANSLATION IN NUMERICAL INTEGRATION                
c         
 2220 continue
      do 200 j = 1, nia
      tempa = z(j) * rcdg
      ax = abs(tempa)
      t = 1.0 / (1.0 + (.2316419 * ax))
      d = .3989422804 * exp(- ((.5 * tempa) * tempa))
c                                                                       
c         
      tail = (d * t) * ((((((((1.330274 * t) - 1.821256) * t) + 1.781478
     &) * t) - 0.3565638) * t) + 0.3193815)
      if (tempa) 213, 214, 214
  213 tail = 1 - tail
  214 qqq = tail
c ----------------------------------------------------------------------
c-        
c GET TRANSITION MATRICES FOR POINTERS                                  
c         
c --------------------------------------------------------------------  
c         
c THE CONDITIONAL PROBABILITY OF CHILD MAJOR GENOTYPE GIVEN MAJOR GENOTY
cPES      
c OF PARENTS, STORED IN TM( , ,1), AS IT IS THE T.M. FOR A CHILD POINTER
c         
c POINTING TO BOTH PARENTS.                                             
c         
c GET THIS FIRST AS UNIDIM. ARRAY FOR EXEC ROUTINE                      
c         
  200 fzdqz(j) = (.3989422804 * exp(- ((.5 * tempa) * tempa))) / qqq
      ntrel = 7
      ndrel = 5
c SKIP THIS SECTION NOT NEEDED WHEN NO MAJOR GENE OR WHEN Q AND X       
c         
c ARE NOT ITERATED ( IN THE LATTER CASE, CALCULATION DONE ONLY          
c         
c AT THE VERY FIRST CALL TO AUXPAR )                                    
c         
c                                                                       
c         
c IF TAU'S NOT MENDELIAN OR AT LEAST ONE OF THEM ITERATED,              
c         
c CALCULATE TRANSITION MATRICES BY SETTING IGOT TO ZERO AND             
c         
c SETTING U TO 1.E-10 TO HAVE NOX=1                                     
c         
      mxrel = ndrel * ntrel
      itau = 0
      if (((itp(10) .eq. 1) .or. (itp(11) .eq. 1)) .or. (itp(12) .eq. 1)
     &) itau = 1
      if ((((tau(1) .lt. 0.99) .or. (tau(2) .lt. 0.49)) .or. (tau(2)
     & .gt. 0.51)) .or. (tau(3) .gt. 0.01)) itau = 1
      if (itau .eq. 0) goto 5555
      igot = 0
      u = 1.e-10
c                                                                       
c         
 5555 continue
      if (qq .eq. 0.0) goto 100
      if (((itp(5) .eq. 0) .and. (itp(8) .eq. 0)) .and. (igot .eq. 1)) 
     &goto 100
   21 nox = 1
      do 2121 i = 1, 3
      ut(i) = (1. - u) * tau(i)
 2121 utt(i) = ((u * tau(i)) + 1.) - tau(i)
      tr(1) = utt(3) * utt(3)
      tr(2) = (2. * ut(3)) * utt(3)
      tr(3) = ut(3) * ut(3)
      tr(4) = utt(2) * utt(3)
      tr(5) = (ut(2) * utt(3)) + (ut(3) * utt(2))
      tr(6) = ut(2) * ut(3)
      tr(7) = utt(1) * utt(3)
      tr(8) = (ut(1) * utt(3)) + (ut(3) * utt(1))
      tr(9) = ut(1) * ut(3)
      tr(10) = utt(2) * utt(2)
      tr(11) = (2. * ut(2)) * utt(2)
      tr(12) = ut(2) * ut(2)
      tr(13) = utt(1) * utt(2)
      tr(14) = (ut(1) * utt(2)) + (ut(2) * utt(1))
      tr(15) = ut(1) * ut(2)
      tr(16) = utt(1) * utt(1)
      tr(17) = (2. * ut(1)) * utt(1)
c ---------------------                                                 
c         
c GET VECTOR OF RELATIVE FITNESS RF( )                                  
c         
      tr(18) = ut(1) * ut(1)
      rf(1) = (1. - sp) * (1. - sp)
      rf(2) = 1. - sp
c     WRITE(IP,4444)(RF(I),I=1,3)                                       
c         
c4444 FORMAT(/1X,6E11.4)                                                
c         
c GET PARENT TO OFFSPRING TRANSITION MATRIX                             
c         
c     TPO(1,1)=QQ                                                       
c         
c     TPO(1,2)=T1                                                       
c         
c     TPO(1,3)=0.0                                                      
c         
c     TPO(2,1)=QQ*TR(1)                                                 
c         
c     TPO(2,2)=T1*TR(1)+QQ*TR(2)                                        
c         
c     TPO(2,3)=T1*TR(2)                                                 
c         
c     TPO(3,1)=U*QQ                                                     
c         
c     TPO(3,2)=QQ*TR(4)+T1*U                                            
c         
c     TPO(3,3)=T1*TR(4)                                                 
c         
      rf(3) = 1.
      trg(1,1) = 1. - tau(3)
      trg(1,2) = tau(3)
      trg(2,1) = 1. - tau(2)
      trg(2,2) = tau(2)
      trg(3,1) = 1. - tau(1)
      trg(3,2) = tau(1)
      trm(1,1) = 1.0
      trm(1,2) = 0.0
      trm(2,1) = u
      trm(2,2) = 1. - u
      trz(1,1) = qq
      trz(1,2) = 1. - qq
      trz(1,3) = 0.0
      trz(2,1) = 0.0
      trz(2,2) = qq
      trz(2,3) = 1. - qq
      do 27 i = 1, 3
      do 27 m = 1, 3
      tpo(i,m) = 0.0
      do 27 l = 1, 2
      ss = 0.0
      do 26 k = 1, 2
   26 ss = ss + (trg(i,k) * trm(k,l))
c GET TPOR ( AFTER SELECTION )                                          
c         
   27 tpo(i,m) = tpo(i,m) + (ss * trz(l,m))
      do 24 i = 1, ng1
      anorm(i) = 0.0
      do 24 k = 1, ng1
c     WRITE(IP,4444)(ANORM(I),I=1,3)                                    
c         
   24 anorm(i) = anorm(i) + (rf(k) * tpo(i,k))
      do 25 i = 1, ng1
      do 25 j = 1, ng1
c GET OFFSPRING TO PARENT TRANSITION MATRIX                             
c         
c     TOP(1,1)=T3*T3*QQ/T2                                              
c         
c     TOP(1,2)=T3*T1*(1.+U)/T2                                          
c         
c     TOP(1,3)=SP*T1*T1/T2                                              
c         
c     TOP(2,1)=T3*T3*QQ*.5/T2                                           
c         
c     TOP(2,2)=T3*(1.+U-2.*U*QQ)*.5/T2                                  
c         
c     TOP(2,3)=T1*(1.+SP-2.*U)*.5/T2                                    
c         
c     TOP(3,1)=0.0                                                      
c         
c     TOP(3,2)=(QQ-U)/TR(4)                                             
c         
c     TOP(3,3)=T1/TR(4)                                                 
c         
   25 tpor(i,j) = (rf(j) * tpo(i,j)) / anorm(i)
      do 31 i = 1, ng1
      do 31 j = 1, ng1
c GET SIBLING TRANSITION MATRIX                                         
c         
   31 top(i,j) = (tpo(j,i) * p(j)) / pgc(i)
c     S(1,1)=(T3+QQ+U)*(T3+QQ+U)*.25                                    
c         
c     S(1,2)=T1*T4*(T3+QQ+U)                                            
c         
c     S(1,3)=T1*T1*T4*T4                                                
c         
c     S(2,1)=QQ*T4*(1.+QQ-SP+U)*.5                                      
c         
c     S(2,2)=1.-T4*(T1+QQ*QQ-U+QQ*U)                                    
c         
c     S(2,3)=(2.-QQ-U)*T1*T4*.5                                         
c         
c     S(3,1)=QQ*QQ*T4*T4                                                
c         
c     S(3,2)=QQ*T4*(2.-QQ-U)                                            
c         
c     S(3,3)=(2.-QQ-U)*(2.-QQ-U)*.25                                    
c         
      t4 = (1. + sp) * .5
      do 33 i = 1, ng1
      do 33 j = 1, ng1
      k = 0
      ss = 0.0
      do 333 l = 1, ng1
      do 333 m = l, ng1
      k = k + 1
      ilm = ((k - 1) * ng1) + i
      jlm = ((k - 1) * ng1) + j
      cons = 1.
      if (l .ne. m) cons = 2.
  333 ss = ss + ((((tr(ilm) * tr(jlm)) * cons) * p(l)) * p(m))
c GET SR ( AFTER SELECTION )                                            
c         
   33 s(i,j) = ss / pgc(i)
      do 28 i = 1, ng1
      anorm(i) = 0.0
      do 28 k = 1, ng1
c     WRITE(IP,4444)(ANORM(I),I=1,3)                                    
c         
   28 anorm(i) = anorm(i) + (rf(k) * s(i,k))
      do 29 i = 1, ng1
      do 29 j = 1, ng1
c --------------------                                                  
c         
c NTREL= NUMBER OF TYPE OF RELATIONSHIP ALLOWED                         
c         
c NDREL= NUMBER OF DEGREE OF RELATIONSHIP ALLOWED                       
c         
   29 sr(i,j) = (rf(j) * s(i,j)) / anorm(i)
      ng1s = (ng1 * (ng1 + 1)) / 2
      mxrel2 = 2 * mxrel
      do 39 i = 1, ng1
      do 39 j = 1, ng1s
      do 39 krel = 1, mxrel2
c --------------------                                                  
c         
c GET NOW TM( , ,1), P(G(CHILD)/G(FATHER),G(MOTHER))                    
c         
   39 tm(j,i,krel) = 0.0
      tm(1,1,1) = tr(1)
      tm(1,2,1) = tr(2)
      tm(1,3,1) = tr(3)
      tm(2,1,1) = tr(4)
      tm(2,2,1) = tr(5)
      tm(2,3,1) = tr(6)
      tm(3,1,1) = tr(7)
      tm(3,2,1) = tr(8)
      tm(3,3,1) = tr(9)
      tm(4,1,1) = tr(10)
      tm(4,2,1) = tr(11)
      tm(4,3,1) = tr(12)
      tm(5,1,1) = tr(13)
      tm(5,2,1) = tr(14)
      tm(5,3,1) = tr(15)
      tm(6,1,1) = tr(16)
      tm(6,2,1) = tr(17)
c GET TMR ( AFTER SELECTION )                                           
c         
      tm(6,3,1) = tr(18)
      do 34 i = 1, ng1s
      anorm(i) = 0.0
      do 34 k = 1, ng1
c     WRITE(IP,4444)(ANORM(I),I=1,6)                                    
c         
   34 anorm(i) = anorm(i) + (rf(k) * tm(i,k,1))
      do 35 i = 1, ng1s
      do 35 j = 1, ng1
c ----------------------------------------------------------------------
c         
c GET 1B = 1C , 1D AND 1E TRANSITION MATRICES                           
c         
   35 tmr(i,j) = (rf(j) * tm(i,j,1)) / anorm(i)
      do 40 i = 1, ng1
      do 40 j = 1, ng1
      tm(i,j,2) = s(i,j)
      tm(i,j,3) = s(i,j)
      tm(i,j,4) = tpo(i,j)
c --------------------                                                  
c         
c GET OTHER B , C , D , E TYPES OF T.M.                                 
c         
   40 tm(i,j,5) = top(i,j)
      kkb = 2
      kkc = kkb + 1
      kkd = kkc + 1
      kke = kkd + 1
      do 45 l = 1, 4
      kb = kkb + ntrel
      kc = kb + 1
      kd = kc + 1
      ke = kd + 1
      do 42 i = 1, ng1
      do 42 j = 1, ng1
      do 42 k = 1, ng1
      tm(i,j,kb) = tm(i,j,kb) + (top(i,k) * tm(k,j,kkb))
      tm(i,j,kc) = tm(i,j,kc) + (sr(i,k) * tm(k,j,kkd))
      tm(i,j,kd) = tm(i,j,kd) + (tpor(i,k) * tm(k,j,kkd))
   42 tm(i,j,ke) = tm(i,j,ke) + (top(i,k) * tm(k,j,kke))
      kkb = kkb + ntrel
      kkc = kkb + 1
      kkd = kkc + 1
c --------------------                                                  
c         
c GET A-TYPE TRANSITION MATRICES                                        
c         
   45 kke = kkd + 1
      ia2 = ntrel + 1
      ia3 = ia2 + ntrel
      ia4 = ia3 + ntrel
      ia5 = ia4 + ntrel
      i2c = ntrel + 3
      i2d = ntrel + 4
      i2e = ntrel + 5
      i3c = (2 * ntrel) + 3
      do 50 i = 1, ng1
      do 50 j = 1, ng1
      do 50 k = 1, ng1
      tm(i,j,ia2) = tm(i,j,ia2) + (top(i,k) * tpo(k,j))
      tm(i,j,ia3) = tm(i,j,ia3) + (top(i,k) * tm(k,j,i2c))
      tm(i,j,ia4) = tm(i,j,ia4) + (tm(i,k,i2e) * tm(k,j,i2d))
c --------------------                                                  
c         
c GET F- AND G-TYPE TRANSITION MATRICES                                 
c         
   50 tm(i,j,ia5) = tm(i,j,ia5) + (tm(i,k,i2e) * tm(k,j,i3c))
      ia2 = ntrel + 1
      id2 = ia2 + 3
      ie2 = id2 + 1
      ia3 = ntrel + ia2
      ic3 = ia3 + 2
      id3 = ic3 + 1
      ia4 = ntrel + ia3
      if3 = ia3 + 5
      if4 = ia4 + 5
      ig4 = if4 + 1
      if5 = ntrel + if4
      ig5 = if5 + 1
      do 58 i = 1, ng1
      do 58 j = 1, ng1
      do 58 k = 1, ng1
      tm(i,j,if3) = tm(i,j,if3) + (top(i,k) * tm(k,j,id2))
      tm(i,j,if4) = tm(i,j,if4) + (top(i,k) * tm(k,j,ic3))
      tm(i,j,ig4) = tm(i,j,ig4) + (top(i,k) * tm(k,j,ia3))
      tm(i,j,if5) = tm(i,j,if5) + (tm(i,k,ie2) * tm(k,j,id3))
c --------------------                                                  
c         
c FOR DESCENDENT POINTERS, GET P(G(POINTER)/MATING TYPE)                
c         
   58 tm(i,j,ig5) = tm(i,j,ig5) + (top(i,k) * tm(j,k,ia4))
      do 60 l = 2, mxrel
      kch = mxrel + l
      do 60 k = 1, ng1
      do 60 ij = 1, ng1s
c USE POST-SELECTION TRANSITION MATRIX FROM PARENTS TO CHILD            
c         
      do 60 ku = 1, ng1
c SINGLE SELECTION OF MULTIPLEX FAMILIES : PROBAND TREATED AS POINTER 1A
c         
c IN ADDITION GET APPROPRIATE 1C MATRIX FOR POINTER TO OFFSPRING        
c         
   60 tm(ij,k,kch) = tm(ij,k,kch) + (tmr(ij,ku) * tm(ku,k,l))
      mxr1 = mxrel + 1
      mxr3 = mxr1 + 2
      do 62 k = 1, ng1
      do 62 ij = 1, ng1s
      tm(ij,k,mxr1) = tm(ij,k,1)
      tm(ij,k,mxr3) = tm(ij,k,1)
   62 continue
      if (u .lt. 1.e-10) nox = 0
  100 continue
c --------------------------------------------------------------------  
c         
      igot = 1
      return 
      end
c                                                                       
c         
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         
c                                                                       
c         
      subroutine bldlt(ff, x, n, xall, nall, itp, ihess, icall, ivar, se
     &, ip, h)
c ----------------------------------------------------------------------
c         
c                                                                       
c         
c OBTAINS APPROXIMATE HESSIAN BY FINITE DIFFERENCE, FACTORIZES, AND FIND
c         
c VARIANCE-COVARIANCE MATRIX WHEN DESIRED                               
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c         
      implicit double precision (o-z, a-h)
      common /mini/ tl(15, 15), d(15), g(15), gsave(15), y(15), p(15), 
     &tmin, fsmf, idif, isw, iret
      common /secc/ seca, seck, b(15, 15), saveb(15, 15), ifail
      dimension x(15), xall(15), itp(15)
c ----------------------------------------------------------------------
c-        
c GET B-MATRIX BY FINITE DIFFERENCE IF ICALL=0 AND IHESS=1 OR IF ICALL=1
c         
      dimension s(15, 15), tu(15), se(15), xsave(15)
      if (icall) 1, 1, 3
    1 if (ihess - 1) 3, 3, 30
    3 continue
      do 5 i = 1, nall
      se(i) = 0.0
    5 xsave(i) = xall(i)
      hh = 10. * h
      do 10 i = 1, n
      temp = x(i)
      x(i) = x(i) + hh
      call fun(f, x, n, xall, nall, itp)
      tu(i) = f
      do 11 j = 1, i
      temp1 = x(j)
      x(j) = x(j) + hh
      call fun(f, x, n, xall, nall, itp)
      b(i,j) = f
      b(j,i) = b(i,j)
   11 x(j) = temp1
   10 x(i) = temp
      do 15 i = 1, nall
   15 xall(i) = xsave(i)
      do 20 i = 1, n
      do 20 j = 1, i
      b(i,j) = (((ff + b(i,j)) - tu(i)) - tu(j)) / (hh * hh)
      b(j,i) = b(i,j)
      saveb(i,j) = b(i,j)
   20 saveb(j,i) = b(i,j)
      write(unit=ip, fmt=2005) 
 2005 format(//10x,14hHESSIAN MATRIX)
      do 2006 i = 1, n
c ----------------------------------------------------------------------
c         
c CARRIES OUT CHOLESKY FACTORIZATION ON B-MATRIX                        
c         
 2006 write(unit=ip, fmt=1006) (b(i,j),j = 1, i)
   30 ifail = 0
      do 100 ic = 1, n
      if (b(ic,ic)) 200, 200, 22
   22 temp = b(ic,ic)
      d(ic) = b(ic,ic)
      tl(ic,ic) = 1.
      if (ic - n) 31, 150, 150
   31 ic1 = ic + 1
      do 50 k = ic1, n
   50 tl(k,ic) = b(k,ic) / temp
      do 100 i = ic1, n
      tlic = tl(i,ic)
      do 100 k = i, n
c ----------------------------------------------------------------------
c---      
c INDICATE WHETHER OR NOT B-MATRIX WAS FOUND TO BE POSITIVE DEFINITE    
c         
  100 b(k,i) = b(k,i) - ((tlic * tl(k,ic)) * temp)
  150 icall = icall - 1
      write(unit=ip, fmt=2000) 
 2000 format(//20x,39hFACTORIZATION OF B-MATRIX HAS SUCCEEDED//)
      goto 300
  200 icall = icall + 1
      write(unit=ip, fmt=2001) 
 2001 format(//20x,36hFACTORIZATION OF B-MATRIX HAS FAILED//)
      ifail = 1
  300 continue
c ----------------------------------------------------------------------
c---      
c GET TL-INVERSE IN S                                                   
c         
      if (icall) 500, 301, 500
  301 do 309 i = 1, n
      if (i - 1) 308, 308, 305
  305 im1 = i - 1
      do 307 k = 1, im1
      su = 0.0
      do 306 j = k, im1
  306 su = su + (tl(i,j) * s(j,k))
  307 s(i,k) = - su
  308 s(i,i) = 1.0
c FIND B-INVERSE                                                        
c         
  309 continue
      do 316 i = 1, n
      do 316 j = i, n
      su = 0.0
      do 315 k = j, n
  315 su = su + ((s(k,i) * s(k,j)) / d(k))
      b(i,j) = su
  316 b(j,i) = su
      write(unit=ip, fmt=1005) 
 1005 format(//10x,17hCOVARIANCE MATRIX)
      do 317 i = 1, n
  317 write(unit=ip, fmt=1006) (b(i,j),j = 1, i)
c ----------------------------------------------------------------------
c-        
 1006 format(1x,13e10.3)
c GET STANDARD ERRORS IF REQUIRED                                       
c         
      if (ivar - 1) 500, 350, 500
  350 k = 0
      do 351 i = 1, nall
      se(i) = 0.0
      if (itp(i) .ne. 1) goto 351
      k = k + 1
      se(i) = sqrt(b(k,k))
  351 continue
  500 continue
c GET CORRELATION MATRIX                                                
c         
      if ((icall .ne. 0) .or. (ivar .ne. 1)) goto 600
      do 501 i = 1, n
      do 501 j = i, n
      s(i,j) = b(i,j) / sqrt(b(i,i) * b(j,j))
  501 s(j,i) = s(i,j)
      write(unit=ip, fmt=1007) 
 1007 format(//10x,18hCORRELATION MATRIX)
      do 502 i = 1, n
  502 write(unit=ip, fmt=1006) (s(i,j),j = 1, i)
  600 continue
      return 
      end
c                                                                       
c         
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         
c                                                                       
c         
      subroutine chkbnd(t, tbnd, n, x, bnd, ip, h, xall, nall, itp, 
     &ihess)
c ----------------------------------------------------------------------
c--       
c                                                                       
c         
c CHECK BOUNDARY CONDITIONS ON PARAMETERS                               
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c--       
      implicit double precision (o-z, a-h)
      common /mini/ tl(15, 15), d(15), g(15), gsave(15), y(15), p(15), 
     &tmin, fsmf, idif, isw, iret
      common /bound/ ibndv
      dimension x(15), bnd(15, 2)
      dimension xall(15), itp(15)
      character*3 kname(13)
c ----------------------------------------------------------------------
c-----    
c CHECK OF CLOSENESS TO A BOUND IS DONE PRIOR TO DETERMINE TBOUND       
c         
      data kname / 'V  ', 'U  ', 'D  ', 'T  ', 'Q  ', 'H  ', 'Z  ', 
     &'X  ', 'W  ', 'T1 ', 'T2 ', 'T3 ', 'R  ' /
      data ik, ii, jk /3*0/
      clb = 1.e10
      do 20 j = 1, 2
      k = 0
      do 20 i = 1, nall
      if (itp(i) - 1) 20, 11, 20
   11 k = k + 1
      check = abs(x(k) - bnd(i,j))
      if (check - clb) 15, 15, 20
   15 clb = check
      ik = k
      ii = i
      jk = j
c IF ONE DESIRES TO AVOID FIXING A PARAMETER TO ITS BOUND, ONE MAY      
c         
c REPLACE 2.1*H BY H AT NEXT STATEMNET                                  
c         
   20 continue
c X(IK), CLOSE TO A BOUND, IS SET TO ITS BOUND                          
c         
c ITERATION WILL RESTART AT BEGINNING, WITH THAT PARAMETER DROPPED      
c         
c FROM THE ITERATED SET                                                 
c         
      if (clb - (.5 * h)) 21, 21, 30
   21 ihess = 0
      if (jk - 1) 210, 210, 211
  210 eh = (- h) - h
      goto 212
  211 eh = h + h
  212 x(ik) = bnd(ii,jk) + eh
      xall(ii) = x(ik)
      itp(ii) = 0
      if (((ii .eq. 4) .or. (ii .eq. 5)) .or. (ii .eq. 6)) ibndv = 1
      k = 0
      do 25 i = 1, nall
      if (itp(i) - 1) 25, 22, 25
   22 k = k + 1
      x(k) = xall(i)
   25 continue
      n = k
      write(unit=ip, fmt=1001) kname(ii)
 1001 format(/1x,10(1h*),1x,a3,1x,20hWAS SET TO A BOUND  ,10(1h*))
c ----------------------------------------------------------------------
c--       
c CHECK THAT INITIAL T-VALUE KEEPS SEARCH WITHIN PERMISSIBLE DOMAIN     
c         
c IF NOT, FIND MAXIMUM STEPSIZE T ADMISSIBLE                            
c         
      return 
   30 continue
      tbnd = 1.e05
      do 10 j = 1, 2
      k = 0
      do 10 i = 1, nall
      if (itp(i) - 1) 10, 12, 10
   12 k = k + 1
      teq = (bnd(i,j) - x(k)) / p(k)
      if (teq - tbnd) 1, 10, 10
    1 if (teq) 10, 10, 2
    2 tbnd = teq
c     WRITE(IP,3000)TBND,T                                              
c         
c3000 FORMAT(1X,'IN CHKBND AT STATEMENT 10 + 1, TBND,T=',2E10.3)        
c         
   10 continue
      if ((t * (2. + h)) - tbnd) 14, 13, 13
   13 t = tbnd * (.5 - h)
   14 call ineq(t, tbnd, n, x, h, xall, nall, itp)
      write(unit=ip, fmt=1000) tbnd, t
 1000 format(10x,5hTBND=,e12.5,10h  RESET T=,e12.5)
      return 
      end
c                                                                       
c         
c **********************************************************************
c******** 
c                                                                       
c         
      subroutine entry(ic, ictrl, nall, xall, itp, hdif, nn, tol, trupb)
      implicit double precision (o-z, a-h)
      integer reclen
      parameter (reclen = 80, nulls = -8388607)
      dimension xall(15), itp(15)
      dimension ia(270), aa(135), mc(9), jname(13), tt(3)
      dimension mf(30)
      character a*1, mc*1, tt*1, mf*1
      character iblank*3, itemp*3, jname*3
      character title*5
      integer idebug
c     PROGRAM TO INTERPRET CONTROL CARD FOR POINTER                     
c         
c ----------------------------------------------------------------------
c-------  
      equivalence (aa, ia), (itemp, tt), (a, title)
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, nev, idebug, np
      common /save/ prep(15), trat
      common /odds/ pis(10), fi(10), ititle(27), itypgr, irap, npi, nfi
      common /sy1a/ ibegin, i1, i2, icol, index, nchar
      common /sy1b/ a(80)
c                                                                       
c        
c ----------------------------                                          
c         
      data mc / 'P', 'I', 'R', 'C', 'A', 'T', 'A', 'C', ' ' /
      data jname / 'V  ', 'U  ', 'D  ', 'T  ', 'Q  ', 'H  ', 'Z  ', 
     &'X  ', 'W  ', 'T1 ', 'T2 ', 'T3 ', 'R  ' /
      data mf / 'C', 'N', 'F', 'T', 2*'S', 'A', 'P', 'S', 'P', 'L', 'P'
     &, 'L', 'P', 'F', 'C', 'A', 'M', 'R', 'I', 'D', 'F', 'O', 'X', 'R'
     &, 'I', 'T', 'Z', 2*'I' /
      write(unit=61, fmt=200) 
  200 format(//1x,130(1h=))
      write(unit=62, fmt=200) 
      iblank = '   '
   20 i1 = 1
      nchar = reclen
   25 read(unit=ic, fmt=26, end=140) (a(i),i = i1, nchar)
   26 format(80a1)
   27 if (nchar .ge. (10 * reclen)) goto 40
      do 404 i = nchar, i1, -1
      if (a(i) .eq. ' ') goto 404
      if (a(i) .eq. '*') goto 35
      goto 40
  404 continue
      goto 55
   35 a(i) = ' '
      i1 = nchar + 1
      nchar = nchar + reclen
      goto 25
   40 icol = 1
      call nsp1
      i1 = icol
c     SKIP ANY CONTROL CARD FROM NUCFAMA                                
c         
      icol = icol + 1
      if (title .eq. 'DEBUG') then
      idebug = 1
      goto 20
      end if
      if (title .eq. 'TITLE') goto 20
      do 44 i = 1, 15
      if ((a(i1) .eq. mf(i)) .and. (a(icol) .eq. mf(i + 15))) goto 20
   44 continue
      write(unit=61, fmt=41) (a(i),i = 1, nchar)
      write(unit=62, fmt=41) (a(i),i = 1, nchar)
c     IDENTIFY CONTROL CARD  --  1=PA,2=IT,3=RA,4=CC                    
c         
   41 format(1x,80a1)
      do 50 i = 1, 4
      if ((a(i1) .ne. mc(i)) .or. (a(icol) .ne. mc(i + 4))) goto 50
   45 goto (95, 105, 130, 140), i
   50 continue
c     ERROR IN CONTROL CARD                                             
c         
      if ((a(i1) .eq. 'R') .and. (a(icol) .eq. 'P')) goto 130
   55 write(unit=61, fmt=56) icol, (a(i),i = 1, nchar)
      write(unit=62, fmt=56) icol, (a(i),i = 1, nchar)
   56 format(/1x,5hERROR,5x,i5,5x,80a1/(21x,80a1))
   57 read(unit=ic, fmt=26, end=140) (a(i),i = 1, reclen)
      write(unit=61, fmt=41) (a(i),i = 1, reclen)
      write(unit=62, fmt=41) (a(i),i = 1, reclen)
      icol = 1
      call nsp1
      icol = icol + 1
      if ((a(i1) .eq. 'C') .and. (a(icol) .eq. 'C')) goto 20
c                                                                       
c         
c     PA(V= ,U= ,D= ,T= ,Q= ,H= ,Z= ,X= ,W= ,T1= ,T2= ,T3= ,R= )(N= )   
c         
c                                                                       
c         
c DEFAULT VALUES                                                        
c         
      goto 57
   95 do 96 i = 1, nall
      prep(i) = 0.
   96 xall(i) = 0.
      prep(1) = 1.
      prep(7) = 1.
      xall(1) = 1.
      xall(7) = 1.
      prep(10) = 1.0
      prep(11) = 0.5
      prep(12) = 0.0
      xall(10) = prep(10)
      xall(11) = prep(11)
      xall(12) = prep(12)
      prep(13) = .5
      xall(13) = .5
      nn = 5
      call par
      if (index) 97, 55, 97
   97 itemp = iblank
      i = 0
   98 i = i + 1
      tt(i) = a(icol)
      icol = icol + 1
      if (((i .ne. 3) .and. (a(icol) .ne. ' ')) .and. (a(icol) .ne. '=')
     &) goto 98
      call nsp1
      do 100 i = 1, nall
      if (itemp .eq. jname(i)) goto 102
  100 continue
      goto 55
  102 xall(i) = xno1(icol)
      prep(i) = xall(i)
  104 if (a(icol) .eq. ')') goto 2000
      call nspace
      goto 97
 2000 call par
      if (index) 2001, 1599, 2001
 2001 call nsp1
      if (a(icol) .ne. 'N') goto 55
      call nspace
      if (a(icol) .ne. '=') goto 55
      call nsp1
      nn = xno1(icol)
      if (a(icol) .eq. ')') goto 1599
      goto 55
 1599 if (xall(7) .eq. 0.) then
      xall(7) = .00001
      prep(7) = .00001
      end if
      if (ibegin) 20, 1600, 20
 1600 ibegin = 1
 1601 read(unit=48) inulls, (ia(i),i = 1, 114)
      if (inulls .ne. nulls) goto 1601
 1611 rewind(unit=48) 
      if (idebug .eq. 1) then
      do 1613 i = 1, 6
 1613 write(unit=61, fmt=1614) (ia(j),j = 1 + (12 * (i - 1)), 12 * i)
 1614 format(1h ,12i10)
      end if
      nread = ia(1)
      nia = ia(2)
      npi = ia(3)
      nfi = ia(4)
      do 1612 i = 1, 10
 1612 yi(i) = aa(i + 2)
      jli = ia(25)
      iiaf = ia(26)
      do 1620 i = 1, 27
 1620 ititle(i) = ia(i + 26)
      itypgr = ia(54)
      do 1621 i = 1, 10
 1621 pis(i) = aa(i + 27)
      do 1622 i = 1, 10
 1622 fi(i) = aa(i + 37)
      do 1623 i = 1, 10
c     WRITE(6,9898)NREAD,NIA,NPI,NFI,JLI,IIAF,ITYPGR,YI,ZI,PIS          
c         
c9898 FORMAT(1X,7I4/1X,10F7.3/1X,10F7.3/1X,10F7.3)                      
c         
 1623 zi(i) = aa(i + 47)
c                                                                       
c         
c     IT(V,U,D,T,Q,H,Z,X,W,T1,T2,T3,R)  (H=,B=,T=)                      
c         
c                                                                       
c         
      goto 20
  105 ictrl = 1
      hdif = 1.e-03
      trupb = 10. * sqrt(hdif)
      tol = 1.e-04
      do 106 i = 1, nall
  106 itp(i) = 0
      call par
      if (index) 107, 55, 107
  107 itemp = iblank
      i = 0
  108 i = i + 1
      tt(i) = a(icol)
      icol = icol + 1
      if ((((i .ne. 3) .and. (a(icol) .ne. ' ')) .and. (a(icol) .ne. ','
     &)) .and. (a(icol) .ne. ')')) goto 108
      call nsp1
      do 110 i = 1, nall
      if (itemp .eq. jname(i)) goto 112
  110 continue
      goto 55
  112 itp(i) = 1
      if (a(icol) .eq. ')') goto 170
      call nspace
      goto 107
  170 ichou = 0
      call par
      if (index) 171, 476, 171
  176 call nspace
  171 if (a(icol) .ne. 'H') goto 173
      call nspace
      if (a(icol) .ne. '=') goto 55
      hdif = xno1(icol)
      ichou = ichou + 1
  172 if (a(icol) .eq. ',') goto 176
      if (a(icol) .ne. ')') goto 55
      if (ichou .eq. 1) trupb = 10. * sqrt(hdif)
      goto 476
  173 if (a(icol) .ne. 'B') goto 174
      call nspace
      if (a(icol) .ne. '=') goto 55
      trupb = xno1(icol)
      ichou = ichou + 2
      goto 172
  174 if (a(icol) .ne. 'T') goto 55
      call nspace
      if (a(icol) .ne. '=') goto 55
      tol = xno1(icol)
c                                                                       
c         
c     RA(V= ,U= ,D= ,T= ,Q= ,H= ,Z= ,X= ,W= ,T1= ,T2= ,T3= ,R= )(TRAT)  
c         
c     RP(V= ,U= ,D= ,T= ,Q= ,H= ,Z= ,X= ,W= ,T1= ,T2= ,T3= ,R= )(TRAT)  
c         
c                                                                       
c         
      goto 172
  130 ictrl = 2
      irap = 0
      if (a(icol) .eq. 'P') then
      irap = 1
      end if
      trat = 1000.
      do 131 i = 1, nall
  131 rxall(i) = 0.
      rxall(1) = 1.
      rxall(7) = 1.
      rxall(10) = 1.0
      rxall(11) = 0.5
      rxall(12) = 0.0
      rxall(13) = .5
      call par
      if (index) 132, 55, 132
  132 itemp = iblank
      i = 0
  133 i = i + 1
      tt(i) = a(icol)
      icol = icol + 1
      if (((i .ne. 3) .and. (a(icol) .ne. ' ')) .and. (a(icol) .ne. '=')
     &) goto 133
      call nsp1
      do 134 i = 1, nall
      if (itemp .eq. jname(i)) goto 136
  134 continue
      goto 55
  136 rxall(i) = xno1(icol)
  138 if (a(icol) .eq. ')') goto 300
      call nspace
      goto 132
  300 call par
      if (index) 301, 476, 301
  301 trat = xno(icol)
c                                                                       
c         
c     CC                                                                
c         
c                                                                       
c         
      goto 476
  140 ictrl = 100
  476 continue
      return 
      end
c                                                                       
c         
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         
c                                                                       
c         
      subroutine fun(slik, x, n, xall, nall, itp)
c                OPTIONAL DIAGNOSTIC PRINTOUT IF DEBUG=1                
c         
c                                                                       
c         
      implicit double precision (o-z, a-h)
      integer idebug
      dimension b(3751)
      dimension eva(100, 14), beva(13, 13)
      equivalence (ib, b)
      common /auxp1/ p(3), g(3), clm2, clm1, cxp2, cxp1(3), rcdg, 
     &cdp, z(10), vv, uu, gg, ca, ck, ea, ek, scals(4), sclq, rsqw
     &, notrs, ng1, nox
      common /auxp2/ ntrel, mxrel, clm2s(3), clm1s(3), rrrs(2, 3, 5), 
     &cxc2s(2, 2), cxc1s(2, 2), cdcs(2, 2), cxd2s(2, 2, 5), cxd1s(2, 2, 
     &5), cdds(2, 2, 5), cxa2s(2, 5), cxa1s(2, 5), cdas(2, 5)
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, neva, idebug, np
      common /odds/ pis(10), fi(10), ititle(27), itypgr, irap, npi, nfi
      common /intg/ xuu(20), auu(20), xdum(20), a(20), fzdqz(10), 
     &tvuu, fcaea, xdev, igot, nuu, nl, nm
      common /data/ tmf(9), tmm(9), qt(18), tmk(18), 
     &zfa, zmo, qtfp, qtf, qtm, pi, 
     &qtmp, qtcp, zphf, zphm, zphc, iph(18), iphz(18), irpf, irpm, irpk, 
     &nc, isel, iphf, iphm, iphzf, iphzm, iphfp, iphmp, iphcp
      common /data2/ cxc2, cxc1, cdc, rrc, cxpf2, 
     &cxpf1, cdpf, cxpm2, cxpm1, cdpm, rrpf, rrpm, rrpc, 
     &cxd2, cxd1, cdd, alik, iquad, iccc, irrpf, irrpm, irrpc
      common /mini/ tl(15, 15), d(15), gdum(15), gsave(15), y(15), pdum(
     &15), tmin, fsmf, idif, isw, iret
      common /trsm/ tr(18), tm(6, 3, 70)
      common ib(7502)
      dimension xall(15), itp(15), x(15)
      data nfam /0/
      slik = 0.0
      neva = neva + 1
      jeva = 0
c ------------------------------------------------------------------    
c         
c FOR EACH ENTRY TO FUN WITH NEW PARAMETERS                             
c         
      nnfam = 0
      k = 0
      do 5 i = 1, nall
      if (itp(i) .ne. 1) goto 5
      k = k + 1
      xall(i) = x(k)
c -------------------------------------------------------------------   
c         
c GET AUXILIARY PARAMETERS                                              
c         
c     CALL BTIME                                                        
c         
    5 continue
c ------------------------------------------------------------------    
c         
c READ BLOCK OF DATA IF REQUIRED                                        
c         
      call auxpar(xall, nall, itp)
      if (nread - 1) 10, 12, 10
   10 read(unit=48, end=13) (ib(i),i = 1, 7502)
   12 nfam = ib(7501)
      nnfam = nnfam + nfam
c -------------------------------------------------------------------   
c         
c FOR EACH FAMILY, GET -LN(LIKELIHOOD)                                  
c         
      iword2 = 0
      do 310 i = 1, nfam
      iword = iword2 + iword2
      k2 = 0
      qtf = 0.
      qtm = 0.
      iphm = 2
      nc = ib(iword + 3)
c ************ TEMP ******************                                  
c         
      iquad = ib(iword + 8)
      if (iquad .eq. 2) iquad = 1
      itype = 2
      if (iquad .eq. 2) itype = 1
      clm2 = clm2s(iquad)
      clm1 = clm1s(iquad)
      isel = ib(iword + 14)
      pi = 0.
c     IQUAD=2, ONE PARENT MISSING AND NO POINTER TO THIS PARENT         
c         
      if (isel .gt. 0) pi = pis(isel)
      if (iquad - 2) 450, 459, 480
  459 if (ib(iword + 4) .eq. 0) goto 18
      iphf = 2
      if (ib(iword + 4) - 2) 463, 463, 460
  460 if (ib(iword + 21) .ne. 2) goto 461
      if (ib(iword + 27) .eq. 2) goto 462
      qtm = b(iword2 + 13)
      iphm = ib(iword + 27)
      iphzm = ib(iword + 28)
      zmo = z(iphzm)
  462 iphf = 2
      k2 = k2 + 16
      goto 25
  461 k2 = 3
      goto 22
  463 if (ib(iword + 21) .eq. 2) goto 470
      goto 22
  480 if (ib(iword + 4) .eq. 0) goto 18
      iphf = 2
      if (ib(iword + 4) - 2) 470, 470, 471
  470 k2 = 13
      goto 25
  471 k2 = 16
c --------------------                                                  
c         
c     FATHER,MOTHER DATA                                                
c         
      goto 25
  450 if (ib(iword + 4) - 2) 19, 22, 23
   19 if (ib(iword + 4)) 21, 18, 21
   18 iphf = 2
      k2 = 10
      goto 25
   23 qtm = b(iword2 + 13)
      iphm = ib(iword + 27)
      iphzm = ib(iword + 28)
      zmo = z(iphzm)
      k2 = 3
   21 qtf = b(iword2 + 10)
      iphf = ib(iword + 21)
      iphzf = ib(iword + 22)
      zfa = z(iphzf)
      k2 = k2 + 13
      goto 25
   22 qtm = b(iword2 + 10)
      iphm = ib(iword + 21)
      iphzm = ib(iword + 22)
      zmo = z(iphzm)
      iphf = 2
c --------------------                                                  
c         
c     CHILDREN DATA                                                     
c         
      k2 = k2 + 13
   25 k = ((k2 + k2) + 1) + iword
      k2 = k2 + iword2
      igec = ib(iword + 9)
      cxc2 = cxc2s(igec,itype)
      cxc1 = cxc1s(igec,itype)
      cdc = cdcs(igec,itype)
      rrc = rrrs(igec,iquad,1)
      iccc = ib(iword + 13)
      do 27 j = 1, nc
      qt(j) = b(k2)
      iph(j) = ib(k)
      iphz(j) = ib(k + 1)
      k2 = k2 + 3
   27 k = k + 6
      k = k - 3
c --------------------                                                  
c         
c     FATHER,MOTHER,CHILDREN'S POINTER                                  
c         
      k2 = k2 - 1
      if (iquad .ne. 2) goto 451
      iphfp = 2
      if (ib(iword + 5)) 29, 29, 428
  428 irpm = ib(iword + 5)
      irrpm = 1 + (irpm / (ntrel + 1))
      igepm = ib(iword + 10)
      rrpm = rrrs(igepm,iquad,irrpm)
      goto 452
  451 if (ib(iword + 5)) 40, 40, 28
   40 iphfp = 2
      goto 29
   28 irpf = ib(iword + 5)
      irrpf = 1 + (irpf / (ntrel + 1))
      igepf = ib(iword + 10)
      rrpf = rrrs(igepf,iquad,irrpf)
      cxpf2 = cxa2s(igepf,irrpf)
      cxpf1 = cxa1s(igepf,irrpf)
      cdpf = cdas(igepf,irrpf)
      qtfp = b(k2 + 1)
      iphfp = ib(k + 3)
      zphf = z(ib(k + 4))
      k = k + 6
      k2 = k2 + 3
   29 if (ib(iword + 6)) 41, 41, 30
   41 iphmp = 2
      goto 31
   30 irpm = ib(iword + 6)
      irrpm = 1 + (irpm / (ntrel + 1))
      igepm = ib(iword + 11)
      rrpm = rrrs(igepm,iquad,irrpm)
  452 cxpm2 = cxa2s(igepm,irrpm)
      cxpm1 = cxa1s(igepm,irrpm)
      cdpm = cdas(igepm,irrpm)
      qtmp = b(k2 + 1)
      iphmp = ib(k + 3)
      zphm = z(ib(k + 4))
      k = k + 6
      k2 = k2 + 3
   31 if (ib(iword + 7)) 42, 42, 32
   42 iphcp = 2
      goto 35
c SINGLE SELECTION OF MULTIPLEX FAMILIES IF IRPK=MXREL+1 (CODE=1A)      
c         
c IF SO, SET ISEL=NPI AND PI=1.0                                        
c         
   32 irpk = ib(iword + 7)
      if ((mxrel + 1) - irpk) 34, 33, 34
   33 isel = npi
      pi = 1.0
   34 continue
      irrpc = irpk - mxrel
      irrpc = 1 + (irrpc / (ntrel + 1))
      igepc = ib(iword + 12)
      rrpc = rrrs(igepc,iquad,irrpc)
      cxd2 = cxd2s(igepc,itype,irrpc)
      cxd1 = cxd1s(igepc,itype,irrpc)
      cdd = cdds(igepc,itype,irrpc)
      qtcp = b(k2 + 1)
      iphcp = ib(k + 3)
c --------------------                                                  
c         
c     GET READY TO CALL EXEC                                            
c         
c SKIP THIS IF Q=0.0 AND CALL LIKEQ0                                    
c         
      zphc = z(ib(k + 4))
   35 if (xall(5)) 36, 36, 350
   36 call likeq0
      goto 830
  350 if (iphmp - 2) 822, 823, 822
  822 tmm(1) = tm(1,1,irpm)
      tmm(2) = tm(1,2,irpm)
      tmm(3) = tm(1,3,irpm)
      tmm(4) = tm(2,1,irpm)
      tmm(5) = tm(2,2,irpm)
      tmm(6) = tm(2,3,irpm)
      tmm(7) = tm(3,1,irpm)
      tmm(8) = tm(3,2,irpm)
      tmm(9) = tm(3,3,irpm)
  823 if (iphfp - 2) 824, 825, 824
  824 tmf(1) = tm(1,1,irpf)
      tmf(2) = tm(1,2,irpf)
      tmf(3) = tm(1,3,irpf)
      tmf(4) = tm(2,1,irpf)
      tmf(5) = tm(2,2,irpf)
      tmf(6) = tm(2,3,irpf)
      tmf(7) = tm(3,1,irpf)
      tmf(8) = tm(3,2,irpf)
      tmf(9) = tm(3,3,irpf)
  825 if (iphcp - 2) 826, 827, 826
  826 tmk(1) = tm(1,1,irpk)
      tmk(2) = tm(1,2,irpk)
      tmk(3) = tm(1,3,irpk)
      tmk(4) = tm(2,1,irpk)
      tmk(5) = tm(2,2,irpk)
      tmk(6) = tm(2,3,irpk)
      tmk(7) = tm(3,1,irpk)
      tmk(8) = tm(3,2,irpk)
      tmk(9) = tm(3,3,irpk)
      tmk(10) = tm(4,1,irpk)
      tmk(11) = tm(4,2,irpk)
      tmk(12) = tm(4,3,irpk)
      tmk(13) = tm(5,1,irpk)
      tmk(14) = tm(5,2,irpk)
      tmk(15) = tm(5,3,irpk)
      tmk(16) = tm(6,1,irpk)
      tmk(17) = tm(6,2,irpk)
c --------------------------------------------------------------        
c         
c GET LIKELIHOOD FOR THAT FAMILY                                        
c         
      tmk(18) = tm(6,3,irpk)
  827 continue
      call like
  830 slik = slik + alik
c     TO GET BETTER K-MATRIX                                            
c         
      if ((neva - n) - 1) 600, 600, 630
  600 if (jeva .eq. 100) goto 630
      jeva = jeva + 1
      eva(jeva,neva) = alik
  630 iword2 = (iword2 + 9) + (ib(iword + 2) * 3)
c --------------------------------------------------------------        
c         
c END OF LOOP OVER FAMILIES                                             
c         
  310 continue
      if (nread .eq. 1) goto 500
      goto 10
   13 rewind(unit=48) 
  500 if (neva .ne. (n + 1)) goto 510
      h = .001
      do 631 i = 1, n
      hi = h * max(abs(x(i)),.1d0)
      do 631 j = i, n
      hj = h * max(abs(x(j)),.1d0)
      temb = 0.
      do 632 k = 1, jeva
  632 temb = temb + (((eva(k,i + 1) - eva(k,1)) * (eva(k,j + 1) - eva(k,
     &1))) / (hi * hj))
      beva(i,j) = temb
  631 beva(j,i) = temb
      if (jeva .ge. nfam) goto 634
      cons = float(nnfam) / float(jeva)
      do 633 i = 1, n
      do 633 j = 1, n
  633 beva(i,j) = cons * beva(i,j)
  634 continue
      if (idebug .ne. 1) goto 520
      do 521 i = 1, n
  521 write(unit=61, fmt=522) (beva(i,j),j = 1, n)
c ----------------------------------------------------------------------
c         
c CARRIES OUT CHOLESKY FACTORIZATION ON B-MATRIX                        
c         
  522 format(1x,4hBEVA,9e13.6)
  520 do 100 ic = 1, n
      if (beva(ic,ic)) 502, 502, 220
  220 temb = beva(ic,ic)
      d(ic) = beva(ic,ic)
      tl(ic,ic) = 1.
      if (ic - n) 231, 510, 510
  231 ic1 = ic + 1
      do 50 k = ic1, n
   50 tl(k,ic) = beva(k,ic) / temb
      do 100 i = ic1, n
      tlic = tl(i,ic)
      do 100 k = i, n
  100 beva(k,i) = beva(k,i) - ((tlic * tl(k,ic)) * temb)
  502 do 504 i = 1, n
      do 503 j = 1, i
  503 tl(i,j) = 0.
      tl(i,i) = 1.
  504 d(i) = 1. / sqrt(x(i))
c     CALL ETIME                                                        
c         
  510 continue
      return 
      end
      function gam(y)
c     IMPLICIT REAL*8(A-H,O-Z)                                          
c         
      implicit double precision (o-z, a-h)
      gam = (((((((((((((((.035868343 * y) - .193527818) * y) + 
     &.482199394) * y) - .756704078) * y) + .918206857) * y) - 
     &.897056937) * y) + .988205891) * y) - .577191652) * y) + 1.0
      return 
      end
      function gamma(x)
c     IMPLICIT REAL*8(A-H,O-Z)                                          
c         
      implicit double precision (o-z, a-h)
      external gam
      z = x
      if (z) 1, 1, 4
    1 gamma = 0.0
      write(unit=61, fmt=2) z
    2 format(/2x,19hARG ERROR FOR GAMMA,e16.8)
      goto 14
    4 if (z - 70.) 6, 1, 1
    6 if (z - 1.) 8, 7, 9
    7 gamma = 1.0
      goto 14
    8 gamma = gam(z) / z
      goto 14
    9 za = 1.0
   10 z = z - 1.
      if (z - 1.) 13, 11, 12
   11 gamma = za
      goto 14
   12 za = za * z
      goto 10
   13 gamma = za * gam(z)
   14 return 
      end
c                                                                       
c         
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       
c                                                                       
c         
      subroutine gemini(h, trupb, nall, ibnd, ihess, xall, itp, n, x, 
     &bnd, f, tol, gtg, ptg, nit, nfe, ip, ivar, se, ic, maxit, ihx, idg
     &)
c ----------------------------------------------------------------------
c-------- 
c                                                                       
c         
c ROUTINE FOR UNCONSTRAINED OPTIMIZATION OF A NON-LINEAR FUNCTION OF N V
cARIABLES 
c                                                                       
c         
c ----------------------------------------------------------------------
c-------- 
      implicit double precision (o-z, a-h)
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, neva, idebug, np
      common /mini/ tl(15, 15), d(15), g(15), gsave(15), y(15), p(15), 
     &tmin, fsmf, idif, isw, iret
      common /bound/ ibndv
c                                                                       
c         
c-----------------------------------------------------------------------
c-        
c                                                                       
c         
c GEMINI, VERSION 1, FEBRUARY 79, J.M.LALOUEL.                          
c         
c                                                                       
c         
c A VARIABLE METRIC ( OR QUASI-NEWTON ) ALGORITHM EVALUATING GRADIENT BY
c         
c FINITE DIFFERENCE APPROXIMATION, WHILE AN APPROXIMATION TO THE HESSIAN
c         
c IS OBTAINED BY RECURRENCE, IMPROVED AT EACH ITERATION BY AN UPDATING  
c         
c PROCEDURE.                                                            
c         
c THE PRESENT ALGORITHM USE RESULTS DUE TO DAVIDON(1959), FLETCHER AND  
c         
c POWELL(1963),STEWART(1967), GILL AND MURRAY(1972), DAVIDON(1975) AND  
c         
c GOLDFARB(1976).                                                       
c         
c THE MOST ESSENTIAL FEATURE OF THIS ALGORITHM CONCERNS THE UPDATING .  
c         
c THE APPROXIMATION TO THE HESSIAN IS RECURRED IN FACTORIZED FORM, AND  
c         
c A MODIFICATION OF AT MOST RANK TWO IS PERFORMED, AT EACH ITERATION,   
c         
c IN REAL PRODUCT FORM ACCORDING TO THE METHOD OF GOLDFARB(1976). IN    
c         
c ADDITION, THE UPDATE IS OPTIMALLY CONDITIONED IN THE SENSE OF DAVIDON.
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c--       
c                                                                       
c         
      common /hang/ iwrong
      dimension xall(15), itp(15), x(15), bnd(15, 2)
c ----------------------------------------------------------------------
c-------- 
c PRECIS,H AND TRUPB ARE DEFINED IN MAIN                                
c         
c ( OR MAY BE READ IN INPUT ROUTINE ).                                  
c         
      dimension v(15), se(15)
      write(unit=ip, fmt=3006) h, trupb
c ----------------------------------------------------------------------
c-------- 
c INITIALIZATION                                                        
c         
 3006 format(/3x,25hDIFFERENTIATION INTERVAL=,e12.5,
     &26hTRUNCATION UPPER    BOUND=,e12.5//)
  700 neva = 0
      nit = 0
      idg = 0
      tolg = 1.e-38
      idif = 1
      t = 0.1
      tmin = 0.0
      fsmf = 0.0
      iwrong = 0
      call fun(f, x, n, xall, nall, itp)
	call flush(61)
	call flush(62)
      if (iwrong .ne. 0) goto 668
c ----------------------------------------------------------------------
c-------- 
c INITIAL B-MATRIX IS TAKEN HERE AS IDENTITY MATRIX, AND STORED IN CHOLE
cSKY      
c     FACTORIZED FORM LDLT. TL IS LOWER TRIANGULAR UNIT MATRIX, WHILE DI
cAGONAL   
c     OF D IS STORED IN VECTOR D.                                       
c         
      nfe = 1
      do 20 i = 1, n
      do 10 j = 1, i
   10 tl(i,j) = 0.0
      tl(i,i) = 1.0
   20 d(i) = 1.0
      if (ihess) 2, 2, 1
c READ OR CALCULATE INITIAL HESSIAN IF DESIRED                          
c         
    1 icall = 0
      call inib(f, x, n, xall, nall, itp, h, ihess, ip, ivar, icall, se
     &, ic)
	call flush(61)
	call flush(62)
      t = 1.0
c ----------------------------------------------------------------------
c-------- 
c START ITERATION WITH CHOICE OF FORWARD DIFFENTIATION                  
c         
    2 continue
c ----------------------------------------------------------------------
c-------- 
c GET GRADIENT BY EITHER FORWARD OR CENTRAL DIFFERENCES                 
c         
  150 isw = 1
c -------------------------------------                                 
c         
c FORWARD DIFFERENCE                                                    
c         
  200 if (idif - 1) 201, 201, 250
  201 do 210 i = 1, n
      hx = h
      if (ihx - 1) 203, 202, 203
  202 hx = h * max(abs(x(i)),.1d0)
  203 continue
      xsave = x(i)
      x(i) = x(i) + hx
      call fun(fxph, x, n, xall, nall, itp)
      g(i) = (fxph - f) / hx
  210 x(i) = xsave
      nfe = nfe + n
c -------------------------------------                                 
c         
c CENTRAL DIFFERENCE                                                    
c         
      if (isw - 1) 300, 300, 400
  250 do 260 i = 1, n
      hx = h
      if (ihx - 1) 253, 252, 253
  252 hx = h * max(abs(x(i)),.1d0)
  253 continue
      xsave = x(i)
      x(i) = x(i) + hx
      call fun(fxph, x, n, xall, nall, itp)
      x(i) = xsave - hx
      call fun(fxmh, x, n, xall, nall, itp)
	call flush(61)
	call flush(62)
      g(i) = (fxph - fxmh) / (hx + hx)
  260 x(i) = xsave
      nfe = (nfe + n) + n
c ----------------------------------------------------------------------
c-------- 
c INCREMENT ITERATION COUNT AND GET SEARCH DIRECTION BY FORWARD AND BACK
cWARD     
c     SUBSTITUTION                                                      
c         
      if (isw - 1) 300, 300, 400
  300 nit = nit + 1
      fsav2 = fsave
      fsave = f
      write(unit=ip, fmt=3001) nit, t, nfe, f
 3001 format(1x,130(1h-),/1x,9hITERATION,i4,6h  T = ,e10.3,5x,6hNFE = 
     &,i5,2x,9hF(-LN L)=, e18.11)
      if (nit - maxit) 3003, 3002, 3002
 3002 idg = 6
      goto 600
 3003 continue
      if (n - 10) 301, 301, 302
  301 write(unit=ip, fmt=1333) (xall(i),i = 1, nall)
 1333 format(1x,33h(V,U,D,T,Q,CK,CA,X,W,T1,T2,T3,R=),9f10.5,/34x,4f10.5)
      write(unit=ip, fmt=1001) (x(i),i = 1, n)
 1001 format(1x,2hX=,4e18.11)
      write(unit=ip, fmt=1002) (g(i),i = 1, n)
c --------------------------------------------------------------------- 
c         
c A TERMINATION CRITERION MAY BE SET ON THE NORM SQRT(GTG)              
c         
 1002 format(1x,2hG=,4e18.11)
  302 gtg = 0.0
      do 303 i = 1, n
      gtg = gtg + (g(i) * g(i))
  303 continue
      gtg = sqrt(gtg)
      if (gtg - tolg) 304, 304, 305
  304 idg = 5
      goto 600
c -------------------------------------                                 
c         
c FORWARD SUBSTITUTION TO SOLVE L*V=-G FOR V                            
c         
  305 continue
      if (n - 1) 310, 310, 315
  310 p(1) = - (g(1) / d(1))
      goto 340
  315 v(1) = - g(1)
      do 320 i = 2, n
      v(i) = - g(i)
      im1 = i - 1
      do 320 j = 1, im1
c -------------------------------------                                 
c         
c BACKWARD SUBSTITUTION TO SOLVE LT*P=DINV*V                            
c         
  320 v(i) = v(i) - (tl(i,j) * v(j))
      p(n) = v(n) / d(n)
      nm1 = n - 1
      do 330 j = 1, nm1
      nmj = n - j
      p(nmj) = v(nmj) / d(nmj)
      nmjp1 = nmj + 1
      do 330 i = nmjp1, n
  330 p(nmj) = p(nmj) - (tl(i,nmj) * p(i))
      write(unit=ip, fmt=1003) (p(i),i = 1, n)
c ----------------------------------------------------------------------
c-------- 
c CHECK IF DIRECTION OF SEARCH IS DOWNHILL                              
c         
 1003 format(1x,2hP=,4e18.11)
  340 ptg = 0.0
c TERMINATE IF EXCESSIVE CANCELLATION LEADS TO SOME P(I) TO BE ZERO     
c         
c THIS TEST APPLIED ONLY IF BOUNDARY CONDITIONS SPECIFIED               
c         
      do 350 i = 1, n
      if (ibnd) 350, 350, 344
  344 if (p(i)) 350, 345, 350
  345 idg = 7
      goto 600
  350 ptg = ptg + (p(i) * g(i))
      idg = 1
      if (ptg) 355, 500, 500
c IF PTG HAS MET SOME TOLERANCE TOL, THE PROCEDURE IS ENDED             
c         
  355 idg = 0
      if (abs(ptg) - tol) 356, 356, 357
  356 idg = 4
      fsmf = fsav2 - f
      write(unit=ip, fmt=4002) fsmf, ptg
 4002 format(/3x,5hFSMF=,e12.5,6h  PTG=,e12.5)
      goto 600
c ----------------------------------------------------------------------
c-------- 
c GET READY TO CARRY OUT LINEAR SEARCH. GET MINIMUM ACCEPTABLE STEP, TMI
cN,       
c     AND INITIAL ESTIMATE OF STEP SIZE T                               
c         
  357 continue
      xpm = 1.e30
      if (nit - 1) 361, 369, 361
  361 do 360 i = 1, n
      if (ihx - 1) 363, 362, 363
  362 xp = max(abs(x(i)),.1d0) / abs(p(i))
      xpm = min(xp,xpm)
      goto 360
  363 xpm = min(abs(1.d0 / p(i)),xpm)
  360 continue
      tmin = ((.5 * xpm) * h) / trupb
      t = - ((2. * fsmf) / ptg)
      if (t) 366, 366, 367
  366 t = 1.0
  367 t = min(t,1.0d0)
      if (idif - 2) 369, 368, 369
  368 t = max(.5d0,t)
  369 continue
      write(unit=ip, fmt=4001) fsmf, ptg, xpm, tmin
 4001 format(3x,5hFSMF=,e12.5,6h  PTG=,e12.5,6h  XPM=,e12.5,7h  TMIN=
     &,e12.5)
      write(unit=ip, fmt=9996) t
 9996 format(3x,10hINITIAL T=,e12.5)
      if (iwrong .eq. 0) goto 667
  668 write(unit=ip, fmt=666) iwrong
  666 format(1x,46hARITHMETIC ERROR DETECTED IN AUXPAR. IWRONG = ,i2/1x,
     &51hABORTED CURRENT IT CARD BEFORE I GET INDIGESTION.  ,6hHARRY.)
      write(unit=62, fmt=666) iwrong
c CHECK BOUNDS IF ANY                                                   
c         
      goto 604
  667 if (ibnd) 371, 371, 370
  370 nsave = n
      call chkbnd(t, tbnd, n, x, bnd, ip, h, xall, nall, itp, ihess)
      if ((n .eq. 0) .or. (ibndv .ne. 0)) goto 604
c A RETURN TO INITIALIZATION OCCURS IF A PARAMETER HAS BEEN SET TO A BOU
cND       
      if (nsave - n) 371, 371, 700
c ----------------------------------------------------------------------
c-------- 
c ENTER LINEAR SEARCH                                                   
c         
  371 continue
	call flush(61)
	call flush(62)
      call step(f, x, n, xall, nall, itp, fsave, t, tbnd, nfe, h, trupb
     &, idg, ibnd, ip)
	call flush(61)
	call flush(62)
c ----------------------------------------------------------------------
c-------- 
c PREPARATION TO UPDATE LDLT                                            
c         
c OLD GRADIENT WAS SAVED IN GSAVE, NEW GRADIENT IS IN G                 
c         
      if (iret - 1) 500, 500, 200
  400 ytp = 0.0
      fsmf = fsave - f
      do 410 i = 1, n
      y(i) = g(i) - gsave(i)
c IF YTP NOT POSITIVE, NO UPDATE AT THIS ITERATION                      
c         
  410 ytp = ytp + (y(i) * p(i))
c ----------------------------------------------------------------------
c-------- 
c CALL UPDATING ROUTINE                                                 
c         
      if (ytp) 300, 300, 450
  450 if (n - 1) 460, 460, 470
  460 d(1) = - ((y(1) * d(1)) / (t * gsave(1)))
      goto 300
  470 call update(n, t, ip, idg)
      if (idg - 3) 471, 600, 471
c ----------------------------------------------------------------------
c-------- 
c RESTART WITH IDIF=2 OR EXIT                                           
c         
  471 goto 300
  500 if (idif - 2) 501, 600, 501
  501 idif = 2
      write(unit=ip, fmt=3005) 
 3005 format(/30x,28(1h-),//30x,28hSWITCH TO CENTRAL DIFFERENCE
     &,//30x,28(1h-)/)
c ----------------------------------------------------------------------
c-------- 
c EXIT SECTION                                                          
c         
      goto 150
  600 if (idg) 602, 150, 602
  602 f = fsave
      k = 0
      do 610 i = 1, nall
      if (itp(i) - 1) 610, 609, 610
  609 k = k + 1
      xall(i) = x(k)
c GET STANDARD ERRORS IF REQUIRED                                       
c         
  610 continue
      if (ivar - 1) 604, 603, 604
  603 icall = 1
      call inib(f, x, n, xall, nall, itp, h, ihess, ip, ivar, icall, se
     &, ic)
      call auxpar(xall, nall, itp)
	call flush(61)
	call flush(62)
  604 continue
      return 
      end
c                                                                       
c         
c **********************************************************************
c*        
c                                                                       
c         
c ----------------------------------------------------------------------
c-        
c                                                                       
c         
c A SERIES OF ROUTINES FOR GAUS-HERMITE QUADRATURE                      
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c-        
c                                                                       
c         
      subroutine hermit(nn, x, a, eps)
c     IMPLICIT REAL*8(A-H,O-Z)                                          
c         
      implicit double precision (o-z, a-h)
      dimension x(20), a(20)
      data ni /0/
      fn = nn
      n1 = nn - 1
      n2 = (nn + 1) / 2
      cc = (1.7724538509 * gamma(fn)) / (2. ** n1)
      s = ((2. * fn) + 1.) ** .16667
      do 10 i = 1, n2
      if (i - 1) 10, 1, 2
    1 xt = (s ** 3) - (1.85575 / s)
      goto 9
    2 if (i - 2) 10, 3, 4
    3 xt = xt - ((1.14 * (fn ** .426)) / xt)
      goto 9
    4 if (i - 3) 10, 5, 6
    5 xt = (1.86 * xt) - (.86 * x(1))
      goto 9
    6 if (i - 4) 10, 7, 8
    7 xt = (1.91 * xt) - (.91 * x(2))
      goto 9
    8 xt = (2. * xt) - x(i - 2)
    9 call hroot(xt, nn, dpn, pn1, eps)
      x(i) = xt
      a(i) = (cc / dpn) / pn1
      ni = (nn - i) + 1
      x(ni) = - xt
   10 a(ni) = a(i)
      return 
      end
      subroutine hrecur(pn, dpn, pn1, x, nn)
c     IMPLICIT REAL*8(A-H,O-Z)                                          
c         
      implicit double precision (o-z, a-h)
      p1 = 1.
      p = x
      dp1 = 0.0
      dp = 1.0
      do 1 j = 2, nn
      fj = j
      fj2 = (fj - 1.) / 2.
      q = (x * p) - (fj2 * p1)
      dq = ((x * dp) + p) - (fj2 * dp1)
      p1 = p
      p = q
      dp1 = dp
    1 dp = dq
      pn = p
      dpn = dp
      pn1 = p1
      return 
      end
      subroutine hroot(x, nn, dpn, pn1, eps)
c     IMPLICIT REAL*8(A-H,O-Z)                                          
c         
      implicit double precision (o-z, a-h)
      iter = 0
    1 iter = iter + 1
      call hrecur(p, dp, pn1, x, nn)
      d = p / dp
      x = x - d
      if (abs(d) - eps) 3, 3, 2
    2 if (iter - 10) 1, 3, 3
    3 dpn = dp
      return 
      end
c                                                                       
c         
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         
c                                                                       
c         
      subroutine ineq(t, tbnd, n, x, h, xall, nall, itp)
c -------------------------------------------------------------------   
c         
c                                                                       
c         
c CHECK INEQUALITY CONSTRAINT IN POINTER                                
c         
c                                                                       
c         
c --------------------------------------------------------------------  
c         
      implicit double precision (o-z, a-h)
      integer idebug
      common /mini/ tl(15, 15), d(15), g(15), gsave(15), y(15), p(15), 
     &tmin, fsmf, idif, isw, iret
      dimension x(15), xall(15), itp(15)
      dimension xx(15)
      dimension teq(3), fteq(3)
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, nev, idebug, np
      common /auxp1/ pl(3), gl(3), clm2, clm1, cxp2, cxp1(3), rcdg
     &, cdp, z(10), vv, uu, gg, ca, ck, ea, ek, scals(4), sclq, 
     &rsqw, notrs, ng1, nox
c ------------------------------------------------------------------    
c         
c GET CURRENT PARAMETER VALUES INTO XX                                  
c         
      dimension pgc(3)
      k = 0
      do 10 i = 1, nall
      if (itp(i) - 1) 9, 8, 9
    8 k = k + 1
      xx(i) = x(k)
      goto 10
    9 xx(i) = xall(i)
c                                                                       
c         
   10 continue
      iflag = 0
    1 continue
      k = 1
      kk = 1
      teq(1) = 0.0
c                                                                       
c         
c ---------------------------------------------------------------       
c         
c FIND IF INEQUALITY SATISFIED                                          
c         
      kount = 0
   25 j = 0
c     WRITE(6,6789)T,TEQ(K),XX(1),CKA                                   
c         
c6789 FORMAT(1X,'T,TEQ(K),XX(1),CKA',4F10.5)                            
c         
      kount = kount + 1
      do 20 i = 1, nall
      if (itp(i) - 1) 20, 18, 20
   18 j = j + 1
      xx(i) = x(j) + (teq(k) * p(j))
c     WRITE(6,6336)(XX(I),I=1,NALL)                                     
c         
c *** CA ***                                                            
c         
   20 continue
      if (itp(7) .eq. 0) xx(7) = xx(6) * zinp
c INTERPRET PARAMETERS                                                  
c         
      cka = max(xx(6),xx(7))
      uu = xx(2)
      dd = xx(3)
c DESCALE Q                                                             
c         
      tt = xx(4) * sqrt(sclq)
      qq = xx(5) / sclq
      sp = xx(8)
c                                                                       
c         
c ALLOWING FOR NO TRANSMISSION                                          
c         
      u = sp * qq
      if (notrs .eq. 0) goto 2020
      if (notrs .eq. 1) xx(10) = 1. - qq
      xx(11) = 1. - xx(10)
      xx(12) = 1. - xx(10)
c                                                                       
c         
c ---------------------------------                                     
c         
c CHARACTERISTICS OF MAJOR LOCUS                                        
c         
 2020 continue
c SP IS X PROPORTION OF SPORADICS                                       
c         
c **** CHECK IF NEXT STATEMENT IS APPROPRIATE ****                      
c         
      gg = 0.0
      if (qq .eq. 0.0) goto 24
      t1 = 1.0 - qq
      t2 = (1. - u) * (1. - u)
      t3 = 1. - sp
      pl(1) = (((qq * qq) * t3) * t3) / t2
      pl(2) = (((2. * qq) * t1) * t3) / t2
      pl(3) = (t1 * t1) / t2
      pgc(1) = qq * qq
      pgc(2) = (2. * qq) * t1
      pgc(3) = t1 * t1
      zu = (uu - ((qq * qq) * tt)) - ((((2. * qq) * t1) * tt) * dd)
      gl(1) = zu + tt
      gl(2) = zu + (tt * dd)
      gl(3) = zu
      do 3 i = 1, 3
c     WRITE(6,6790)GG                                                   
c         
c6790 FORMAT(1X,'GG',F10.5)                                             
c         
    3 gg = gg + ((pgc(i) * (gl(i) - uu)) * (gl(i) - uu))
   24 fteq(k) = (gg + cka) - xx(1)
c     WRITE(6,6791)FTEQ(K),TEQ(K)                                       
c         
c6791 FORMAT(1X,'FTEQ(K),TEQ(K)',2E10.3)                                
c         
c ----------------------------------------------------------------      
c         
      if (kount .gt. 100) goto 40
      if (kk .eq. 3) goto 55
      k = k + 1
      kk = k
      if (k - 2) 50, 50, 51
   50 teq(k) = tbnd * .5
      goto 25
   51 teq(k) = tbnd
c                                                                       
c         
c -----------------------                                               
c         
      goto 25
c FTEQ(2) ACTIVE                                                        
c         
   55 if (fteq(2) .lt. 0.0) goto 60
      k = 2
      temp = teq(2)
      fteq(3) = fteq(2)
      teq(2) = (teq(1) + temp) * .5
      teq(3) = temp
c FTEQ(2) PASSIVE                                                       
c         
      goto 25
   60 if (fteq(2) .lt. (-0.01)) goto 65
      tbnd = teq(2)
      goto 75
c AND FTEQ(3) ACTIVE                                                    
c         
   65 if (fteq(3) .lt. 0.0) goto 70
      k = 2
      temp = teq(2)
      fteq(1) = fteq(2)
      teq(2) = (teq(2) + teq(3)) * .5
      teq(1) = temp
c AND FTEQ(3) PASSIVE                                                   
c         
      goto 25
c ----------------------------                                          
c         
   70 tbnd = teq(3)
   75 if ((t * (2. + h)) - tbnd) 40, 26, 26
   26 if (iflag .eq. 1) goto 41
      iflag = 1
      tbnd = tbnd / 2.
      goto 1
   41 t = tbnd * (.5 - h)
c     WRITE(6,6792)T,TBND,FTEQ(K),XX(1),CKA,GG                          
c         
c6792 FORMAT(3X,'T,TBND,FTEQ(K),XX(1),CKA,GG',6E10.3)                   
c         
   40 continue
      return 
      end
c                                                                       
c         
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         
c                                                                       
c         
      subroutine inib(f, x, n, xall, nall, itp, h, ihess, ip, ivar, 
     &icall, se, ic)
c ----------------------------------------------------------------------
c-        
c                                                                       
c         
c OBTAINS INITIAL APPROXIMATION TO HESSIAN WHEN DESIRED                 
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c--       
      implicit double precision (o-z, a-h)
      common /mini/ tl(15, 15), d(15), g(15), gsave(15), y(15), p(15), 
     &tmin, fsmf, idif, isw, iret
      common /secc/ seca, seck, b(15, 15), saveb(15, 15), ifail
      dimension x(15), xall(15), itp(15)
c ----------------------------------------------------------------------
c-------  
      dimension se(15)
c READ B-MATRIX FROM INPUT FILE IF IHESS=2                              
c         
      if (icall - 1) 4, 6, 4

    4 if (ihess - 1) 5, 5, 1
    1 do 2 i = 1, n
    2 read(unit=ic, fmt=1000) (b(i,j),j = 1, i)
 1000 format(6e12.5)
      do 3 i = 1, n
      do 3 j = 1, i
    3 b(j,i) = b(i,j)
c ----------------------------------------------------------------------
c         
c GET FACTORIZATION OF INITIAL B-MATRIX                                 
c         
    5 continue
c IF FACTORIZATION FAILED , B WAS NOT P.D.; START WITH A  DIAGONAL      
c         
c APPROXIMATION TO B                                                    
c         
    6 call bldlt(f, x, n, xall, nall, itp, ihess, icall, ivar, se, ip, h
     &)
      if (icall - 1) 20, 10, 20
   10 do 15 i = 1, n
      do 12 j = 1, i
   12 tl(i,j) = 0.0
      tl(i,i) = 1.0
   15 d(i) = 1.0
c GET S.E. OF H AND Z FOR FINAL OUTPUT                                  
c         
   20 continue
      if (icall .eq. 0) call sehz(xall, nall, itp, se)
      return 
      end
c                                                                       
c         
c **********************************************************************
c         
c                                                                       
c         
      subroutine iter(xall, nall, itp, ip, iq, ictrl, h, asclq, tol, 
     &trupb)
c --------------------------------------------------------------------- 
c         
c                                                                       
c         
c SUBROUTINE CONTROLLING ITERATIVE PROCEDURE IN POINTER                 
c         
c                                                                       
c         
c --------------------------------------------------------------------- 
c         
c                                                                       
c         
      implicit double precision (o-z, a-h)
      common /bound/ ibndv
c -------------------------------------------------------------------   
c         
c -------------------------------------------------------------------   
c         
c GET READY TO CALL GEMINI                                              
c         
c H,TOL,TRUPB DEFINED IN ENTRY ( DEFAULT VALUE 1.E-04 )                 
c         
      dimension xall(15), itp(15), x(15), bnd(15, 2), se(15)
      ihess = 0
      ivar = 1
      ihx = 1
      ibnd = 1
c SET BOUNDARY CONDITIONS                                               
c         
      twoh = h + h
      bnd(1,1) = twoh
      bnd(1,2) = 10000.
      bnd(2,1) = -10000.
      bnd(2,2) = 10000.
      bnd(3,1) = twoh
      bnd(3,2) = 1.0 - twoh
      bnd(4,1) = twoh
      bnd(4,2) = 10000. / sqrt(asclq)
      bnd(5,1) = twoh
      bnd(5,2) = asclq - twoh
      bnd(6,1) = twoh
      bnd(6,2) = 10000.
      bnd(7,1) = .00001 + twoh
      bnd(7,2) = 10000.
      bnd(8,1) = twoh
      bnd(8,2) = 1.0 - twoh
      bnd(9,1) = twoh
      bnd(9,2) = 10000.
      bnd(10,1) = twoh
      bnd(10,2) = 1.0 - twoh
      bnd(11,1) = twoh
      bnd(11,2) = 1.0 - twoh
      bnd(12,1) = twoh
      bnd(12,2) = 1.0 - twoh
      bnd(13,1) = twoh
      bnd(13,2) = 0.707 - twoh
      do 6 i = 1, nall
      if (itp(i) .eq. 0) goto 6
      if (xall(i) .le. bnd(i,1)) xall(i) = bnd(i,1) + (2. * twoh)
      if (xall(i) .ge. bnd(i,2)) xall(i) = bnd(i,2) - (2. * twoh)
c SET VECTOR OF ITERATED PARAMETERS                                     
c         
    6 continue
      n = 0
      do 5 i = 1, nall
      if (itp(i) .ne. 1) goto 5
    1 n = n + 1
      x(n) = xall(i)
c --------------------------------------------------------------------  
c         
c START ITERATIVE PROCESS                                               
c         
    5 continue
      maxit = 20 * n
      ibndv = 0
      call gemini(h, trupb, nall, ibnd, ihess, xall, itp, n, x, bnd, f, 
     &tol, gtg, ptg, nit, nfe, ip, ivar, se, ic, maxit, ihx, idg)
      if (ibndv .eq. 0) goto 100
      write(unit=ip, fmt=1000) 
 1000 format(//10x,37hSTOP. PECULIAR BOUNDARY VALUE REACHED)
      write(unit=iq, fmt=1000) 
      goto 150
  100 call output(nall, xall, nit, nfe, f, ptg, idg, iq, ivar, se, ictrl
     &, h, tol, trupb, n, itp)
  150 continue
      return 
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         
c                                                                       
c         
c                                                                       
c         
c ********************************************************************  
c         
c                                                                       
c         
      subroutine like()
c --------------------------------------------------------------------- 
c         
c                                                                       
c         
c PERFORMS LIKELIHOOD CALCULATION FOR GIVEN FAMILY                      
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c         
      implicit double precision (o-z, a-h)
      integer idebug
      dimension c1(20), c2(20), trnsl(2), scale(2)
      dimension zcdc(10), gccdc(3), gcdd(3), isave(10)
      dimension save1(10), save2(10), save3(10)
      dimension sygc(3), syiljm(6), tempm(20)
      dimension pjm1(20), pjm2(20), pjm3(20), pijlm(6), qijlm(6)
      dimension pil1(20), pil2(20), pil3(20)
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, nev, idebug, np
      common /data/ tmf(9), tmm(9), qt(18), tmk(18),
     &zfa, zmo, qtfp, qtf, qtm, pi,
     &qtmp, qtcp, zphf, zphm, zphc, iph(18), iphz(18), irpf, irpm, irpk,
     &nc, isel, iphf, iphm, iphzf, iphzm, iphfp, iphmp, iphcp
      common /data2/ cxc2, cxc1, cdc, rrc, cxpf2, 
     &cxpf1, cdpf, cxpm2, cxpm1, cdpm, rrpf, rrpm, rrpc, 
     &cxd2, cxd1, cdd, alik, iquad, iccc, irrpf, irrpm, irrpc
      common /auxp1/ p(3), g(3), clm2, clm1, cxp2, cxp1(3), rcdg, 
     &cdp, z(10), vv, uu, gg, ca, ck, ea, ek, scals(4), sclq, rsqw
     &, notrs, ng1, nox
      common /intg/ xuu(20), auu(20), x(20), a(20), fzdqz(10), xdev, 
     &tvuu, fcaea, igot, nuu, nl, nm
c ----------------------------------------------------------------------
c         
c SCALING AND TRANSLATION AND ITR BUSINESS                              
c         
c FOR NOW, TRANSLATION DONE ONLY FOR QT GIVEN                           
c         
      common /trsm/ tr(18), tm(6, 3, 70)
      data sypoi1, sypoi2, chm1, chm2, chm3, chq1, chq2 /7*0./
      data tsy1, tsy2, tsy3, tsy4, tsy5, tsy6 /6*0./
      nn = max0(nm,nl)
      itr = 0
      trnsl(1) = 0.0
c GET APPROX. TO E(C(P)/KIDS)                                           
c         
      trnsl(2) = 0.0
      aaa = 0.0
      nefkd = 0
      iflag = 0
      do 826 kid = 1, nc
      iphzk = iphz(kid)
      if (iph(kid) - 2) 823, 826, 825
  823 if (iph(kid) - 1) 826, 824, 826
  824 aaa = aaa + fzdqz(iphzk)
      nefkd = nefkd + 1
      goto 826
  825 aaa = aaa + qt(kid)
      nefkd = nefkd + 1
      iflag = 1
  826 continue
c BRANCHING DEPENDING ON IQUAD                                          
c         
      if (nefkd .ne. 0) aaa = aaa / float(nefkd + 1)
c UNKNOWN X UNKNOWN MATING                                              
c         
      if (iquad - 2) 821, 821, 822
  822 trnsl(1) = 0.0
      trnsl(2) = fcaea * (aaa - tvuu)
      nn = nuu
      rscal = 1. / scals(3)
      do 8260 i = 1, nn
      c1(i) = 0.
 8260 c2(i) = (xuu(i) + trnsl(2)) * rscal
      goto 827
c OTHER MATING TYPES                                                    
c         
  821 continue
      scale(1) = scals(1)
      scale(2) = scals(1)
      if (iflag .eq. 1) scale(1) = scals(4)
c GET E(C/AFF.STAT.) INTO QTF, QTM                                      
c         
      if (iflag .eq. 1) scale(2) = scals(4)
      if (iphf .eq. 1) qtf = fzdqz(iphzf)
      if (iphm .eq. 1) qtm = fzdqz(iphzm)
      trnsl(1) = fcaea * (qtf - tvuu)
      if (iphf .lt. 3) trnsl(1) = trnsl(1) + (fcaea * (aaa - tvuu))
      if ((iphf .ge. 3) .or. (iphfp .ge. 3)) scale(1) = scals(4)
      trnsl(2) = fcaea * (qtm - tvuu)
      if (iphm .lt. 3) trnsl(2) = trnsl(2) + (fcaea * (aaa - tvuu))
      if ((iphm .ge. 3) .or. (iphmp .ge. 3)) scale(2) = scals(4)
c ONE PARENT ONLY                                                       
c         
      if (iquad - 2) 2400, 2401, 2400
 2401 scale(1) = 1.0
      rscal = 1. / scale(2)
      do 2402 i = 1, nn
      c1(i) = 0.
 2402 c2(i) = (x(i) + trnsl(2)) * rscal
      goto 827
 2400 if (abs(trnsl(1)) - xdev) 25, 25, 26
   25 trnsl(1) = 0.0
   26 if (abs(trnsl(2)) - xdev) 27, 27, 28
   27 trnsl(2) = 0.0
   28 ttrns = trnsl(1) + trnsl(2)
      if (ttrns) 30, 29, 30
   29 itr = 1
   30 continue
      rscal = 1. / (scale(1) * scale(2))
      do 810 i = 1, nn
      c1(i) = (x(i) + trnsl(1)) / scale(1)
c ----------------------------------------------------------------------
c         
  810 c2(i) = (x(i) + trnsl(2)) / scale(2)
  827 fpi = 1.0
      do 820 ikh = 1, nia
c *** DO THIS BETTER LATER ***                                          
c         
c     CALL BTIME                                                        
c         
  820 zcdc(ikh) = cdc * z(ikh)
      nll = nl
      nmm = nm
c     ONE DIMENSIONAL INTEGRATION FOR U X U ONLY                        
c         
      if (iquad .eq. 1) goto 835
      nll = 1
      nmm = nuu
      c1(1) = 0.
c -----------------------------------------------------------------     
c         
c                                                                       
c         
c     SECTION TO CALCULATE QZXW1 AND QZXW2                              
c         
c     USING INDIVIDUALS WITH QT AND AF STATUS                           
c         
c     WORKS ONLY IF PARAMETER W IS NOT EQUAL TO 0                       
c         
c                                                                       
c         
  835 frscal = fpi * rscal
      qzxw1 = 1.
      qzxw2 = 1.
      irtn = 1
c     MOTHER                                                            
c         
      if (rsqw .lt. 1.e-05) goto 4555
      if (iphm - 3) 4006, 4006, 4000
 4000 sy1 = zmo - qtm
      iphnow = iphm
 4001 xx = sy1 * rsqw
      ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 4002, 4003, 4003
 4002 t3 = 1.0 - t3
 4003 if (iphnow - 4) 4005, 4004, 4005
 4004 t3 = 1. - t3
 4005 qzxw2 = qzxw2 * t3
 4006 continue
c POINTER TO MOTHER                                                     
c         
      goto (4010, 4020, 4030, 4040, 4050), irtn
 4010 if (iphmp - 3) 4020, 4020, 4011
 4011 sy1 = zphm - qtmp
      iphnow = iphmp
      irtn = 2
c     FATHER                                                            
c         
      goto 4001
 4020 if (iphf - 3) 4030, 4030, 4021
 4021 sy1 = zfa - qtf
      iphnow = iphf
      irtn = 3
c POINTER TO FATHER                                                     
c         
      goto 4001
 4030 if (iphfp - 3) 4040, 4040, 4031
 4031 sy1 = zphf - qtfp
      iphnow = iphfp
      irtn = 4
c POINTER TO CHILDREN                                                   
c         
      goto 4001
 4040 if (iphcp - 3) 4050, 4050, 4041
 4041 sy1 = zphc - qtcp
      iphnow = iphcp
      irtn = 5
c     CHILDREN                                                          
c         
      goto 4001
 4050 irtn = 6
      qzxw1 = qzxw2
      do 4055 kid = 1, nc
      if (iph(kid) - 3) 4055, 4055, 4051
 4051 ksy = iphz(kid)
      sy1 = z(ksy) - qt(kid)
      iphnow = iph(kid)
      xx = sy1 * rsqw
      ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 4062, 4063, 4063
 4062 t3 = 1.0 - t3
 4063 if (iphnow - 4) 4052, 4064, 4052
 4064 t3 = 1. - t3
 4052 qzxw1 = qzxw1 * t3
 4055 continue
c     END OF QZXW1 AND QZXW2 SECTION                                    
c         
c -------------------------------------------------------------------   
c         
 4555 continue
      gcdd(1) = g(1) * cdd
      gcdd(2) = g(2) * cdd
      gcdd(3) = g(3) * cdd
      do 700 m = 1, nmm
      am = a(m)
      if (iquad .eq. 3) am = auu(m)
c -----------------------------------------------------------------     
c         
c -------------------------------------------------------------------   
c         
c ----------------------------------------------------------------      
c         
cALCULATE PJM(1,2,3=J,M)                                                
c         
      tempm(m) = (clm1 * exp((c2(m) * c2(m)) * clm2)) * am
      if (iphm - 2) 92, 701, 91
  701 pjm1(m) = p(1)
      pjm2(m) = p(2)
      pjm3(m) = p(3)
c     QUANTITATIVE DATA                                                 
c         
      goto 99
   91 sy1 = qtm - c2(m)
      sy2 = sy1 - g(1)
      pjm1(m) = cxp1(1) * exp((cxp2 * sy2) * sy2)
      sy2 = sy1 - g(2)
      pjm2(m) = cxp1(2) * exp((cxp2 * sy2) * sy2)
      sy2 = sy1 - g(3)
      pjm3(m) = cxp1(3) * exp((cxp2 * sy2) * sy2)
c     AFFECTION STATUS                                                  
c         
      goto 99
   92 sy1 = zmo - c2(m)
      xx = (sy1 - g(1)) * cdp
      iqfn = -1
   11 ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 12, 13, 13
   12 t3 = 1.0 - t3
   13 if (iphm) 93, 14, 93
   14 t3 = 1.0 - t3
   93 if (iqfn) 94, 95, 97
   94 pjm1(m) = t3 * p(1)
      xx = (sy1 - g(2)) * cdp
      iqfn = 0
      goto 11
   95 pjm2(m) = t3 * p(2)
      xx = (sy1 - g(3)) * cdp
      iqfn = 1
      goto 11
c END OF PJM                                                            
c         
cALCULATE POINTER TO MOTHER                                             
c         
   97 pjm3(m) = t3 * p(3)
c     QUANTITATIVE DATA                                                 
c         
   99 if (iphmp - 2) 82, 702, 81
   81 sy1 = qtmp - (c2(m) * rrpm)
      sy2 = sy1 - g(1)
      sypoi1 = cxpm1 * exp((cxpm2 * sy2) * sy2)
      sy2 = sy1 - g(2)
      sypoi2 = cxpm1 * exp((cxpm2 * sy2) * sy2)
      sy2 = sy1 - g(3)
      sypoi3 = cxpm1 * exp((cxpm2 * sy2) * sy2)
c     AFFECTION STATUS                                                  
c         
      goto 89
   82 sy1 = zphm - (c2(m) * rrpm)
      xx = (sy1 - g(1)) * cdpm
      iqfn = -1
   83 ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 84, 85, 85
   84 t3 = 1.0 - t3
   85 if (iphmp) 87, 86, 87
   86 t3 = 1.0 - t3
   87 if (iqfn) 77, 78, 79
   77 sypoi1 = t3
      xx = (sy1 - g(2)) * cdpm
      iqfn = 0
      goto 83
   78 sypoi2 = t3
      xx = (sy1 - g(3)) * cdpm
      iqfn = 1
      goto 83
   79 sypoi3 = t3
   89 poim1 = ((tmm(1) * sypoi1) + (tmm(2) * sypoi2)) + (tmm(3) * sypoi3
     &)
      poim2 = ((tmm(4) * sypoi1) + (tmm(5) * sypoi2)) + (tmm(6) * sypoi3
     &)
      poim3 = ((tmm(7) * sypoi1) + (tmm(8) * sypoi2)) + (tmm(9) * sypoi3
     &)
      pjm1(m) = pjm1(m) * poim1
      pjm2(m) = pjm2(m) * poim2
c END OF POINTER TO MOTHER                                              
c         
cALCULATE PIL(1,2,3=I,L)                                                
c         
      pjm3(m) = pjm3(m) * poim3
  702 if (iphf - 2) 102, 110, 101
  110 pil1(m) = p(1)
      pil2(m) = p(2)
      pil3(m) = p(3)
c     QUANTITATIVE DATA                                                 
c         
      goto 116
  101 sy1 = qtf - c1(m)
      sy2 = sy1 - g(1)
      pil1(m) = cxp1(1) * exp((cxp2 * sy2) * sy2)
      sy2 = sy1 - g(2)
      pil2(m) = cxp1(2) * exp((cxp2 * sy2) * sy2)
      sy2 = sy1 - g(3)
      pil3(m) = cxp1(3) * exp((cxp2 * sy2) * sy2)
c     AFFECTION STATUS                                                  
c         
      goto 116
  102 sy1 = zfa - c1(m)
      xx = (sy1 - g(1)) * cdp
      iqfn = -1
    1 ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 2, 3, 3
    2 t3 = 1.0 - t3
    3 if (iphf) 103, 4, 103
    4 t3 = 1.0 - t3
  103 if (iqfn) 104, 105, 107
  104 pil1(m) = t3 * p(1)
      xx = (sy1 - g(2)) * cdp
      iqfn = 0
      goto 1
  105 pil2(m) = t3 * p(2)
      xx = (sy1 - g(3)) * cdp
      iqfn = 1
      goto 1
c END OF PIL                                                            
c         
cALCULATE POINTER FOR FATHER                                            
c         
  107 pil3(m) = t3 * p(3)
c     QUANTITATIVE DATA                                                 
c         
  116 if (iphfp - 2) 402, 700, 401
  401 sy1 = qtfp - (c1(m) * rrpf)
      sy2 = sy1 - g(1)
      sypoi1 = cxpf1 * exp((cxpf2 * sy2) * sy2)
      sy2 = sy1 - g(2)
      sypoi2 = cxpf1 * exp((cxpf2 * sy2) * sy2)
      sy2 = sy1 - g(3)
      sypoi3 = cxpf1 * exp((cxpf2 * sy2) * sy2)
c     AFFECTION STATUS                                                  
c         
      goto 409
  402 sy1 = zphf - (c1(m) * rrpf)
      xx = (sy1 - g(1)) * cdpf
      iqfn = -1
  403 ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 404, 405, 405
  404 t3 = 1.0 - t3
  405 if (iphfp) 407, 406, 407
  406 t3 = 1.0 - t3
  407 if (iqfn) 417, 418, 419
  417 sypoi1 = t3
      xx = (sy1 - g(2)) * cdpf
      iqfn = 0
      goto 403
  418 sypoi2 = t3
      xx = (sy1 - g(3)) * cdpf
      iqfn = 1
      goto 403
  419 sypoi3 = t3
  409 poif1 = ((tmf(1) * sypoi1) + (tmf(2) * sypoi2)) + (tmf(3) * sypoi3
     &)
      poif2 = ((tmf(4) * sypoi1) + (tmf(5) * sypoi2)) + (tmf(6) * sypoi3
     &)
      poif3 = ((tmf(7) * sypoi1) + (tmf(8) * sypoi2)) + (tmf(9) * sypoi3
     &)
      pil1(m) = pil1(m) * poif1
      pil2(m) = pil2(m) * poif2
c END OF POINTER TO FATHER                                              
c         
      pil3(m) = pil3(m) * poif3
  700 continue
      slm1 = 0.
c                                                                       
c         
c     M A I N   L O O P  STARTS HERE                                    
c         
c                                                                       
c         
c **** TEMP ****                                                        
c         
c     NTIME=0                                                           
c         
      slm2 = 0.
      do 100 l = 1, nll
      if (iquad - 2) 830, 829, 829
  829 templ = 1.0
      goto 831
  830 templ = (clm1 * exp((c1(l) * c1(l)) * clm2)) * a(l)
  831 m1 = ((l - 1) * itr) + 1
      do 90 m = m1, nmm
      sij1 = 0.
      sij2 = 0.
      isave(1) = 0
      isave(2) = 0
      isave(3) = 0
      isave(4) = 0
      isave(5) = 0
      isave(6) = 0
      isave(7) = 0
      isave(8) = 0
      isave(9) = 0
      isave(10) = 0
cALCULATE POINTER TO CHILDREN                                           
c         
      c1l2m = c1(l) + c2(m)
      if (iphcp - 2) 422, 432, 421
  432 qijlm(1) = 1.
      qijlm(2) = 1.
      qijlm(3) = 1.
      qijlm(4) = 1.
      qijlm(5) = 1.
      qijlm(6) = 1.
      poic1 = 1.
      poic2 = 1.
      poic3 = 1.
      poic4 = 1.
      poic5 = 1.
      poic6 = 1.
c     QUANTITATIVE DATA                                                 
c         
      goto 431
  421 t2 = qtcp - (c1l2m * rrpc)
      t1 = t2 - g(1)
      sypoi1 = cxd1 * exp((t1 * t1) * cxd2)
      t1 = t2 - g(2)
      sypoi2 = cxd1 * exp((t1 * t1) * cxd2)
      t1 = t2 - g(3)
      sypoi3 = cxd1 * exp((t1 * t1) * cxd2)
c     AFFECTION STATUS                                                  
c         
      goto 429
  422 temp = (zphc - (c1l2m * rrpc)) * cdd
      xx = temp - gcdd(1)
      iqfn = -1
  423 ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 424, 425, 425
  424 t3 = 1.0 - t3
  425 if (iphcp) 427, 426, 427
  426 t3 = 1.0 - t3
  427 if (iqfn) 437, 438, 439
  437 sypoi1 = t3
      xx = temp - gcdd(2)
      iqfn = 0
      goto 423
  438 sypoi2 = t3
      xx = temp - gcdd(3)
      iqfn = 1
      goto 423
  439 sypoi3 = t3
  429 poic1 = ((tmk(1) * sypoi1) + (tmk(2) * sypoi2)) + (tmk(3) * sypoi3
     &)
      poic2 = ((tmk(4) * sypoi1) + (tmk(5) * sypoi2)) + (tmk(6) * sypoi3
     &)
      poic3 = ((tmk(7) * sypoi1) + (tmk(8) * sypoi2)) + (tmk(9) * sypoi3
     &)
      poic4 = ((tmk(10) * sypoi1) + (tmk(11) * sypoi2)) + (tmk(12) * 
     &sypoi3)
      poic5 = ((tmk(13) * sypoi1) + (tmk(14) * sypoi2)) + (tmk(15) * 
     &sypoi3)
c END OF POINTER TO CHILDREN                                            
c         
cALCULATE QIJLM(1,2,3,4,5,6=IJ)                                         
c         
      poic6 = ((tmk(16) * sypoi1) + (tmk(17) * sypoi2)) + (tmk(18) * 
     &sypoi3)
      qijlm(1) = poic1
      qijlm(2) = poic2
      qijlm(3) = poic3
      qijlm(4) = poic4
      qijlm(5) = poic5
c CACULATE PLM                                                          
c         
      qijlm(6) = poic6
  431 plm = templ * tempm(m)
      clm = c1l2m * rrc
      sygc(1) = g(1) + clm
      sygc(2) = g(2) + clm
      sygc(3) = g(3) + clm
      gccdc(1) = cdc * sygc(1)
      gccdc(2) = cdc * sygc(2)
cALCULATE PIJLM(1,2,3,4,5,6=IJ)                                         
c         
      gccdc(3) = cdc * sygc(3)
      pijlm(1) = qijlm(1)
      pijlm(2) = qijlm(2)
      pijlm(3) = qijlm(3)
      pijlm(4) = qijlm(4)
      pijlm(5) = qijlm(5)
c     KIDS  LOOP                                                        
c         
      pijlm(6) = qijlm(6)
      do 200 kid = 1, nc
      ksy = iphz(kid)
c     QUANTITATIVE DATA                                                 
c         
      if (iph(kid) - 2) 302, 200, 201
  201 t1 = qt(kid) - sygc(1)
      chm1 = cxc1 * exp((t1 * t1) * cxc2)
      t1 = qt(kid) - sygc(2)
      chm2 = cxc1 * exp((t1 * t1) * cxc2)
      t1 = qt(kid) - sygc(3)
      chm3 = cxc1 * exp((t1 * t1) * cxc2)
      if (isel) 2500, 2500, 301
cCCC  IF(KSY.EQ.1.AND.ISAVE(1).EQ.1)GO TO 252                           
c         
 2500 if (nox) 2250, 2250, 250
  301 if (ksy - isave(1)) 305, 252, 305
  305 if (isave(ksy)) 306, 202, 306
  306 chq1 = pi * save1(ksy)
      chq2 = pi * save2(ksy)
      chq3 = pi * save3(ksy)
c     AFFECTION STATUS                                                  
c         
      if (nox) 251, 2251, 251
  302 if (isave(ksy)) 307, 202, 307
  307 if (iph(kid)) 308, 303, 308
  308 chm1 = save1(ksy)
      chm2 = save2(ksy)
      chm3 = save3(ksy)
      if (isel) 2002, 2002, 309
  309 chq1 = chm1 * pi
      chq2 = chm2 * pi
      chq3 = chm3 * pi
      if (nox) 251, 2251, 251
  303 chm1 = 1. - save1(ksy)
      chm2 = 1. - save2(ksy)
      chm3 = 1. - save3(ksy)
      if (isel) 2002, 2002, 311
  311 chq1 = save1(ksy) * pi
      chq2 = save2(ksy) * pi
      chq3 = save3(ksy) * pi
      if (nox) 251, 2251, 251
 2002 if (nox) 250, 2250, 250
  202 xx = zcdc(ksy) - gccdc(1)
      iqfn = -1
  211 ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 212, 213, 213
  212 t3 = 1.0 - t3
  213 t4 = t3 * pi
      t33 = t3
      if (iph(kid) - 1) 214, 215, 2115
 2115 if (iqfn) 263, 264, 265
  214 t3 = 1.0 - t3
  215 if (iqfn) 203, 204, 205
  203 chm1 = t3
  263 chq1 = t4
      save1(ksy) = t33
      xx = zcdc(ksy) - gccdc(2)
      iqfn = 0
      goto 211
  204 chm2 = t3
  264 chq2 = t4
      save2(ksy) = t33
      xx = zcdc(ksy) - gccdc(3)
      iqfn = 1
      goto 211
  205 chm3 = t3
  265 chq3 = t4
      save3(ksy) = t33
c     PARAMETERS X OR Q = 0 , DON'T NEED TR VECTOR.                     
c         
      if (nox) 2001, 2000, 2001
 2000 if (isel) 2266, 2266, 2267
 2266 isave(ksy) = 1
      goto 2250
 2267 if (ksy - 1) 2255, 2310, 2255
 2310 isave(1) = 1
      tsy1 = 1. - chq1
      tsy2 = 1. - (.5 * (chq1 + chq2))
      tsy3 = 1. - chq2
      tsy4 = 1. - ((.25 * (chq1 + chq3)) + (.5 * chq2))
      tsy5 = 1. - (.5 * (chq2 + chq3))
      tsy6 = 1. - chq3
      goto 252
 2255 isave(ksy) = 1
 2251 qijlm(1) = qijlm(1) * (1. - chq1)
      qijlm(2) = qijlm(2) * (1. - (.5 * (chq1 + chq2)))
      qijlm(3) = qijlm(3) * (1. - chq2)
      qijlm(4) = qijlm(4) * (1. - ((.25 * (chq1 + chq3)) + (.5 * chq2)))
      qijlm(5) = qijlm(5) * (1. - (.5 * (chq2 + chq3)))
      qijlm(6) = qijlm(6) * (1. - chq3)
 2250 pijlm(1) = pijlm(1) * chm1
      pijlm(2) = (pijlm(2) * .5) * (chm1 + chm2)
      pijlm(3) = pijlm(3) * chm2
      pijlm(4) = pijlm(4) * ((.25 * (chm1 + chm3)) + (.5 * chm2))
      pijlm(5) = (pijlm(5) * .5) * (chm2 + chm3)
      pijlm(6) = pijlm(6) * chm3
c     X AND Q NOT EQUAL TO 0, NEED TR VECTOR.                           
c         
      goto 200
 2001 if (isel) 266, 266, 267
  266 isave(ksy) = 1
c     ONE LIABILITY CLASS SPEED ENHANCER                                
c         
      goto 250
  267 if (ksy - 1) 255, 310, 255
  310 isave(1) = 1
      tsy1 = ((1. - (tr(1) * chq1)) - (tr(2) * chq2)) - (tr(3) * chq3)
      tsy2 = ((1. - (tr(4) * chq1)) - (tr(5) * chq2)) - (tr(6) * chq3)
      tsy3 = ((1. - (tr(7) * chq1)) - (tr(8) * chq2)) - (tr(9) * chq3)
      tsy4 = ((1. - (tr(10) * chq1)) - (tr(11) * chq2)) - (tr(12) * chq3
     &)
      tsy5 = ((1. - (tr(13) * chq1)) - (tr(14) * chq2)) - (tr(15) * chq3
     &)
      tsy6 = ((1. - (tr(16) * chq1)) - (tr(17) * chq2)) - (tr(18) * chq3
     &)
  252 qijlm(1) = qijlm(1) * tsy1
      qijlm(2) = qijlm(2) * tsy2
      qijlm(3) = qijlm(3) * tsy3
      qijlm(4) = qijlm(4) * tsy4
      qijlm(5) = qijlm(5) * tsy5
      qijlm(6) = qijlm(6) * tsy6
      if (nox) 250, 2250, 250
  255 isave(ksy) = 1
  251 qijlm(1) = qijlm(1) * (((1. - (tr(1) * chq1)) - (tr(2) * chq2)) - 
     &(tr(3) * chq3))
      qijlm(2) = qijlm(2) * (((1. - (tr(4) * chq1)) - (tr(5) * chq2)) - 
     &(tr(6) * chq3))
      qijlm(3) = qijlm(3) * (((1. - (tr(7) * chq1)) - (tr(8) * chq2)) - 
     &(tr(9) * chq3))
      qijlm(4) = qijlm(4) * (((1. - (tr(10) * chq1)) - (tr(11) * chq2))
     & - (tr(12) * chq3))
      qijlm(5) = qijlm(5) * (((1. - (tr(13) * chq1)) - (tr(14) * chq2))
     & - (tr(15) * chq3))
      qijlm(6) = qijlm(6) * (((1. - (tr(16) * chq1)) - (tr(17) * chq2))
     & - (tr(18) * chq3))
  250 pijlm(1) = pijlm(1) * (((tr(1) * chm1) + (tr(2) * chm2)) + (tr(3)
     & * chm3))
      pijlm(2) = pijlm(2) * (((tr(4) * chm1) + (tr(5) * chm2)) + (tr(6)
     & * chm3))
      pijlm(3) = pijlm(3) * (((tr(7) * chm1) + (tr(8) * chm2)) + (tr(9)
     & * chm3))
      pijlm(4) = pijlm(4) * (((tr(10) * chm1) + (tr(11) * chm2)) + (tr(
     &12) * chm3))
      pijlm(5) = pijlm(5) * (((tr(13) * chm1) + (tr(14) * chm2)) + (tr(
     &15) * chm3))
      pijlm(6) = pijlm(6) * (((tr(16) * chm1) + (tr(17) * chm2)) + (tr(
     &18) * chm3))
c END OF PILJM                                                          
c         
  200 continue
      if (isel) 274, 274, 273
  273 qijlm(1) = poic1 - qijlm(1)
      qijlm(2) = poic2 - qijlm(2)
      qijlm(3) = poic3 - qijlm(3)
      qijlm(4) = poic4 - qijlm(4)
      qijlm(5) = poic5 - qijlm(5)
c END OF QIJLM                                                          
c         
      qijlm(6) = poic6 - qijlm(6)
  274 syiljm(1) = pil1(l) * pjm1(m)
      syiljm(2) = (pil1(l) * pjm2(m)) + (pil2(l) * pjm1(m))
      syiljm(3) = (pil1(l) * pjm3(m)) + (pil3(l) * pjm1(m))
      syiljm(4) = pil2(l) * pjm2(m)
      syiljm(5) = (pil2(l) * pjm3(m)) + (pil3(l) * pjm2(m))
      syiljm(6) = pil3(l) * pjm3(m)
      if (itr * (m - l)) 272, 272, 271
  271 syiljm(1) = syiljm(1) + (pil1(m) * pjm1(l))
      syiljm(2) = (syiljm(2) + (pil1(m) * pjm2(l))) + (pil2(m) * pjm1(l)
     &)
      syiljm(3) = (syiljm(3) + (pil1(m) * pjm3(l))) + (pil3(m) * pjm1(l)
     &)
      syiljm(4) = syiljm(4) + (pil2(m) * pjm2(l))
      syiljm(5) = (syiljm(5) + (pil2(m) * pjm3(l))) + (pil3(m) * pjm2(l)
     &)
      syiljm(6) = syiljm(6) + (pil3(m) * pjm3(l))
  272 sij1 = (((((pijlm(1) * syiljm(1)) + (pijlm(2) * syiljm(2))) + (
     &pijlm(3) * syiljm(3))) + (pijlm(4) * syiljm(4))) + (pijlm(5) * 
     &syiljm(5))) + (pijlm(6) * syiljm(6))
      if (iccc) 2721, 2722, 2722
 2721 sij2 = ((((((qijlm(1) * p(1)) * p(1)) + (((qijlm(2) * p(1)) * p(2)
     &) * 2.)) + (((qijlm(3) * p(1)) * p(3)) * 2.)) + ((qijlm(4) * p(2))
     & * p(2))) + (((qijlm(5) * p(2)) * p(3)) * 2.)) + ((qijlm(6) * p(3)
     &) * p(3))
      goto 2723
 2722 sij2 = (((((qijlm(1) * syiljm(1)) + (qijlm(2) * syiljm(2))) + (
     &qijlm(3) * syiljm(3))) + (qijlm(4) * syiljm(4))) + (qijlm(5) * 
     &syiljm(5))) + (qijlm(6) * syiljm(6))
 2723 continue
      slm1 = slm1 + (sij1 * plm)
c     WRITE(61,1701)SIJ1,PLM,PIJLM,SYILJM                               
c         
c1701 FORMAT(1X,'HELP',3E18.11)                                         
c         
c **** TEMP ****                                                        
c         
c     NTIME=NTIME+1                                                     
c         
      slm2 = slm2 + (sij2 * plm)
   90 continue
c **** TEMP ****                                                        
c         
c     WRITE(6,3939)IQUAD,NTIME,NLL,NM,NL,NN                             
c         
c3939 FORMAT(1X,'IQUAD,NTIME',2I5,'NLL,NM',2I5,'NL,NN',2I5)             
c         
  100 continue
      slm1 = slm1 * frscal
c     WRITE(6,3939)IQUAD,NL,NM,NN,SLM1,SLM2,QZXW1,QZXW2                 
c         
c3939 FORMAT(1X,4I3,4(1X,E15.8))                                        
c         
      slm2 = slm2 * rscal
      slm1 = slm1 * qzxw1
      slm2 = slm2 * qzxw2
      if (isel) 812, 813, 813
  812 pfam = slm1
      goto 814
  813 pfam = slm1 / slm2
  814 continue
c     CALL ETIME                                                        
c         
c PRINT OUTPUT                                                          
c         
c     WRITE(61,7001)ALIK                                                
c         
c7001 FORMAT(1X,'ALOG',E18.11)                                          
c         
      alik = - log(pfam)
      return 
      end
c                                                                       
c         
c **********************************************************************
c**       
c                                                                       
c         
c ********************************************************************  
c         
c                                                                       
c         
c NOT LISTED HERE ARE FUNCTIONS QFN,FN AND FUNCTIONS REQUIRED FOR       
c         
c OBTENTION OF WEIGHTS AND ABSCISSAS FOR NUMERICAL INTEGRATION.         
c         
c THIS ROUTINES ARE IN LIBRARY JMLIB.                                   
c         
c                                                                       
c         
c ********************************************************************* 
c         
c                                                                       
c         
c **********************************************************************
c**       
c                                                                       
c         
      subroutine likeq0()
c --------------------------------------------------------------------- 
c         
c                                                                       
c         
c PERFORMS LIKELIHOOD CALCULATION FOR GIVEN FAMILY WHEN Q=0.0           
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c         
      implicit double precision (o-z, a-h)
      integer idebug
      dimension c1(20), c2(20), trnsl(2), scale(2)
      dimension zcdc(10), isave(10)
      dimension save1(10)
      dimension tempm(20)
      dimension pjm1(20)
      dimension pil1(20)
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, nev, idebug, np
      common /data/ tmf(9), tmm(9), qt(18), tmk(18),
     &zfa, zmo, qtfp, qtf, qtm, pi,
     &qtmp, qtcp, zphf, zphm, zphc, iph(18), iphz(18), irpf, irpm, irpk,
     &nc, isel, iphf, iphm, iphzf, iphzm, iphfp, iphmp, iphcp
      common /data2/ cxc2, cxc1, cdc, rrc, cxpf2, 
     &cxpf1, cdpf, cxpm2, cxpm1, cdpm, rrpf, rrpm, rrpc, 
     &cxd2, cxd1, cdd, alik, iquad, iccc, irrpf, irrpm, irrpc
      common /auxp1/ p(3), g(3), clm2, clm1, cxp2, cxp1(3), rcdg, 
     &cdp, z(10), vv, uu, gg, ca, ck, ea, ek, scals(4), sclq, rsqw
     &, notrs, ng1, nox
      common /intg/ xuu(20), auu(20), x(20), a(20), fzdqz(10), xdev, 
     &tvuu, fcaea, igot, nuu, nl, nm
c ----------------------------------------------------------------------
c         
c SCALING AND TRANSLATION AND ITR BUSINESS                              
c         
c FOR NOW, TRANSLATION DONE ONLY FOR QT GIVEN                           
c         
      common /trsm/ tr(18), tm(6, 3, 70)
      data chm1, tsy1 /2*0./
      nn = max0(nm,nl)
      itr = 0
      trnsl(1) = 0.0
c GET APPROX. TO E(C(P)/KIDS)                                           
c         
      trnsl(2) = 0.0
      aaa = 0.0
      nefkd = 0
      iflag = 0
      do 826 kid = 1, nc
      iphzk = iphz(kid)
      if (iph(kid) - 2) 823, 826, 825
  823 if (iph(kid) - 1) 826, 824, 826
  824 aaa = aaa + fzdqz(iphzk)
      nefkd = nefkd + 1
      goto 826
  825 aaa = aaa + qt(kid)
      nefkd = nefkd + 1
      iflag = 1
  826 continue
c BRANCHING DEPENDING ON IQUAD                                          
c         
      if (nefkd .ne. 0) aaa = aaa / float(nefkd + 1)
c UNKNOWN X UNKNOWN MATING                                              
c         
      if (iquad - 2) 821, 821, 822
  822 trnsl(1) = 0.0
      trnsl(2) = fcaea * (aaa - tvuu)
      nn = nuu
      rscal = 1. / scals(3)
      do 8260 i = 1, nn
      c1(i) = 0.
 8260 c2(i) = (xuu(i) + trnsl(2)) * rscal
      goto 827
c OTHER MATING TYPES                                                    
c         
  821 continue
      scale(1) = scals(1)
      scale(2) = scals(1)
      if (iflag .eq. 1) scale(1) = scals(4)
c GET E(C/AFF.STAT.) INTO QTF, QTM                                      
c         
      if (iflag .eq. 1) scale(2) = scals(4)
      if (iphf .eq. 1) qtf = fzdqz(iphzf)
      if (iphm .eq. 1) qtm = fzdqz(iphzm)
      trnsl(1) = fcaea * (qtf - tvuu)
      if (iphf .lt. 3) trnsl(1) = trnsl(1) + (fcaea * (aaa - tvuu))
      if ((iphf .ge. 3) .or. (iphfp .ge. 3)) scale(1) = scals(4)
      trnsl(2) = fcaea * (qtm - tvuu)
      if (iphm .lt. 3) trnsl(2) = trnsl(2) + (fcaea * (aaa - tvuu))
      if ((iphm .ge. 3) .or. (iphmp .ge. 3)) scale(2) = scals(4)
c ONE PARENT ONLY                                                       
c         
      if (iquad - 2) 2400, 2401, 2400
 2401 scale(1) = 1.0
      rscal = 1. / scale(2)
      do 2402 i = 1, nn
      c1(i) = 0.
 2402 c2(i) = (x(i) + trnsl(2)) * rscal
      goto 827
 2400 if (abs(trnsl(1)) - xdev) 25, 25, 26
   25 trnsl(1) = 0.0
   26 if (abs(trnsl(2)) - xdev) 27, 27, 28
   27 trnsl(2) = 0.0
   28 ttrns = trnsl(1) + trnsl(2)
      if (ttrns) 30, 29, 30
   29 itr = 1
   30 continue
      rscal = 1. / (scale(1) * scale(2))
      do 810 i = 1, nn
      c1(i) = (x(i) + trnsl(1)) / scale(1)
c ----------------------------------------------------------------------
c         
  810 c2(i) = (x(i) + trnsl(2)) / scale(2)
  827 fpi = 1.0
      do 820 ikh = 1, nia
  820 zcdc(ikh) = cdc * z(ikh)
      nll = nl
      nmm = nm
c     ONE DIMENSIONAL INTEGRATION FOR U X U ONLY                        
c         
      if (iquad .eq. 1) goto 835
      nll = 1
      nmm = nuu
c *** DO THIS BETTER LATER ***                                          
c         
      c1(1) = 0.
c     CALL BTIME                                                        
c         
c -----------------------------------------------------------------     
c         
c                                                                       
c         
c     SECTION TO CALCULATE QZXW1 AND QZXW2                              
c         
c     USING INDIVIDUALS WITH QT AND AF STATUS                           
c         
c     WORKS ONLY IF PARAMETER W IS NOT EQUAL TO 0                       
c         
c                                                                       
c         
  835 continue
      qzxw1 = 1.
      qzxw2 = 1.
      irtn = 1
c     MOTHER                                                            
c         
      if (rsqw .lt. 1.e-05) goto 4555
      if (iphm - 3) 4006, 4006, 4000
 4000 sy1 = zmo - qtm
      iphnow = iphm
 4001 xx = sy1 * rsqw
      ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 4002, 4003, 4003
 4002 t3 = 1.0 - t3
 4003 if (iphnow - 4) 4005, 4004, 4005
 4004 t3 = 1. - t3
 4005 qzxw2 = qzxw2 * t3
 4006 continue
c POINTER TO MOTHER                                                     
c         
      goto (4010, 4020, 4030, 4040, 4050), irtn
 4010 if (iphmp - 3) 4020, 4020, 4011
 4011 sy1 = zphm - qtmp
      iphnow = iphmp
      irtn = 2
c     FATHER                                                            
c         
      goto 4001
 4020 if (iphf - 3) 4030, 4030, 4021
 4021 sy1 = zfa - qtf
      iphnow = iphf
      irtn = 3
c POINTER TO FATHER                                                     
c         
      goto 4001
 4030 if (iphfp - 3) 4040, 4040, 4031
 4031 sy1 = zphf - qtfp
      iphnow = iphfp
      irtn = 4
c POINTER TO CHILDREN                                                   
c         
      goto 4001
 4040 if (iphcp - 3) 4050, 4050, 4041
 4041 sy1 = zphc - qtcp
      iphnow = iphcp
      irtn = 5
c     CHILDREN                                                          
c         
      goto 4001
 4050 irtn = 6
      qzxw1 = qzxw2
      do 4055 kid = 1, nc
      if (iph(kid) - 3) 4055, 4055, 4051
 4051 ksy = iphz(kid)
      sy1 = z(ksy) - qt(kid)
      iphnow = iph(kid)
      xx = sy1 * rsqw
      ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 4062, 4063, 4063
 4062 t3 = 1.0 - t3
 4063 if (iphnow - 4) 4052, 4064, 4052
 4064 t3 = 1. - t3
 4052 qzxw1 = qzxw1 * t3
 4055 continue
c     END OF QZXW1 AND QZXW2 SECTION                                    
c         
c -------------------------------------------------------------------   
c         
 4555 continue
      gcdd = g(1) * cdd
      do 700 m = 1, nmm
      am = a(m)
      if (iquad .eq. 3) am = auu(m)
c -----------------------------------------------------------------     
c         
c -------------------------------------------------------------------   
c         
cALCULATE PJM(1,2,3=J,M)                                                
c         
      tempm(m) = (clm1 * exp((c2(m) * c2(m)) * clm2)) * am
      if (iphm - 2) 92, 701, 91
  701 pjm1(m) = p(1)
c     QUANTITATIVE DATA                                                 
c         
      goto 99
   91 sy2 = (qtm - g(1)) - c2(m)
      pjm1(m) = cxp1(1) * exp((cxp2 * sy2) * sy2)
c     AFFECTION STATUS                                                  
c         
      goto 99
   92 xx = ((zmo - g(1)) - c2(m)) * cdp
      ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 12, 13, 13
   12 t3 = 1.0 - t3
   13 if (iphm) 93, 14, 93
   14 t3 = 1.0 - t3
c END OF PJM                                                            
c         
cALCULATE POINTER TO MOTHER                                             
c         
   93 pjm1(m) = t3 * p(1)
c     QUANTITATIVE DATA                                                 
c         
   99 if (iphmp - 2) 82, 702, 81
   81 sy2 = (qtmp - g(1)) - (c2(m) * rrpm)
      sypoi1 = cxpm1 * exp((cxpm2 * sy2) * sy2)
c     AFFECTION STATUS                                                  
c         
      goto 89
   82 xx = ((zphm - g(1)) - (c2(m) * rrpm)) * cdpm
      ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 84, 85, 85
   84 t3 = 1.0 - t3
   85 if (iphmp) 87, 86, 87
   86 t3 = 1.0 - t3
   87 sypoi1 = t3
c END OF POINTER TO MOTHER                                              
c         
cALCULATE PIL(1,2,3=I,L)                                                
c         
   89 pjm1(m) = pjm1(m) * sypoi1
  702 if (iphf - 2) 102, 110, 101
  110 pil1(m) = p(1)
c     QUANTITATIVE DATA                                                 
c         
      goto 116
  101 sy2 = (qtf - g(1)) - c1(m)
      pil1(m) = cxp1(1) * exp((cxp2 * sy2) * sy2)
c     AFFECTION STATUS                                                  
c         
      goto 116
  102 xx = ((zfa - g(1)) - c1(m)) * cdp
      ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 2, 3, 3
    2 t3 = 1.0 - t3
    3 if (iphf) 103, 4, 103
    4 t3 = 1.0 - t3
c END OF PIL                                                            
c         
cALCULATE POINTER FOR FATHER                                            
c         
  103 pil1(m) = t3 * p(1)
c     QUANTITATIVE DATA                                                 
c         
  116 if (iphfp - 2) 402, 700, 401
  401 sy2 = (qtfp - g(1)) - (c1(m) * rrpf)
      sypoi1 = cxpf1 * exp((cxpf2 * sy2) * sy2)
c     AFFECTION STATUS                                                  
c         
      goto 409
  402 xx = ((zphf - g(1)) - (c1(m) * rrpf)) * cdpf
      ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 404, 405, 405
  404 t3 = 1.0 - t3
  405 if (iphfp) 407, 406, 407
  406 t3 = 1.0 - t3
  407 sypoi1 = t3
c END OF POINTER TO FATHER                                              
c         
  409 pil1(m) = pil1(m) * sypoi1
  700 continue
      slm1 = 0.
c                                                                       
c         
c     M A I N   L O O P  STARTS HERE                                    
c         
c                                                                       
c         
      slm2 = 0.
      do 100 l = 1, nll
      if (iquad - 2) 830, 829, 829
  829 templ = 1.0
      goto 831
  830 templ = (clm1 * exp((c1(l) * c1(l)) * clm2)) * a(l)
  831 m1 = ((l - 1) * itr) + 1
      do 90 m = m1, nmm
      sij1 = 0.
      sij2 = 0.
      isave(1) = 0
      isave(2) = 0
      isave(3) = 0
      isave(4) = 0
      isave(5) = 0
      isave(6) = 0
      isave(7) = 0
      isave(8) = 0
      isave(9) = 0
      isave(10) = 0
cALCULATE POINTER TO CHILDREN                                           
c         
      c1l2m = c1(l) + c2(m)
      if (iphcp - 2) 422, 432, 421
  432 qijlm = 1.
      poic1 = 1.
c     QUANTITATIVE DATA                                                 
c         
      goto 431
  421 t1 = (qtcp - g(1)) - (c1l2m * rrpc)
      sypoi1 = cxd1 * exp((t1 * t1) * cxd2)
c     AFFECTION STATUS                                                  
c         
      goto 429
  422 xx = ((zphc - (c1l2m * rrpc)) * cdd) - gcdd
      ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 424, 425, 425
  424 t3 = 1.0 - t3
  425 if (iphcp) 427, 426, 427
  426 t3 = 1.0 - t3
  427 sypoi1 = t3
c END OF POINTER TO CHILDREN                                            
c         
cALCULATE QIJLM(1,2,3,4,5,6=IJ)                                         
c         
  429 poic1 = sypoi1
c CACULATE PLM                                                          
c         
      qijlm = poic1
  431 plm = templ * tempm(m)
      clm = c1l2m * rrc
      sygc = g(1) + clm
cALCULATE PIJLM(1,2,3,4,5,6=IJ)                                         
c         
      gccdc = cdc * sygc
c     KIDS  LOOP                                                        
c         
      pijlm = qijlm
      do 200 kid = 1, nc
      ksy = iphz(kid)
c     QUANTITATIVE DATA                                                 
c         
      if (iph(kid) - 2) 302, 200, 201
  201 t1 = qt(kid) - sygc
      chm1 = cxc1 * exp((t1 * t1) * cxc2)
cCCC  IF(KSY.EQ.1.AND.ISAVE(1).EQ.1)GO TO 252                           
c         
      if (isel) 250, 250, 301
  301 if (ksy - isave(1)) 305, 252, 305
  305 if (isave(ksy)) 306, 202, 306
  306 chq1 = pi * save1(ksy)
c     AFFECTION STATUS                                                  
c         
      goto 251
  302 if (isave(ksy)) 307, 202, 307
  307 if (iph(kid)) 308, 303, 308
  308 chm1 = save1(ksy)
      if (isel) 250, 250, 309
  309 chq1 = chm1 * pi
      goto 251
  303 chm1 = 1. - save1(ksy)
      if (isel) 250, 250, 311
  311 chq1 = save1(ksy) * pi
      goto 251
  202 xx = zcdc(ksy) - gccdc
      ax = xx
      if (xx .lt. 0.) ax = - xx
      t1 = 1.0 / (1.0 + (.2316419 * ax))
      t2 = 0.3989422804 * exp(- ((.5 * xx) * xx))
      t3 = (t2 * t1) * ((((((((1.330274 * t1) - 1.821256) * t1) + 
     &1.781478) * t1) - 0.3565638) * t1) + 0.3193815)
      if (xx) 212, 213, 213
  212 t3 = 1.0 - t3
  213 t4 = t3 * pi
      t33 = t3
      if (iph(kid) - 1) 214, 215, 263
  214 t3 = 1.0 - t3
  215 chm1 = t3
  263 chq1 = t4
      save1(ksy) = t33
      if (isel) 266, 266, 267
  266 isave(ksy) = 1
c     ONE LIABILITY CLASS SPEED ENHANCER                                
c         
      goto 250
  267 if (ksy - 1) 255, 310, 255
  310 isave(1) = 1
      tsy1 = 1. - chq1
  252 qijlm = qijlm * tsy1
      goto 250
  255 isave(ksy) = 1
  251 qijlm = qijlm * (1. - chq1)
  250 pijlm = pijlm * chm1
c END OF PILJM                                                          
c         
  200 continue
      if (isel) 274, 274, 273
c END OF QIJLM                                                          
c         
  273 qijlm = poic1 - qijlm
  274 syiljm = pil1(l) * pjm1(m)
      if (itr * (m - l)) 272, 272, 271
  271 syiljm = syiljm + (pil1(m) * pjm1(l))
  272 sij1 = pijlm * syiljm
      if (iccc) 2721, 2722, 2722
 2721 sij2 = qijlm
      goto 2723
 2722 sij2 = qijlm * syiljm
 2723 continue
      slm1 = slm1 + (sij1 * plm)
c     WRITE(61,810)SIJ1,PLM,PIJLM,SYILJM                                
c         
c 810 FORMAT(1X,'HELP',3E18.11)                                         
c         
      slm2 = slm2 + (sij2 * plm)
   90 continue
  100 continue
      slm1 = ((slm1 * fpi) * rscal) * qzxw1
      slm2 = (slm2 * rscal) * qzxw2
      if (isel) 812, 813, 813
  812 pfam = slm1
      goto 814
  813 pfam = slm1 / slm2
  814 continue
c     CALL ETIME                                                        
c         
c PRINT OUTPUT                                                          
c         
c     WRITE(61,7001)ALIK                                                
c         
c7001 FORMAT(1X,'ALOG',E18.11)                                          
c         
      alik = - log(pfam)
      return 
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         
c                                                                       
c         
      subroutine lratio(xall, nall, itp, iq, ictrl, nn)
      implicit double precision (o-z, a-h)
c               OPTIONAL PRINTOUT IF DEBUG=1                            
c         
c                                                                       
c         
      parameter (nulls = -8388607)
      dimension b(3751), aout(100), jname(13), title(7)
      character jname*3
      character oldped*6, cped*6, cfam*6
      character title*12
      logical start, newped
      integer idebug
      equivalence (ib, b), (iaout, aout)
      common /auxp1/ p(3), g(3), clm2, clm1, cxp2, cxp1(3), rcdg, 
     &cdp, z(10), vv, uu, gg, ca, ck, ea, ek, scals(4), sclq, rsqw
     &, notrs, ng1, nox
      common /auxp2/ ntrel, mxrel, clm2s(3), clm1s(3), rrrs(2, 3, 5), 
     &cxc2s(2, 2), cxc1s(2, 2), cdcs(2, 2), cxd2s(2, 2, 5), cxd1s(2, 2, 
     &5), cdds(2, 2, 5), cxa2s(2, 5), cxa1s(2, 5), cdas(2, 5)
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, nev, idebug, np
      common /save/ prep(15), trat
      common /odds/ pis(10), fi(10), ititle(27), itypgr, irap, npi, nfi
      common /intg/ xuu(20), auu(20), x(20), a(20), fzdqz(10), xdev, 
     &tvuu, fcaea, igot, nuu, nl, nm
      common /data/ tmf(9), tmm(9), qt(18), tmk(18),
     &zfa, zmo, qtfp, qtf, qtm, pi,
     &qtmp, qtcp, zphf, zphm, zphc, iph(18), iphz(18), irpf, irpm, irpk,
     &nc, isel, iphf, iphm, iphzf, iphzm, iphfp, iphmp, iphcp
      common /data2/ cxc2, cxc1, cdc, rrc, cxpf2, 
     &cxpf1, cdpf, cxpm2, cxpm1, cdpm, rrpf, rrpm, rrpc, 
     &cxd2, cxd1, cdd, alik, iquad, iccc, irrpf, irrpm, irrpc
      common /trsm/ tr(18), tm(6, 3, 70)
      common ib(7502), iaout(200)
      dimension xall(15), itp(15)
c ----------------------------------------------------------------------
c-------  
c FOR EACH OF TWO HYPOTHESES                                            
c         
      data jname / 'V  ', 'U  ', 'D  ', 'T  ', 'Q  ', 'H  ', 'Z  ', 
     &'X  ', 'W  ', 'T1 ', 'T2 ', 'T3 ', 'R  ' /
      data rped /0./
      write(unit=iq, fmt=1003) 
 1003 format(1h1///30x,28hASSESSING LIKELIHOOD OF TWO ,
     &20hCOMPETING HYPOTHESES//30x,48(1h*))
      schi2 = 0.
      sl1 = 0.
      sl2 = 0.
      nfamj = 0
      do 60 ih = 1, 2
      rewind(unit=52) 
      read(unit=52, fmt=1) inum, nfam1
      read(unit=52, fmt=2) title
    1 format(2i6)
    2 format(a12)
    5 format(2a6)
      if (idebug .eq. 1) write(unit=61, fmt=901) inum, nfam1, title
  901 format(6h INUM=,i3,7h NFAM1=,i3/1h ,7a12)
      start = .true.
      if (ih .eq. 1) goto 20
      do 74 i = 1, nall
c -----------------------------------                                   
c         
c ALLOWING FOR NO TRANSMISSION                                          
c         
c IF T1,T2,T3 ALL EQUAL , NO TRANSMISSION MODEL                         
c         
   74 xall(i) = rxall(i)
c NOTRS = 0 : MENDELIAN                                                 
c         
c NOTRS = 1 : T1=T2=T3=1.-QQ                                            
c         
c NOTRS = 2 : T1=T2=T3=T                                                
c         
      notrs = 0
      testt = abs((1. - xall(5)) - xall(10))
      if ((xall(10) .eq. xall(11)) .and. (xall(11) .eq. xall(12))) notrs
     & = 2
      if ((notrs .eq. 2) .and. (testt .lt. 1.e-05)) notrs = 1
c *** CA ***                                                            
c         
      xall(6) = xall(6) * xall(1)
      zinp = xall(7)
      xall(7) = xall(6) * xall(7)
      nl = 1
      nm = 1
      a(1) = 1.0
      x(1) = 0.0
      auu(1) = 1.0
      xuu(1) = 0.
      nuu = 1
      if (xall(6) .eq. 0.) goto 20
c ----------------------------------------------------------------      
c         
c GET WEIGHTS AND ABSCISSAS FOR QUADRATURE                              
c         
      eps = 1.e-10
      call hermit(nn, x, a, eps)
      do 4 i = 1, nn
    4 a(i) = a(i) * exp(x(i) * x(i))
      nl = nn
      nm = nn
      nuu = min(2 * nn,20)
      call hermit(nuu, xuu, auu, eps)
      do 44 i = 1, nuu
c *** CA ***                                                            
c         
   44 auu(i) = auu(i) * exp(xuu(i) * xuu(i))
      xdev = (.5 * xall(7)) / xall(1)
   20 hloc = xall(6) / xall(1)
      if (xall(6) .ne. 0.) zloc = xall(7) / xall(6)
      write(unit=iq, fmt=1006) ih, jname, (xall(i),i = 1, 5), hloc, zloc
     &, (xall(j),j = 8, 13)
 1006 format(///1x,10hHYPOTHESIS,i3//1x,13(3x,a3,3x)//1x,13(f8.5,1x)/)
      if (ih .eq. 1) goto 137
      if ((irap .eq. 0) .and. (itypgr .ne. 0)) then
      write(unit=iq, fmt=1007) 
 1007 format(//47x,16hLIKELIHOOD-RATIO//1x,10hPEDIGREE  ,6hFAMILY,7x,
     &8h-2LN(PA),4x,7x,8h-2LN(RA),4x,5x,12h2(LN(PA/RA)),4x,7x,8h-2LN(PA)
     &,4x,7x,8h-2LN(RA),4x,5x,12h2(LN(PA/RA))/)
      else
      if ((irap .eq. 1) .and. (itypgr .ne. 0)) then
      write(unit=iq, fmt=1008) 
 1008 format(//47x,16hLIKELIHOOD-RATIO//1x,10hPEDIGREE  ,7x,8h-2LN(PA)
     &,4x,7x,8h-2LN(RA),4x,3x,12h2(LN(PA/RA))/)
      else
      write(unit=iq, fmt=1009) 
 1009 format(//47x,16hLIKELIHOOD-RATIO//11x,7hFAMILY ,7x,8h-2LN(PA)
     &,4x,9x,8h-2LN(RA),4x,3x,12h2(LN(PA/RA))/)
      end if
      end if
  137 nout = 0
c GET AUXILIARY PARAMETERS                                              
c         
      nfami = 0
      igot = 0
      nox = 0
c ------------------------------------------------------------------    
c         
c READ BLOCK OF DATA IF REQUIRED                                        
c         
      call auxpar(xall, nall, itp)
      if (nread - 1) 10, 12, 10
   10 read(unit=48, end=13) (ib(i),i = 1, 7502)
      if (ib(1) .eq. nulls) goto 13
   12 nfam = ib(7501)
c -------------------------------------------------------------------   
c         
c FOR EACH FAMILY, GET -LN(LIKELIHOOD)                                  
c         
      iword2 = 0
      do 310 i = 1, nfam
      iword = iword2 + iword2
      k2 = 0
      qtf = 0.
      qtm = 0.
      iphm = 2
      nc = ib(iword + 3)
c ********************* TEMP *************************                  
c         
      iquad = ib(iword + 8)
      if (iquad .eq. 2) iquad = 1
      itype = 2
      if (iquad .eq. 2) itype = 1
      clm2 = clm2s(iquad)
      clm1 = clm1s(iquad)
      isel = ib(iword + 14)
      if ((itypgr .eq. 1) .or. (iiaf .eq. 1)) read(unit=52, fmt=5) cfam
     &, cped
      if (itypgr .ne. 1) rped = b(iword2 + 9)
      pi = 0.
c     IQUAD=2, ONE PARENT MISSING AND NO POINTER TO THIS PARENT.        
c         
      if (isel .gt. 0) pi = pis(isel)
      if (iquad - 2) 450, 459, 480
  459 if (ib(iword + 4) .eq. 0) goto 18
      iphf = 2
      if (ib(iword + 4) - 2) 463, 463, 460
  460 if (ib(iword + 21) .ne. 2) goto 461
      if (ib(iword + 27) .eq. 2) goto 462
      qtm = b(iword2 + 13)
      iphm = ib(iword + 27)
      iphzm = ib(iword + 28)
      zmo = z(iphzm)
  462 iphf = 2
      k2 = k2 + 16
      goto 25
  461 k2 = 3
      goto 22
  463 if (ib(iword + 21) .eq. 2) goto 470
      goto 22
  480 if (ib(iword + 4) .eq. 0) goto 18
      iphf = 2
      if (ib(iword + 4) - 2) 470, 470, 471
  470 k2 = 13
      goto 25
  471 k2 = 16
c --------------------                                                  
c         
c     FATHER,MOTHER DATA                                                
c         
      goto 25
  450 if (ib(iword + 4) - 2) 19, 22, 23
   19 if (ib(iword + 4)) 21, 18, 21
   18 iphf = 2
      k2 = 10
      goto 25
   23 qtm = b(iword2 + 13)
      iphm = ib(iword + 27)
      iphzm = ib(iword + 28)
      zmo = z(iphzm)
      k2 = 3
   21 qtf = b(iword2 + 10)
      iphf = ib(iword + 21)
      iphzf = ib(iword + 22)
      zfa = z(iphzf)
      k2 = k2 + 13
      goto 25
   22 qtm = b(iword2 + 10)
      iphm = ib(iword + 21)
      iphzm = ib(iword + 22)
      zmo = z(iphzm)
      iphf = 2
c --------------------                                                  
c         
c     CHILDREN DATA                                                     
c         
      k2 = k2 + 13
   25 k = ((k2 + k2) + 1) + iword
      k2 = k2 + iword2
      igec = ib(iword + 9)
      cxc2 = cxc2s(igec,itype)
      cxc1 = cxc1s(igec,itype)
      cdc = cdcs(igec,itype)
      rrc = rrrs(igec,iquad,1)
      iccc = ib(iword + 13)
      do 27 j = 1, nc
      qt(j) = b(k2)
      iph(j) = ib(k)
      iphz(j) = ib(k + 1)
      k2 = k2 + 3
   27 k = k + 6
      k = k - 3
c --------------------                                                  
c         
c     FATHER,MOTHER,CHILDREN'S POINTER                                  
c         
      k2 = k2 - 1
      if (iquad .ne. 2) goto 451
      iphfp = 2
      if (ib(iword + 5)) 29, 29, 428
  428 irpm = ib(iword + 5)
      irrpm = 1 + (irpm / (ntrel + 1))
      igepm = ib(iword + 10)
      rrpm = rrrs(igepm,iquad,irrpm)
      goto 452
  451 if (ib(iword + 5)) 40, 40, 28
   40 iphfp = 2
      goto 29
   28 irpf = ib(iword + 5)
      irrpf = 1 + (irpf / (ntrel + 1))
      igepf = ib(iword + 10)
      rrpf = rrrs(igepf,iquad,irrpf)
      cxpf2 = cxa2s(igepf,irrpf)
      cxpf1 = cxa1s(igepf,irrpf)
      cdpf = cdas(igepf,irrpf)
      qtfp = b(k2 + 1)
      iphfp = ib(k + 3)
      zphf = z(ib(k + 4))
      k = k + 6
      k2 = k2 + 3
   29 if (ib(iword + 6)) 41, 41, 30
   41 iphmp = 2
      goto 31
   30 irpm = ib(iword + 6)
      irrpm = 1 + (irpm / (ntrel + 1))
      igepm = ib(iword + 11)
      rrpm = rrrs(igepm,iquad,irrpm)
  452 cxpm2 = cxa2s(igepm,irrpm)
      cxpm1 = cxa1s(igepm,irrpm)
      cdpm = cdas(igepm,irrpm)
      qtmp = b(k2 + 1)
      iphmp = ib(k + 3)
      zphm = z(ib(k + 4))
      k = k + 6
      k2 = k2 + 3
   31 if (ib(iword + 7)) 42, 42, 32
   42 iphcp = 2
      goto 35
c SINGLE SELECTION OF MULTIPLEX FAMILIES IF IRPK=MXREL+1 (CODE=1A)      
c         
c IF SO,SET ISEL=NPI AND PI=1.0                                         
c         
   32 irpk = ib(iword + 7)
      if ((mxrel + 1) - irpk) 34, 33, 34
   33 isel = npi
      pi = 1.0
   34 continue
      irrpc = irpk - mxrel
      irrpc = 1 + (irrpc / (ntrel + 1))
      igepc = ib(iword + 12)
      rrpc = rrrs(igepc,iquad,irrpc)
      cxd2 = cxd2s(igepc,itype,irrpc)
      cxd1 = cxd1s(igepc,itype,irrpc)
      cdd = cdds(igepc,itype,irrpc)
      qtcp = b(k2 + 1)
      iphcp = ib(k + 3)
c --------------------                                                  
c         
c     GET READY TO CALL EXEC                                            
c         
c SKIP THIS IF Q=0.0 AND CALL LIKEQ0                                    
c         
      zphc = z(ib(k + 4))
   35 if (xall(5)) 36, 36, 350
   36 call likeq0
      goto 830
  350 if (iphmp - 2) 822, 823, 822
  822 tmm(1) = tm(1,1,irpm)
      tmm(2) = tm(1,2,irpm)
      tmm(3) = tm(1,3,irpm)
      tmm(4) = tm(2,1,irpm)
      tmm(5) = tm(2,2,irpm)
      tmm(6) = tm(2,3,irpm)
      tmm(7) = tm(3,1,irpm)
      tmm(8) = tm(3,2,irpm)
      tmm(9) = tm(3,3,irpm)
  823 if (iphfp - 2) 824, 825, 824
  824 tmf(1) = tm(1,1,irpf)
      tmf(2) = tm(1,2,irpf)
      tmf(3) = tm(1,3,irpf)
      tmf(4) = tm(2,1,irpf)
      tmf(5) = tm(2,2,irpf)
      tmf(6) = tm(2,3,irpf)
      tmf(7) = tm(3,1,irpf)
      tmf(8) = tm(3,2,irpf)
      tmf(9) = tm(3,3,irpf)
  825 if (iphcp - 2) 826, 827, 826
  826 tmk(1) = tm(1,1,irpk)
      tmk(2) = tm(1,2,irpk)
      tmk(3) = tm(1,3,irpk)
      tmk(4) = tm(2,1,irpk)
      tmk(5) = tm(2,2,irpk)
      tmk(6) = tm(2,3,irpk)
      tmk(7) = tm(3,1,irpk)
      tmk(8) = tm(3,2,irpk)
      tmk(9) = tm(3,3,irpk)
      tmk(10) = tm(4,1,irpk)
      tmk(11) = tm(4,2,irpk)
      tmk(12) = tm(4,3,irpk)
      tmk(13) = tm(5,1,irpk)
      tmk(14) = tm(5,2,irpk)
      tmk(15) = tm(5,3,irpk)
      tmk(16) = tm(6,1,irpk)
      tmk(17) = tm(6,2,irpk)
c --------------------------------------------------------------        
c         
c GET LIKELIHOOD FOR THAT FAMILY                                        
c         
      tmk(18) = tm(6,3,irpk)
  827 continue
      call like
  830 alik = 2. * alik
      if (ih .eq. 2) goto 8030
      nfamj = nfamj + 1
      sl1 = sl1 + alik
      nout = nout + 1
      aout(nout) = alik
      if (nout .ne. 100) goto 8099
      write(unit=47) (iaout(j),j = 1, 200)
 8002 nout = 0
      goto 8099
 8030 nfami = nfami + 1
      sl2 = sl2 + alik
      nout = nout + 1
      if (itypgr .eq. 1) then
      newped = cped .ne. oldped
      else
      newped = rped .ne. roldpd
      end if
      if (.not. newped) then
      xih1 = xih1 + aout(nout)
      xih2 = xih2 + alik
      else
      if (start) then
      start = .false.
      else
      chi2 = xih2 - xih1
      if ((irap .eq. 0) .and. (itypgr .eq. 1)) write(unit=62, fmt=8041) 
     &oldped, xih1, xih2, chi2
      if ((irap .eq. 0) .and. (itypgr .eq. 2)) write(unit=62, fmt=8042) 
     &roldpd, xih1, xih2, chi2
      if ((irap .eq. 1) .and. (itypgr .eq. 1)) write(unit=62, fmt=8043) 
     &oldped, xih1, xih2, chi2
      if ((irap .eq. 1) .and. (itypgr .eq. 2)) write(unit=62, fmt=8044) 
     &roldpd, xih1, xih2, chi2
 8041 format(1x,a6,4x,6x,3(19x),2x,3(4x,e15.8))
 8042 format(1x,f10.0,6x,3(19x),2x,3(4x,e15.8))
 8043 format(1x,a6,4x,3(4x,e15.8))
 8044 format(1x,f10.0,3(4x,e15.8))
      end if
      if (itypgr .eq. 1) then
      oldped = cped
      else
      roldpd = rped
      end if
      xih1 = aout(nout)
      xih2 = alik
      end if
      chi2 = alik - aout(nout)
      schi2 = schi2 + chi2
      if (trat .eq. 1000.) goto 608
      if (abs(chi2) .lt. trat) goto 63
      goto 618
  608 if ((irap .eq. 1) .and. (itypgr .ne. 0)) goto 63
  618 if (iiaf - 1) 61, 61, 609
  609 write(unit=iq, fmt=62) b(iword2 + 8), aout(nout), alik, chi2
   62 format(1x,6x,f10.0,3(4x,e15.8))
      goto 63
   61 write(unit=iq, fmt=1002) cfam, aout(nout), alik, chi2
 1002 format(1x,10x,a6,3(4x,e15.8))
   63 if (nout .ne. 100) goto 8099
      if (nfami .eq. nfamj) goto 13
      read(unit=47) (iaout(j),j = 1, 200)
 8032 nout = 0
 8099 iword2 = (iword2 + 9) + (ib(iword + 2) * 3)
      if (idebug .ne. 1) goto 310
      write(unit=61, fmt=6161) i, iquad, isel, iccc, igec, nc, iphf, 
     &iphm, (iph(kkk),kkk = 1, nc), (iphz(kkz),kkz = 1, nc)
 6161 format(1x,30i3)
      write(unit=61, fmt=6162) qtf, qtm, (qt(kkk),kkk = 1, nc)
 6162 format(5x,10e10.3)
c --------------------------------------------------------------        
c         
c END OF LOOP OVER FAMILIES                                             
c         
  310 continue
      if (nread .eq. 1) goto 605
      goto 10
   13 rewind(unit=48) 
  605 if (ih .eq. 2) goto 60
      if (nout .eq. 0) goto 600
      write(unit=47) (iaout(i),i = 1, 200)
  600 rewind(unit=47) 
      read(unit=47) (iaout(i),i = 1, 200)
   60 continue
      rewind(unit=47) 
      chi2 = xih2 - xih1
      if ((irap .eq. 0) .and. (itypgr .eq. 1)) write(unit=62, fmt=8041) 
     &oldped, xih1, xih2, chi2
      if ((irap .eq. 0) .and. (itypgr .eq. 2)) write(unit=62, fmt=8042) 
     &roldpd, xih1, xih2, chi2
      if ((irap .eq. 1) .and. (itypgr .eq. 1)) write(unit=62, fmt=8043) 
     &oldped, xih1, xih2, chi2
      if ((irap .eq. 1) .and. (itypgr .eq. 2)) write(unit=62, fmt=8044) 
     &roldpd, xih1, xih2, chi2
      do 66 i = 1, nall
   66 xall(i) = prep(i)
      if ((irap .eq. 0) .and. (itypgr .ne. 0)) then
      write(unit=62, fmt=1005) sl1, sl2, schi2
      else
      write(unit=62, fmt=1004) sl1, sl2, schi2
      end if
 1004 format(/1x,3hALL,7x,6x,3(4x,e15.8))
 1005 format(/1x,3hALL,7x,6x,3(19x),2x,3(4x,e15.8))
 8098 return 
      end
c                                                                       
c         
c **********************************************************************
c******** 
c                                                                       
c         
      subroutine nspace()
      implicit double precision (o-z, a-h)
      character a*1
      common /sy1a/ ibegin, i1, i2, icol, index, nchar
      common /sy1b/ a(80)
c                                                                       
c         
c     FINDS NEXT NON-BLANK COLUMN                                       
c         
c                                                                       
c         
      icol = icol + 1
      entry nsp1()
    5 index = 1
    6 if (a(icol) .eq. ' ') then
      icol = icol + 1
      if (icol .gt. nchar) goto 55
      else
      goto 7
      end if
      goto 6
    7 return 
c                                                                       
c         
c     CHECKS FOR RIGHT AND LEFT PARENTHESIS STARTING WITH ICOL          
c         
c     I1 ON RETURN = COL. AFTER (   ;   I2 = COL. BEFORE )              
c         
c     ICOL ON RETURN = 1ST NON-BLANK COLUMN AFTER (                     
c         
c                                                                       
c         
      entry par()
    8 if (a(icol) .ne. '(') then
      icol = icol + 1
      if (icol .gt. nchar) goto 55
      else
      goto 9
      end if
      goto 8
    9 icol = icol + 1
      i1 = icol
   10 if (a(icol) .ne. ')') then
      icol = icol + 1
      if ((icol .gt. nchar) .or. (a(icol) .eq. '(')) goto 55
      else
      goto 11
      end if
      goto 10
   11 i2 = icol - 1
      icol = i1
      goto 5
   55 index = 0
      end
c                                                                       
c         
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         
c                                                                       
c         
      subroutine output(nall, xall, nit, nfe, f, ptg, idg, iq, ivar, se
     &, ictrl, h, tol, trupb, n, itp)
c ----------------------------------------------------------------------
c-        
c                                                                       
c         
c WRITE FINAL OUTPUT ON LOGICAL UNIT IQ                                 
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c-        
      implicit double precision (o-z, a-h)
      character conjnt*50
      integer idebug
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, nev, idebug, np
      common /auxp1/ pee(3), gee(3), clm2, clm1, cxp2, cxp1(3), 
     &rcdg, cdp, z(10), vv, uu, gg, ca, ck, ea, ek, scals(4), sclq
     &, rsqw, notrs, ng1, nox
      common /mini/ tl(15, 15), d(15), g(15), gsave(15), y(15), p(15), 
     &tmin, fsmf, idif, isw, iret
      common /secc/ seca, seck, b(15, 15), saveb(15, 15), ifail
      common /save/ prep(15), trat
      common /odds/ pis(10), fi(10), ititle(27), itypgr, irap, npi, nfi
      common /data/ tmf(9), tmm(9), qt(18), tmk(18),
     &zfa, zmo, qtfp, qtf, qtm, pi,
     &qtmp, qtcp, zphf, zphm, zphc, iph(18), iphz(18), irpf, irpm, irpk,
     &nc, isel, iphf, iphm, iphzf, iphzm, iphfp, iphmp, iphcp
      dimension se(15), itp(15)
      dimension xall(15), ggall(15)
      igg = 0
      do 11 i = 1, nall
      ggall(i) = 0.
      if (itp(i) .eq. 0) goto 11
      igg = igg + 1
      ggall(i) = g(igg)
c DESCALE T AND Q IF THEY WERE SCALED                                   
c         
   11 continue
      xall(4) = xall(4) * sqrt(sclq)
      se(4) = se(4) * sqrt(sclq)
      xall(5) = xall(5) / sclq
c                                                                       
c         
c ALLOWING FOR NO TRANSMISSION                                          
c         
      se(5) = se(5) / sclq
      if (notrs .eq. 0) goto 7070
      if (notrs .eq. 1) xall(10) = 1. - xall(5)
      xall(11) = xall(10)
      xall(12) = xall(10)
 7070 continue
      do 71 i = 1, nall
   71 prep(i) = xall(i)
      prep(6) = xall(6) / xall(1)
      if (xall(6)) 72, 73, 72
   72 prep(7) = xall(7) / xall(6)
      goto 74
   73 prep(7) = 0.
c NOTE THAT FINAL OUTPUT IS SAVED IN DIFFERENT OUTPUT AREA              
c         
   74 continue
      write(unit=iq, fmt=1006) (ititle(i),i = 1, 27)
 1006 format(/1x,27a3)
      write(unit=iq, fmt=400) nit, nfe
  400 format(//1x,14hITERATION NO. ,i4,10x,6hNFE = ,i3/1x,18(1h*)//)
      f2 = 2. * f
      if (isel .ge. 0) then
      conjnt = '(CHILDREN CONDITIONAL ON PARENTS)'
      else
      conjnt = '(JOINT PARENTS AND CHILDREN)'
      end if
      write(unit=iq, fmt=401) f2, conjnt, fsmf, ptg
  401 format(1x,11h-2 LN(L) = ,e18.11,5x,a50//1x,6h# F = ,e18.11//1x,
     &15hGOODNESS-OF-FIT,5x,8hUK-1U = ,f10.5/)
      write(unit=iq, fmt=402) (prep(i),i = 1, nall)
  402 format(1x,11hPARAMETERS:,4x,2hV ,6x,4x,2hU ,6x,4x,2hD ,6x,4x,2hT 
     &,6x,4x,2hQ ,6x,4x,2hH ,6x,4x,2hZ ,6x,4x,2hX ,6x,4x,2hW ,/,1x,
     &11hESTIMATES: ,9f12.7,//3x,11x,4x,2hT1,6x,4x,2hT2,6x,4x,2hT3
     &,6x,4x,2hR ,/1x,11x,4f12.7/)
      write(unit=iq, fmt=404) (ggall(i),i = 1, nall)
  404 format(1x,11hU- SCORES: ,9f10.5,/1x,11x,4f10.5)
      write(unit=iq, fmt=405) (se(i),i = 1, nall)
  405 format(1x,11hS. ERRORS: ,9f10.5,/1x,11x,4f10.5)
      write(unit=iq, fmt=814) 
  814 format(1x,130(1h*))
cc JH Zhao change from 1x130(1h*)
      write(unit=iq, fmt=406) 
  406 format(1x,10hB - MATRIX)
      do 711 i = 1, n
  711 write(unit=iq, fmt=407) (saveb(i,j),j = 1, n)
  407 format(1x,13f10.2)
      write(unit=iq, fmt=408) 
  408 format(1x)
      write(unit=iq, fmt=409) 
  409 format(1x,11hB - INVERSE)
      do 712 i = 1, n
  712 write(unit=iq, fmt=410) (b(i,j),j = 1, n)
  410 format(1x,13f10.5)
      write(unit=iq, fmt=408) 
      if (ifail .eq. 0) then
      write(unit=iq, fmt=413) 
  413 format(1x,36hFACTORIZATION OF B-MATRIX SUCCEEDED.)
      else
      write(unit=iq, fmt=412) 
  412 format(1x,33hFACTORIZATION OF B-MATRIX FAILED.)
      end if
      write(unit=iq, fmt=814) 
      write(unit=iq, fmt=1004) gg, ck, seck, ca, seca
 1004 format(1x,21hAUXILIARY PARAMETERS:,4x,3hGG=,e11.4,5h  CK=,e11.4,
     &9h  SE(CK)=,e11.4,5h  CA=,e11.4,9h  SE(CA)=,e11.4)
c          ADDITIONAL OUTPUT GOTTEN IN AUXOUT                           
c         
      write(unit=iq, fmt=814) 
      if (nia .eq. 0) goto 1999
      call auxout(iq)
c --------------------------------------------------------------------  
c         
c *** KEEP *** USER OUGHT TO INSERT THIS SECTION IN HIS OUTPUT OUTINE   
c         
c INDICATES REASON FOR TERMINATION                                      
c         
      write(unit=iq, fmt=814) 
 1999 write(unit=iq, fmt=2000) 
 2000 format(1x,37hITERATIVE PROCESS TERMINATED BECAUSE:)
      idg = idg + 1
      goto (1, 2, 3, 4, 5, 6, 7, 8), idg
    1 write(unit=iq, fmt=2001) 
 2001 format(25x,41h*** MAXIMUM POSSIBLE ACCURACY REACHED ***)
      goto 9
    2 write(unit=iq, fmt=2002) 
 2002 format(25x,45h*** SEARCH DIRECTION NOT DOWNHILL ANYMORE ***)
      goto 9
    3 write(unit=iq, fmt=2003) 
 2003 format(25x,45h*** ACCUMULATION OF ROUNDING ERRORS PREVENTS ,
     &20hFURTHER PROGRESS ***)
      goto 9
    4 write(unit=iq, fmt=2004) 
 2004 format(25x,40h*** ALL SIGNIFICANT DIGITS LOST THROUGH ,
     &40hCANCELLATION IN OPTIMAL CONDITIONING ***)
      goto 9
    5 write(unit=iq, fmt=2005) 
 2005 format(25x,46h*** SPECIFIED TOLERANCE ON NORMALIZED GRADIENT,
     &12h WAS MET ***)
      goto 9
    6 write(unit=iq, fmt=2006) 
 2006 format(25x,47h*** SPECIFIED TOLERANCE ON GRADIENT WAS MET ***)
      goto 9
    7 write(unit=iq, fmt=2007) 
 2007 format(25x,44h*** MAXIMUM NUMBER OF ITERATIONS REACHED ***)
      goto 9
    8 write(unit=iq, fmt=2008) 
 2008 format(25x,39h*** EXCESSIVE CANCELLATION IN GRADIENT ,
     &15hCALCULATION ***)
    9 write(unit=iq, fmt=814) 
      write(unit=iq, fmt=420) h, trupb, tol
c *** KEEP *** END OF SECTION TO SAVE                                   
c         
  420 format(42h NUMERICAL DIFFERENTIATION INCREMENT(H) = ,e14.7,1h;,4x,
     &32hTRUNCATION UPPER BOUND(TRUPB) = ,e14.7,1h;,4x,6hTOL = ,e14.7)
      return 
      end
c **********************************************************************
c**       
c                                                                       
c         
c                          P O I N T E R                                
c         
c                          -------------                                
c         
c COMPLEX SEGREGATION ANALYSIS OF NUCLEAR FAMILIES WITH/WITHOUT POINTERS
c         
c                     --- AUTOSOMAL CASE ---                            
c         
c 532                                                                   
c            
c **********************************************************************
c***      
c                                                                       
c         
c     program pointr
      subroutine pointr(jobfile,profile,terfile)
      implicit double precision (o-z, a-h)
      character ctime*24
      character cfam*6, cgr*6
      character title*12
      character ia*1, jobfile*30,profile*30,terfile*30
      integer idebug
      dimension cfam(250), cgr(250), title(7)
      dimension xall(15), itp(15)
      common /para/ nread, nia, yi(10), zi(10), rxall(15), zinp, jli
     &, iiaf, nev, idebug, np
      common /save/ prep(15), trat
      common /sy1a/ ibegin, i1, i2, icol, index, nchar
      common /sy1b/ ia(80)
      common /intg/ xuu(20), auu(20), x(20), a(20), fzdqz(10), xdev, 
     &tvuu, fcaea, igot, nuu, nl, nm
      common /auxp1/ p(3), g(3), clm2, clm1, cxp2, cxp1(3), rcdg, 
     &cdp, z(10), vv, uu, gg, ca, ck, ea, ek, scals(4), sclq, rsqw
     &, notrs, ng1, nox
      common ib(7502)
c     call getarg(1,jobfile)
c     call getarg(2,profile)
c     call getarg(3,terfile)
      open(unit=52, file='INTERMED.ALP', status='OLD') 
      open(unit=47, form='unformatted',status='SCRATCH') 
      open(unit=60, file=jobfile, status='OLD') 
      open(unit=61, file=profile, status='OLD',access='append')
cc JH Zhao 6/6/2001 from fileopt='EOF' 
      open(unit=62, file=terfile, status='OLD',access='append')
cc JH Zhao 6/6/2001 from fileopt='EOF' 
      open(unit=48,file='INTERMEDFILE',status='OLD',
     & form='unformatted') 

      call fdate(ctime)
c                                                                       
c         
      write(unit=62, fmt=1361) ctime
      write(unit=61, fmt=1362) ctime
 1361 format(//1x,a24,9x,41h RPOINTA CREATED 08-30-85.  TERSE OUTPUT.)
 1362 format(//1x,a24,9x,42h RPOINTA CREATED 08-30-85.  PROLIX OUTPUT.)
      ibegin = 0
      ictrl = 100
      ic = 60
      ip = 61
      iq = 62
      id = 49
      id1 = 48
c ------------------------------------------------------------------    
c         
c INTERPRET CONTROL 'CARD'                                              
c         
      nall = 13
c -----------------------------------                                   
c         
c ALLOWING FOR NO TRANSMISSION                                          
c         
c IF T1,T2,T3 ALL EQUAL , NO TRANSMISSION MODEL                         
c         
10      call entry(ic, ictrl, nall, xall, itp, h, nn, tol, trupb)
	call flush(61)
	call flush(62)
c NOTRS = 0 : MENDELIAN                                                 
c         
c NOTRS = 1 : T1=T2=T3=1.-QQ                                            
c         
c NOTRS = 2 : T1=T2=T3=T                                                
c         
      notrs = 0
      testt = abs((1. - xall(5)) - xall(10))
      if ((xall(10) .eq. xall(11)) .and. (xall(11) .eq. xall(12))) notrs
     & = 2
      if ((notrs .eq. 2) .and. (testt .lt. 1.e-05)) notrs = 1
      if (notrs .eq. 2) itp(10) = 1
      if (notrs .ne. 0) write(unit=ip, fmt=1000) 
c -----------------------------------                                   
c         
c THIS NEXT VARIABLE CONTROLS CALCULATION OF TVUU AND FCAEA             
c         
c IN AUXPAR, CALCULATION BEING DONE ONLY AT VERY FIRST CALL             
c         
c ALSO, CONTROLS CALCULATION OF TRANSITION MATRICES IN AUXPAR,          
c         
c REQUIRED AT VERY FIRST CALL TO AUXPAR ONLY IN CASE Q AND X            
c         
c ARE NOT ITERATED                                                      
c         
 1000 format(/10x,31hNO TRANSMISSION OF MAJOR EFFECT)
c     WRITE(6,666)(XALL(I),I=1,NALL),(ITP(J),J=1,NALL)                  
c         
c 666 FORMAT(1X,'XALL',13F8.5,/1X,'ITP ',13I3)                          
c         
c                                                                       
c         
      igot = 0
c                                                                       
c         
c -------------------------------------------------------------------   
c         
c Q WILL BE RESCALED IF FOUND SMALL RELATIVE TO DIFFERENTIATION INTERVAL
c         
c THIS SCALING NOT POSSIBLE FOR LIKELIHOOD-RATIOS                       
c         
      if (ictrl .eq. 100) goto 500
      sclq = 1.
      goto (11, 12), ictrl
   11 continue
      if ((xall(5) .gt. 0.0) .and. (xall(5) .lt. (4. * h))) sclq = 1. / 
     &xall(5)
      xall(5) = xall(5) * sclq
      xall(4) = xall(4) / sqrt(sclq)
c     CHECK IF H=0.0                                                    
c         
c IN FACT H REPLACED BY CK HERE                                         
c         
   12 continue
c *** CA ***                                                            
c         
      xall(6) = prep(6) * xall(1)
      zinp = prep(7)
c NN DEFINED IN ENTRY ( DEFAULT VALUE NN=5 )                            
c         
      xall(7) = xall(6) * prep(7)
      write(unit=ip, fmt=3000) nn
 3000 format(10x,15h********** NN =,i3,11h **********/)
      nl = 1
      nm = 1
      a(1) = 1.0
      x(1) = 0.0
      auu(1) = 1.0
      xuu(1) = 0.0
      nuu = 1
      if (xall(6) .eq. 0.) goto 5
c ----------------------------------------------------------------      
c         
c GET WEIGHTS AND ABSCISSAS FOR QUADRATURE                              
c         
      eps = 1.e-10
      call hermit(nn, x, a, eps)
      do 4 i = 1, nn
    4 a(i) = a(i) * exp(x(i) * x(i))
      nl = nn
      nm = nn
      nuu = min(2 * nn,20)
      call hermit(nuu, xuu, auu, eps)
      do 40 i = 1, nuu
c *** CA ***                                                            
c         
   40 auu(i) = auu(i) * exp(xuu(i) * xuu(i))
c ----------------------------------------------------------------      
c         
c     IF INPUT = ONLY ONE BLOCK, READ IN DATA                           
c         
      xdev = (.5 * xall(7)) / xall(1)
    5 if (nread .ne. 1) goto 7
      rewind(unit=id1) 
c -----------------------------------------------------------------     
c         
c BRANCH TO ITERATION OR LIKELIHOOD-RATIOS                              
c         
c                                                                       
c         
      read(unit=id1,err=3260) (ib(i),i = 1, 7502)
    7 goto (1, 2), ictrl
    1 nox = 0
      nev = 0
      asclq = sclq
	call flush(61)
	call flush(62)
      call iter(xall, nall, itp, ip, iq, ictrl, h, asclq, tol, trupb)
      goto 10   
2     call lratio(xall, nall, itp, iq, ictrl, nn)
      goto 10
  500 rewind(unit=52) 
      read(unit=52, fmt=3) inum, nfam1
    3 format(2i6)
      if (inum .eq. 2) then
      inum = 1
      read(unit=52, fmt=13) title
   13 format(a12)
      do 17 i = 1, nfam1
   17 read(unit=52, fmt=15) cfam(1), cgr(1)
      i = 1
   14 read(unit=52, fmt=15, end=16) cfam(i), cgr(i)
   15 format(2a6)
      i = i + 1
      goto 14
   16 i = i - 1
      rewind(unit=52) 
      write(unit=52, fmt=3) inum - 1, nfam1
      write(unit=52, fmt=13) title
      do 18 j = 1, i
   18 write(unit=52, fmt=15) cfam(i), cgr(i)
      endfile(unit=52) 
      end if
      call fdate(ctime)
c                                                                       
c         
      write(unit=62, fmt=1361) ctime
      write(unit=61, fmt=1362) ctime
      write(unit=61, fmt=1363) 
      write(unit=62, fmt=1363) 
 1363 format(1h0,130(1h=))
      return 
3260  write(unit=61,fmt=2000)
      write(unit=62,fmt=2000)
2000  format(30(1h*),' W A R N I N G ',30(1h*),/,30x,
     &     'No records have been read from the intermediate files',/,
     &     30x,'This may not be an error but may be due to selection',
     &     ' via EMX')
      call fdate(ctime)
      write(unit=62, fmt=1361) ctime
      write(unit=61, fmt=1362) ctime
      write(unit=61, fmt=1363) 
      write(unit=62, fmt=1363) 

      end
c                                                                       
c         
c **********************************************************************
c**       
c                                                                       
c         
      function qf(x)
c --------------------------------------------------------------------  
c         
c                                                                       
c         
c FINDS UPPER/LOWER TAIL AREA OF A NORMAL CURVE WITH ADAMS METHOD       
c         
c AS PRESENTED BY HILL(73) APPLIED STATISTICS ALGORITHM 66.             
c         
c                                                                       
c         
c --------------------------------------------------------------------- 
c         
      implicit double precision (o-z, a-h)
      logical up
c                                                                       
c         
c CONSTANTS C1,C2 AND C3 MUST ADJUSTED TO COMPUTER CHARACTERISTICS      
c         
c C1:=(N+9)/3 WHERE N=NUMBER OF DECIMAL DIGITS AVAILABLE                
c         
c C2 IS SUCH THAT EXP(-0.5*SQR(C2))/(C2*SQRT(2*PI))                     
c         
c    IS JUST GREATER THAN THE SMALLEST ALLOWABLE REAL NUMBER.           
c         
c                                                                       
c         
c HARRIS                                                                
c         
c     C1=7.                                                             
c         
c     C2=13.1                                                           
c         
c                                                                       
c         
c VAX                                                                   
c         
      up = .true.
      c1 = 8
c                                                                       
c         
      c2 = 13.01
      c3 = 1.28
      z = x
      if (z .ge. 0.0) goto 10
      up = .not. up
      z = - z
   10 if ((z .le. c1) .or. (up .and. (z .le. c2))) goto 20
      qf = 0.0
      goto 40
   20 y = (.5 * z) * z
c                                                                       
c         
      if (z .gt. c3) goto 30
      qf = .5 - (z * (0.398942280444 - ((0.399903438504 * y) / ((y + 
     &5.75885480458) - (29.8213557808 / ((y + 2.62433121679) + (
     &48.6959930692 / (y + 5.92885724438))))))))
c                                                                       
c         
      goto 40
c                                                                       
c         
   30 qf = (0.398942280385 * exp(- y)) / ((z - 3.8052e-8) + (
     &1.00000615302 / ((z + 3.98064794e-04) + (1.98615381364 / ((z - 
     &0.151679116635) + (5.29330324926 / ((z + 4.8385912808) - (
     &15.1508972451 / ((z + 0.742380924027) + (30.789933034 / (z + 
     &3.99019417011)))))))))))
   40 if (.not. up) qf = 1. - qf
      return 
      end
c                                                                       
c         
c *****************************************************************     
c         
c                                                                       
c         
      subroutine sehz(xall, nall, itp, se)
c                                                                       
c         
c ------------------------------------------------------------------    
c         
c                                                                       
c         
c GET S.E. OF H AND Z FOR FINAL OUTPUT                                  
c         
c                                                                       
c         
c ------------------------------------------------------------------    
c         
      implicit double precision (o-z, a-h)
      common /auxp1/ p(3), g(3), clm2, clm1, cxp2, cxp1(3), rcdg, 
     &cdp, z(10), vv, uu, gg, ca, ck, ea, ek, scals(4), sclq, rsqw
     &, notrs, ng1, nox
      common /secc/ seca, seck, b(15, 15), saveb(15, 15), ifail
      dimension xall(15), itp(15), se(15)
      seca = 0.0
      seck = 0.0
      ick = 0
      ica = 0
      k = 0
      do 10 i = 1, nall
      if (itp(i) .ne. 1) goto 10
      k = k + 1
      if (i .eq. 6) ick = k
      if (i .eq. 7) ica = k
   10 continue
c CK ITERATED                                                           
c         
      if (ick .eq. 0) goto 50
c V NOT ITERATED                                                        
c         
      if (itp(1)) 20, 20, 30
   20 seck = se(6)
      se(6) = sqrt(b(ick,ick) / (xall(1) * xall(1)))
c V ITERATED                                                            
c         
      goto 50
   30 seck = se(6)
      v2 = xall(1) * xall(1)
      se(6) = ((b(ick,ick) / v2) - (((2. * xall(6)) * b(1,ick)) / (v2 * 
     &xall(1)))) + (((xall(6) * xall(6)) * b(1,1)) / (v2 * v2))
      se(6) = sqrt(se(6))
c CA ITERATED. IT IS ASSUMED CK WAS ALSO ITERATED                       
c         
   50 if (ica .eq. 0) goto 100
      seca = se(7)
      ck2 = xall(6) * xall(6)
      se(7) = ((b(ica,ica) / ck2) - (((2. * xall(7)) * b(ick,ica)) / (
     &ck2 * xall(6)))) + (((xall(7) * xall(7)) * b(ick,ick)) / (ck2 * 
     &ck2))
      se(7) = sqrt(se(7))
  100 continue
      return 
      end
c                                                                       
c         
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+        
c                                                                       
c         
      subroutine step(f, x, n, xall, nall, itp, fsave, t, tbnd, nfe, h, 
     &trupb, idg, ibnd, ip)
c-----------------------------------------------------------------------
c         
c                                                                       
c         
c PERFORMS LINE SEARCH ALONG DIRECTION P                                
c         
c                                                                       
c         
c-----------------------------------------------------------------------
c         
      implicit double precision (o-z, a-h)
      common /mini/ tl(15, 15), d(15), g(15), gsave(15), y(15), p(15), 
     &tmin, fsmf, idif, isw, iret
c-----------------------------------------------------------------------
c         
c INITIALIZATION                                                        
c         
      dimension x(15), xall(15), itp(15), xt(15)
      sumt = 0.0
      idg = 0
c-----------------------------------------------------------------------
c         
c OBTAIN FT=F(X+T*P)                                                    
c         
      curv = 0.75
      do 10 i = 1, n
   10 xt(i) = x(i) + (t * p(i))
      call fun(ft, xt, n, xall, nall, itp)
c-----------------------------------------------------------------------
c         
c TEST SLOPE OF LINE FROM F TO FT                                       
c         
      nfe = nfe + 1
c                                                                       
c         
c-----------------------------------------------------------------------
c         
c                                                                       
c         
c DECREASING, HENCE T DOUBLED                                           
c         
c                                                                       
c         
      if (f - ft) 200, 200, 100
  100 sumt = sumt + t
  101 continue
c CHECK BOUNDS IF ANY                                                   
c         
      twot = t + t
      if (ibnd) 110, 110, 102
  102 if (tbnd) 110, 104, 104
c ACTIVE BOUNDARY CONSTRAINT                                            
c         
  104 if (twot - tbnd) 110, 106, 106
  106 f = ft
      do 108 i = 1, n
  108 x(i) = xt(i)
      write(unit=ip, fmt=1001) 
 1001 format(/1x,10(1h*),26hACTIVE BOUNDARY CONSTRAINT,10(1h*))
c-----------------------------------------------------------------------
c         
c NO ACTIVE BOUNDARY CONSTRAINT                                         
c         
      goto 310
  110 continue
      do 112 i = 1, n
      x(i) = xt(i)
  112 xt(i) = x(i) + (t * p(i))
      call fun(f2t, xt, n, xall, nall, itp)
c-----------------------------------------------------------------------
c         
c TEST SLOPE OF LINE FROM FT TO F2T                                     
c         
      nfe = nfe + 1
      if (f2t - ft) 120, 115, 115
  115 f = ft
c TRANSFER TO EXIT SECTION IF MINIMUM BRACKETED                         
c         
c IF NOT, F, FT, F2T MONOTONIC DECREASING                               
c         
c THEN, CHECK CURVATURE                                                 
c         
      goto 300
c WILL GO TO 124 IF CONCAVE, TO 122 IF CONVEX                           
c         
c THEN, CHECK IF CONVEX STEEP OR CONVEX FLAT                            
c         
  120 if (((f2t + f) - ft) - ft) 124, 122, 122
c COMES TO 124 IF CONVEX STEEP OR IF CONCAVE, GO TO 126 IF CONVEX AND FL
cAT       
  122 if ((ft - f2t) - (curv * (f - ft))) 126, 124, 124
  124 sumt = sumt + t
      t = t + t
      ft = f2t
      goto 101
  126 sumt = sumt + t
      f = f2t
      do 128 i = 1, n
c TRANSFER TO EXIT SECTION                                              
c         
  128 x(i) = xt(i)
c                                                                       
c         
c-----------------------------------------------------------------------
c         
c ASCENDING LINE FROM F TO FT                                           
c         
c                                                                       
c         
c T WILL BE HALVED IF SOME CONDITIONS TESTED HERE ARE MET               
c         
      goto 300
  200 if (f - fsave) 208, 202, 208
c COMES TO 206 IF NO PROGRESS YET MADE, SWITCH TO CENTRAL DIFFERENCE HAS
c         
c ALREADY OCCURRED, AND LIMIT ON POSSIBLE ACCURACY MAY HAVE BEEN REACHED
c.        
c IN FACT, WILL EXIT IF T LESS THAN TMIN                                
c         
  202 if (idif - 2) 208, 206, 208
  206 idg = 2
c-----------------------------------------------------------------------
c         
c T IS HALVED                                                           
c         
  208 if (t - tmin) 400, 220, 220
  220 t = .5 * t
      do 222 i = 1, n
  222 xt(i) = x(i) + (t * p(i))
      call fun(ft2, xt, n, xall, nall, itp)
c CHECK IF MINIMUM IS BRACKETED. IF YES, GO TO EXIT SECTION. IF NOT, WIL
cL        
c FIND SMALLER T VALUE IN NEXT BOX                                      
c         
      nfe = nfe + 1
c MINIMUM BRACKETED                                                     
c         
      if (ft2 - f) 240, 243, 243
  240 sumt = sumt + t
      f = ft2
      idg = 0
      do 242 i = 1, n
  242 x(i) = xt(i)
c-----------------------------------------------------------------------
c         
c MINIMUM NOT BRACKETED, F, FT2, FT MONOTONIC INCREASING                
c         
c P BEING A DESCENT DIRECTION, THERE MUST BE A SMALLER T FOR WHICH F IS 
c         
c DECREASED, THAT LIES BETWEEN T AND T/2                                
c         
c CHECK CURVATURE                                                       
c         
      if (sumt - tmin) 400, 300, 300
c CONCAVE                                                               
c         
  243 if (((ft + f) - ft2) - ft2) 244, 244, 246
  244 scal = .1
c CONVEX, PARABOLIC INTERPOLATION PERFORMED                             
c         
      goto 248
  246 scal = 1.0 + ((.5 * (f - ft)) / (((f + ft) - ft2) - ft2))
c-----------------------------------------------------------------------
c         
c RESTART WITH SMALLER T                                                
c         
      scal = max(.1d0,scal)
  248 t = scal * t
      do 250 i = 1, n
  250 xt(i) = x(i) + (t * p(i))
      call fun(ft, xt, n, xall, nall, itp)
      nfe = nfe + 1
      if (f - ft) 200, 200, 252
  252 sumt = sumt + t
      idg = 0
      f = ft
      do 254 i = 1, n
c IF MINIMUM NOT BRACKETED, CONTINUE SEARCH BY RETURNING TO 200         
c         
c IF MINIMUM BRACKETED, GO TO EXIT SECTION                              
c         
  254 x(i) = xt(i)
c                                                                       
c         
c-----------------------------------------------------------------------
c         
c                                                                       
c         
c EXIT SECTION                                                          
c         
      if (t - tmin) 400, 300, 300
  300 t = sumt
  310 continue
      do 312 i = 1, n
  312 gsave(i) = g(i)
      isw = 2
      iret = 2
      goto 500
  400 iret = 1
  500 continue
      return 
      end
c                                                                       
c         
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         
c                                                                       
c         
      subroutine update(n, t, ip, idg)
c ----------------------------------------------------------------------
c-------- 
c                                                                       
c         
c     GOLDFARB METHOD FOR VARIABLE METRIC B-MATRIX UPDATING             
c         
c     UPDATE OF FACTORIZED B-MATRIX BY RANK TWO MODIFICATION            
c         
c     IN REAL PRODUCT FORM WITH FORMULA (3.13) OR (3.16) , DEPENDING WHE
cTHER     
c     NONE OR AT LEAST ONE LAMBA(J)**2 IS GREATER THAN 4.               
c         
c     OPTIMAL CONDITIONING OF DAVIDON ALSO INCORPORATED                 
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c-------- 
      implicit double precision (o-z, a-h)
      common /mini/ tl(15, 15), d(15), g(15), gsave(15), y(15), p(15), 
     &tmin, fsmf, idif, isw, iret
      dimension wtil(15), ztil(15), w(15), z(15), wtjp1(15), ztjp1(15), 
     &s(15), dp(15)
      real nu, muj, lambj, lambj2
c                                                                       
c         
c ****** WARNING ****** NOTE THAT GSTAR IS IN ARRAY G WHILE G IS IN ARRA
cY GSAVE, 
c     GSTAR BEING GRADIENT AT IT (NIT+1) AND G GRADIENT AT IT (NIT)     
c         
c                                                                       
c         
c ----------------------------------------------------------------------
c-------- 
c SOLVES (3.10) AND (3.11) FOR WTIL AND ZTIL                            
c         
      iswup = 1
      wtil(1) = - gsave(1)
      ztil(1) = y(1)
      do 10 k = 2, n
      wtil(k) = - gsave(k)
      ztil(k) = y(k)
      km1 = k - 1
c SAVE TL LOWER TRIANG IN UPPER TRIANG OF ARRAY TL, IN CASE A SWITCH TO 
c         
c     (3.16) SHOULD BE REQUIRED                                         
c         
      do 10 i = 1, km1
      tl(i,k) = tl(k,i)
      wtil(k) = wtil(k) - (tl(k,i) * wtil(i))
c ----------------------------------------------------------------------
c-------- 
c GET SCALARS B,C AND D, DENOTED SB,SC,SD HERE                          
c         
   10 ztil(k) = ztil(k) - (tl(k,i) * ztil(i))
      sb = 0.0
      sc = 0.0
      sd = 0.0
      do 20 i = 1, n
      sb = sb + (((t * ztil(i)) * wtil(i)) / d(i))
      sc = sc + ((((t * t) * wtil(i)) * wtil(i)) / d(i))
c ----------------------------------------------------------------------
c-------- 
c OPTIMAL CONDITIONING OF UPDATE, ACCORDING TO THEOREM 3 OF DAVIDON (197
c5)       
c --------------------                                                  
c         
c OPTIMAL CONDITIONING MAY BE DISABLED BY SETTING IOC TO 1              
c         
c IN SUCH CASE, UPDATE IS DFP-UPDATE FOR B (ALPHA=0)                    
c         
   20 sd = sd + ((ztil(i) * ztil(i)) / d(i))

      ioc = 0
      if (ioc - 1) 205, 24, 205
c PREVENTING UNDERFLOW                                                  
c         
c IF ANY OF SB,SC,SD IS SMALLER THAN 1.E-10, USE SR1-UPDATE FOR B       
c         
  205 continue
      fbcd = min(sb,sc,sd)
      if (fbcd - 1.e-10) 21, 21, 210
c --------------------                                                  
c         
  210 continue
      fbcd = ((2. * sc) * (sd / sb)) / (sc + sd)
c SR1 UPDATE                                                            
c         
      if (fbcd - 1.0) 21, 22, 22
   21 alpha = -1.0
      if (sc - sb) 211, 600, 211
  211 if (sb - sd) 212, 600, 212
  212 if ((sd - (2. * sb)) + sc) 213, 600, 213
  213 sa = 1. / (sqrt(abs(sc - sb)) * sqrt(abs(sb - sd)))
      thet1 = - ((((sd - sb) * sa) + 1.) / ((sd - (2. * sb)) + sc))
      thet2 = sa + (((sa * (sb - sc)) + 1.) / ((sd - (2. * sb)) + sc))
c --------------------                                                  
c         
c RANK TWO UPDATE                                                       
c         
c GET ALPHA                                                             
c         
      goto 29
   22 aa = ((sb / sc) - (2. * (sd / sb))) + (sd / sc)
      bb = (sb / sc) - 1.
      cc = 1. - (sb / sd)
      del2 = (bb * bb) - (aa * cc)
c IF DISCRIMINANT NEGATIVE OR EQUAL ZERO, TAKE ALPHA EQUAL TO ZERO      
c         
c ---------------                                                       
c         
      if (del2 - 1.e-08) 24, 24, 25
c THIS IS DFP-UPDATE FOR B                                              
c         
   24 alpha = 0.0
      sa = 1. / (sqrt(sb) * sqrt(sc))
      thet1 = - (1. / sc)
      thet2 = sa
      goto 29
   25 del = sqrt(del2)
      alph1 = ((- bb) + del) / aa
c FOR NOW, ALWAYS CHOSE ROOT OF SMALLEST MODULUS                        
c         
      alph2 = ((- bb) - del) / aa
      alpha = alph1
      if (abs(alph1) - abs(alph2)) 252, 252, 251
  251 alpha = alph2
c IF ALPHA VERY SMALL, ALPHA TAKEN EQUAL TO ZERO ( DFP-UPDATE FOR B)    
c         
  252 continue
c GET SA                                                                
c         
      if (abs(alpha) - 1.e-05) 24, 24, 240
  240 sa = (((((alpha + 1.) * (alpha + 1.)) + (sc / sb)) - (((alpha * 
     &alpha) * (sc / sb)) * (sd / sb))) - 1.0) + ((alpha * alpha) * (sd
     & / sb))
      if (sa) 26, 26, 27
   26 sa = 0.0
      goto 28
   27 sa = sqrt(sa)
c GET THET1 AND THET2 FOR NON TRIVIAL ALPHA                             
c         
      sa = 1. / (sa * sb)
   28 rdiv = 1. / ((((alpha * alpha) * sd) + ((2. * alpha) * sb)) + sc)
      thet1 = - (((sa * (alpha * ((alpha * sd) + sb))) + 1.) * rdiv)
      thet2 = sa + ((((alpha * sa) * (sc + (alpha * sb))) - alpha) * 
     &rdiv)
c --------------------                                                  
c         
c ----------------------------------------------------------------------
c-------- 
c GET W AND Z AS GIVEN JUST AFTER (3.11)                                
c         
   29 continue
      do 30 i = 1, n
      w(i) = (t * wtil(i)) + (alpha * ztil(i))
c ----------------------------------------------------------------------
c-------- 
c READY TO APPLY GOLDFARB RECURRENCE 3 TO SOLVE CONCURRENTLY            
c         
c     (3.13) , (3.12) BEING ALSO SOLVED SIMULTANEOUSLY                  
c         
c                                                                       
c         
c GET THE S(J) FIRST (P.802,TOP) AND INITIALIZE                         
c         
   30 z(i) = ((t * thet1) * wtil(i)) + (thet2 * ztil(i))
      np1 = n + 1
      nm1 = n - 1
      s(n) = 0.0
      do 35 k = 1, nm1
      j = np1 - k
      jm1 = j - 1
   35 s(jm1) = s(j) + ((w(j) * w(j)) / d(j))
      nu = 1.0
c ----------------------------------------------------------------------
c-------- 
c INITIALIZE                                                            
c         
      eta = 0.0
  350 if (iswup - 1) 351, 351, 352
  351 do 37 i = 1, n
      wtjp1(i) = - gsave(i)
   37 ztjp1(i) = y(i)
      goto 353
  352 do 38 i = 1, n
      wtil(i) = (- (t * g(i))) + (alpha * y(i))
c ----------------------------------------------------------------------
c-------- 
   38 ztil(i) = (- ((t * thet1) * g(i))) + (thet2 * y(i))
  353 do 50 j = 1, nm1
      jp1 = j + 1
      if (iswup - 1) 354, 354, 355
  354 do 39 k = jp1, n
      wtjp1(k) = wtjp1(k) - (wtil(j) * tl(k,j))
   39 ztjp1(k) = ztjp1(k) - (ztil(j) * tl(k,j))
      goto 356
  355 do 40 k = jp1, n
      wtjp1(k) = wtil(k) - (w(j) * tl(k,j))
c RECURRENCE 3 TO GET AJ,BJ,ETC...                                      
c         
   40 ztjp1(k) = ztil(k) - (z(j) * tl(k,j))
  356 aj = (nu * z(j)) - (eta * w(j))
      thj = 1.0 + ((aj * w(j)) / d(j))
      lambj2 = (thj * thj) + (((aj * aj) * s(j)) / d(j))
      if (iswup - 1) 449, 449, 451
c SWITCH TO (3.16) UPDATE                                               
c         
  449 if (lambj2 - 10.) 451, 451, 450
  450 iswup = 2
      do 452 k = 2, n
      km1 = k - 1
      do 452 i = 1, km1
  452 tl(k,i) = tl(i,k)
      goto 350
c UPDATES DPLUS(J)                                                      
c         
  451 continue
      dp(j) = d(j) * lambj2
c TAKES SIGN OF LAMBDA OPPOSITE TO THAT OF THETA                        
c         
      lambj = sqrt(lambj2)
      if (thj) 42, 42, 41
   41 lambj = - lambj
   42 muj = thj - lambj
      bj = (thj * w(j)) + (aj * s(j))
      gamlj = (bj * nu) / (lambj2 * d(j))
c NOTE THAT GAMLJ AND BETLJ STAND FOR GAMMA TILDA AND BETA TILDA RESP. I
cN (3.14) 
      betlj = (aj - (bj * eta)) / (lambj2 * d(j))
      nu = - (nu / lambj)
c UPDATES J-TH COLUMN OF TL                                             
c         
      eta = - ((eta + ((aj * aj) / (muj * d(j)))) / lambj)
      if (iswup - 1) 357, 357, 358
  357 do 45 k = jp1, n
   45 tl(k,j) = (tl(k,j) + ((t * (betlj + (thet1 * gamlj))) * wtjp1(k)))
     & + (((alpha * betlj) + (thet2 * gamlj)) * ztjp1(k))
      goto 50
  358 do 46 k = jp1, n
      tl(k,j) = ((tl(k,j) / lambj2) + (betlj * wtil(k))) + (gamlj * ztil
     &(k))
      wtil(k) = wtjp1(k)
   46 ztil(k) = ztjp1(k)
   50 continue
      aj = (nu * z(n)) - (eta * w(n))
c NOTICE THAT THIS IS ACTUALLY LAMBDA(N), NEEDED TO UPDATE D(N)         
c         
      lambj = 1.0 + ((aj * w(n)) / d(n))
c TRANSFER DP IN D                                                      
c         
      dp(n) = (d(n) * lambj) * lambj
      do 60 i = 1, n
   60 d(i) = dp(i)
c COMES HERE IF DIVISION BY ZERO WOULD OCCUR IN OPTIMAL CONDITIONING    
c         
      return 
  600 idg = 3
      return 
      end
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c:::::::: 
      function xno(icol)
      implicit double precision (o-z, a-h)
      character a*1, d*1
      character finput*30
      character format*7
      dimension d(10)
      common /sy1a/ ibegin, i1, i2, idum1, index, nchar
      common /sy1b/ a(80)
      data d / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
      call nsp1
c                                                                       
c         
      goto 5
c                                                                       
c         
      entry xno1(icol)
c                                                                       
c         
      call nspace
c                                                                       
c         
      entry xno2(icol)
    5 finput = ' '
      format = '(F  . )'
      index = 1
      n = 0
      l = 0
c                                                                       
c         
c     SEARCHES FOR  NUMBER - REAL OR INTEGER, POSITIVE OR NEGATIVE.     
c         
c          ICOL ON RETURN = 1ST NON-BLANK COLUMN AFTER NUMBER.          
c         
c          INDEX ON RETURN = 1 IF OKAY  ;  = 0 IF SOMETHING'S WRONG.    
c         
c                                                                       
c         
      m = 0
      if (a(icol) .eq. '-') then
      n = 1
      finput(1:1) = a(icol)
      icol = icol + 1
      end if
    6 if (((((((((((a(icol) .eq. '0') .or. (a(icol) .eq. '1')) .or. (a(
     &icol) .eq. '2')) .or. (a(icol) .eq. '3')) .or. (a(icol) .eq. '4'))
     & .or. (a(icol) .eq. '5')) .or. (a(icol) .eq. '6')) .or. (a(icol)
     & .eq. '7')) .or. (a(icol) .eq. '8')) .or. (a(icol) .eq. '9'))
     & .or. (a(icol) .eq. '.')) then
      n = n + 1
      finput(n:n) = a(icol)
      icol = icol + 1
      if (icol .gt. nchar) goto 38
      if (l .eq. 1) m = m + 1
      if (a(icol) .eq. '.') l = 1
      else
      goto 7
      end if
      goto 6
    7 if (n .ne. 0) goto 40
   38 xno = 0.
      index = 0
      return 
   40 m1 = 0
      if (n .gt. 9) then
      m1 = n / 10
      n = n - (10 * m1)
      end if
      format(3:3) = d(m1 + 1)
      format(4:4) = d(n + 1)
      format(6:6) = d(m + 1)
c     BRINGS ICOL TO NEXT NON-BLANK CHARACTER                           
c         
      read(unit=finput, fmt=format) xno
      jcol = icol
      call nsp1
      end
