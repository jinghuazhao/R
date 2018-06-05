       subroutine family(famdata,famsize,pobs,p,stat,toenum,tailp,sump,
     & nenum)
       integer famsize, toenum
       integer famdata(famsize,3)
c      based on a version for Sun (without implicit) dated 27/2/2003
c      Exact test of disease clustering for family members


c       DZ   3/96,   7/99, 1/02


      double precision zero,const,p,stat(20)
      integer maxfac,m(20),sib,aff,freq,oldsib,ns,
     &   naff,nsibs,nfam,fm(20,20),j,maxsize
c     logical done
      double precision fac(0:8000),fac0(8001),ptail,psum,dcase
      double precision pobs,tailp,sump,nenum
      equivalence (fac(1),fac0(2))
      common/factab/fac0
      common /jhzhao/psum,ptail,dcase
      data maxfac/8000/, zero/0.000d0/, fm/400*0/
          
c                  build table of log factorials w/ zero subscript
      fac(0)=zero
      do 1 j=1,maxfac
      const=j
   1  fac(j)=fac(j-1)+ dlog(const)

c   read frequency data, build frequency matrix table -fm-
c     open(10,file='family.dat')
      oldsib=-1
      maxsize=1
c     do 20 j=1,10000000
      do 20 i=1,famsize
        sib=famdata(i,1)
        aff=famdata(i,2)
        freq=famdata(i,3)
        fm(aff+1,sib)=freq
        if(sib .GT. maxsize)maxsize=sib
*       if(sib .NE. oldsib)write(6,*)
        oldsib=sib
   20 continue

c Find marginal totals from -fm- and constant part of the likelihood
      call build(fm,m,1,maxsize,nfam,nsibs,naff,const)
c   fm(i,j)= # of families with j sibs, (i-1) of whom 
c               are affected
c   m(j)   = # of families with j sibs
c   nfam   = # number of families
c   nsibs  = # of siblings
c   naff   = # of affected sibs
c   const  = constant part of the log-likelihood

c  Probability of the observed table
      call prob(fm,1,maxsize,const,p)
      pobs=p
c  Generate test statistics, expected values under null hypothesis
      call test(fm,m,1,maxsize,stat,ns,naff,nsibs,.true.)
c  Enumerate all possible tables for exact tail area
      if(toenum.ne.1) return
      call enum(nsibs,naff,nfam,m,maxsize,const,p)
      tailp=ptail
      sump=psum
      nenum=dcase

c     stop 9999
      return
      end

c===========================================================

      subroutine test(fm,m,first,last,stat,ns,naff,nsibs,trace)

c    Compute test statistics for the set of family 
c       frequencies in -fm-. If -trace- then print them out

      double precision fac(0:8000),fac0(8001)
      equivalence (fac(1),fac0(2))
      common/factab/fac0
      integer first,last,m(1),fm(20,20),naff,nsibs,ns,j,i,
     & ia,ib,ic,nssave
      double precision stat(20),zero,one,two,dexp,dlog,be,
     & binp,sbe,he,she,fmij,dble,eps,chi1,chi2,dsqrt
      logical trace
      data zero/0.0d0/, one/1.0d0/, two/2.0d0/, eps/1.0d-9/
      data nssave/5/

*     if(trace)write(6,1000)
      ns=nssave
      do 10 j=1,ns
         stat(j)=zero
   10 continue  
      binp=dble(naff)/dble(nsibs)           
c         j=# of sibs
c         fm(i,j)=# of families with jsibs, i-1 affected
c         be=binomial expectation
c         he=hypergeometric expectation
*     write(6,1001)
      do 100 j=first,last
         if(m(j) .LE. 0)go to 100
         sbe=zero
         she=zero
         do 80 i=1,j+1
            ia=i-1
            be=dexp(fac(j)-fac(ia)-fac(j-ia))*m(j)
            he=be
            if(ia .GT. 0)be=be*binp**ia 
            if(ia .LT. j)be=be*(one-binp)**(j-ia)
            do 30 ib=1,j
               ic=ib-1
               he=he/(nsibs-ic)
               if(ic .LT. ia)then
                  he=he*(naff-ic)
               else
                  he=he*(nsibs-naff+ib-j)
               endif  
   30       continue
            fmij=dble(fm(i,j))
c                      two chi residuals
            chi1=zero
            chi2=zero
            if(be .GT. eps)chi1=(fmij-be)/dsqrt(be)
            if(he .GT. eps)chi2=(fmij-he)/dsqrt(he)
*           if(trace)write(6,1001)j,ia,fm(i,j),be,he,chi1,chi2
            sbe=sbe+be
            she=she+he
c                       test statistics:
c                                             deviance
            if(fmij .GT. eps)then
               if(be .GT. eps)stat(1)=stat(1)+fmij*dlog(fmij/be)
               if(he .GT. eps)stat(2)=stat(2)+fmij*dlog(fmij/he)
            endif
c                                             chi-squared
            if(be .GT. eps)stat(3)=stat(3)+(fmij-be)**2/be
            if(he .GT. eps)stat(4)=stat(4)+(fmij-he)**2/he
c                                    exact log-likelihood
            stat(ns)=stat(ns)+fmij*(fac(ia)+fac(j-ia))
   80    continue
*        if(trace)write(6,1002)m(j),sbe,she
*        if(trace)write(6,1001)
  100 continue
      stat(1)=stat(1)*two
      stat(2)=stat(2)*two
*     if(trace)write(6,1003)(stat(j),j=1,ns)
      return
 
 1000 format(/t20,'Test Statistics for Family Clusters',
     &  //,t31,'Obs''d',t42,'Expectations:',4x,
     &    'Chi Residuals',/ 
     &  t14,'# sibs   # aff    freq     Bin.        Hyp.',
     &  '    Bin.    Hyp.')
 1001 format(7x,3i9,2f11.5,2f8.2)
 1002 format(t14,'Totals:',t26,i9,2f11.5)
 1003 format(t14,'Deviance:',     t35,2f11.5,/,
     &       t14,'Chi-squared:',  t35,2f11.5,/,
     &       t14,'Log-likelihood',t35,f11.5)
      end

c===========================================================

      subroutine enum(nsibs,naff,nfam,m,maxsize,
     &    const,p0)

      logical done(20),alldone
      double precision p,ptail,const,
     &  psum,p0,one,zero,dcase,pasum,proba
      integer freq(20,20),m(20),maxsize,level,j,naff,nsibs,
     &  nfam,alloc(20),per,mod,ncase
      common /jhzhao/psum,ptail,dcase
c                    periodicity of output tables
      data per/1000000/,  one/1.00d0/, zero/0.00d0/
 
c     OUTER LOOP: allocate cases to various sized families
      alldone=.true.
      psum=zero
      pasum=zero
      ncase=0
      dcase=zero
   10 call cmulte(alloc,naff,maxsize,alldone)
      if(alldone .OR. alloc(maxsize) .GT. maxsize*m(maxsize))then
*        write(6,1002)p0,psum,ptail,dcase
         return
      endif
c       alloc(i) cases in families of size i:  Is this valid?
      do 20 j=1,maxsize
         if(alloc(j) .GT. j*m(j))go to 10
   20 continue
      call pralloc(alloc,m,maxsize,naff,nsibs,proba)
      pasum=pasum+proba
c     Shortcut: If proba < p0 then skip inner loop
      if(proba .LE. p0)then
         psum=psum+proba
         ptail=ptail+proba
         go to 10
      endif

c INNER LOOP: generate family frequencies with these allocations
      level=maxsize
   40 done(level)=.true.
   50 call cfe(m(level),alloc(level),level,freq(1,level),done(level))
c                            go up one level when done at this level
      if(.NOT. done(level))go to 70
      if(level .GE. maxsize)go to 10
      level=level+1
      go to 50
c                                 evaluate and go down to next level
   70 level=level-1
      if(level .GE. 1)go to 40
      ncase=ncase+1
      dcase=dcase+one
      call prob(freq,1,maxsize,const,p)
      psum=psum+p
c                                           accumulate tail area
      if(p .LE. p0)ptail=ptail+p
c                                 periodically print out tables
      if(mod(ncase,per) .EQ. 0)then
c        write(6,1001)ncase,(alloc(j),j=1,maxsize)
         call out(freq,m,nsibs,naff,nfam,1,maxsize)
c        write(6,1002)p,psum,ptail,dcase
      endif
      level=1
      go to 50

cc 1001 format(/' ENUM:ncase: ',i14,' alloc: ',15i3)
 1002 format(/' Prob of this table:',e16.8,  
     &   4x, ' Cumulative prob: ',e16.8, /
     &       ' Tail probability:  ',e16.8,
     &   4x, ' Case number:     ',e16.8)
      end

c  ========================================================

      subroutine build(fm,m,first,last,nfam,nsibs,naff,const)

c    Find marginal totals from -fm- and constant part 
c      of the likelihood
      
c   fm(i,j)= # of families with j sibs, i-1 are affected
c   m(j)   = # of families with j sibs
c   nfam   = # number of families
c   nsibs  = # of siblings
c   naff   = # of affected sibs
c   const  = constant part of the log-likelihood

      integer fm(20,20),naff,nsibs,first,last,i,j,
     &    m(20),fmij,nfam
      double precision const,zero
      double precision fac(0:8000),fac0(8001)
      equivalence (fac(1),fac0(2))
      common/factab/fac0
      data zero/0.0000d0/

      nfam=0
      nsibs=0
      naff=0
      const=zero
      do 20 j=first,last
         m(j)=0
         do 10 i=1,j+1
            fmij=fm(i,j)
            naff=naff+(i-1)*fmij
            m(j)=m(j)+fmij
   10    continue
         const=const+fac(j)*m(j)+fac(m(j))
         nfam=nfam+m(j)
         nsibs=nsibs+j*m(j)
   20 continue       
      const=const-fac(nsibs)
      const=const+fac(naff)+fac(nsibs-naff)
      return
      end

c ====================================================

      subroutine prob(fm,first,last,const,p)

c  Find the probability of the table -fm-

      integer fm(20,20),first,last,i,j
      double precision zero,const,p,dexp,lminf,slf
      double precision fac(0:8000),fac0(8001)
      equivalence (fac(1),fac0(2))
      common/factab/fac0
c       lminf is log of smallest positive dp number
      data zero/0.00d0/, lminf/-708.75d0/

      slf=const
      do 20 j=first,last
         do 10 i=1,j+1
            slf=slf - (fac(i-1)+fac(j-i+1))*fm(i,j)
            slf=slf - fac(fm(i,j))
   10    continue
   20 continue
c  Is the probability positive or zero in double precision?
      p=zero
      if(slf.gt.lminf)p=dexp(slf)
      return
      end

c =======================================================

      subroutine out(freq,m,nsibs,naff,nfam,first,last)

c Print a table of frequencies and check for inconsistencies

      logical error
      integer freq(20,20),m(20),nsibs,naff,nfam,i,j,
     &    first,last,csibs,cfam,caff,cm(20)

      csibs=0
      cfam=0
      caff=0
      error=.false.
c     write(6,1001)nsibs,naff,nfam
      error = error .OR. (nsibs .LT. 0)
      error = error .OR. (naff .LT. 0)
      error = error .OR. (nfam .LT. 0)
      do 30 j=first,last
         cm(j)=0
c        write(6,1006)m(j),(freq(i,j),i=1,j+1)
         error = error .OR. (m(j) .LT. 0)
         do 20 i=1,j+1
            cfam=cfam+freq(i,j)
            cm(j)=cm(j)+freq(i,j)
            caff=caff+(i-1)*freq(i,j)
            error=error .OR. (freq(i,j) .LT. 0)
   20    continue
   30 continue
c                    Check for inconsistencies
      if(error)go to 900
      if(caff .NE. naff)go to 900
      if(cfam .NE. nfam)go to 900
      do 50 j=first,last
         if(cm(j) .NE. m(j))go to 900
   50 continue
      return
c                        ERRORS !
c 900 write(6,1000)
  900 continue
c     write(6,1001)csibs,caff,cfam
c     write(6,1006)(cm(j),j=first,last)
      call rexit("8888")

cc 1000 format(' OUT: error detected',3i5)
cc 1001 format(/1x,i5,' sibs',i8,' affected',i8,' families')
cc 1006 format(1x,i4,5x,20i4)
      end

c =====================================================

       subroutine cfe(n,m,i,x,done)

c   Enumerate all possible cluster frequencies for n families, 
c   all of size i with m infected individuals.  x(j) is the 
c   number of families with j-1 infected individuals.  
c   Done signifies the start and end of the sequence.

      logical done
      integer x(20),m,n,i,j,k,sn,sm,ip1

      ip1=i+1
c                     test for valid parameter values
      if(i .LT. 1) call rexit("440")
      if(m .LT. 0) call rexit("441")
      if(m .GT. i*n) call rexit("442")
c      Special cases where only one outcome is possible:
c                                                    m=i*n
      if(m .EQ. i*n)then
         done=.NOT. done
         do 2 j=1,i
    2    x(j)=0
         x(ip1)=n
         return
      endif
c                                                    n=0 or 1
      if(n .LE. 1)then
         done=.NOT. done
         do 4 j=1,ip1
    4    x(j)=0
         if(m .GT. i) call rexit("443")
         x(m+1)=n
         return
      endif
c                                               m=0 or 1; or i=1
      if(i .EQ. 1 .OR. m .LE. 1)then
         done=.NOT. done
         do 6 j=1,ip1
    6    x(j)=0
         x(1)=n-m
         x(2)=m
         return
      endif
c                                Initialize the general case
      if(done)then
         j=m/n
         j=j+1
         if(j .GT. i) call rexit("444")
         do 20 k=1,ip1
   20    x(k)=0
c                     two smallest possible frequencies
         x(j+1)=m-(j-1)*n
         x(j)=j*n-m
         done=.false.
         return
      endif
c        Otherwise update existing frequencies in -x-
   25 j=3
   30 x(j)=x(j)+1
c               How many frequencies are aleady used ?
      sn=n
      sm=m
      do 40 k=3,ip1
         sn=sn-x(k)
         sm=sm-(k-1)*x(k)
   40 continue
c              Are a valid number of frequencies left to us?
      if(0 .LE. sm  .AND.  sm .LE. sn)then
         x(2)=sm
         x(1)=sn-sm
         return
      endif
      if(0 .LE. sn .AND. sn .LT. sm)go to 25
c      Next column of 'odometer'; done when we run off the end
      x(j)=0
      j=j+1
      if(j .LE. ip1)go to 30
      done=.true.
      return
      end

c ========================================================

      subroutine pralloc(alloc,m,maxsize,naff,nsibs,proba)

c  Multivariate hypergeometric probability of an 
c  allocation -alloc- ('m' in paper) of cases

      integer m(1),alloc(1),maxsize,naff,nsibs,j
      double precision proba,dexp,lminf
c   table of log factorials, with zero subscript
      double precision fac(0:8000),fac0(8001)
      equivalence (fac(1),fac0(2))
      common/factab/fac0
c    log of smallest positive double precision number
      data lminf/-708.75d0/


c                  remove these comments to bypass the shortcut
c                  i=0
c                  proba=1.00d0
c                  if(i.eq.i)return


      proba=fac(naff)+fac(nsibs-naff)-fac(nsibs)
      do 10 j=1,maxsize
         proba=proba+fac(j*m(j))
         proba=proba-fac(alloc(j))
         proba=proba-fac(j*m(j)-alloc(j))
   10 continue
      if(proba .LT. lminf)proba=lminf
      proba=dexp(proba)
      return
      end

c ======================================================

      subroutine cmulte(n,m,k,done)

c  On successive calls, generate the complete multinomial outcomes
c  in k categories with sample size m into vector n().  
c  DONE signifies completion of the task or initialization on input.
       
      integer i,j,k,m,n(1),sum
      logical done

      if(k .EQ. 1)then
         n(1)=m
         done=.not.done
         return
      endif
      if(m .EQ. 0)then
         done=.NOT. done
         do 10 j=1,k
   10    n(j)=0
         return
      endif

      if(done)go to 500

      j=2
c                        generate next vector in the sequence
  100 n(j)=n(j)+1
c                                             find cumulative sum
      sum=0
      do 200 i=j,k
  200 sum=sum+n(i)
      if(sum .GT. m)go to 300
c                     n(1) is what ever is left over
      n(1)=m-sum
      return
c                    clear this column and carry over to the next
  300 n(j)=0
      j=j+1
      if(j .LE. k)go to 100
C                                       done when we run off the end
      done=.TRUE.
      return
c                        initialize on the first call to cmulte
  500 do 600 i=1,k
  600 n(I)=0
      n(1)=m
      done=.FALSE.
      return
      end


c--------------------- end of this file ----------------------
