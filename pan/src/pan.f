C***********************************************************************
C The following subroutines are used by both pan() and ecme().
C***********************************************************************
      subroutine istfin(ntot,subj,m,ist,ifin)
C creates vectors of starting and finishing positions for subjects
      integer ntot,subj(ntot),m,ist(m),ifin(m),scur
      scur=-999
      icur=0
      do 100 i=1,ntot
         if(subj(i).ne.scur) then
            icur=icur+1
            ist(icur)=i
            scur=subj(i)
         endif
 100  continue
      do 200 i=1,m-1
         ifin(i)=ist(i+1)-1
 200  continue
      ifin(m)=ntot
      return
      end
C***********************************************************************
      subroutine chfc(p,pw,s)
C overwrites s (upper tri) with cholesky factor
      integer p,pw
      double precision s(p,p),sum
      do 50 i=1,pw
         sum=dble(0.)
         do 10 k=1,i-1
            sum=sum+s(k,i)**2
 10      continue
         s(i,i)=dsqrt(s(i,i)-sum)
         do 40 j=i+1,pw
            sum=dble(0.)
            do 30 k=1,i-1
               sum=sum+s(k,i)*s(k,j)
 30         continue
            s(i,j)=(s(i,j)-sum)/s(i,i)
 40      continue
 50   continue
      return
      end
C***********************************************************************
      subroutine chl(p,pw,m,s,l)
C overwrites lth layer of s (upper tri) with cholesky factor
      integer p,pw,m,l
      double precision s(p,p,m),sum
      do 50 i=1,pw
         sum=dble(0.)
         do 10 k=1,i-1
            sum=sum+s(k,i,l)**2
 10      continue
         s(i,i,l)=dsqrt(s(i,i,l)-sum)
         do 40 j=i+1,pw
            sum=dble(0.)
            do 30 k=1,i-1
               sum=sum+s(k,i,l)*s(k,j,l)
 30         continue
            s(i,j,l)=(s(i,j,l)-sum)/s(i,i,l)
 40      continue
 50   continue
      return
      end
C***********************************************************************
      subroutine bkslv(p,pw,s)
C inverts an upper triangular matrix
      integer p,pw
      double precision s(p,p),sum
      s(1,1)=dble(1.)/s(1,1)
      do 10 k=2,pw
         s(k,k)=dble(1.)/s(k,k)
         do 5 j=1,k-1
            sum=dble(0.)
            do 3 i=j,k-1
               sum=sum+s(j,i)*s(i,k)
 3          continue
            s(j,k)=-sum*s(k,k)
 5       continue
 10   continue
      return
      end
C***********************************************************************
      subroutine mm(p,pw,wm,cm)
C calculates upper tri part of cm = wm%*%t(wm), where wm is upper tri
      integer p,pw
      double precision wm(p,p),cm(p,p),sum
      do 10 i=1,pw
         do 5 j=i,pw
            sum=dble(0.)
            do 2 k=max(i,j),pw
               sum=sum+wm(i,k)*wm(j,k)
 2          continue
            cm(i,j)=sum
 5       continue
 10   continue
      return
      end
C***********************************************************************
      subroutine bkslvl(p,pw,m,s,l)
C inverts an upper triangular matrix in layer l of s
      integer p,pw,m,l
      double precision s(p,p,m),sum
      s(1,1,l)=dble(1.)/s(1,1,l)
      do 10 k=2,pw
         s(k,k,l)=dble(1.)/s(k,k,l)
         do 5 j=1,k-1
            sum=dble(0.)
            do 3 i=j,k-1
               sum=sum+s(j,i,l)*s(i,k,l)
 3          continue
            s(j,k,l)=-sum*s(k,k,l)
 5       continue
 10   continue
      return
      end
C***********************************************************************
      subroutine mmul(p,pw,m,wm,l,cm)
C calculates upper tri part of cm = wm%*%t(wm), where wm is upper tri
C(works in layer l)
      integer p,pw,m,l
      double precision wm(p,p,m),cm(p,p),sum
      do 10 i=1,pw
         do 5 j=i,pw
            sum=dble(0.)
            do 2 k=max(i,j),pw
               sum=sum+wm(i,k,l)*wm(j,k,l)
 2          continue
            cm(i,j)=sum
 5       continue
 10   continue
      return
      end
C***********************************************************************
C The following are used only by pan().
C***********************************************************************
      subroutine prelimm(ntot,subj,m,ist,ifin,pcol,pred,q,
     /     zcol,ztz,patt,nstar,p,xcol,xtxinv,wkpp)
C Preliminary manipulations for pan. After execution, 
C     ist  = vector of starting positions for subject i, i=1,...,m
C     ifin = vector of finishing positions for subject i, i=1,...,m
C     ztz  = t(z_i)%*%z, i=1,...,m
C     nstar = total number of rows in y containing data
C     xtxinv = inv( sum of t(X_i)%*%X_i )
      integer ntot,subj(ntot),m,ist(m),ifin(m),pcol,
     /     q,zcol(q),nstar,s,patt(ntot),p,xcol(p)
      double precision pred(ntot,pcol),ztz(q,q,m),sum,xtxinv(p,p),
     /     wkpp(p,p)
      call istfin(ntot,subj,m,ist,ifin)
      nstar=0
      do 10 i=1,ntot
         if(patt(i).ne.0) nstar=nstar+1
 10   continue
      do 100 s=1,m
         do 90 i=1,q
            do 80 j=i,q
               sum=dble(0.)
               do 60 k=ist(s),ifin(s)
                  if(patt(k).ne.0) then
                     sum=sum+pred(k,zcol(i))*pred(k,zcol(j))
                  endif
 60            continue
               ztz(i,j,s)=sum
               if(i.ne.j) ztz(j,i,s)=sum
 80         continue
 90      continue
 100  continue
      do 150 i=1,p
         do 140 j=i,p
            sum=dble(0.)
            do 130 k=1,ntot
               if(patt(k).ne.0) then
                  sum=sum+pred(k,xcol(i))*pred(k,xcol(j))
               endif
 130        continue
            wkpp(i,j)=sum
 140     continue
 150  continue
      call chfc(p,p,wkpp)
      call bkslv(p,p,wkpp)
      call mm(p,p,wkpp,xtxinv)
      do 160 i=1,p
         do 155 j=i,p
            xtxinv(j,i)=xtxinv(i,j)
 155     continue
 160  continue
      return
      end
C***********************************************************************
      subroutine mgibbsbd(ntot,subj,m,ist,ifin,pcol,pred,q,zcol,ztz,
     /     patt,nstar,r,y,p,xcol,npatt,rmat,sflag,beta,sigma,psi,b,
     /     xtxinv,wkpp,wkpr,eps,wkrqrq1,wkrqrq2,sig,wkrr1,wkrr2,iter,
     /     wkqr1,wkqr2,wkqrv,nhyp,hyp,delta,iposn,pstfin,betas,
     /     sigmas,psis,wkqq1,wkqq2)
C conventional Gibbs sampler for the block-diagonal version of pan
      integer ntot,subj(ntot),m,ist(m),ifin(m),pcol,sflag,
     /     q,zcol(q),nstar,patt(ntot),r,p,xcol(p),npatt,rmat(npatt,r),
     /     iter,nhyp,iposn(ntot),pstfin(npatt,2),it
      integer oc(100),mc(100)
      double precision pred(ntot,pcol),ztz(q,q,m),y(ntot,r),
     /     beta(p,r),sigma(r,r),psi(q,q,r),b(q,r,m),xtxinv(p,p),
     /     wkpp(p,p),wkpr(p,r),eps(ntot,r),wkrqrq1(r*q,r*q),
     /     wkrqrq2(r*q,r*q),sig(r*q,r*q,m),wkrr1(r,r),wkrr2(r,r),
     /     wkqr1(q,r),wkqr2(q,r),wkqrv(q*r),hyp(nhyp),
     /     delta(ntot,r),betas(p,r,iter),sigmas(r,r,iter),
     /     psis(q,q,r,iter),wkqq1(q,q),wkqq2(q,q)
      double precision wkr(100)
      real junk
      junk=gauss()
      call prelimm(ntot,subj,m,ist,ifin,pcol,pred,q,zcol,ztz,patt,
     /     nstar,p,xcol,xtxinv,wkpp)
      if(sflag.ne.1) then
C get starting values of parameters
         call mimpy(ntot,r,y,patt,npatt,rmat)
         call mkxty(ntot,r,y,pcol,pred,p,xcol,patt,wkpr)
         call mkbeta(p,r,xtxinv,wkpr,beta)
         call mkeps1(ntot,r,y,pcol,pred,p,xcol,beta,eps,patt)
         call mksigma(ntot,r,eps,nstar,sigma,patt)
         call mksigbd(r,q,m,psi,sigma,ztz,sig,wkrr1,wkrr2,
     /     wkrqrq2,1,nhyp,hyp,wkqq1,wkqq2)
         call mkpsi0bd(r,q,m,psi,sig,wkrqrq1)
      else
C set starting value of y to eps
         do 50 i=1,ntot
            do 40 j=1,r
               y(i,j)=eps(i,j)
 40         continue
 50      continue
      endif
      do 100 it=1,iter
         call mkeps1(ntot,r,y,pcol,pred,p,xcol,beta,eps,patt)
         call mksigbd(r,q,m,psi,sigma,ztz,sig,wkrr1,wkrr2,
     /     wkrqrq2,0,nhyp,hyp,wkqq1,wkqq2)
         call drb(r,q,m,b,wkrr2,eps,pcol,pred,zcol,wkqr1,
     /     wkqr2,ist,ifin,patt,ntot,sig,wkrqrq2,wkqrv)
         call drpsibd(r,q,m,psi,wkqq1,wkqq2,wkqrv,nhyp,hyp,b)
         call mkeps2(ntot,m,r,y,pcol,pred,q,zcol,b,eps,patt,ist,ifin)
         call mkxty(ntot,r,eps,pcol,pred,p,xcol,patt,wkpr)
         call mkbeta(p,r,xtxinv,wkpr,beta)
         call mkeps1(ntot,r,eps,pcol,pred,p,xcol,beta,delta,patt)
         call drsigma(ntot,r,delta,nstar,sigma,patt,nhyp,hyp,wkrr1,
     /     wkrr2,p)
         call drbeta(r,sigma,xtxinv,p,beta,wkrr1,wkpp,wkpr)
         call mkeps1(ntot,r,eps,pcol,pred,p,xcol,beta,delta,patt)
         call dreps(100,100,oc,mc,100,wkr,ntot,iposn,npatt,pstfin,
     /     r,rmat,delta,sigma,wkrr1,wkrr2)
         call mky(ntot,r,pcol,pred,delta,y,p,xcol,q,zcol,beta,m,b,
     /     ist,ifin,npatt,rmat,patt)
         call storebd(iter,it,p,r,beta,q,psi,sigma,betas,sigmas,psis)
 100  continue
      return
      end
C***********************************************************************
      subroutine mgibbs(ntot,subj,m,ist,ifin,pcol,pred,q,zcol,ztz,patt,
     /     nstar,r,y,p,xcol,npatt,rmat,sflag,beta,sigma,psi,b,
     /     xtxinv,wkpp,wkpr,eps,wkrqrq1,wkrqrq2,sig,wkrr1,wkrr2,iter,
     /     wkqr1,wkqr2,wkqrv,nhyp,hyp,delta,iposn,pstfin,betas,
     /     sigmas,psis)
C conventional Gibbs sampler.
      integer ntot,subj(ntot),m,ist(m),ifin(m),pcol,sflag,
     /     q,zcol(q),nstar,patt(ntot),r,p,xcol(p),npatt,rmat(npatt,r),
     /     iter,nhyp,iposn(ntot),pstfin(npatt,2),it
      integer oc(100),mc(100)
      double precision pred(ntot,pcol),ztz(q,q,m),y(ntot,r),
     /     beta(p,r),sigma(r,r),psi(r*q,r*q),b(q,r,m),xtxinv(p,p),
     /     wkpp(p,p),wkpr(p,r),eps(ntot,r),wkrqrq1(r*q,r*q),
     /     wkrqrq2(r*q,r*q),sig(r*q,r*q,m),wkrr1(r,r),wkrr2(r,r),
     /     wkqr1(q,r),wkqr2(q,r),wkqrv(q*r),hyp(nhyp),
     /     delta(ntot,r),betas(p,r,iter),sigmas(r,r,iter),
     /     psis(r*q,r*q,iter)
      double precision wkr(100)
      real junk
      junk=gauss()
      call prelimm(ntot,subj,m,ist,ifin,pcol,pred,q,zcol,ztz,patt,
     /     nstar,p,xcol,xtxinv,wkpp)
      if(sflag.ne.1) then
C get starting values of parameters
         call mimpy(ntot,r,y,patt,npatt,rmat)
         call mkxty(ntot,r,y,pcol,pred,p,xcol,patt,wkpr)
         call mkbeta(p,r,xtxinv,wkpr,beta)
         call mkeps1(ntot,r,y,pcol,pred,p,xcol,beta,eps,patt)
         call mksigma(ntot,r,eps,nstar,sigma,patt)
         call mksig(r,q,m,psi,sigma,ztz,sig,wkrr1,wkrr2,wkrqrq1,
     /     wkrqrq2,1,nhyp,hyp)
         call mkpsi0(r,q,m,psi,sig,wkrqrq1)
      else
C set starting value of y to eps
         do 50 i=1,ntot
            do 40 j=1,r
               y(i,j)=eps(i,j)
 40         continue
 50      continue
      endif
      do 100 it=1,iter
         call mkeps1(ntot,r,y,pcol,pred,p,xcol,beta,eps,patt)
         call mksig(r,q,m,psi,sigma,ztz,sig,wkrr1,wkrr2,wkrqrq1,
     /     wkrqrq2,0,nhyp,hyp)
         call drb(r,q,m,b,wkrr2,eps,pcol,pred,zcol,wkqr1,
     /     wkqr2,ist,ifin,patt,ntot,sig,wkrqrq2,wkqrv)
         call drpsi(r,q,m,psi,wkrqrq1,wkrqrq2,wkqrv,nhyp,hyp,b)
         call mkeps2(ntot,m,r,y,pcol,pred,q,zcol,b,eps,patt,ist,ifin)
         call mkxty(ntot,r,eps,pcol,pred,p,xcol,patt,wkpr)
         call mkbeta(p,r,xtxinv,wkpr,beta)
         call mkeps1(ntot,r,eps,pcol,pred,p,xcol,beta,delta,patt)
         call drsigma(ntot,r,delta,nstar,sigma,patt,nhyp,hyp,wkrr1,
     /     wkrr2,p)
         call drbeta(r,sigma,xtxinv,p,beta,wkrr1,wkpp,wkpr)
         call mkeps1(ntot,r,eps,pcol,pred,p,xcol,beta,delta,patt)
         call dreps(100,100,oc,mc,100,wkr,ntot,iposn,npatt,pstfin,
     /     r,rmat,delta,sigma,wkrr1,wkrr2)
         call mky(ntot,r,pcol,pred,delta,y,p,xcol,q,zcol,beta,m,b,
     /     ist,ifin,npatt,rmat,patt)
         call store(iter,it,p,r,beta,q,psi,sigma,betas,sigmas,psis)
 100  continue
      return
      end
C***********************************************************************
      subroutine storebd(iter,it,p,r,beta,q,psi,sigma,betas,sigmas,psis)
      integer iter,it,p,r,q
      double precision beta(p,r),betas(p,r,iter),sigma(r,r),
     /     sigmas(r,r,iter),psi(q,q,r),psis(q,q,r,iter)
      do 10 j=1,r
         do 5 i=1,p
            betas(i,j,it)=beta(i,j)
 5       continue
 10   continue
      do 21 l=1,r
         do 20 j=1,q
            do 15 i=1,q
               psis(i,j,l,it)=psi(i,j,l)
 15         continue
 20      continue
 21   continue
      do 30 j=1,r
         do 25 i=1,r
            sigmas(i,j,it)=sigma(i,j)
 25      continue
 30   continue
      return
      end
C***********************************************************************
      subroutine store(iter,it,p,r,beta,q,psi,sigma,betas,sigmas,psis)
      integer iter,it,p,r,q
      double precision beta(p,r),betas(p,r,iter),sigma(r,r),
     /     sigmas(r,r,iter),psi(r*q,r*q),psis(r*q,r*q,iter)
      do 10 j=1,r
         do 5 i=1,p
            betas(i,j,it)=beta(i,j)
 5       continue
 10   continue
      do 20 j=1,r*q
         do 15 i=1,r*q
            psis(i,j,it)=psi(i,j)
 15      continue
 20   continue
      do 30 j=1,r
         do 25 i=1,r
            sigmas(i,j,it)=sigma(i,j)
 25      continue
 30   continue
      return
      end
C***********************************************************************
      subroutine mky(ntot,r,pcol,pred,delta,y,p,xcol,q,zcol,beta,m,b,
     /     ist,ifin,npatt,rmat,patt)
C calculates y = delta + X%*%beta  + Z%*%b for missing observations
      integer ntot,r,pcol,p,xcol(p),q,zcol(q),m,ist(m),ifin(m),npatt,
     /     rmat(npatt,r),patt(ntot),s,st,fin
      double precision pred(ntot,pcol),delta(ntot,r),y(ntot,r),
     /     beta(p,r),b(q,r,m),sum
      do 100 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 90 i=st,fin
            if(patt(i).eq.0) then
               do 20 j=1,r
                  sum=dble(0.)
                  do 10 k=1,p
                     sum=sum+pred(i,xcol(k))*beta(k,j)
 10               continue
                  do 15 k=1,q
                     sum=sum+pred(i,zcol(k))*b(k,j,s)
 15               continue
                  y(i,j)=sum+delta(i,j)
 20            continue
            else
               do 85 j=1,r
                  if(rmat(patt(i),j).eq.0) then
                     sum=dble(0.)
                     do 40 k=1,p
                        sum=sum+pred(i,xcol(k))*beta(k,j)
 40                  continue
                     do 50 k=1,q
                        sum=sum+pred(i,zcol(k))*b(k,j,s)
 50                  continue
                     y(i,j)=sum+delta(i,j)
                  endif
 85            continue
            endif
 90      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine dreps(loc,lmc,oc,mc,lwkr,wkr,ntot,iposn,npatt,pstfin,
     /     r,rmat,delta,sigma,wkrr1,wkrr2)
C draws residuals given parameters and random effects
      integer loc,lmc,oc(loc),mc(lmc),lwkr,nmc,noc,ntot,iposn(ntot),
     /     npatt,pstfin(npatt,2),r,rmat(npatt,r),patt(ntot),pt
      double precision wkr(lwkr),delta(ntot,r),sigma(r,r),wkrr1(r,r),
     /     wkrr2(r,r),sum
C take care of completely missing cases
      if(pstfin(1,1).ne.1) then
         do 5 i=1,r
            do 4 j=i,r
               wkrr1(i,j)=sigma(i,j)
 4          continue
 5       continue
         call chfc(r,r,wkrr1)
         do 100 i=1,(pstfin(1,1)-1)
            do 10 j=1,r
               wkr(j)=dble(gauss())
 10         continue
            do 50 j=r,1,-1
               sum=dble(0.)
               do 40 k=1,j
                  sum=sum+wkrr1(k,j)*wkr(k)
 40            continue
               wkr(j)=sum
 50         continue
            do 70 j=1,r
               delta(iposn(i),j)=wkr(j)
 70         continue
 100     continue
      endif
C now take care of partially observed cases
      do 110 i=1,r
         do 105 j=i,r
            wkrr1(i,j)=sigma(i,j)
 105     continue
 110  continue
      do 500 pt=1,npatt
         call swpobs(r,wkrr1,npatt,rmat,pt)
         do 120 i=1,r
            do 119 j=i+1,r
               wkrr1(j,i)=wkrr1(i,j)
 119        continue
 120     continue
         call getmc(r,npatt,rmat,pt,lmc,mc,nmc)
         call getoc(r,npatt,rmat,pt,loc,oc,noc)
         call chsub(r,wkrr1,lmc,mc,nmc,wkrr2)
         do 450 i=pstfin(pt,1),pstfin(pt,2)
            do 150 j=1,nmc
               wkr(j)=dble(gauss())
 150        continue
            do 200 j=nmc,1,-1
               sum=dble(0.)
               do 190 k=1,j
                  sum=sum+wkrr2(k,j)*wkr(k)
 190           continue
               wkr(j)=sum
 200        continue
            do 300 j=1,nmc
               sum=dble(0.)
               do 210 k=1,noc
                  sum=sum+wkrr1(oc(k),mc(j))*delta(iposn(i),oc(k))
 210           continue
               delta(iposn(i),mc(j))=sum+wkr(j)
 300        continue
 450     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine chsub(r,sigma,lmc,mc,nmc,wkrr)
C takes cholesky factor of submatrix of sigma corresponding to mc,
C storing it in the upper-left (nmc x nmc) corner of wkrr
      integer r,lmc,mc(lmc),nmc
      double precision sigma(r,r),wkrr(r,r)
      do 10 i=1,nmc
         do 5 j=i,nmc
            wkrr(i,j)=sigma(mc(i),mc(j))
 5       continue
 10   continue
      call chfc(r,nmc,wkrr)
      return
      end
C***********************************************************************
      subroutine getmc(r,npatt,rmat,pt,lmc,mc,nmc)
      integer r,npatt,rmat(npatt,r),pt,lmc,mc(lmc),nmc,posn
      nmc=0
      posn=0
      do 10 j=1,r
         if(rmat(pt,j).eq.0) then
            nmc=nmc+1
            posn=posn+1
            mc(posn)=j
         endif
 10   continue
      return
      end
C***********************************************************************
      subroutine getoc(r,npatt,rmat,pt,loc,oc,noc)
      integer r,npatt,rmat(npatt,r),pt,loc,oc(loc),noc,posn
      noc=0
      posn=0
      do 10 j=1,r
         if(rmat(pt,j).eq.1) then
            noc=noc+1
            posn=posn+1
            oc(posn)=j
         endif
 10   continue
      return
      end
C***********************************************************************
      subroutine swpobs(r,sigma,npatt,rmat,pt)
C sweeps sigma on positions of rmat(pt,*) containing 1's
      integer r,npatt,rmat(npatt,r),pt,col
      double precision sigma(r,r)
      do 100 col=1,r
         if((rmat(pt,col).eq.1).and.(sigma(col,col).gt.dble(0.))) 
     /        call swp(r,sigma,col)
         if((rmat(pt,col).eq.0).and.(sigma(col,col).le.dble(0.)))
     /        call rsw(r,sigma,col)
 100  continue
      return
      end
C***********************************************************************
      subroutine swp(n,mat,p)
C sweep upper-triangular portion of mat on position p
      integer n,p
      double precision mat(n,n)
      mat(p,p)=-dble(1.)/mat(p,p)
      do 1 i=1,p-1
         mat(i,p)=-mat(i,p)*mat(p,p)
 1    continue
      do 2 j=p+1,n
         mat(p,j)=-mat(p,j)*mat(p,p)
 2    continue
      do 6 i=1,p-1
         do 4 j=i,p-1
            mat(i,j)=mat(i,j)+mat(i,p)*mat(j,p)/mat(p,p)
 4       continue
         do 5 j=p+1,n
            mat(i,j)=mat(i,j)+mat(i,p)*mat(p,j)/mat(p,p)
 5       continue
 6    continue
      do 8 i=p+1,n
         do 7 j=i,n
            mat(i,j)=mat(i,j)+mat(p,i)*mat(p,j)/mat(p,p)
 7       continue
 8    continue
      return
      end
C***********************************************************************
      subroutine rsw(n,mat,p)
C reverse-sweep upper-triangular portion of mat on position p
      integer n,p
      double precision mat(n,n)
      mat(p,p)=-dble(1.)/mat(p,p)
      do 1 i=1,p-1
         mat(i,p)=mat(i,p)*mat(p,p)
 1    continue
      do 2 j=p+1,n
         mat(p,j)=mat(p,j)*mat(p,p)
 2    continue
      do 6 i=1,p-1
         do 4 j=i,p-1
            mat(i,j)=mat(i,j)+mat(i,p)*mat(j,p)/mat(p,p)
 4       continue
         do 5 j=p+1,n
            mat(i,j)=mat(i,j)+mat(i,p)*mat(p,j)/mat(p,p)
 5       continue
 6    continue
      do 8 i=p+1,n
         do 7 j=i,n
            mat(i,j)=mat(i,j)+mat(p,i)*mat(p,j)/mat(p,p)
 7       continue
 8    continue
      return
      end
C***********************************************************************
      subroutine drbeta(r,sigma,xtxinv,p,beta,wkrr1,wkpp,wkpr)
C draws adds normal(0, sigma Otimes xtxinv) noise to beta
      integer r,p
      double precision sigma(r,r),xtxinv(p,p),beta(p,r),wkrr1(r,r),
     /     wkpp(p,p),wkpr(p,r),sum
      do 2 i=1,r
         do 1 j=i,r
            wkrr1(i,j)=sigma(i,j)
 1       continue
 2    continue
      do 5 i=1,p
         do 4 j=i,p
            wkpp(i,j)=xtxinv(i,j)
 4       continue
 5    continue
      call chfc(r,r,wkrr1)
      call chfc(p,p,wkpp)
      do 10 i=1,p
         do 9 j=1,r
            wkpr(i,j)=dble(gauss())
 9       continue
 10   continue
C multiply each column of wkpr by t(wkpp)
      do 100 j=1,r
         do 80 i=p,1,-1
            sum=dble(0.)
            do 70 k=1,i
               sum=sum+wkpp(k,i)*wkpr(k,j)
 70         continue
            wkpr(i,j)=sum
 80      continue
 100  continue
C add noise
      do 200 j=1,r
         do 180 i=1,j
            do 170 k=1,p
               beta(k,j)=beta(k,j)+wkrr1(i,j)*wkpr(k,i)
 170        continue
 180     continue
 200  continue
      return
      end
C***********************************************************************
      subroutine drsigma(ntot,r,delta,nstar,sigma,patt,nhyp,hyp,wkrr1,
     /     wkrr2,p)
C draws sigma
      integer ntot,r,patt(ntot),nstar,nhyp,p
      double precision delta(ntot,r),sigma(r,r),hyp(nhyp),a,
     /     wkrr1(r,r),wkrr2(r,r)
C initialize wkrr1 to Binv
      k=1
      a=hyp(k)
      do 5 i=1,r
         do 4 j=i,r
            k=k+1
            wkrr1(i,j)=hyp(k)
 4       continue
 5    continue
C accumulate sum_i t(delta_i)%*%delta_i
      do 100 i=1,ntot
         if(patt(i).ne.0) then
            do 80 j=1,r
               do 70 k=j,r
                  wkrr1(j,k)=wkrr1(j,k)+delta(i,j)*delta(i,k)
 70            continue
 80         continue
         endif
 100  continue
C draw inverse Wishart
      call chfc(r,r,wkrr1)
      call bfac(r,sngl(a)+float(nstar-p),sigma)
      call bkslv(r,r,sigma)
      do 200 i=1,r
         do 190 j=1,r
            sum=dble(0.)
            do 180 k=1,min(i,j)
               sum=sum+wkrr1(k,i)*sigma(k,j)
 180        continue
            wkrr2(i,j)=sum
 190     continue
 200  continue
      do 300 i=1,r
         do 280 j=i,r
            sum=dble(0.)
            do 250 k=1,r
               sum=sum+wkrr2(i,k)*wkrr2(j,k)
 250        continue
            sigma(i,j)=sum
            if(i.ne.j) sigma(j,i)=sigma(i,j)
 280        continue
 300     continue
      return
      end
C***********************************************************************
      subroutine mkeps2(ntot,m,r,y,pcol,pred,q,zcol,b,eps,patt,ist,ifin)
C calculates eps = y - Z%*% b
      integer ntot,r,patt(ntot),pcol,q,zcol(q),ist(m),ifin(m),s,
     /     st,fin,m
      double precision y(ntot,r),pred(ntot,pcol),b(q,r,m),eps(ntot,r),
     /     sum
      do 100 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 90 i=st,fin
            if(patt(i).ne.0) then
            do 80 j=1,r
               sum=dble(0.)
               do 70 k=1,q
                  sum=sum+pred(i,zcol(k))*b(k,j,s)
 70            continue
               eps(i,j)=y(i,j)-sum   
 80            continue
            endif
 90      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine drpsibd(r,q,m,psi,wkqq1,wkqq2,wkqrv,nhyp,hyp,b)
C Version of drpsi for block-diagonal pan
      integer r,q,m,s,nhyp
      double precision psi(q,q,r),sig(r*q,r*q,m),wkqq1(q,q),
     /     hyp(nhyp),wkqrv(q*r),c,b(q,r,m),sum,wkqq2(q,q)
      do 400 l=1,r
C initialize wkqq1 to appropriate portion of Dinv
         k=1 + r*(r+1)/2 + l
         c=hyp(k)
         k=1 + r*(r+1)/2 + r + (l-1)*q*(q+1)/2
         do 5 i=1,q
            do 4 j=i,q
               k=k+1
               wkqq1(i,j)=hyp(k)
 4          continue
 5       continue
C accumulate BtB into wkrr1
         do 100 s=1,m
            do 15 i=1,q
               wkqrv(i)=b(i,l,s)
 15         continue
            do 40 i=1,q
               do 30 j=i,q
                  wkqq1(i,j)=wkqq1(i,j)+wkqrv(i)*wkqrv(j)
 30            continue
 40         continue
 100     continue
C draw inverse Wishart
         call chfc(q,q,wkqq1)
         call bfac(q,sngl(c)+float(m),wkqq2)
         call bkslv(q,q,wkqq2)
         do 200 i=1,q
            do 190 j=1,q
               sum=dble(0.)
               do 180 k=1,min(i,j)
                  sum=sum+wkqq1(k,i)*wkqq2(k,j)
 180           continue
               psi(i,j,l)=sum
 190        continue
 200     continue
         do 205 i=1,q
            do 204 j=1,q
               wkqq2(i,j)=psi(i,j,l)
 204        continue
 205     continue
         do 300 i=1,q
            do 280 j=i,q
               sum=dble(0.)
               do 250 k=1,q
                  sum=sum+wkqq2(i,k)*wkqq2(j,k)
 250           continue
               psi(i,j,l)=sum
               if(i.ne.j) psi(j,i,l)=psi(i,j,l)
 280        continue
 300     continue
 400  continue
      return
      end
C***********************************************************************
      subroutine drpsi(r,q,m,psi,wkrqrq1,wkrqrq2,wkqrv,nhyp,hyp,b)
C draws new value of psi
      integer r,q,m,s,nhyp
      double precision psi(r*q,r*q),sig(r*q,r*q,m),wkrqrq1(r*q,r*q),
     /     hyp(nhyp),wkqrv(q*r),c,b(q,r,m),sum,wkrqrq2(r*q,r*q)
C initialize wkrqrq1 to Dinv
      k=r*(r+1)/2 + 2
      c=hyp(k)
      do 5 i=1,r*q
         do 4 j=i,r*q
            k=k+1
            wkrqrq1(i,j)=hyp(k)
 4       continue
 5    continue
C put vec(b_i) into wkqrv
      do 100 s=1,m
         ia=0
         do 20 j=1,r
            do 15 i=1,q
               ia=ia+1
               wkqrv(ia)=b(i,j,s)
 15         continue
 20      continue
         do 40 i=1,r*q
            do 30 j=i,r*q
               wkrqrq1(i,j)=wkrqrq1(i,j)+wkqrv(i)*wkqrv(j)
 30         continue
 40      continue
 100  continue
C draw inverse Wishart
      call chfc(r*q,r*q,wkrqrq1)
      call bfac(r*q,sngl(c)+float(m),psi)
      call bkslv(r*q,r*q,psi)
      do 200 i=1,r*q
         do 190 j=1,r*q
            sum=dble(0.)
            do 180 k=1,min(i,j)
               sum=sum+wkrqrq1(k,i)*psi(k,j)
 180        continue
            wkrqrq2(i,j)=sum
 190     continue
 200  continue
      do 300 i=1,r*q
         do 280 j=i,r*q
            sum=dble(0.)
            do 250 k=1,r*q
               sum=sum+wkrqrq2(i,k)*wkrqrq2(j,k)
 250        continue
            psi(i,j)=sum
            if(i.ne.j) psi(j,i)=psi(i,j)
 280        continue
 300     continue
      return
      end
C***********************************************************************
      subroutine drb(r,q,m,b,sigmainv,eps,pcol,pred,zcol,wkqr1,
     /     wkqr2,ist,ifin,patt,ntot,sig,wkrqrq2,wkqrv)
C draws b_i, i=1,...,m. 
      integer r,q,m,s,pcol,zcol(q),ist(m),ifin(m),st,fin,patt(ntot),
     /     ntot
      double precision sig(r*q,r*q,m),eps(ntot,r),pred(ntot,pcol),
     /     sigmainv(r,r),b(q,r,m),wkqr1(q,r),wkqr2(q,r),
     /     wkrqrq1(r*q,r*q),wkrqrq2(r*q,r*q),sum,wkqrv(q*r)
      do 500 s=1,m
         st=ist(s)
         fin=ifin(s)
C put t(Z_i)%*%(jth column of eps_i) into jth column of wkqr1, j=1..r
         do 100 j=1,r
            do 90 i=1,q
               sum=dble(0.)
               do 80 k=st,fin
                  if(patt(k).ne.0) sum=sum+pred(k,zcol(i))*eps(k,j)
 80            continue
               wkqr1(i,j)=sum
 90         continue
 100     continue
C put sum_k=1^r sigmainv(j,k) *t(Z_i)%*% (col k of eps_i) into wkqr2
         do 200 j=1,r
            do 190 i=1,q
               sum=dble(0.)
               do 180 k=1,r
                  sum=sum+sigmainv(j,k)*wkqr1(i,k)
 180           continue
               wkqr2(i,j)=sum
 190        continue
 200     continue
C put sig_i into wkrqrq2 and calculate b_i
         call mmul(r*q,r*q,m,sig,s,wkrqrq2)
         do 220 i=1,r*q
            do 210 j=i+1,r*q
               wkrqrq2(j,i)=wkrqrq2(i,j)
 210        continue
 220     continue
         ia=0
         do 400 i=1,r
            do 380 ii=1,q
               ia=ia+1
               ja=0
               sum=dble(0.)
               do 370 j=1,r
                  do 360 jj=1,q
                     ja=ja+1
                     sum=sum+wkrqrq2(ia,ja)*wkqr2(jj,j)
 360              continue
 370           continue
               b(ii,i,s)=sum
 380        continue
 400     continue
         do 410 i=1,r*q
            wkqrv(i)=dble(gauss())
 410     continue
         ia=0
         do 450 i=1,r
            do 440 ii=1,q
               ia=ia+1
               sum=dble(0.)
               do 430 ja=ia,r*q
                  sum=sum+sig(ia,ja,s)*wkqrv(ja)
 430           continue
               b(ii,i,s)=b(ii,i,s)+sum
 440        continue
 450     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine mkpsi0bd(r,q,m,psi,sig,wkrqrq1)
C calculates initial estimate of psi given sig.
      integer r,q,m,s
      double precision psi(q,q,r),sig(r*q,r*q,m),wkrqrq1(r*q,r*q)
      do 6 l=1,r
         do 5 i=1,q
            do 4 j=i,q
               psi(i,j,l)=dble(0.)
 4          continue
 5       continue
 6    continue
      do 100 s=1,m
         call mmul(r*q,r*q,m,sig,s,wkrqrq1)
         do 30 l=1,r
            do 20 i=1,q
               do 10 j=i,q
                  ii=(l-1)*q + i
                  jj=(l-1)*q + j
                  psi(i,j,l)=psi(i,j,l)+wkrqrq1(ii,jj)
 10            continue
 20         continue
 30      continue
 100  continue
      do 120 l=1,r
         do 110 i=1,q
            do 105 j=i,q
               psi(i,j,l)=psi(i,j,l)/dble(m)
               if(i.ne.j) psi(j,i,l)=psi(i,j,l)
 105        continue
 110     continue
 120  continue
      return
      end
C***********************************************************************
      subroutine mkpsi0(r,q,m,psi,sig,wkrqrq1)
C calculates initial estimate of psi given sig.
      integer r,q,m,s
      double precision psi(r*q,r*q),sig(r*q,r*q,m),wkrqrq1(r*q,r*q)
      do 5 i=1,r*q
         do 4 j=i,r*q
            psi(i,j)=dble(0.)
 4       continue
 5    continue
      do 100 s=1,m
         call mmul(r*q,r*q,m,sig,s,wkrqrq1)
         do 20 i=1,r*q
            do 10 j=i,r*q
               psi(i,j)=psi(i,j)+wkrqrq1(i,j)
 10         continue
 20      continue
 100  continue
      do 110 i=1,r*q
         do 105 j=i,r*q
            psi(i,j)=psi(i,j)/dble(m)
            if(i.ne.j) psi(j,i)=psi(i,j)
 105     continue
 110  continue
      return
      end
C***********************************************************************
      subroutine mksigbd(r,q,m,psi,sigma,ztz,sig,wkrr1,wkrr2,
     /     wkrqrq2,zflag,nhyp,hyp,wkqq1,wkqq2)
C Version of mksig for the block-diagonal pan.
C calculates the square root of
C       sig = inv( inv(psi) + ( inv(sigma) Otimes ztz ) ), i=1,..m
C Use zflag = 1 to set psi equal to its prior guess.
C When finished, wkrr2 contains inv(sigma).
      integer r,q,m,s,zflag,nhyp
      double precision psi(q,q,r),sigma(r,r),ztz(q,q,m),
     /     sig(r*q,r*q,m),wkrr1(r,r),wkrr2(r,r),wkrqrq1(r*q,r*q),
     /     wkrqrq2(r*q,r*q),hyp(nhyp),c,wkqq1(q,q),wkqq2(q,q)
      if(zflag.eq.1) then
         k1= 1 + r*(r+1)/2
         k2= 1 + r*(r+1)/2 + r
         do 5 l=1,r
            k1=k1+1
            c=hyp(k1)
            do 4 i=1,q
               do 3 j=i,q
                  k2=k2+1
                  psi(i,j,l)=hyp(k2)/c
 3             continue
 4          continue
 5       continue
      endif
C invert psi and put into wkrqrq2, using the fact that it's block-diagonal
      do 102 i=1,r*q
         do 101 j=i,r*q
            wkrqrq2(i,j)=dble(0.)
 101     continue
 102  continue
      do 110 l=1,r
         do 105 i=1,q
            do 104 j=i,q
               wkqq1(i,j)=psi(i,j,l)
 104        continue
 105     continue
         call chfc(q,q,wkqq1)
         call bkslv(q,q,wkqq1)
         call mm(q,q,wkqq1,wkqq2)
         do 108 i=1,q
            do 107 j=i,q
               ii=(l-1)*q + i
               jj=(l-1)*q + j
               wkrqrq2(ii,jj)=wkqq2(i,j)
 107        continue
 108     continue
 110  continue
C invert sigma and put into wkrr2
      do 7 i=1,r
         do 6 j=i,r
            wkrr1(i,j)=sigma(i,j)
 6       continue
 7    continue
      call chfc(r,r,wkrr1)
      call bkslv(r,r,wkrr1)
      call mm(r,r,wkrr1,wkrr2)
      do 9 i=1,r
         do 8 j=i+1,r
            wkrr2(j,i)=wkrr2(i,j)
 8       continue
 9    continue
      do 100 s=1,m
C initialize sig(,,s) to inv(sigma) Otimes t(Z_i)%*%Z_i
         do 30 i=1,r
            do 20 j=i,r
               do 15 ii=1,q
                  do 10 jj=1,q
                     ia=(i-1)*q+ii
                     ja=(j-1)*q+jj
                     sig(ia,ja,s)=wkrr2(i,j)*ztz(ii,jj,s)
 10               continue
 15            continue
 20         continue
 30      continue
C add in inv(psi) and take inverse square-root
         do 40 i=1,r*q
            do 35 j=i,r*q
               sig(i,j,s)=sig(i,j,s)+wkrqrq2(i,j)
 35         continue
 40      continue
         call chl(r*q,r*q,m,sig,s)
         call bkslvl(r*q,r*q,m,sig,s)
 100  continue
      return
      end
C***********************************************************************
      subroutine mksig(r,q,m,psi,sigma,ztz,sig,wkrr1,wkrr2,wkrqrq1,
     /     wkrqrq2,zflag,nhyp,hyp)
C calculates the square root of
C       sig = inv( inv(psi) + ( inv(sigma) Otimes ztz ) ), i=1,..m
C Use zflag = 1 to set psi equal to its prior guess.
C When finished, wkrr2 contains inv(sigma).
      integer r,q,m,s,zflag,nhyp
      double precision psi(r*q,r*q),sigma(r,r),ztz(q,q,m),
     /     sig(r*q,r*q,m),wkrr1(r,r),wkrr2(r,r),wkrqrq1(r*q,r*q),
     /     wkrqrq2(r*q,r*q),hyp(nhyp),c
      if(zflag.eq.1) then
         k=r*(r+1)/2 + 2
         c=hyp(k)
         do 4 i=1,r*q
            do 3 j=i,r*q
               k=k+1
               psi(i,j)=hyp(k)/c
 3          continue
 4       continue
      endif
C invert psi and put into wkrqrq2
      do 2 i=1,r*q
         do 1 j=i,r*q
            wkrqrq1(i,j)=psi(i,j)
 1       continue
 2    continue
      call chfc(r*q,r*q,wkrqrq1)
      call bkslv(r*q,r*q,wkrqrq1)
      call mm(r*q,r*q,wkrqrq1,wkrqrq2)
C invert sigma and put into wkrr2
      do 6 i=1,r
         do 5 j=i,r
            wkrr1(i,j)=sigma(i,j)
 5       continue
 6    continue
      call chfc(r,r,wkrr1)
      call bkslv(r,r,wkrr1)
      call mm(r,r,wkrr1,wkrr2)
      do 9 i=1,r
         do 8 j=i+1,r
            wkrr2(j,i)=wkrr2(i,j)
 8       continue
 9    continue
      do 100 s=1,m
C initialize sig(,,s) to inv(sigma) Otimes t(Z_i)%*%Z_i
         do 30 i=1,r
            do 20 j=i,r
               do 15 ii=1,q
                  do 10 jj=1,q
                     ia=(i-1)*q+ii
                     ja=(j-1)*q+jj
                     sig(ia,ja,s)=wkrr2(i,j)*ztz(ii,jj,s)
 10               continue
 15            continue
 20         continue
 30      continue
C add in inv(psi) and take inverse square-root
         do 40 i=1,r*q
            do 35 j=i,r*q
               sig(i,j,s)=sig(i,j,s)+wkrqrq2(i,j)
 35         continue
 40      continue
         call chl(r*q,r*q,m,sig,s)
         call bkslvl(r*q,r*q,m,sig,s)
 100  continue
      return
      end
C***********************************************************************
      subroutine mksigma(ntot,r,eps,nstar,sigma,patt)
C calculates (1/nstar)*sum of t(eps_i)%*%eps_i
      integer ntot,r,patt(ntot),nstar
      double precision eps(ntot,r),sigma(r,r)
      do 2 i=1,r
         do 1 j=i,r
            sigma(i,j)=dble(0.)
 1       continue
 2    continue
      do 100 i=1,ntot
         if(patt(i).ne.0) then
            do 90 j=1,r
               do 80 k=j,r
                  sigma(j,k)=sigma(j,k)+eps(i,j)*eps(i,k)
 80            continue
 90         continue
         endif
 100  continue
      do 110 i=1,r
         do 108 j=i,r
            sigma(i,j)=sigma(i,j)/dble(nstar)
            if(i.ne.j) sigma(j,i)=sigma(i,j)
 108     continue
 110  continue
      return
      end
C***********************************************************************
      subroutine mkeps1(ntot,r,y,pcol,pred,p,xcol,beta,eps,patt)
C calculates eps = y - X%*% beta
      integer ntot,r,patt(ntot),pcol,p,xcol(p)
      double precision y(ntot,r),pred(ntot,pcol),beta(p,r),eps(ntot,r),
     /     sum
      do 100 i=1,ntot
         if(patt(i).ne.0) then
            do 90 j=1,r
               sum=dble(0.)
               do 80 k=1,p
                  sum=sum+pred(i,xcol(k))*beta(k,j)
 80            continue
               eps(i,j)=y(i,j)-sum   
 90         continue
         endif
 100  continue
      continue
      return
      end
C***********************************************************************
      subroutine mkbeta(p,r,xtxinv,xty,beta)
C calculates betahat
      integer p,r
      double precision xtxinv(p,p),xty(p,r),beta(p,r),sum
      do 100 i=1,p
         do 90 j=1,r
            sum=dble(0.)
            do 50 k=1,p
               sum=sum+xtxinv(i,k)*xty(k,j)
 50         continue
            beta(i,j)=sum
 90      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine mkxty(ntot,r,y,pcol,pred,p,xcol,patt,xty)
      integer ntot,r,pcol,p,xcol(p),patt(ntot)
      double precision y(ntot,r),pred(ntot,pcol),xty(p,r),sum
      do 100 i=1,p
         do 80 j=1,r
            sum=dble(0.)
            do 50 k=1,ntot
               if(patt(k).ne.0) then
                  sum=sum+pred(k,xcol(i))*y(k,j)
               endif
 50         continue
            xty(i,j)=sum
 80      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine mimpy(ntot,r,y,patt,npatt,rmat)
C unconditional mean imputation for variables in y
      integer ntot,r,patt(ntot),npatt,rmat(npatt,r),denom,rij
      double precision y(ntot,r),sum,mean
      do 200 j=1,r
         sum=dble(0.)
         denom=0
         do 100 i=1,ntot
            if(patt(i).ne.0) then
               rij=rmat(patt(i),j)
               sum=sum+dble(rij)*y(i,j)
               denom=denom+rij
            endif
 100     continue
         mean=sum/dble(denom)
         do 150 i=1,ntot
            if(patt(i).ne.0) then
               if(rmat(patt(i),j).eq.0) y(i,j)=mean
            endif
 150     continue
 200  continue
      return
      end
C***********************************************************************
        real function rangen(init)
        integer a,p,ix,b15,b16,xhi,xalo,leftflo,fhi,k,init
        save ix
        data a/16807/,b15/32768/,b16/65536/,p/2147483647/
        if(init.ne.0) ix=init
        xhi=ix/b16
        xalo=(ix-xhi*b16)*a
        leftflo=xalo/b16
        fhi=xhi*a+leftflo
        k=fhi/b15
        ix=(((xalo-leftflo*b16)-p)+(fhi-k*b15)*b16)+k
        if (ix.lt.0)ix=ix+p
        rangen=float(ix)*4.656612875E-10
        return
        end
C***********************************************************************
        real function rngs(seed)
C initializes rangen with seed
        integer seed
        rngs=rangen(seed)
        return
        end
C***********************************************************************
        function gamm(a)
C Generates a random gamma(a) deviate. If a>=1, uses the method of 
C Fishman (1976); if 0<a<1, the method of Ahrens (1974)
        real a,u,y,q,e,b,p,u1,lq
        data e/2.718282/
        if(a.ge.1)then
1          continue
           u=rangen(0)
           y=-log(rangen(0))
C           q=(y/exp(y-1))**(a-1)
           lq=(a-1.)*(log(y)-(y-1.))
           q=exp(lq)
           if(u.le.q)then
              gamm=a*y
           else
              goto 1
           endif
        else
2          continue
           u=rangen(0)
           b=(e+a)/e
           p=b*u
           if(p.gt.1) goto 4
           continue
           x=p**(1/a)
           u1=rangen(0)
           if(u1.gt.(e**(-x)))then
              goto 2
           else
              gamm=x
              goto 10
           endif
4          continue
           x=-log((b-p)/a)
           u1=rangen(0)
           if(u1.gt.(x**(a-1)))then
              goto 2
           else
              gamm=x
              goto 10
           endif
        endif
10      continue
        return
        end
C***********************************************************************
	function gauss()
	integer alt
	real next
        save alt,next
	data pi/3.141593/
        if((alt.ne.0).and.(alt.ne.1)) alt=0
	if(alt.eq.0)then
	  u1=rangen(0)
	  u2=rangen(0)
	  gauss=sqrt(-2*log(u1))*cos(2*pi*u2)
	  next=sqrt(-2*log(u1))*sin(2*pi*u2)
	  alt=1
	else
	  gauss=next
	  alt=0
	endif
	return
	end
C***********************************************************************
      subroutine bfac(p,m,b)
C draws upper triangular square-root of a Wishart(m,I) using Bartlett
C decomposition.
      integer p
      real m,jnk
      double precision b(p,p)
      jnk=gauss()
      do 10 j=1,p
         b(j,j)=dble(sqrt(2.*gamm((m-float(j)+1.)/2.)))
 10   continue
      do 30 j=1,p-1
         do 20 k=j+1,p
            b(j,k)=gauss()
 20      continue
 30   continue
      return
      end
C***********************************************************************
C The following are used only by ecme().
C***********************************************************************
      subroutine prelim(ntot,subj,m,ist,ifin,occ,nmax,vmax,vh,vi,pcol,
     /     pred,q,zcol,ztv,sig0,iflag)
C Preliminary manipulations for ECME. After execution, 
C     ist  = vector of starting positions for subject i, i=1,...,m
C     ifin = vector of finishing positions for subject i, i=1,...,m
C     vh   = inverse sqrt of V_i, i=1,...,m
C     vi   = inverse of V_i, i=1,...,m
C     ztv  = t(z_i)%*% inverse of V_i, i=1,...,m
C     sig0 = t(z_i)%*%inv(V_i)%*%(z_i), i=1,...,m
C Note: if V_i's are all identity, then let vmax, vh, vi be blank arrays
C and set iflag = 1.
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,pcol,
     /     q,zcol(q),iflag
      double precision vmax(nmax,nmax),vh(nmax,nmax,m),vi(nmax,nmax,m),
     /     pred(ntot,pcol),ztv(q,nmax,m),sig0(q,q,m)
      call istfin(ntot,subj,m,ist,ifin)
      if(iflag.ne.1) then
         call mkv(m,nmax,vmax,ntot,occ,ist,ifin,vh)
         call chv(nmax,m,vh,ntot,occ,ist,ifin)
         call bkv(nmax,m,vh,ntot,occ,ist,ifin)
         call mmulv(nmax,m,vh,vi,ntot,occ,ist,ifin)
      endif
      call mmu(ntot,pcol,pred,q,zcol,nmax,m,vh,occ,ist,ifin,ztv,iflag)
      call mmtm(q,nmax,m,ztv,ntot,occ,ist,ifin,sig0)
      if(iflag.ne.1) then
         call mml(ntot,q,nmax,m,vh,occ,ist,ifin,ztv)
      endif
      return
      end
C***********************************************************************
      subroutine stval(ntot,m,ist,ifin,occ,nmax,vi,vh,pcol,
     /     pred,q,zcol,ztv,sig0,iflag,sig,psi,sigma2,p,xcol,beta,wkq1,
     /     wkq2,wkq3,y,delta,b,wk,w,xtw,xtwx,xtwy,xtwxinv)
C Get default starting values for ECME
      integer ntot,m,ist(m),ifin(m),occ(ntot),nmax,pcol,
     /     q,zcol(q),iflag,p,xcol(p),s
      double precision vi(nmax,nmax,m),pred(ntot,pcol),
     /     ztv(q,nmax,m),sig0(q,q,m),sig(q,q,m),psi(q,q),sigma2,
     /     beta(p),wkq1(q,q),wkq2(q,q),y(ntot),delta(ntot),b(q,m),
     /     wk(q,nmax,m),vh(nmax,nmax,m),wkq3(q,q),w(nmax,nmax,m),
     /     xtw(p,nmax),xtwx(p,p),xtwy(p),xtwxinv(p,p)
      call gls(ntot,m,ist,ifin,occ,nmax,vi,vh,pcol,
     /     pred,q,iflag,sig,sigma2,p,xcol,beta,
     /     y,delta,w,xtw,xtwx,xtwy,xtwxinv)
C get psi
      do 310 i=1,q
         do 305 j=i,q
            wkq1(i,j)=dble(0.)
 305     continue
 310  continue
      do 500 s=1,m
         do 400 i=1,q
            do 350 j=i,q
               wkq1(i,j)=wkq1(i,j)+sig0(i,j,s)
 350        continue
 400     continue
 500  continue
      call chfc(q,q,wkq1)
      call bkslv(q,q,wkq1)
      call mm(q,q,wkq1,psi)
      do 600 i=1,q
         do 550 j=i,q
            psi(i,j)=psi(i,j)*sigma2*dble(m)
            if(i.ne.j) then
               psi(j,i)=psi(i,j)
            endif
 550     continue
 600  continue
      return
      end
C***********************************************************************
      subroutine nopsi(ntot,m,ist,ifin,occ,nmax,vi,vh,pcol,
     /     pred,q,zcol,ztv,sig0,iflag,sig,psi,sigma2,p,xcol,beta,wkq1,
     /     wkq2,wkq3,y,delta,b,wk,w,xtw,xtwx,xtwy,xtwxinv,ll)
C gets GLS estimates for beta and sigma2, assuming psi=0, and evaluates
C loglikelihood
      integer ntot,m,ist(m),ifin(m),occ(ntot),nmax,pcol,
     /     q,zcol(q),iflag,p,xcol(p)
      double precision vi(nmax,nmax,m),pred(ntot,pcol),
     /     ztv(q,nmax,m),sig0(q,q,m),sig(q,q,m),psi(q,q),sigma2,
     /     beta(p),wkq1(q,q),wkq2(q,q),y(ntot),delta(ntot),b(q,m),
     /     wk(q,nmax,m),vh(nmax,nmax,m),wkq3(q,q),w(nmax,nmax,m),
     /     xtw(p,nmax),xtwx(p,p),xtwy(p),xtwxinv(p,p),ll
      call gls(ntot,m,ist,ifin,occ,nmax,vi,vh,pcol,
     /     pred,q,iflag,sig,sigma2,p,xcol,beta,
     /     y,delta,w,xtw,xtwx,xtwy,xtwxinv)
      call mkll2(nmax,m,w,ntot,delta,occ,ist,ifin,ll)
      return
      end
C***********************************************************************
      subroutine mkll2(nmax,m,w,ntot,delta,occ,ist,ifin,ll)
C evaluates loglikelihood
      integer nmax,m,ntot,occ(ntot),ist(m),ifin(m),s,st,fin
      double precision w(nmax,nmax,m),delta(ntot),ll,sum
      ll=dble(0.)
      call chv(nmax,m,w,ntot,occ,ist,ifin)
      do 900 s=1,m
         st=ist(s)
         fin=ifin(s)
         sum=dble(0.)
         do 100 i=st,fin
            sum=sum+dlog(w(occ(i),occ(i),s))
 100     continue
         ll=ll+dble(2.)*sum
         do 300 i=st,fin
            sum=dble(0.)
            do 200 j=i,fin
               sum=sum+w(occ(i),occ(j),s)*delta(j)
 200        continue
            ll=ll-sum**2
 300     continue
 900  continue
      ll=ll/dble(2.)
      return
      end
C***********************************************************************
      subroutine gls(ntot,m,ist,ifin,occ,nmax,vi,vh,pcol,
     /     pred,q,iflag,sig,sigma2,p,xcol,beta,
     /     y,delta,w,xtw,xtwx,xtwy,xtwxinv)
C gets GLS estimates for beta and sigma2, assuming psi=0
      integer ntot,m,ist(m),ifin(m),occ(ntot),nmax,pcol,
     /     q,zcol(q),iflag,p,xcol(p),st,fin,s
      double precision vi(nmax,nmax,m),pred(ntot,pcol),
     /     ztv(q,nmax,m),sig0(q,q,m),sig(q,q,m),psi(q,q),sigma2,
     /     beta(p),wkq1(q,q),wkq2(q,q),y(ntot),delta(ntot),b(q,m),
     /     wk(q,nmax,m),vh(nmax,nmax,m),wkq3(q,q),w(nmax,nmax,m),
     /     xtw(p,nmax),xtwx(p,p),xtwy(p),xtwxinv(p,p)
C initialize W_i to inv(V_i) and  get beta
      do 10 i=1,p
         xtwy(i)=dble(0.)
         do 5 j=i,p
            xtwx(i,j)=dble(0.)
 5       continue
 10   continue
      do 100 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 12 i=st,fin
            do 11 j=i,fin
               if(iflag.ne.1) then
                  w(occ(i),occ(j),s)=vi(occ(i),occ(j),s)
               else
                  if(i.eq.j) then
                     w(occ(i),occ(j),s)=dble(1.)
                  else
                     w(occ(i),occ(j),s)=dble(0.)
                  endif
               endif
 11         continue
 12      continue
         call mkxtw(ntot,pcol,pred,p,xcol,occ,st,fin,nmax,w,xtw,s,m)
         call mkxtwx(ntot,pcol,pred,p,xcol,occ,st,fin,nmax,xtw,xtwx)
         call mkxtwy(ntot,p,occ,st,fin,nmax,xtw,y,xtwy)
 100  continue
      call chfc(p,p,xtwx)
      call bkslv(p,p,xtwx)
      call mm(p,p,xtwx,xtwxinv)
      do 130 i=1,p
         sum=dble(0.)
         do 110 j=1,i
            sum=sum+xtwxinv(j,i)*xtwy(j)
 110     continue
         do 120 j=i+1,p
            sum=sum+xtwxinv(i,j)*xtwy(j)
 120     continue
         beta(i)=sum
 130  continue
C get sigma2
      call mkdel(ntot,pcol,pred,p,xcol,y,beta,delta)
      sigma2=dble(0.)
      do 300 s=1,m
         st=ist(s)
         fin=ifin(s)
         if(iflag.ne.1) then
            do 200 i=st,fin
               sum1=dble(0.)
               do 190 j=st,i
                  sum1=sum1+delta(j)*vh(occ(j),occ(i),s)
 190            continue
               sigma2=sigma2+sum1**2
 200        continue
         else
            do 210 i=st,fin
               sigma2=sigma2+delta(i)**2
 210        continue
         endif
 300  continue
      sigma2=sigma2/dble(ntot)
      return
      end
C***********************************************************************
      subroutine mkb(q,nmax,m,wk,ntot,delta,b,occ,ist,ifin)
C calculates bwig from wk and delta
      integer q,nmax,m,ntot,occ(ntot),ist(m),ifin(m),s,st,fin
      double precision wk(q,nmax,m),delta(ntot),b(q,m),sum
      do 500 s=1,m
         do 400 i=1,q
            st=ist(s)
            fin=ifin(s)
            sum=dble(0.)
            do 300 j=st,fin
               sum=sum+wk(i,occ(j),s)*delta(j)
 300           continue
            b(i,s)=sum
 400     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine ecme3(ntot,subj,m,ist,ifin,occ,nmax,vi,vh,pcol,
     /     pred,q,zcol,ztv,sig0,iflag,sig,psi,sigma2,p,xcol,beta,wkq1,
     /     wkq2,wkq3,y,delta,b,wk,w,xtw,xtwx,xtwy,xtwxinv,llk,vmax,
     /     sflag,eps,obeta,opsi,maxits,iter,cvgd)
C set sflag=1 if starting values supplied
      integer ntot,m,ist(m),ifin(m),occ(ntot),nmax,pcol,
     /     q,zcol(q),iflag,p,xcol(p),sflag,c1,c2,c3,maxits,
     /     iter,cvgd,subj(ntot)
      double precision vi(nmax,nmax,m),pred(ntot,pcol),
     /     ztv(q,nmax,m),sig0(q,q,m),sig(q,q,m),psi(q,q),sigma2,
     /     beta(p),wkq1(q,q),wkq2(q,q),y(ntot),delta(ntot),b(q,m),
     /     wk(q,nmax,m),vh(nmax,nmax,m),wkq3(q,q),w(nmax,nmax,m),
     /     xtw(p,nmax),xtwx(p,p),xtwy(p),xtwxinv(p,p),ll,llk(maxits),
     /     ldxi,ldsig,vmax(nmax,nmax),eps,obeta(p),opsi(q,q),osigma2
      call prelim(ntot,subj,m,ist,ifin,occ,nmax,vmax,vh,vi,pcol,
     /     pred,q,zcol,ztv,sig0,iflag)
      if(sflag.ne.1) then
         call stval(ntot,m,ist,ifin,occ,nmax,vi,vh,pcol,
     /        pred,q,zcol,ztv,sig0,iflag,sig,psi,sigma2,p,xcol,beta,
     /        wkq1,wkq2,wkq3,y,delta,b,wk,w,xtw,xtwx,xtwy,xtwxinv)
      else
         call mkdel(ntot,pcol,pred,p,xcol,y,beta,delta)
      endif
      iter=0
      cvgd=0
 1    continue
         iter=iter+1
         do 2 i=1,p
            obeta(i)=beta(i)
 2       continue
         do 10 i=1,q
            do 5 j=i,q
               opsi(i,j)=psi(i,j)
               wkq3(i,j)=psi(i,j)/sigma2
 5          continue
 10      continue
         osigma2=sigma2
         call mksig3(q,wkq3,m,sig0,sig,wkq1,wkq2,ldxi,ldsig)
         call mkwk3(q,m,sig,nmax,ztv,wk,ntot,occ,ist,ifin)
         call mkb(q,nmax,m,wk,ntot,delta,b,occ,ist,ifin)
         call mkxi(q,m,b,sig,wkq3,sigma2)
         do 20 i=1,q
            do 15 j=i,q
               psi(i,j)=wkq3(i,j)*sigma2
               if(i.ne.j)  psi(j,i)=psi(i,j)
 15         continue
 20      continue
         call mkbeta3(q,nmax,m,wk,ztv,vi,w,ntot,occ,ist,ifin,
     /        pcol,pred,p,xcol,y,xtw,xtwx,xtwy,xtwxinv,beta,iflag)
         call mkdel(ntot,pcol,pred,p,xcol,y,beta,delta)
         call mksig23(ntot,delta,m,sigma2,nmax,occ,ist,ifin,w)
         ll=-dble(.5)*dble(ntot)*dlog(sigma2)+dble(m)*ldxi+ldsig
         ll=ll-dble(.5)*dble(ntot)
         llk(iter)=ll
         c1=0
         do 30 i=1,p
            if(dabs(beta(i)-obeta(i)).gt.(eps*dabs(obeta(i)))) c1=1
 30      continue
         c2=0
         do 40 i=1,q
            do 35 j=i,q
               if(dabs(psi(i,j)-opsi(i,j)).gt.(eps*dabs(opsi(i,j))))
     /              c2=1
 35         continue
 40      continue
         c3=0
         if(dabs(sigma2-osigma2).gt.(eps*dabs(osigma2))) c3=1
         if((c1.eq.0).and.(c2.eq.0).and.(c3.eq.0)) cvgd=1
         if((cvgd.eq.0).and.(iter.lt.maxits)) goto 1
      return
      end
C***********************************************************************
      subroutine mkw3(q,nmax,m,wk,ztv,vi,s,w,ntot,occ,st,fin,iflag)
C for ECME-3
C makes (upper triangular part of) weight matrix w for subject s
      integer q,nmax,m,s,ntot,occ(ntot),st,fin,iflag
      double precision wk(q,nmax,m),ztv(q,nmax,m),vi(nmax,nmax,m),
     /     w(nmax,nmax,m),sum
      if(iflag.ne.1) then
         do 300 i=st,fin
            do 200 j=i,fin
               sum=dble(0.)
               do 100 k=1,q
                  sum=sum+ztv(k,occ(i),s)*wk(k,occ(j),s)
 100           continue
               w(occ(i),occ(j),s)=(vi(occ(i),occ(j),s)-sum)
 200        continue
 300     continue
      else
         do 600 i=st,fin
            do 500 j=i,fin
               sum=dble(0.)
               do 400 k=1,q
                  sum=sum+ztv(k,occ(i),s)*wk(k,occ(j),s)
 400           continue
               if(i.eq.j) then
                  w(occ(i),occ(j),s)=(dble(1.)-sum)
               else
                  w(occ(i),occ(j),s)=-sum
               endif
 500        continue
 600     continue
      endif
      return
      end
C***********************************************************************
      subroutine mkll(nmax,m,w,ntot,delta,occ,ist,ifin,ll,ldpsi,ldsig,
     /     sigma2)
C Calculates loglikelihood value given w, delta, sigma2. Assumes that 
C  -.5*logdet(psi) is in ldpsi, and +.5*sum(logdet(sig_i)) is in ldsig.
      integer nmax,m,ntot,occ(ntot),ist(m),ifin(m),s,st,fin
      double precision w(nmax,nmax,m),delta(ntot),ll,sum,ldpsi,ldsig,
     /     sigma2
      ll=dble(0.)
      do 900 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 300 i=st,fin
            sum=dble(0.)
            do 100 j=st,i
               sum=sum+delta(j)*w(occ(j),occ(i),s)
 100        continue
            do 200 j=i+1,fin
               sum=sum+delta(j)*w(occ(i),occ(j),s)
 200        continue
            ll=ll+sum*delta(i)
 300     continue
 900  continue
      ll=-dble(.5)*ll-dble(.5)*dble(ntot)*dlog(sigma2)
      ll=ll+dble(m)*ldpsi+ldsig
      return
      end
C***********************************************************************
      subroutine mkbeta3(q,nmax,m,wk,ztv,vi,w,ntot,occ,ist,ifin,
     /     pcol,pred,p,xcol,y,xtw,xtwx,xtwy,xtwxinv,beta,iflag)
C for ECME-3
      integer q,nmax,m,s,ntot,occ(ntot),ist(m),ifin(m),st,fin,p,pcol,
     /     xcol(p),iflag
      double precision wk(q,nmax,m),ztv(q,nmax,m),vi(nmax,nmax,m),
     /     w(nmax,nmax,m),pred(ntot,pcol),xtw(p,nmax),
     /     xtwx(p,p),xtwy(p),beta(p),y(ntot),xtwxinv(p,p),sum
      do 10 i=1,p
         xtwy(i)=dble(0.)
         do 5 j=i,p
            xtwx(i,j)=dble(0.)
 5       continue
 10   continue
      do 100 s=1,m
         st=ist(s)
         fin=ifin(s)
         call mkw3(q,nmax,m,wk,ztv,vi,s,w,ntot,occ,st,fin,iflag)
         call mkxtw(ntot,pcol,pred,p,xcol,occ,st,fin,nmax,w,xtw,s,m)
         call mkxtwx(ntot,pcol,pred,p,xcol,occ,st,fin,nmax,xtw,xtwx)
         call mkxtwy(ntot,p,occ,st,fin,nmax,xtw,y,xtwy)
 100  continue
      call chfc(p,p,xtwx)
      call bkslv(p,p,xtwx)
      call mm(p,p,xtwx,xtwxinv)
      do 130 i=1,p
         sum=dble(0.)
         do 110 j=1,i
            sum=sum+xtwxinv(j,i)*xtwy(j)
 110     continue
         do 120 j=i+1,p
            sum=sum+xtwxinv(i,j)*xtwy(j)
 120     continue
         beta(i)=sum
 130  continue
      return
      end
C***********************************************************************
      subroutine mkxtwy(ntot,p,occ,st,fin,nmax,xtw,y,xtwy)
C increments xtwy for subject s
      integer ntot,p,xcol(p),occ(ntot),st,fin,nmax
      double precision xtw(p,nmax),y(ntot),xtwy(p),sum
      do 200 i=1,p
         sum=dble(0.)
         do 100 j=st,fin
            sum=sum+xtw(i,occ(j))*y(j)
 100     continue
         xtwy(i)=xtwy(i)+sum
 200  continue
      return
      end
C***********************************************************************
      subroutine mkxtwx(ntot,pcol,pred,p,xcol,occ,st,fin,nmax,xtw,xtwx)
C increments xtwx for subject s
      integer ntot,pcol,p,xcol(p),occ(ntot),st,fin,nmax
      double precision pred(ntot,pcol),xtw(p,nmax),xtwx(p,p),sum
      do 300 i=1,p
         do 200 j=i,p
            sum=dble(0.)
            do 100 k=st,fin
               sum=sum+xtw(i,occ(k))*pred(k,xcol(j))
 100        continue
            xtwx(i,j)=xtwx(i,j)+sum
 200     continue
 300  continue
      return
      end
C***********************************************************************
      subroutine mkxtw(ntot,pcol,pred,p,xcol,occ,st,fin,nmax,w,xtw,s,m)
      integer ntot,pcol,p,xcol(p),occ(ntot),st,fin,nmax,s,m
      double precision pred(ntot,pcol),w(nmax,nmax,m),xtw(p,nmax),sum
      do 500 i=1,p
         do 400 j=st,fin
            sum=dble(0.)
            do 200 k=st,j
               sum=sum+pred(k,xcol(i))*w(occ(k),occ(j),s)
 200        continue
            do 300 k=j+1,fin
               sum=sum+pred(k,xcol(i))*w(occ(j),occ(k),s)
 300        continue
            xtw(i,occ(j))=sum
 400     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine mkxi(q,m,b,sig,xi,sigma2)
C calculates new estimate of xi from b, sig, and sigma2
      integer q,m,s
      double precision b(q,m),sig(q,q,m),xi(q,q),sigma2
      do 10 i=1,q
         do 5 j=i,q
            xi(i,j)=dble(0.)
 5       continue
 10   continue
      do 100 s=1,m
         do 90 i=1,q
            do 80 j=i,q
               xi(i,j)=xi(i,j)+sig(i,j,s)+b(i,s)*b(j,s)/sigma2
 80         continue
 90      continue
 100  continue
      do 110 i=1,q
         do 105 j=i,q
            xi(i,j)=xi(i,j)/dble(m)
 105     continue
 110  continue
      return
      end
C***********************************************************************
      subroutine mkwk3(q,m,sig,nmax,ztv,wk,ntot,occ,ist,ifin)
C Version of mkwk for ECME-3
      integer q,m,nmax,occ(ntot),ist(m),ifin(m),st,fin,s
      double precision sig(q,q,m),ztv(q,nmax,m),wk(q,nmax,m),sum
      do 400 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 300 i=1,q
            do 200 j=st,fin
               sum=dble(0.)
               do 50 k=1,i-1
                  sum=sum+sig(k,i,s)*ztv(k,occ(j),s)
 50            continue
               do 100 k=i,q
                  sum=sum+sig(i,k,s)*ztv(k,occ(j),s)
 100           continue
               wk(i,occ(j),s)=sum
 200        continue
 300     continue
 400  continue
      return
      end
C***********************************************************************
      subroutine mkv(m,nmax,vmax,ntot,occ,ist,ifin,v)
C makes array of submatrices of vmax
      integer m,nmax,ntot,occ(ntot),ist(m),ifin(m),s
      double precision vmax(nmax,nmax),v(nmax,nmax,m)
      do 100 s=1,m
         do 90 i=ist(s),ifin(s)
            do 80 j=i,ifin(s)
               v(occ(i),occ(j),s)=vmax(occ(i),occ(j))
 80         continue
 90      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine chv(nmax,m,v,ntot,occ,ist,ifin)
C overwrites v with its cholesky factors
      integer nmax,m,ntot,occ(ntot),ist(m),ifin(m),s
      double precision v(nmax,nmax,m),sum
      do 100 s=1,m
         do 50 i=ist(s),ifin(s)
            sum=dble(0.)
            do 10 k=ist(s),i-1
               sum=sum+v(occ(k),occ(i),s)**2
 10         continue
            v(occ(i),occ(i),s)=dsqrt(v(occ(i),occ(i),s)-sum)
            do 40 j=i+1,ifin(s)
               sum=dble(0.)
               do 30 k=ist(s),i-1
                  sum=sum+v(occ(k),occ(i),s)*v(occ(k),occ(j),s)
 30            continue
               v(occ(i),occ(j),s)=(v(occ(i),occ(j),s)-sum)/
     /              v(occ(i),occ(i),s)
 40         continue
 50      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine bkv(nmax,m,v,ntot,occ,ist,ifin)
C inverts upper-triangular portions of v by backsolve
      integer nmax,m,ntot,occ(ntot),ist(m),ifin(m),s,st,fin
      double precision v(nmax,nmax,m),sum
      do 100 s=1,m
         st=ist(s)
         fin=ifin(s)
         v(occ(st),occ(st),s)=dble(1.)/v(occ(st),occ(st),s)
         do 10 k=st+1,fin
            v(occ(k),occ(k),s)=dble(1.)/v(occ(k),occ(k),s)
            do 5 j=st,k-1
               sum=dble(0.)
               do 3 i=j,k-1
                  sum=sum+v(occ(j),occ(i),s)*v(occ(i),occ(k),s)
 3             continue
               v(occ(j),occ(k),s)=-sum*v(occ(k),occ(k),s)
 5          continue
 10      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine mmu(ntot,pcol,pred,q,zcol,nmax,m,v,occ,ist,ifin,ztv,
     /     iflag)
C multiplies z^t by upper-triangular v, stores result in ztv
      integer ntot,pcol,q,zcol(q),nmax,m,occ(ntot),ist(m),ifin(m),s,
     /     st,fin,iflag
      double precision pred(ntot,pcol),v(nmax,nmax,m),ztv(q,nmax,m),
     /     sum
      if(iflag.ne.1) then
         do 300 s=1,m
            st=ist(s)
            fin=ifin(s)
            do 200 i=1,q
               do 100 j=st,fin
                  sum=dble(0.)
                  do 50 k=st,j
                     sum=sum+pred(k,zcol(i))*v(occ(k),occ(j),s)
 50               continue
                  ztv(i,occ(j),s)=sum
 100           continue
 200        continue
 300     continue
      else
         do 600 s=1,m
            st=ist(s)
            fin=ifin(s)
            do 500 i=1,q
               do 400 j=st,fin
                  ztv(i,occ(j),s)=pred(j,zcol(i))
 400           continue
 500        continue
 600     continue
      endif
      return
      end
C***********************************************************************
      subroutine mml(ntot,q,nmax,m,v,occ,ist,ifin,ztv)
C multiplies ztv by lower-triangular v, stores result in ztv
      integer ntot,q,nmax,m,occ(ntot),ist(m),ifin(m),s,st,fin
      double precision v(nmax,nmax,m),ztv(q,nmax,m),sum
      do 300 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 200 i=1,q
            do 100 j=st,fin
               sum=dble(0.)
               do 50 k=j,fin
                  sum=sum+ztv(i,occ(k),s)*v(occ(j),occ(k),s)
 50            continue
               ztv(i,occ(j),s)=sum
 100        continue
 200     continue
 300  continue
      return
      end
C***********************************************************************
      subroutine mmulv(nmax,m,vh,vi,ntot,occ,ist,ifin)
C calculates upper tri part of vi = vh%*%t(vh), where vh is upper tri
C(works in layers s=1,...,m)
      integer nmax,m,ntot,occ(ntot),ist(m),ifin(m),s,st,fin
      double precision vh(nmax,nmax,m),vi(nmax,nmax,m),sum
      do 300 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 200 i=st,fin
            do 100 j=i,fin
               sum=dble(0.)
               do 50 k=max(i,j),fin
                  sum=sum+vh(occ(i),occ(k),s)*vh(occ(j),occ(k),s)
 50            continue
               vi(occ(i),occ(j),s)=sum
 100        continue
 200        continue
 300  continue
      return
      end
C***********************************************************************
      subroutine mmtm(q,nmax,m,ztv,ntot,occ,ist,ifin,sig0)
C multiply ztv by its transpose, store result (upper-tri part) in sig0
      integer q,nmax,m,ntot,occ(ntot),ist(m),ifin(m),s
      double precision ztv(q,nmax,m),sig0(q,q,m),sum
      do 500 s=1,m
         do 400 i=1,q
            do 300 j=i,q
               sum=dble(0.)
               do 200 k=ist(s),ifin(s)
                  sum=sum+ztv(i,occ(k),s)*ztv(j,occ(k),s)
 200           continue
               sig0(i,j,s)=sum
 300        continue
 400     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine mksig23(ntot,delta,m,sigma2,nmax,occ,ist,ifin,w)
C for ECME-3
      integer ntot,m,s,nmax,occ(ntot),ist(m),ifin(m),st,fin
      double precision delta(ntot),sigma2,sum,vh(nmax,nmax,m),
     /     w(nmax,nmax,m)
      sigma2=dble(0.)
      do 900 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 300 i=st,fin
            sum=dble(0.)
            do 100 j=st,i
               sum=sum+delta(j)*w(occ(j),occ(i),s)
 100        continue
            do 200 j=i+1,fin
               sum=sum+delta(j)*w(occ(i),occ(j),s)
 200        continue
            sigma2=sigma2+sum*delta(i)
 300     continue
 900  continue
      sigma2=sigma2/dble(ntot)
      return
      end
C***********************************************************************
      subroutine mkdel(ntot,pcol,pred,p,xcol,y,beta,delta)
      integer ntot,p,xcol(p),pcol
      double precision pred(ntot,pcol),y(ntot),beta(p),delta(ntot),sum
      do 100 i=1,ntot
         sum=dble(0.)
         do 70 j=1,p
            sum=sum+pred(i,xcol(j))*beta(j)
 70      continue
         delta(i)=y(i)-sum
 100  continue
      return
      end
C***********************************************************************
      subroutine mksig3(q,xi,m,sig0,sig,wkq1,wkq2,ldxi,ldsig)
C Version of mksig for ECME-3. Given xi=psi/sigma2.
C After execution, xi contains the inverse cholesky factor
C of the original xi, and wkq1 contains the inverse of the original
C xi (upper-tri part only).
C Also calculates ldxi = -.5*logdet(xi) and
C                 ldsig = +.5*sum( logdet(sig_i) ), i=1,...,m
      integer q,m,s
      double precision xi(q,q),sig0(q,q,m),sig(q,q,m),
     /     wkq1(q,q),wkq2(q,q),ldxi,ldsig
C invert xi 
      call chfc(q,q,xi)
      call bkslv(q,q,xi)
      ldxi=dble(0.)
      do 1 i=1,q
         ldxi=ldxi+dlog(xi(i,i))
 1    continue
      call mm(q,q,xi,wkq1)
      ldsig=dble(0.)
      do 500 s=1,m
         do 100 i=1,q
            do 80 j=i,q
               sig(i,j,s)=wkq1(i,j)+sig0(i,j,s)
 80         continue
 100     continue
         call chl(q,q,m,sig,s)
         call bkslvl(q,q,m,sig,s)
         do 110 i=1,q
            ldsig=ldsig+dlog(sig(i,i,s))
 110     continue
         call mmul(q,q,m,sig,s,wkq2)
         do 200 i=1,q
            do 150 j=i,q
               sig(i,j,s)=wkq2(i,j)
 150        continue
 200     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine bdiag(q,m,sig)
C fills in elements of each layer of sig below the diagonal
      integer q,m,s
      double precision sig(q,q,m)
      do 10 s=1,m
         do 5 i=1,q
            do 1 j=1,i-1
               sig(i,j,s)=sig(j,i,s)
 1          continue
 5       continue
 10   continue
      return
      end
C***********************************************************************
