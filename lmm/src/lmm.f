C***********************************************************************
C The following subroutines are used by both pan() and ecme().
C***********************************************************************
      subroutine istfin(ntot,subj,m,ist,ifin)
C creates vectors of starting and finishing positions for subjects
      integer ntot,subj(ntot),m,ist(m),ifin(m),scur,i,icur
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
      subroutine bkslv(p,pw,s)
C inverts an upper triangular matrix
      integer p,pw,j,k,i
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
      integer p,pw,i,j,k
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
      integer p,pw,m,l,i,j,k
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
      integer p,pw,m,l,i,j,k
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
C The following subroutines are used by mgibbs().
C***********************************************************************
        real function rangen(init)
        integer a,p,ix,b15,b16,xhi,xalo,leftflo,fhi,k,init
        data a/16807/,b15/32768/,b16/65536/,p/2147483647/
        save ix
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
        real tmp,rangen
        rngs=rangen(seed)
        return
        end
C***********************************************************************
        function gamm(a)
C Generates a random gamma(a) deviate. If a>=1, uses the method of 
C Fishman (1976); if 0<a<1, the method of Ahrens (1974)
        real a,u,y,q,e,b,p,u1,lq,x
        real rangen
        data e/2.718282/
        if(a.ge.1)then
1          continue
           u=rangen(0)
           y=-log(rangen(0))
           lq=(a-1.)*(log(y)-(y-1.))
           q=exp(lq)
C           q=(y/exp(y-1))**(a-1)
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
3          continue
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
	real next,pi,u1,u2,gauss,rangen
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
      subroutine mku2(q,xi,m,ztvinvz,u,wkqq1,wkqq2,ldxi,ldu,err,sqrtu)
C Like mku, except that it also saves the square roots of U_i, 
C i=1,...,m in sqrtu.
C Calculates U_i, i=1,...,m from xi.
C After execution, wkqq1 contains the inverse of the original
C xi (upper-tri part only).
C Also calculates ldxi = -.5*logdet(xi) and
C                 ldu = +.5*sum( logdet(U_i) ), i=1,...,m
C Sets err=1 if the input value of xi is not positive definite, or if
C   the value of any (t(Zi)%*%inv(Vi)%*%Zi + inv(xi)) is not pos. def.
      integer q,m,s,err,i,j
      double precision xi(q,q),ztvinvz(q,q,m),u(q,q,m),
     /     wkqq1(q,q),wkqq2(q,q),ldxi,ldu,sqrtu(q,q,m)
      err=0
C *** put the inverse of xi into wkqq1
      do 2 i=1,q
         do 1 j=i,q
            wkqq2(i,j)=xi(i,j)
 1       continue
 2    continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) goto 999
      call bkslv(q,q,wkqq2)
      ldxi=dble(0.)
      do 10 i=1,q
         ldxi=ldxi+dlog(wkqq2(i,i))
 10   continue
      call mm(q,q,wkqq2,wkqq1)
      ldu=dble(0.)
      do 500 s=1,m
         do 100 i=1,q
            do 80 j=i,q
               sqrtu(i,j,s)=wkqq1(i,j)+ztvinvz(i,j,s)
 80         continue
 100     continue
         call chle(q,q,m,sqrtu,s,err)
         call bkslvl(q,q,m,sqrtu,s)
         do 110 i=1,q
            ldu=ldu+dlog(sqrtu(i,i,s))
 110     continue
         call mmul(q,q,m,sqrtu,s,wkqq2)
         do 200 i=1,q
            do 150 j=i,q
               u(i,j,s)=wkqq2(i,j)
 150        continue
 200     continue
 500  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine drbeta(p,beta,wkpp,sigma2)
C Adds N(0,sigma2*Gamma) noise to beta. Assumes that wkpp contains
C the upper-tri square root of Gamma, wkpp%*%t(wkpp)=Gamma.
      integer p,i,j
      real gauss
      double precision beta(p),wkpp(p,p),sigma2,z,sqrtsig2
      sqrtsig2=dsqrt(sigma2)
      do 20 i=1,p
         z=dble(gauss())
         do 10 j=1,i
            beta(j)=beta(j)+wkpp(j,i)*z*sqrtsig2
 10      continue
 20   continue
      return
      end
C***********************************************************************
      subroutine drb(m,q,sqrtu,sigma2,b)
C Adds N(0,sigma2*U_i) noise to b_i, i=1,...,m
      integer m,q,i,j,s
      double precision sqrtu(q,q,m),b(q,m),sigma2,z,sqrtsig2
      real gauss
      sqrtsig2=dsqrt(sigma2)
      do 100 s=1,m
         do 20 i=1,q
            z=dble(gauss())
            do 15 j=1,i
               b(j,s)=b(j,s)+sqrtu(j,i,s)*z*sqrtsig2
 15         continue
 20      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine bfac(p,m,b)
C draws upper triangular square-root of a Wishart(m,I) using Bartlett
C decomposition.
      integer p,j,k
      real m,jnk,gauss,gamm
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
      subroutine drxi(m,q,b,xi,wkqq1,wkqq2,dinv,sigma2,abc)
C Draws xi from its augmented-data posterior
      integer m,q,i,j,k,s,err
      double precision b(q,m),xi(q,q),wkqq1(q,q),dinv(q,q),sigma2,
     /     abc(3),sum,wkqq2(q,q)
C     *** put inverse square root of scale matrix into wkqq1 *****
      do 5 i=1,q
         do 4 j=i,q
            wkqq1(i,j)=dinv(i,j)
 4       continue
 5    continue
      do 20 s=1,m
         do 10 i=1,q
            do 9 j=i,q
               wkqq1(i,j)=wkqq1(i,j)+b(i,s)*b(j,s)
 9          continue
 10      continue
 20   continue
      do 25 i=1,q
         do 24 j=i,q
            wkqq1(i,j)=wkqq1(i,j)/sigma2
 24      continue
 25   continue
      call chfce(q,q,wkqq1,err)
C     *** draw xi ****************************
      call bfac(q,float(m)+sngl(abc(3)),xi)
      call bkslv(q,q,xi)
      do 40 i=1,q
         do 39 j=1,q
            sum=dble(0.)
            do 38 k=1,min(i,j)
               sum=sum+wkqq1(k,i)*xi(k,j)
 38         continue
            wkqq2(i,j)=sum
 39      continue
 40   continue
      do 50 i=1,q
         do 49 j=i,q
            sum=dble(0.)
             do 48 k=1,q
                sum=sum+wkqq2(i,k)*wkqq2(j,k)
 48          continue
             xi(i,j)=sum
             if(i.ne.j) xi(j,i)=xi(i,j)
 49       continue
 50    continue
       return
       end
C***********************************************************************
C Algorithm MGIBBS: modified Gibbs sampler
      subroutine mgibbs(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,vinv,
     /     pcol,pred,q,zcol,ztvinv,ztvinvz,isflag,err,msg,u,sigma2,
     /     p,xcol,beta,y,delta,xtw,xtwx,xtwy,xtwxinv,wkqq1,
     /     wkqq2,xi,wkqnm,b,maxits,abc,dinv,sqrtu,sigma2s,psis)
C set sflag=1 if starting values supplied.
C Returned error codes in msg:
C     0 = no errors
C     1 = supplied V_i matrix is not positive definite
C     2 = GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank
C     3 = inadequate information to obtain starting value of xi
C     4 = value of xi or inv(Ui) became non-pos.def. during iterations
C     5 = t(X)%*%W%*%X became non-pos.def. during iterations
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,
     /     pcol,q,zcol(q),iflag,err,msg,iter,sflag,p,xcol(p),i,j,
     /     maxits,isflag(2)
      double precision vmax(nmax,nmax),w(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),
     /     ztvinv(q,nmax,m),ztvinvz(q,q,m),ldv,u(q,q,m),sigma2,
     /     beta(p),y(ntot),delta(ntot),xtw(p,nmax),xtwx(p,p),xtwy(p),
     /     xtwxinv(p,p),wkqq1(q,q),wkqq2(q,q),xi(q,q),
     /     wkqnm(q,nmax,m),b(q,m),osigma2,ldu,ldxi,abc(3),dinv(q,q),
     /     trdx,aprime,sqrtu(q,q,m),sigma2s(maxits),psis(q,q,maxits)
      real bprime,jnk,gauss,gamm
      iflag=isflag(1)
      sflag=isflag(2)
      jnk=gauss()
      msg=0
      iter=0
      call preecme1(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,ldv,err)
      if(err.eq.1) then
         msg=1
         goto 999
      endif
      if(sflag.ne.1) then
         call stval1(ntot,m,ist,ifin,occ,nmax,vinv,pcol,pred,
     /        q,ztvinv,ztvinvz,iflag,err,msg,sigma2,p,xcol,
     /        beta,y,delta,xtw,xtwx,xtwy,xtwxinv,wkqq1,wkqq2,xi,
     /        wkqnm,b)
      endif
C********* start of main iteration *****************
 1    continue
         iter=iter+1
         osigma2=sigma2
         call mku2(q,xi,m,ztvinvz,u,wkqq1,wkqq2,ldxi,ldu,err,sqrtu)
         if(err.eq.1) then
            msg=4
            goto 999
         endif
         call trdixi(q,trdx,wkqq1,dinv)
         call mkwkqnm(q,m,u,nmax,ztvinv,wkqnm,ntot,occ,ist,ifin)
         call mkw(q,nmax,m,ist,ifin,wkqnm,ztvinv,vinv,w,ntot,occ,
     /        iflag)
         call gls(ntot,m,ist,ifin,occ,nmax,pcol,pred,0,sigma2,
     /        p,xcol,beta,y,delta,w,xtw,xtwx,xtwy,xtwxinv,err)
         if(err.eq.1) then
            msg=5
            goto 999
         endif
         aprime=abc(1)+trdx+dfloat(ntot)*sigma2
         bprime=sngl(abc(2)+dfloat(ntot-p)+abc(3)*dfloat(q))
         sigma2=aprime/dble(2.*gamm(bprime/2.))
         call drbeta(p,beta,xtwx,osigma2)
         call mkdel(ntot,pcol,pred,p,xcol,y,beta,delta)
         call mkb(q,nmax,m,wkqnm,ntot,delta,b,occ,ist,ifin)
         call drb(m,q,sqrtu,osigma2,b)
         call drxi(m,q,b,xi,wkqq1,wkqq2,dinv,osigma2,abc)
         sigma2s(iter)=sigma2
         do 10 i=1,q
            do 9 j=1,q
               psis(i,j,iter)=xi(i,j)*sigma2
 9          continue
 10      continue
         if(iter.lt.maxits) goto 1
C********* end of main iteration *****************
      call bdiag(q,m,u)
      do 70 i=1,p
         do 60 j=i+1,p
            xtwxinv(j,i)=xtwxinv(i,j)
 60      continue
 70   continue
 999  continue
      return
      end
C***********************************************************************
      subroutine drcand(g,wkg,gmax,estarhat,wkgg,q,xi,sigma2,df,
     /     logdens,wkqq1,ntries)
C Draws a candidate for Metropolis-Hastings, storing it in sigma2
C and xi, and calculates its approximate density on the eta scale.
C   ntries = no. of attempts required to get a positive definite
C   candidate.
      integer q,g,gmax,i,j,gi,err,ntries
      real gauss,gamm
      double precision wkg(0:g),estarhat(0:gmax),wkgg(0:g,0:g),
     /     xi(q,q),sigma2,df,sum,logdens,chsq,tmp,wkqq1(q,q)
      ntries=0
 1    continue
      ntries=ntries+1
C     *** draw the chisquare and normal variates, and calculate ***
C     *** the log-density on the eta scale ************************
      chsq=dble(2.*gamm(sngl(df)/2.))
      sum=dble(0.)
      do 5 i=0,g
         wkg(i)=dble(gauss())
         sum=sum+wkg(i)**2
 5    continue
      tmp = dble(1.)+sum/df
      logdens= (-(df+dfloat(g+1))/dble(2.))*dlog(tmp)
C     **** premultiply wkg by wkgg ***********************
      do 15 i=0,g
         sum=dble(0.)
         do 14 j=i,g
            sum=sum+wkgg(i,j)*wkg(j)
 14      continue
         wkg(i)=sum
 15   continue
C     *** finish off the candidate value in wkg ***********
      do 20 i=0,g
         wkg(i)=wkg(i)*dsqrt((df+dfloat(g+1))/chsq) + estarhat(i)
 20   continue
C     *** store the candidate values in sigma2 and xi while ***
C     ***** calculating the Jacobian **************************
C     ***** and reject if outside the parameter space *********
      sum=wkg(0)
      sigma2=dexp(-wkg(0))
      gi=0
      do 30 i=1,q
         do 29 j=i,q
            gi=gi+1
            if(i.eq.j) then
               sum=sum+wkg(gi)
               wkqq1(i,j)=dexp(wkg(gi))
            else
               wkqq1(i,j)=wkg(gi)
            endif
 29      continue
 30   continue
      logdens=logdens-sum
      call chfce(q,q,wkqq1,err)
      if(err.eq.1) goto 1
      call bkslv(q,q,wkqq1)
      call mm(q,q,wkqq1,xi)
      do 35 i=1,q
         do 34 j=i+1,q
            xi(j,i)=xi(i,j)
 34      continue
 35   continue
      return
      end
C***********************************************************************
      subroutine appxdens(q,xi,sigma2,g,wkg,wkgg2,df,gmax,estarhat,
     /     wkqq1,wkqq2,logdens)
C Calculates the approximation density at (xi,sigma2) 
      integer q,g,gmax,err,gi,i,j
      double precision xi(q,q),sigma2,wkg(0:g),wkgg2(0:g,0:g),df,
     /     estarhat(0:gmax),wkqq1(q,q),wkqq2(q,q),sum,logdens,tmp
C     *** store etastar-estarhat in wkg ****************
      wkg(0)=-dlog(sigma2)-estarhat(0)
      do 5 i=1,q
         do 4 j=1,q
            wkqq1(i,j)=xi(i,j)
 4       continue
 5    continue
      call chfce(q,q,wkqq1,err)
      call bkslv(q,q,wkqq1)
      call mm(q,q,wkqq1,wkqq2)
      gi=0
      do 10 i=1,q
         do 9 j=i,q
            gi=gi+1
            if(i.eq.j) then
              wkg(gi)=dlog(wkqq2(i,j))-estarhat(gi)
            else
              wkg(gi)=wkqq2(i,j)-estarhat(gi)
            endif
 9       continue
 10   continue
C     **** premultiply wkg by wkgg2 ***********************
      do 15 i=0,g
         sum=dble(0.)
         do 14 j=i,g
            sum=sum+wkgg2(i,j)*wkg(j)
 14      continue
         wkg(i)=sum
 15   continue
C     **** calculate quadratic form *****************
      sum=dble(0.)
      do 20 i=0,g
         sum=sum+wkg(i)**2
 20   continue
C     ********* calculate density on etastar scale **********
      tmp = dble(1.)+sum/(df+dfloat(g+1))
      logdens= (-(df+dfloat(g+1))/dble(2.))*dlog(tmp)
C     ******** add in the log-Jacobian **********************
      sum=-dlog(sigma2)
      do 25 i=1,q
         sum=sum+dlog(wkqq2(i,i))
 25   continue
      logdens=logdens-sum
      return
      end
C***********************************************************************
C Algorithm FAST-MCMC
      subroutine fastmcmc(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,vinv,
     /     pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,err,msg,u,sigma2,
     /     p,xcol,beta,y,delta,xtw,xtwx,xtwy,xtwxinv,wkqq1,
     /     wkqq2,xi,wkqnm,b,maxits,abc,dinv,sqrtu,sigma2s,psis,
     /     g,wkgg,wkgg2,wkg,sig2hat,xihat,xigibbs,reject,ratios,df)
C Returned error codes in msg:
C     0 = no errors
C     1 = supplied V_i matrix is not positive definite
C     2 = GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank
C     3 = inadequate information to obtain starting value of xi
C     4 = value of xi or inv(Ui) became non-pos.def. during iterations
C     5 = t(X)%*%W%*%X became non-pos.def. during iterations
C     6 = supplied xihat is non-pos.def.
      integer gmax
      data gmax/50/
      double precision estarhat(0:50)
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,
     /     pcol,q,zcol(q),iflag,err,msg,iter,p,xcol(p),i,j,
     /     maxits,g,gi,reject(maxits),accept,ntries
      double precision vmax(nmax,nmax),w(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),
     /     ztvinv(q,nmax,m),ztvinvz(q,q,m),ldv,u(q,q,m),sigma2,
     /     beta(p),y(ntot),delta(ntot),xtw(p,nmax),xtwx(p,p),xtwy(p),
     /     xtwxinv(p,p),wkqq1(q,q),wkqq2(q,q),xi(q,q),
     /     wkqnm(q,nmax,m),b(q,m),ldu,ldxi,abc(3),dinv(q,q),
     /     trdx,aprime,sqrtu(q,q,m),sigma2s(maxits),psis(q,q,maxits),
     /     wkgg(0:g,0:g),wkgg2(0:g,0:g),wkg(0:g),
     /     xihat(q,q),sig2hat,df,lh,s2hat,ldxtwx,lp,s2gibbs,
     /     xigibbs(q,q),lhcand,ratios(maxits),olp,olh,lratio
      real bprime,jnk,gauss,gamm,rangen
      jnk=gauss()
      msg=0
      iter=0
      call preecme1(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,ldv,err)
      if(err.eq.1) then
         msg=1
         goto 999
      endif
C     *** store modal value of etastar in estarhat ***
      estarhat(0)=-dlog(sig2hat)
      do 5 i=1,q
         do 4 j=1,q
            wkqq1(i,j)=xihat(i,j)
 4       continue
 5    continue
      call chfce(q,q,wkqq1,err)
      if(err.eq.1) then
         msg=6
         goto 999
      endif
      call bkslv(q,q,wkqq1)
      call mm(q,q,wkqq1,wkqq2)
      gi=0
      do 10 i=1,q
         do 9 j=i,q
            gi=gi+1
            if(i.eq.j) then
              estarhat(gi)=dlog(wkqq2(i,j))
            else
              estarhat(gi)=wkqq2(i,j)
            endif
 9       continue
 10   continue
C     *** store inverse of wkgg in wkgg2 *************
      do 15 i=0,g
         do 14 j=i,g
            wkgg2(i,j)=wkgg(i,j)
 14      continue
 15   continue
      call bkslv(g+1,g+1,wkgg2)
C********* start of main iteration *****************
 50   continue
         iter=iter+1
         reject(iter)=0
 53      continue
         olp=lp
         olh=lh
         call mku2(q,xi,m,ztvinvz,u,wkqq1,wkqq2,ldxi,ldu,err,sqrtu)
         if(err.eq.1) then
            msg=4
            goto 999
         endif
         call trdixi(q,trdx,wkqq1,dinv)
         call mkwkqnm(q,m,u,nmax,ztvinv,wkqnm,ntot,occ,ist,ifin)
         call mkw(q,nmax,m,ist,ifin,wkqnm,ztvinv,vinv,w,ntot,occ,
     /        iflag)
         call gls(ntot,m,ist,ifin,occ,nmax,pcol,pred,0,s2hat,
     /        p,xcol,beta,y,delta,w,xtw,xtwx,xtwy,xtwxinv,err)
         if(err.eq.1) then
            msg=5
            goto 999
         endif
         s2hat= (abc(1) + s2hat*dfloat(ntot) + trdx)/
     /        (abc(2)+dfloat(ntot-p-2)+abc(3)*dfloat(q))
C        *** evaluate log-posterior at current parameters *********
C        *** note that xtwx contains the square root of xtwxinv ***
C        ***** first we get its log-determinant *******************
         ldxtwx=dble(0.)
         do 65 i=1,p
            ldxtwx=ldxtwx+dlog(xtwx(i,i))
 65      continue
         lp = -dble(.5)*(dfloat(ntot-p-2)+abc(2)+abc(3)*dfloat(q))
     /        *(dlog(sigma2) + s2hat/sigma2)
     /        + (dfloat(m-q-1)+abc(3))*ldxi + ldu 
     /        + ldxtwx
C        *** evaluate approximate density at current parameters *****
         call appxdens(q,xi,sigma2,g,wkg,wkgg2,df,gmax,estarhat,
     /     wkqq1,wkqq2,lh)
C        *** calculate log-acceptance ratio and decide whether to ***
C        ***** keep the Metropolis candidate ************************
         if(iter.gt.1) then
            if(reject(iter-1).eq.0) then
               lratio=lp-olp-lh+olh
               if(lratio.gt.dble(0.)) then
                  accept=1
               else if(lratio.lt.dble(-30.)) then
                  accept=0
               else
                  if(dble(rangen(0)).le.dexp(lratio)) then
                     accept=1
                  else
                     accept=0
                  endif
               endif
               ratios(iter-1)=dexp(lratio)
               if(accept.eq.1) then
                  sigma2s(iter-1)=sigma2
                  do 70 i=1,q
                     do 69 j=1,q
                        psis(i,j,iter-1)=xi(i,j)*sigma2
 69                  continue
 70               continue
               else
                  reject(iter-1)=1
                  sigma2=s2gibbs
                  sigma2s(iter-1)=sigma2
                  do 80 i=1,q
                     do 79 j=1,q
                        xi(i,j)=xigibbs(i,j)
                        psis(i,j,iter-1)=xi(i,j)*sigma2
 79                  continue
 80               continue
                  goto 53
               endif
            endif
         endif
C        *** finish out the Gibbs cycle *****************************
         aprime=s2hat*(abc(2)+dfloat(ntot-p-2)+abc(3)*dfloat(q))
         bprime=sngl(abc(2)+dfloat(ntot-p)+abc(3)*dfloat(q))
         s2gibbs=aprime/dble(2.*gamm(bprime/2.))
         call drbeta(p,beta,xtwx,sigma2)
         call mkdel(ntot,pcol,pred,p,xcol,y,beta,delta)
         call mkb(q,nmax,m,wkqnm,ntot,delta,b,occ,ist,ifin)
         call drb(m,q,sqrtu,sigma2,b)
         call drxi(m,q,b,xigibbs,wkqq1,wkqq2,dinv,sigma2,abc)
C        *** draw candidate for Metropolis-Hastings ****************
C        **** storing the values in sigma2 and xi ******************
         call drcand(g,wkg,gmax,estarhat,wkgg,q,xi,sigma2,df,
     /     lhcand,wkqq1,ntries)
         if(iter.le.maxits) goto 50
C********* end of main iteration *****************
      call bdiag(q,m,u)
      do 170 i=1,p
         do 160 j=i+1,p
            xtwxinv(j,i)=xtwxinv(i,j)
 160     continue
 170  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine trdixi(q,trdx,xiinv,dinv)
C Finds trace of inv(D)%*%inv(xi)
      integer i,j,q
      double precision trdx,xiinv(q,q),dinv(q,q)
C     *** symmetrize xiinv ************
      do 2 i=1,q
         do 1 j=i+1,q
            xiinv(j,i)=xiinv(i,j)
 1       continue
 2    continue
C     *** calculate trace *************
      trdx=dble(0.)
      do 10 i=1,q
         do 5 j=1,q
            trdx=trdx+dinv(i,j)*xiinv(j,i)
 5       continue
 10   continue
      return
      end
C***********************************************************************
      subroutine mkocc(ntot,occ,m,ist,ifin)
C Fills in occ with 1,...,ni for subjects i=1,...,m
      integer ntot,occ(ntot),m,ist(m),ifin(m),s,i,j
      do 20 s=1,m
         j=0
         do 10 i=ist(s),ifin(s)
            j=j+1
            occ(i)=j
 10      continue
 20   continue
      return
      end
C***********************************************************************
      subroutine preecme1(ntot,subj,m,ist,ifin,occ,nmax,vmax,wknnm,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,ldv,err)
C Preliminary manipulations for ECME-ML. After execution, 
C     ist  = vector of starting positions for subject i, i=1,...,m
C     ifin = vector of finishing positions for subject i, i=1,...,m
C     vinv   = inverse of V_i, i=1,...,m
C     ztvinv  = t(z_i)%*% inverse of V_i, i=1,...,m
C     ztvinvz = t(z_i)%*%inv(V_i)%*%(z_i), i=1,...,m
C     ldv = .5*sum of log det(V_i)
C Needs a workspace wknnm.
C Note: if V_i's are all identity (indicated by iflag=1) then 
C vmax, wknnm, vinv, are ignored and occ is filled in with 1,...,ni.
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,pcol,
     /     q,zcol(q),iflag,err
      double precision vmax(nmax,nmax),wknnm(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),ztvinv(q,nmax,m),
     /     ztvinvz(q,q,m),ldv
      call istfin(ntot,subj,m,ist,ifin)
      if(iflag.ne.1) then
         call mkv(m,nmax,vmax,ntot,occ,ist,ifin,wknnm)
         call chv(nmax,m,wknnm,ntot,occ,ist,ifin,ldv,err)
         if(err.eq.1) goto 999
         call bkv(nmax,m,wknnm,ntot,occ,ist,ifin)
         call mmulv(nmax,m,wknnm,vinv,ntot,occ,ist,ifin)
      else
         call mkocc(ntot,occ,m,ist,ifin)
         ldv=dble(0.)
      endif
      call mmu(ntot,pcol,pred,q,zcol,nmax,m,wknnm,occ,ist,ifin,ztvinv,
     /     iflag)
      call mmtm(q,nmax,m,ztvinv,ntot,occ,ist,ifin,ztvinvz)
      if(iflag.ne.1) then
         call mml(ntot,q,nmax,m,wknnm,occ,ist,ifin,ztvinv)
      endif
 999  continue
      return
      end
C***********************************************************************
C The following subroutines are used by fastrml().
C***********************************************************************
C Algorithm CEBAYES: To be run after FAST-RML. Provides corrected
C empirical Bayes estimates of V(beta), V(b_i), and Cov(b_i,beta)
C for i=1,...,m.
      subroutine cebayes(m,p,q,b,u,a,xtwxinv,sigma2,ztvinvx,
     /     g,wkgg,wkqp,wkpg,wkqg,wkp,varbeta,varb,covbbeta,
     /     wkpg2,wkqg2,xi,wkqq1,ntot,err)
      integer m,p,q,g,ntot
      double precision b(q,m),u(q,q,m),a(q,q,m),xtwxinv(p,p),
     /     sigma2,ztvinvx(q,p,m),wkgg(0:g,0:g),wkqp(q,p),
     /     wkpg(p,g),wkqg(q,g),wkp(p),varbeta(p,p),varb(q,q,m),
     /     covbbeta(q,p,m),wkpg2(p,0:g),wkqg2(q,0:g),
     /     xi(q,q),wkqq1(q,q)
      integer i,j,k,s,gi,ii,gj,jj,jjmin,err
      double precision sum,trahah,trahaj,trajaj
C     ******************************************************
C     ******************************************************
C     *** this section was added to make sure that wkgg ****
C     *** contains the second derivatives with respect to **
C     *** eta rather than etastar **************************
      do 508 i=0,g
         do 507 j=0,g
            wkgg(i,j)=dble(0.)
 507     continue
 508  continue
      do 800 s=1,m
C        *** store (xi-U_i) in wkqq1 ******************
         do 585 i=1,q
            do 584 j=i,q
                wkqq1(i,j)=xi(i,j)-u(i,j,s)
                if(i.ne.j) wkqq1(j,i)=wkqq1(i,j)
 584        continue
 585     continue
C        *** accumulate zeroeth row of wkgg ************
         gi=0
         do 600 i=1,q
            do 590 j=i,q
               gi=gi+1
               if(i.eq.j) then
                  wkgg(0,gi)=wkgg(0,gi)+wkqq1(i,i)
               else
                  wkgg(0,gi)=wkgg(0,gi)+dble(2.)*wkqq1(i,j)
               endif
 590        continue
 600     continue
C        *** accumulate rest of wkgg ********************
         gi=0
         do 700 i=1,q
            do 690 j=i,q
               gi=gi+1
               gj=gi-1
               do 680 ii=i,q
                  if(ii.eq.i) then
                     jjmin=j
                  else
                     jjmin=ii
                  endif
                  do 670 jj=jjmin,q
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahah(q,wkqq1,i,ii)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,i,ii,jj)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,ii,i,j)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trajaj(q,wkqq1,i,j,ii,jj)
                        endif
                     endif
 670              continue
 680           continue
 690        continue
 700     continue
 799     continue
 800  continue
      wkgg(0,0)=dfloat(ntot-p)*(sigma2**2)/dble(2.)
      gi=0
      do 825 i=1,q
         do 824 j=i,q
            gi=gi+1
            wkgg(0,gi)=wkgg(0,gi)*sigma2/dble(2.)
            wkgg(gi,0)=wkgg(0,gi)
            do 822 gj=gi,g
               wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
               if(gi.ne.gj) wkgg(gj,gi)=wkgg(gi,gj)
 822        continue
 824     continue
 825  continue
C     **** end of computation of wkgg **********************
C     ******************************************************
C     *** now find its inverse square-root *****************
      err=0
      call chfce(g+1,g+1,wkgg,err)
      if(err.eq.1) goto 999
      call bkslv(g+1,g+1,wkgg)
C     ******************************************************
C     *** initialize wkpg **********************************
      do 3 gi=1,g
         do 2 i=1,p
               wkpg(i,gi)=dble(0.)
 2       continue
 3    continue
C     *** first pass ******************************************
      do 100 s=1,m
C        *** find U_i%*%gamma_i and store it in covbbeta(*,*,s) ****
         do 20 i=1,q
            do 19 j=1,p
               sum=dble(0.)
               do 18 k=1,q
                  sum=sum+u(i,k,s)*ztvinvx(k,j,s)
 18            continue
               covbbeta(i,j,s)=sum
 19         continue
 20      continue
C        *** accumulate t(gamma_i)%*%U_i%*%G_j%*%b_i into wkpg ***
C        ****** for j=1,...,g ************************************
         gi=0
         do 50 i=1,q
            do 40 j=i,q
               gi=gi+1
               if(i.eq.j) then
                  do 35 ii=1,p
                     wkpg(ii,gi)=wkpg(ii,gi)+covbbeta(i,ii,s)*b(i,s)
 35               continue
               else
                  do 30 ii=1,p
                     wkpg(ii,gi)=wkpg(ii,gi)+covbbeta(j,ii,s)*b(i,s)
     /                    + covbbeta(i,ii,s)*b(j,s)
 30               continue
               endif
 40         continue
 50      continue
 100  continue
C     *** symmetrize xtwxinv **********************
      do 110 i=1,p
         do 109 j=i+1,p
            xtwxinv(j,i)=xtwxinv(i,j)
 109     continue
 110  continue
C     *** premultiply columns of wkpg by xtwxinv to get **********
C     ***** the derivatives of betawig with respect to omega_j ***
      do 150 gi=1,g
         do 140 i=1,p
            sum=dble(0.)
            do 130 j=1,p
               sum=sum+xtwxinv(i,j)*wkpg(j,gi)
 130        continue
            wkp(i)=sum
 140     continue
         do 145 i=1,p
            wkpg(i,gi)=wkp(i)
 145     continue
 150  continue
C     *** store (d beta/d eta)%*%wkgg in wkpg2 *****************
C     ***** Notice that wkgg contains the square root of Cinv **
      do 160 i=1,p
         wkpg2(i,0)=dble(0.)
         do 159 j=1,g
            sum=dble(0.)
            do 158 k=1,j
               sum=sum+wkpg(i,k)*wkgg(k,j)
 158        continue
            wkpg2(i,j)=sum
 159     continue
 160  continue
C     *** calculate varbeta **********************************
      do 170 i=1,p
         do 169 j=i,p
            sum=dble(0.)
            do 168 k=0,g
               sum=sum+wkpg2(i,k)*wkpg2(j,k)
 168        continue
            varbeta(i,j)=sigma2*xtwxinv(i,j)+sum
            if(i.ne.j) varbeta(j,i)=varbeta(i,j)
 169     continue
 170  continue
C     *** second pass ****************************************
      do 400 s=1,m
C        *** calculate the derivatives of b_i with respect to ***
C        ****  omega_j, j=1,...,g and store in wkqg *************
         gi=0
         do 200 i=1,q
            do 190 j=i,q
               gi=gi+1
               if(i.eq.j) then
                  do 175 ii=1,q
                     wkqg(ii,gi)=-u(ii,i,s)*b(i,s)
 175              continue
               else
                  do 180 ii=1,q
                     wkqg(ii,gi)=-u(ii,j,s)*b(i,s)-u(ii,i,s)*b(j,s)
 180              continue
               endif
 190           continue
 200     continue
         do 250 gi=1,g
            do 240 i=1,q
               sum=dble(0.)
               do 235 j=1,p
                  sum=sum+covbbeta(i,j,s)*wkpg(j,gi)
 235           continue
               wkqg(i,gi)=wkqg(i,gi)-sum
 240        continue
 250     continue
C     *** store (d b_i/d eta)%*%wkgg in wkqg2 *****************
         do 260 i=1,q
            wkqg2(i,0)=dble(0.)
            do 259 j=1,g
               sum=dble(0.)
               do 258 k=1,j
                  sum=sum+wkqg(i,k)*wkgg(k,j)
 258           continue
               wkqg2(i,j)=sum
 259        continue
 260     continue
C     *** calculate varb(*,*,s) **********************************
         do 270 i=1,q
            do 269 j=i,q
               sum=dble(0.)
               do 268 k=0,g
                  sum=sum+wkqg2(i,k)*wkqg2(j,k)
 268           continue
               varb(i,j,s)=sigma2*(u(i,j,s)+a(i,j,s))+sum
               if(i.ne.j) varb(j,i,s)=varb(i,j,s)
 269        continue
 270     continue
C        *** postmultiply covbbeta(*,*,s) by -sigma2*xtwxinv *****
C        *** and store result in wkqp ****************************
         do 300 i=1,q
            do 290 j=1,p
               sum=dble(0.)
               do 285 k=1,p
                  sum=sum+covbbeta(i,k,s)*xtwxinv(k,j)
 285           continue
               wkqp(i,j)=-sigma2*sum
 290        continue
 300     continue
C        *** calculate covbbeta **********************************
         do 320 i=1,q
            do 315 j=1,p
               sum=dble(0.)
               do 310 k=0,g
                  sum=sum+wkqg2(i,k)*wkpg2(j,k)
 310           continue
               covbbeta(i,j,s)=wkqp(i,j)+sum
 315        continue
 320     continue
 400  continue
 999  continue
      return
      end
C***********************************************************************
C Algorithm FAST-RML: Finds RML estimates when V_i are known.
      subroutine fastrml(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,vinv,
     /     pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,err,msg,u,iter,
     /     sflag,sigma2,p,xcol,beta,y,delta,xtw,xtwx,xtwy,
     /     xtwxinv,wkqq1,wkqq2,xi,wkqnm,b,cvgd,oxi,maxits,llvec,
     /     eps,xiecme,g,reject,ztvinvx,a,wkqp,wkg,wkgg)
C set sflag=1 if starting values supplied.
C Returned error codes in msg:
C     0 = no errors
C     1 = supplied V_i matrix is not positive definite
C     2 = GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank
C     3 = inadequate information to obtain starting value of xi
C     4 = value of xi or inv(Ui) became non-pos.def. during iterations
C     5 = t(X)%*%W%*%X became non-pos.def. during iterations
C    10 = non-positive definite xi at input to scoring step
C    11 = non-positive definite wkgg in scoring step
C
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,
     /     pcol,q,zcol(q),iflag,err,msg,iter,sflag,p,xcol(p),cvgd,
     /     maxits,c2,c3,i,j,g,reject(maxits),notconcave
      double precision vmax(nmax,nmax),w(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),
     /     ztvinv(q,nmax,m),ztvinvz(q,q,m),ldv,u(q,q,m),sigma2,
     /     beta(p),y(ntot),delta(ntot),xtw(p,nmax),xtwx(p,p),xtwy(p),
     /     xtwxinv(p,p),wkqq1(q,q),wkqq2(q,q),xi(q,q),
     /     wkqnm(q,nmax,m),b(q,m),oxi(q,q),osigma2,ldu,ldxi,
     /     ll,llvec(maxits),eps,xiecme(q,q),
     /     ztvinvx(q,p,m),a(q,q,m),wkqp(q,p),ldxtwx,
     /     sig2ecme,osig2ecme,wkg(0:g),wkgg(0:g,0:g)
      msg=0
      iter=0
      notconcave=0
      call prefstrm(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,p,xcol,ztvinvx,
     /     iflag,ldv,err)
      if(err.eq.1) then
         msg=1
         goto 999
      endif
      if(sflag.ne.1) then
         call stval1(ntot,m,ist,ifin,occ,nmax,vinv,pcol,pred,
     /        q,ztvinv,ztvinvz,iflag,err,msg,sigma2,p,xcol,
     /        beta,y,delta,xtw,xtwx,xtwy,xtwxinv,wkqq1,wkqq2,xi,
     /        wkqnm,b)
      endif
      cvgd=0
C********* start of main iteration *****************
 1    continue
         iter=iter+1
         reject(iter)=0
 3       continue
         osigma2=sigma2
C        *** calculate U_i, U_i%*%t(Z_i)%*%inv(V_i), and W_i ****
         call mku(q,xi,m,ztvinvz,u,wkqq1,wkqq2,ldxi,ldu,err)
         if(err.eq.1) then
            msg=4
            goto 999
         endif
         call mkwkqnm(q,m,u,nmax,ztvinv,wkqnm,ntot,occ,ist,ifin)
         call mkw(q,nmax,m,ist,ifin,wkqnm,ztvinv,vinv,w,ntot,occ,
     /        iflag)
C        *** save current value of sigma2ecme ******************
         osig2ecme=sig2ecme
C        *** calculate beta, delta_i, and sigma2ecme ***********
         call gls(ntot,m,ist,ifin,occ,nmax,pcol,pred,0,sig2ecme,
     /        p,xcol,beta,y,delta,w,xtw,xtwx,xtwy,xtwxinv,err)
         if(err.eq.1) then
            msg=5
            goto 999
         endif
         sig2ecme=sig2ecme*dfloat(ntot)/dfloat(ntot-p)
C        *** evaluate loglikelihood ************************
C        *** note that xtwx contains the square root of xtwxinv ***
C        ***** first we get its log-determinant *******************
         ldxtwx=dble(0.)
         do 15 i=1,p
            ldxtwx=ldxtwx+dlog(xtwx(i,i))
 15      continue
         ll = -.5*dfloat(ntot-p)*dlog(osigma2)
     /        + dfloat(m)*ldxi + ldu + ldxtwx
     /        -.5*dfloat(ntot-p)*sig2ecme/osigma2
         llvec(iter)=ll
         if(iter.gt.1) then
C           *** if loglikelihood has gone down, replace the scoring ****
C           ******  estimate of xi by the ecme estimate ****************
            if(reject(iter-1).eq.0) then
               if(llvec(iter).lt.llvec(iter-1)) then
                  sigma2=osig2ecme
                  do 17 i=1,q
                     do 16 j=1,q
                        xi(i,j)=xiecme(i,j)
 16                  continue
 17               continue
                  reject(iter-1)=1
                  goto 3
               endif
            endif
         endif
C        *** check for convergence **********************************
         if(iter.gt.1) then
            c2=0
            do 20 i=1,q
               do 19 j=i,q
                  if(dabs(xi(i,j)-oxi(i,j)).gt.(eps*dabs(oxi(i,j))))
     /                 c2=1
 19            continue
 20         continue
            c3=0
            if(dabs(sigma2-osigma2).gt.(eps*dabs(osigma2))) c3=1
            if((c2.eq.0).and.(c3.eq.0)) then
               cvgd=1
               goto 50
            endif
         endif
C        *** calculate b_i ********************************
         call mkb(q,nmax,m,wkqnm,ntot,delta,b,occ,ist,ifin)
C        *** save current parameters **********************
         osigma2=sigma2
         do 27 i=1,q
            do 26 j=1,q
               oxi(i,j)=xi(i,j)
 26         continue
 27      continue
C        *** perform scoring step, putting result in xi and sigma2,**
C        ****** and put ECME estimate of xi in xiecme ***************
         call fscovr2(m,q,b,u,a,xi,oxi,xiecme,wkqq1,wkqq2,
     /        p,xtwxinv,xtwx,wkqp,g,wkg,wkgg,sigma2,
     /        msg,osigma2,ntot,sig2ecme,ztvinvx)
         if(msg.eq.10) goto 999
         if(msg.eq.11) then
C        *** wkgg not pos.def., substitute ECME value instead ****
            notconcave=1
            do 29 i=1,q
               do 28 j=1,q
                  xi(i,j)=xiecme(i,j)
 28            continue
 29         continue
            sigma2=sig2ecme
            reject(iter)=1
         endif
         if(iter.lt.maxits) goto 1
C********* end of main iteration *****************
 50   continue
      if(notconcave.eq.1) msg=11
      call bdiag(q,m,u)
      do 70 i=1,p
         do 60 j=i+1,p
            xtwxinv(j,i)=xtwxinv(i,j)
 60      continue
 70   continue
 999  continue
      return
      end
C***********************************************************************
      subroutine fscovr2(m,q,b,u,a,xi,oxi,xiecme,wkqq1,wkqq2,
     /     p,xtwxinv,wkpp,wkqp,g,wkg,wkgg,sigma2,msg,
     /     osigma2,ntot,sig2ecme,ztvinvx)
C Fisher scoring algorithm to perform RML estimation for xi and sigma2.
C Like fscovr, but it performs scoring on eta^* rather than eta.
C       oxi = old value of xi
C        xi = new value
C    xiecme = ecme version of xi
C Error messages returned through msg:
C     0 = no error
C    10 = non-positive definite xi
C    11 = non-positive definite wkgg
      integer m,q,p,g,msg,ntot
      double precision b(q,m),u(q,q,m),a(q,q,m),
     /     xi(q,q),oxi(q,q),xiecme(q,q),wkqq1(q,q),wkqq2(q,q),
     /     xtwxinv(p,p),wkpp(p,p),wkqp(q,p),wkg(0:g),
     /     wkgg(0:g,0:g),sigma2,sum,deflate,
     /     osigma2,sig2ecme,ztvinvx(q,p,m),oldtau,tau
      integer s,i,j,k,ii,jj,jjmin,gi,gj,err
      double precision trahah,trahaj,trajaj
      msg=0
C     *** initialize xiecme *****
      do 2 i=1,q
         do 1 j=i,q
            xiecme(i,j)=dble(0.)
 1       continue
 2    continue
C     *** store cholesky of xtwxinv in wkpp ******
      do 5 i=1,p
         do 4 j=i,p
            wkpp(i,j)=xtwxinv(i,j)
 4       continue
 5    continue
      call chfce(p,p,wkpp,err)
C     *** initialize the workspaces wkg and wkgg *******
      do 8 i=0,g
         wkg(i)=dble(0.)
         do 7 j=0,g
            wkgg(i,j)=dble(0.)
 7       continue
 8    continue
C     *** main loop to accumulate wkgg and xiecme ******
      do 300 s=1,m
C        *** symmetrize Ui ******************************************
         do 25 i=1,q
            do 24 j=i+1,q
               u(j,i,s)=u(i,j,s)
 24         continue
 25      continue
C        *** store U_i%*%t(Z_i)%*%inv(V_i)%*%X_i in wkqp ***********
         do 30 i=1,q
            do 29 j=1,p
               sum=dble(0.)
               do 27 k=1,q
                  sum=sum+u(i,k,s)*ztvinvx(k,j,s)
 27            continue
               wkqp(i,j)=sum
 29         continue
 30      continue
C        *** multiply wkqp by lower cholesky of xtwxinv *****
         do 60 i=1,q
            do 55 j=1,p
               sum=dble(0.)
               do 52 k=j,p
                  sum=sum+wkqp(i,k)*wkpp(j,k)
 52            continue
               wkqp(i,j)=sum
 55         continue
 60      continue
C        *** multiply wkqp by its transpose, using symmetry *****
C        ***** to get A_i ***************************************
         do 80 i=1,q
            do 70 j=i,q
               sum=dble(0.)
               do 68 k=1,p
                  sum=sum+wkqp(i,k)*wkqp(j,k)
 68            continue
               a(i,j,s)=sum
               if(i.ne.j) a(j,i,s)=sum
 70         continue
 80      continue
C        *** store (xi-U_i) in wkqq1 and ******************
C        ***** accumulate xiecme **********************
         do 85 i=1,q
            do 84 j=i,q
               sum=oxi(i,j)-u(i,j,s)
               wkqq1(i,j)=sum
               if(i.ne.j) wkqq1(j,i)=wkqq1(i,j)
               xiecme(i,j)=xiecme(i,j)+a(i,j,s)+u(i,j,s)
     /              +b(i,s)*b(j,s)/osigma2
 84         continue
 85      continue
C        *** accumulate zeroeth row of wkgg ************
         gi=0
         do 100 i=1,q
            do 90 j=i,q
               gi=gi+1
               if(i.eq.j) then
                  wkgg(0,gi)=wkgg(0,gi)+wkqq1(i,i)
               else
                  wkgg(0,gi)=wkgg(0,gi)+dble(2.)*wkqq1(i,j)
               endif
 90         continue
 100     continue
C        *** accumulate rest of wkgg ********************
         gi=0
         do 200 i=1,q
            do 190 j=i,q
               gi=gi+1
               gj=gi-1
               do 180 ii=i,q
                  if(ii.eq.i) then
                     jjmin=j
                  else
                     jjmin=ii
                  endif
                  do 170 jj=jjmin,q
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahah(q,wkqq1,i,ii)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,i,ii,jj)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,ii,i,j)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trajaj(q,wkqq1,i,j,ii,jj)
                        endif
                     endif
 170              continue
 180           continue
 190        continue
 200     continue
 300  continue
C     *** finish computing xiecme **************
      do 302 i=1,q
         do 301 j=i,q
            xiecme(i,j)=xiecme(i,j)/dfloat(m)
            if(i.ne.j) xiecme(j,i)=xiecme(i,j)
 301     continue
 302  continue
C     *** put the inverse of oxi into wkqq1 ****
      do 308 i=1,q
         do 307 j=i,q
            wkqq2(i,j)=oxi(i,j)
 307     continue
 308  continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) then
         msg=10
         goto 999
      endif
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,wkqq1)
C     ************ finish off wkgg ***********************
      wkgg(0,0)=dfloat(ntot-p)/dble(2.)
      gi=0
      do 325 i=1,q
         do 324 j=i,q
            gi=gi+1
            wkgg(0,gi)=wkgg(0,gi)/dble(2.)
            if(i.eq.j) wkgg(0,gi)=wkgg(0,gi)*wkqq1(i,i)
            wkgg(gi,0)=wkgg(0,gi)
            gj=gi-1
            do 320 ii=i,q
               if(ii.eq.i) then
                  jjmin=j
               else
                  jjmin=ii
               endif
               do 318 jj=jjmin,q
                  gj=gj+1
                  wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                  if(i.eq.j) then
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(i,i)*
     /                       wkqq1(ii,ii)
                     else
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(i,i)
                     endif
                  else
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(ii,ii)
                     endif
                  endif
                  wkgg(gj,gi)=wkgg(gi,gj)
 318           continue
 320        continue
 324     continue
 325  continue
C     **** compute wkg *******************************
      sum=dble(0.)
      gi=0
      do 344 i=1,q
         do 343 j=i,q
            gi=gi+1
            if(i.eq.j) then
               sum=sum+wkgg(0,gi)*dlog(wkqq1(i,i))
            else
               sum=sum+wkgg(0,gi)*wkqq1(i,j)
            endif
 343     continue
 344  continue
      wkg(0)=dfloat(ntot-p)*dble(.5)*(dble(1.)-sig2ecme/osigma2)
     /     - wkgg(0,0)*dlog(osigma2) + sum
      gi=0
      do 352 i=1,q
         do 351 j=i,q
            gi=gi+1
            if(i.eq.j) then
               wkg(gi)=dfloat(m)*dble(.5)*(oxi(i,i)-xiecme(i,i))
     /              *wkqq1(i,i) - wkgg(0,gi)*dlog(osigma2)
            else
               wkg(gi)=dfloat(m)*(oxi(i,j)-xiecme(i,j))
     /              - wkgg(0,gi)*dlog(osigma2)
            endif
            gj=0
            sum=dble(0.)
            do 350 ii=1,q
               do 349 jj=ii,q
                  gj=gj+1
                  if(ii.eq.jj) then
                     sum=sum+wkgg(gi,gj)*dlog(wkqq1(ii,ii))
                  else
                     sum=sum+wkgg(gi,gj)*wkqq1(ii,jj)
                  endif
 349           continue
 350        continue
            wkg(gi)=wkg(gi)+sum
 351     continue
 352  continue
C     *** solve the linear system, storing the result in wkg ******
C     **** wkgg is converted to its inverse square root ***********
      call chfce(g+1,g+1,wkgg,err)
      if(err.eq.1) then
         msg=11
         goto 999
      endif
      call bkslv(g+1,g+1,wkgg)
      do 360 i=g,0,-1
         sum=dble(0.)
         do 359 j=0,i
            sum=sum+wkgg(j,i)*wkg(j)
 359     continue
         wkg(i)=sum
 360  continue
      do 370 i=0,g
         sum=dble(0.)
         do 369 j=i,g
            sum=sum+wkgg(i,j)*wkg(j)
 369     continue
         wkg(i)=sum
 370  continue
C     *** invert wkg, putting the result into xi **************
C     ****** and get new sigma2 *******************************
C     *** step-halving is used here if xi is not pos.def. ****
      deflate=dble(1.)
      oldtau=-dlog(osigma2)
      do 420 i=1,q
         wkqq1(i,i)=dlog(wkqq1(i,i))
 420  continue
 425  continue
      tau=oldtau+deflate*(wkg(0)-oldtau)
      gi=0
      do 430 i=1,q
         do 429 j=i,q
            gi=gi+1
            wkqq2(i,j)=wkqq1(i,j)+deflate*(wkg(gi)-wkqq1(i,j))
 429     continue
         wkqq2(i,i)=dexp(wkqq2(i,i))
 430  continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) then
         deflate=deflate/dfloat(2)
         goto 425
      endif
      sigma2=dexp(-tau)
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,xi)
      do 440 i=1,q
         do 439 j=i+1,q
            xi(j,i)=xi(i,j)
 439     continue
 440  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine fscovr(m,q,b,u,a,xi,oxi,xiecme,wkqq1,wkqq2,
     /     p,xtwxinv,wkpp,wkqp,g,wkg,wkgg,sigma2,msg,
     /     osigma2,ntot,sig2ecme,ztvinvx)
C Fisher scoring algorithm to perform RML estimation for xi and sigma2.
C       oxi = old value of xi
C        xi = new value
C    xiecme = ecme version of xi
C Error messages returned through msg:
C     0 = no error
C    10 = non-positive definite xi
C    11 = non-positive definite wkgg
      integer m,q,p,g,msg,ntot
      double precision b(q,m),u(q,q,m),a(q,q,m),
     /     xi(q,q),oxi(q,q),xiecme(q,q),wkqq1(q,q),wkqq2(q,q),
     /     xtwxinv(p,p),wkpp(p,p),wkqp(q,p),wkg(0:g),
     /     wkgg(0:g,0:g),sigma2,sum,deflate,
     /     osigma2,sig2ecme,ztvinvx(q,p,m),oldtau,tau
      integer s,i,j,k,ii,jj,jjmin,gi,gj,err
      double precision trahah,trahaj,trajaj
      msg=0
C     *** initialize xiecme *****
      do 2 i=1,q
         do 1 j=i,q
            xiecme(i,j)=dble(0.)
 1       continue
 2    continue
C     *** store cholesky of xtwxinv in wkpp ******
      do 5 i=1,p
         do 4 j=i,p
            wkpp(i,j)=xtwxinv(i,j)
 4       continue
 5    continue
      call chfce(p,p,wkpp,err)
C     *** initialize the workspaces wkg and wkgg *******
      do 8 i=0,g
         wkg(i)=dble(0.)
         do 7 j=0,g
            wkgg(i,j)=dble(0.)
 7       continue
 8    continue
C     *** main loop to accumulate wkgg and xiecme ******
      do 300 s=1,m
C        *** symmetrize Ui ******************************************
         do 25 i=1,q
            do 24 j=i+1,q
               u(j,i,s)=u(i,j,s)
 24         continue
 25      continue
C        *** store U_i%*%t(Z_i)%*%inv(V_i)%*%X_i in wkqp ***********
         do 30 i=1,q
            do 29 j=1,p
               sum=dble(0.)
               do 27 k=1,q
                  sum=sum+u(i,k,s)*ztvinvx(k,j,s)
 27            continue
               wkqp(i,j)=sum
 29         continue
 30      continue
C        *** multiply wkqp by lower cholesky of xtwxinv *****
         do 60 i=1,q
            do 55 j=1,p
               sum=dble(0.)
               do 52 k=j,p
                  sum=sum+wkqp(i,k)*wkpp(j,k)
 52            continue
               wkqp(i,j)=sum
 55         continue
 60      continue
C        *** multiply wkqp by its transpose, using symmetry *****
C        ***** to get A_i ***************************************
         do 80 i=1,q
            do 70 j=i,q
               sum=dble(0.)
               do 68 k=1,p
                  sum=sum+wkqp(i,k)*wkqp(j,k)
 68            continue
               a(i,j,s)=sum
               if(i.ne.j) a(j,i,s)=sum
 70         continue
 80      continue
C        *** store (xi-U_i) in wkqq1 and ******************
C        ***** accumulate xiecme **********************
         do 85 i=1,q
            do 84 j=i,q
               sum=oxi(i,j)-u(i,j,s)
               wkqq1(i,j)=sum
               if(i.ne.j) wkqq1(j,i)=wkqq1(i,j)
               xiecme(i,j)=xiecme(i,j)+a(i,j,s)+u(i,j,s)
     /              +b(i,s)*b(j,s)/osigma2
 84         continue
 85      continue
C        *** accumulate zeroeth row of wkgg ************
         gi=0
         do 100 i=1,q
            do 90 j=i,q
               gi=gi+1
               if(i.eq.j) then
                  wkgg(0,gi)=wkgg(0,gi)+wkqq1(i,i)
               else
                  wkgg(0,gi)=wkgg(0,gi)+dble(2.)*wkqq1(i,j)
               endif
 90         continue
 100     continue
C        *** accumulate rest of wkgg ********************
         gi=0
         do 200 i=1,q
            do 190 j=i,q
               gi=gi+1
               gj=gi-1
               do 180 ii=i,q
                  if(ii.eq.i) then
                     jjmin=j
                  else
                     jjmin=ii
                  endif
                  do 170 jj=jjmin,q
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahah(q,wkqq1,i,ii)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,i,ii,jj)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,ii,i,j)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trajaj(q,wkqq1,i,j,ii,jj)
                        endif
                     endif
 170              continue
 180           continue
 190        continue
 200     continue
 300  continue
C     *** finish computing xiecme **************
      do 302 i=1,q
         do 301 j=i,q
            xiecme(i,j)=xiecme(i,j)/dfloat(m)
            if(i.ne.j) xiecme(j,i)=xiecme(i,j)
 301     continue
 302  continue
C     ************ finish off wkgg ***********************
      wkgg(0,0)=dfloat(ntot-p)*(osigma2**2)/dble(2.)
      gi=0
      do 325 i=1,q
         do 324 j=i,q
            gi=gi+1
            wkgg(0,gi)=wkgg(0,gi)*osigma2/dble(2.)
            wkgg(gi,0)=wkgg(0,gi)
            do 322 gj=gi,g
                wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                if(gi.ne.gj) wkgg(gj,gi)=wkgg(gi,gj)
 322        continue
 324      continue
 325   continue
C     *** put the inverse of oxi into wkqq1 ****
      do 335 i=1,q
         do 334 j=i,q
            wkqq2(i,j)=oxi(i,j)
 334     continue
 335  continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) then
         msg=10
         goto 999
      endif
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,wkqq1)
C     **** compute wkg *******************************
      sum=dble(0.)
      gi=0
      do 344 i=1,q
         do 343 j=i,q
            gi=gi+1
            sum=sum+wkgg(0,gi)*wkqq1(i,j)
 343     continue
 344  continue
      wkg(0)=dfloat(ntot-p)*(osigma2-sig2ecme/dble(2.))+sum
      gi=0
      do 352 i=1,q
         do 351 j=i,q
            gi=gi+1
            if(i.eq.j) then
               wkg(gi)=dfloat(m)*(oxi(i,i)-xiecme(i,i))/dfloat(2)
     /              + wkgg(0,gi)/osigma2
            else
               wkg(gi)=dfloat(m)*(oxi(i,j)-xiecme(i,j))
     /              + wkgg(0,gi)/osigma2
            endif
            gj=0
            sum=dble(0.)
            do 350 ii=1,q
               do 349 jj=ii,q
                  gj=gj+1
                  sum=sum+wkgg(gi,gj)*wkqq1(ii,jj)
 349           continue
 350        continue
            wkg(gi)=wkg(gi)+sum
 351     continue
 352  continue
C     *** solve the linear system, storing the result in wkg ******
C     **** wkgg is converted to its inverse square root ***********
      call chfce(g+1,g+1,wkgg,err)
      if(err.eq.1) then
         msg=11
         goto 999
      endif
      call bkslv(g+1,g+1,wkgg)
      do 360 i=g,0,-1
         sum=dble(0.)
         do 359 j=0,i
            sum=sum+wkgg(j,i)*wkg(j)
 359     continue
         wkg(i)=sum
 360  continue
      do 370 i=0,g
         sum=dble(0.)
         do 369 j=i,g
            sum=sum+wkgg(i,j)*wkg(j)
 369     continue
         wkg(i)=sum
 370  continue
C     *** invert wkg, putting the result into xi **************
C     ****** and get new sigma2 *******************************
C     *** step-halving is used here if wkg is not pos.def. ****
C     ****** or if sigma2 is negative *************************
      deflate=dble(1.)
      oldtau=dble(1.)/osigma2
 425  continue
      tau=oldtau+deflate*(wkg(0)-oldtau)
      gi=0
      do 430 i=1,q
         do 429 j=i,q
            gi=gi+1
            wkqq2(i,j)=wkqq1(i,j)+deflate*(wkg(gi)-wkqq1(i,j))
 429     continue
 430  continue
      call chfce(q,q,wkqq2,err)
      if((err.eq.1).or.(tau.le.dble(0))) then
         deflate=deflate/dfloat(2)
         goto 425
      endif
      sigma2=dble(1.)/tau
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,xi)
      do 440 i=1,q
         do 439 j=i+1,q
            xi(j,i)=xi(i,j)
 439     continue
 440  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine prefstrm(ntot,subj,m,ist,ifin,occ,nmax,vmax,wknnm,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,p,xcol,ztvinvx,
     /     iflag,ldv,err)
C Preliminary manipulations for FAST-RML. After execution, 
C     ist  = vector of starting positions for subject i, i=1,...,m
C     ifin = vector of finishing positions for subject i, i=1,...,m
C     vinv   = inverse of V_i, i=1,...,m
C     ztvinv  = t(Z_i)%*% inverse of V_i, i=1,...,m
C     ztvinvz = t(Z_i)%*%inv(V_i)%*%(Z_i), i=1,...,m
C     ztvinvx = t(Z_i)%*%inv(V_i)%*%(X_i), i=1,...,m
C     ldv = .5*sum of log det(V_i)
C Needs a workspace wknnm.
C Note: if V_i's are all identity (indicated by iflag=1) then 
C vmax, wknnm, vinv, are ignored and occ is filled in with 1,...,ni.
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,pcol,
     /     q,zcol(q),p,xcol(p),iflag,err
      double precision vmax(nmax,nmax),wknnm(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),ztvinv(q,nmax,m),
     /     ztvinvz(q,q,m),ztvinvx(q,p,m),ldv
      call istfin(ntot,subj,m,ist,ifin)
      if(iflag.ne.1) then
         call mkv(m,nmax,vmax,ntot,occ,ist,ifin,wknnm)
         call chv(nmax,m,wknnm,ntot,occ,ist,ifin,ldv,err)
         if(err.eq.1) goto 999
         call bkv(nmax,m,wknnm,ntot,occ,ist,ifin)
         call mmulv(nmax,m,wknnm,vinv,ntot,occ,ist,ifin)
      else
         call mkocc(ntot,occ,m,ist,ifin)
         ldv=dble(0.)
      endif
      call mmu(ntot,pcol,pred,q,zcol,nmax,m,wknnm,occ,ist,ifin,ztvinv,
     /     iflag)
      call mmtm(q,nmax,m,ztvinv,ntot,occ,ist,ifin,ztvinvz)
      if(iflag.ne.1) then
         call mml(ntot,q,nmax,m,wknnm,occ,ist,ifin,ztvinv)
      endif
      call mkztvix(ntot,q,nmax,m,occ,ist,ifin,ztvinv,pcol,pred,
     /     p,xcol,ztvinvx)
C     *** fill in elements of ztvinvz below the diagonal ***
      call bdiag(q,m,ztvinvz)
 999  continue
      return
      end
C***********************************************************************
      function trajaj(q,a,j,k,l,m)
C Calculates trace of A%*%J_jk%*%A%*%J_lm, where A is a symmetric matrix
C and J_jk is the matrix with ones in positions (j,k) and (k,j) 
C and zeroes elsewhere.
C Note that A must be filled in above and below the diagonal.
      integer q,j,k,l,m
      double precision a(q,q),trajaj
      trajaj=dble(2.)*(a(j,l)*a(k,m)+a(k,l)*a(j,m))
      return
      end
C***********************************************************************
      function trahaj(q,a,i,j,k)
C Calculates trace of A%*%H_i%*%A%*%J_jk, where A is a symmetric matrix,
C H_i is the matrix with a one in position (i,i) and zeroes elsewhere, 
C and J_jk is the matrix with  ones in positions (j,k) and (k,j) and 
C zeroes elsewhere.
C Note that A must be filled in above and below the diagonal.
      integer q,i,j,k
      double precision a(q,q),trahaj
      trahaj=dble(2.)*a(i,j)*a(i,k)
      return
      end
C***********************************************************************
      function trahah(q,a,i,j)
C Calculates trace of A%*%H_i%*%A%*%H_j, where A is a symmetric matrix,
C and H_i is the matrix with a one in position (i,i) and zeroes 
C elsewhere.
C Note that A must be filled in above and below the diagonal.
      integer q,i,j
      double precision a(q,q),trahah
      trahah=a(i,j)*a(i,j)
      return
      end
C***********************************************************************
C The following subroutines are used by ecmerml().
C***********************************************************************
C Algorithm ECME-RML: Finds RML estimates when V_i are known.
      subroutine ecmerml(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,vinv,
     /     pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,err,msg,u,iter,
     /     sflag,sigma2,p,xcol,beta,y,delta,xtw,xtwx,xtwy,
     /     xtwxinv,wkqq1,wkqq2,xi,wkqnm,b,cvgd,obeta,oxi,maxits,
     /     llvec,eps,ztvinvx,a,wkqp)
C set sflag=1 if starting values supplied.
C Returned error codes in msg:
C     0 = no errors
C     1 = supplied V_i matrix is not positive definite
C     2 = GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank
C     3 = inadequate information to obtain starting value of xi
C     4 = value of xi or inv(Ui) became non-pos.def. during iterations
C     5 = t(X)%*%W%*%X became non-pos.def. during iterations
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,
     /     pcol,q,zcol(q),iflag,err,msg,iter,sflag,p,xcol(p),cvgd,
     /     maxits,c1,c2,c3,i,j
      double precision vmax(nmax,nmax),w(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),
     /     ztvinv(q,nmax,m),ztvinvz(q,q,m),ldv,u(q,q,m),sigma2,
     /     beta(p),y(ntot),delta(ntot),xtw(p,nmax),xtwx(p,p),xtwy(p),
     /     xtwxinv(p,p),wkqq1(q,q),wkqq2(q,q),xi(q,q),
     /     wkqnm(q,nmax,m),b(q,m),obeta(p),oxi(q,q),osigma2,ldu,ldxi,
     /     ll,llvec(maxits),eps,ztvinvx(q,p,m),a(q,q,m),wkqp(q,p),
     /     ldxtwx
      msg=0
      iter=0
      call preecme2(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,p,xcol,ztvinvx,
     /     iflag,ldv,err)
      if(err.eq.1) then
         msg=1
         goto 999
      endif
      if(sflag.ne.1) then
         call stval1(ntot,m,ist,ifin,occ,nmax,vinv,pcol,pred,
     /        q,ztvinv,ztvinvz,iflag,err,msg,sigma2,p,xcol,
     /        beta,y,delta,xtw,xtwx,xtwy,xtwxinv,wkqq1,wkqq2,xi,
     /        wkqnm,b)
      endif
      cvgd=0
C********* start of main iteration *****************
 1    continue
         iter=iter+1
C        *** save current estimates ********************
         osigma2=sigma2
         do 2 i=1,p
            obeta(i)=beta(i)
 2       continue
         do 10 i=1,q
            do 5 j=i,q
               oxi(i,j)=xi(i,j)
 5          continue
 10      continue
C        *** calculate U_i and W_i *********************
         call mku(q,xi,m,ztvinvz,u,wkqq1,wkqq2,ldxi,ldu,err)
         if(err.eq.1) then
            msg=4
            goto 999
         endif
         call mkwkqnm(q,m,u,nmax,ztvinv,wkqnm,ntot,occ,ist,ifin)
         call mkw(q,nmax,m,ist,ifin,wkqnm,ztvinv,vinv,w,ntot,occ,
     /        iflag)
C        *** calculate beta and sigma2 *******************
         call gls(ntot,m,ist,ifin,occ,nmax,pcol,pred,0,sigma2,
     /        p,xcol,beta,y,delta,w,xtw,xtwx,xtwy,xtwxinv,err)
         if(err.eq.1) then
            msg=5
            goto 999
         endif
         sigma2=sigma2*dfloat(ntot)/dfloat(ntot-p)
C        *** evaluate loglikelihood ************************
C        *** note that xtwx contains the square root of xtwxinv ***
         ldxtwx=dble(0.)
         do 15 i=1,p
            ldxtwx=ldxtwx+dlog(xtwx(i,i))
 15      continue
         ll = -dble(.5)*dfloat(ntot-p)*dlog(osigma2)
     /        + dfloat(m)*ldxi + ldu + ldxtwx
     /        -.5*dfloat(ntot-p)*sigma2/osigma2
         llvec(iter)=ll
C        *** calculate b_i, A_i, and xi ********************
         call mkb(q,nmax,m,wkqnm,ntot,delta,b,occ,ist,ifin)
         call bdiag(q,m,u)
         call mka(q,m,p,u,ztvinvx,xtwxinv,xtwx,wkqp,a)
         call mkxi2(q,m,b,u,a,xi,osigma2)
         c1=0
         do 30 i=1,p
            if(dabs(beta(i)-obeta(i)).gt.(eps*dabs(obeta(i)))) c1=1
 30      continue
         c2=0
         do 40 i=1,q
            do 35 j=i,q
               if(dabs(xi(i,j)-oxi(i,j)).gt.(eps*dabs(oxi(i,j))))
     /              c2=1
 35         continue
 40      continue
         c3=0
         if(dabs(sigma2-osigma2).gt.(eps*dabs(osigma2))) c3=1
         if((c1.eq.0).and.(c2.eq.0).and.(c3.eq.0)) cvgd=1
         if((cvgd.eq.0).and.(iter.lt.maxits)) goto 1
C********* end of main iteration *****************
      do 70 i=1,p
         do 60 j=i+1,p
            xtwxinv(j,i)=xtwxinv(i,j)
 60      continue
 70   continue
 999  continue
      return
      end
C***********************************************************************
      subroutine mka(q,m,p,u,ztvinvx,xtwxinv,wkpp,wkqp,a)
C calculates A_i, i=1,...m from u, ztvinvx and xtwxinv.
C Assumes that elements of u below the diagonal have been filled in.
      integer q,m,p,s,i,j,k,err
      double precision u(q,q,m),ztvinvx(q,p,m),xtwxinv(p,p),
     /     wkpp(p,p),wkqp(q,p),a(q,q,m),sum
C     *** store cholesky of xtwxinv in wkpp ***
      do 5 i=1,p
         do 4 j=i,p
            wkpp(i,j)=xtwxinv(i,j)
 4       continue
 5    continue
      call chfce(p,p,wkpp,err)
      do 300 s=1,m
C        *** store U_i%*%t(Z_i)%*%inv(V_i)%*%X_i in wkqp ***
         do 40 i=1,q
            do 35 j=1,p
               sum=dble(0.)
               do 32 k=1,q
                  sum=sum+u(i,k,s)*ztvinvx(k,j,s)
 32            continue
               wkqp(i,j)=sum
 35         continue
 40      continue
C        *** multiply wkqp by lower cholesky of xtwxinv *****
         do 60 i=1,q
            do 55 j=1,p
               sum=dble(0.)
               do 52 k=j,p
                  sum=sum+wkqp(i,k)*wkpp(j,k)
 52            continue
               wkqp(i,j)=sum
 55         continue
 60      continue
C        *** multiply wkqp by its transpose, using symmetry *****
         do 80 i=1,q
            do 70 j=i,q
               sum=dble(0.)
               do 68 k=1,p
                  sum=sum+wkqp(i,k)*wkqp(j,k)
 68            continue
               a(i,j,s)=sum
               if(i.ne.j) a(j,i,s)=sum
 70         continue
 80      continue
 300  continue
      return
      end
C***********************************************************************
      subroutine mkxi2(q,m,b,u,a,xi,sigma2)
C calculates new estimate of xi from b, u, a, and sigma2
      integer q,m,s,i,j
      double precision b(q,m),u(q,q,m),a(q,q,m),xi(q,q),sigma2
      do 10 i=1,q
         do 5 j=i,q
            xi(i,j)=dble(0.)
 5       continue
 10   continue
      do 100 s=1,m
         do 90 i=1,q
            do 80 j=i,q
               xi(i,j)=xi(i,j)+u(i,j,s)+a(i,j,s)+
     /              b(i,s)*b(j,s)/sigma2
 80         continue
 90      continue
 100  continue
      do 110 i=1,q
         do 105 j=i,q
            xi(i,j)=xi(i,j)/dfloat(m)
            if(i.ne.j) xi(j,i)=xi(i,j)
 105     continue
 110  continue
      return
      end
C***********************************************************************
      subroutine mkb(q,nmax,m,wkqnm,ntot,delta,b,occ,ist,ifin)
C calculates bwig from wkqnm and delta
      integer q,nmax,m,ntot,occ(ntot),ist(m),ifin(m),s,st,fin,i,j
      double precision wkqnm(q,nmax,m),delta(ntot),b(q,m),sum
      do 500 s=1,m
         do 400 i=1,q
            st=ist(s)
            fin=ifin(s)
            sum=dble(0.)
            do 300 j=st,fin
               sum=sum+wkqnm(i,occ(j),s)*delta(j)
 300           continue
            b(i,s)=sum
 400     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine mkw(q,nmax,m,ist,ifin,wkqnm,ztvinv,vinv,w,ntot,occ,
     /     iflag)
C Given wkqnm, calculates (upper-tri part of) weight matrices W_i
      integer q,nmax,m,s,ntot,occ(ntot),st,fin,iflag,ist(m),ifin(m),
     /     i,j,k
      double precision wkqnm(q,nmax,m),ztvinv(q,nmax,m),
     /     vinv(nmax,nmax,m),w(nmax,nmax,m),sum
      if(iflag.ne.1) then
         do 350 s=1,m
            st=ist(s)
            fin=ifin(s)            
            do 300 i=st,fin
               do 200 j=i,fin
                  sum=dble(0.)
                  do 100 k=1,q
                     sum=sum+ztvinv(k,occ(i),s)*wkqnm(k,occ(j),s)
 100              continue
                  w(occ(i),occ(j),s)=(vinv(occ(i),occ(j),s)-sum)
 200           continue
 300        continue
 350     continue
      else
         do 650 s=1,m
            st=ist(s)
            fin=ifin(s)            
            do 600 i=st,fin
               do 500 j=i,fin
                  sum=dble(0.)
                  do 400 k=1,q
                     sum=sum+ztvinv(k,occ(i),s)*wkqnm(k,occ(j),s)
 400              continue
                  if(i.eq.j) then
                     w(occ(i),occ(j),s)=(dble(1.)-sum)
                  else
                     w(occ(i),occ(j),s)=-sum
                  endif
 500           continue
 600        continue
 650     continue
      endif
      return
      end
C***********************************************************************
      subroutine mkwkqnm(q,m,u,nmax,ztvinv,wkqnm,ntot,occ,ist,ifin)
C Given Ui, i=1,...,m calculates Ui%*%t(Zi)%*%inv(Vi)
      integer q,m,nmax,ntot,occ(ntot),ist(m),ifin(m),st,fin,s,i,j,k
      double precision u(q,q,m),ztvinv(q,nmax,m),wkqnm(q,nmax,m),sum
      do 400 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 300 i=1,q
            do 200 j=st,fin
               sum=dble(0.)
               do 50 k=1,i-1
                  sum=sum+u(k,i,s)*ztvinv(k,occ(j),s)
 50            continue
               do 100 k=i,q
                  sum=sum+u(i,k,s)*ztvinv(k,occ(j),s)
 100           continue
               wkqnm(i,occ(j),s)=sum
 200        continue
 300     continue
 400  continue
      return
      end
C***********************************************************************
      subroutine mku(q,xi,m,ztvinvz,u,wkqq1,wkqq2,ldxi,ldu,err)
C Calculates U_i, i=1,...,m from xi.
C After execution, wkqq1 contains the inverse of the original
C xi (upper-tri part only).
C Also calculates ldxi = -.5*logdet(xi) and
C                 ldu = +.5*sum( logdet(U_i) ), i=1,...,m
C Sets err=1 if the input value of xi is not positive definite, or if
C   the value of any (t(Zi)%*%inv(Vi)%*%Zi + inv(xi)) is not pos. def.
      integer q,m,s,err,i,j
      double precision xi(q,q),ztvinvz(q,q,m),u(q,q,m),
     /     wkqq1(q,q),wkqq2(q,q),ldxi,ldu
      err=0
C *** put the inverse of xi into wkqq1
      do 2 i=1,q
         do 1 j=i,q
            wkqq2(i,j)=xi(i,j)
 1       continue
 2    continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) goto 999
      call bkslv(q,q,wkqq2)
      ldxi=dble(0.)
      do 10 i=1,q
         ldxi=ldxi+dlog(wkqq2(i,i))
 10   continue
      call mm(q,q,wkqq2,wkqq1)
      ldu=dble(0.)
      do 500 s=1,m
         do 100 i=1,q
            do 80 j=i,q
               u(i,j,s)=wkqq1(i,j)+ztvinvz(i,j,s)
 80         continue
 100     continue
         call chle(q,q,m,u,s,err)
         call bkslvl(q,q,m,u,s)
         do 110 i=1,q
            ldu=ldu+dlog(u(i,i,s))
 110     continue
         call mmul(q,q,m,u,s,wkqq2)
         do 200 i=1,q
            do 150 j=i,q
               u(i,j,s)=wkqq2(i,j)
 150        continue
 200     continue
 500  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine stval1(ntot,m,ist,ifin,occ,nmax,vinv,pcol,
     /     pred,q,ztvinv,ztvinvz,iflag,err,msg,sigma2,p,xcol,
     /     beta,y,delta,xtw,xtwx,xtwy,xtwxinv,wkqq1,wkqq2,xi,
     /     wkqnm,b)
C Starting values for lmm.
C Sets msg=2 if t(X)%*%W%*%X is not full rank.
C Sets msg=3 if there are no individual subject-level regressions from
C which to estimate starting value for xi.
      integer ntot,m,ist(m),ifin(m),occ(ntot),nmax,pcol,q,
     /     iflag,err,msg,p,xcol(p),s,st,fin,mstar,err1,i,j,k
      double precision vinv(nmax,nmax,m),pred(ntot,pcol),
     /     ztvinv(q,nmax,m),ztvinvz(q,q,m),sigma2,beta(p),
     /     y(ntot),delta(ntot),xtw(p,nmax),xtwx(p,p),xtwy(p),
     /     xtwxinv(p,p),wkqq1(q,q),wkqq2(q,q),xi(q,q),sum,
     /     wkqnm(q,nmax,m),b(q,m)
      err=0
C get GLS estimates with weights w=vinv
      call gls(ntot,m,ist,ifin,occ,nmax,pcol,pred,iflag,sigma2,
     /     p,xcol,beta,y,delta,vinv,xtw,xtwx,xtwy,xtwxinv,err)
      if(err.eq.1) then
         msg=2
         goto 999
      endif
C now perform unit-level regressions where possible
      mstar=0
      do 10 i=1,q
         do 5 j=i,q
            xi(i,j)=dble(0.)
 5       continue
 10   continue
      do 300 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 20 i=1,q
            do 15 j=i,q
               wkqq1(i,j)=ztvinvz(i,j,s)
 15         continue
 20      continue
         call chfce(q,q,wkqq1,err1)
         if(err1.eq.1) goto 299
         mstar=mstar+1
         call bkslv(q,q,wkqq1)
         call mm(q,q,wkqq1,wkqq2)
         do 100 i=1,q
            do 80 j=st,fin
               sum=dble(0.)
               do 40 k=1,i-1
                  sum=sum+wkqq2(k,i)*ztvinv(k,occ(j),s)
 40            continue
               do 50 k=i,q
                  sum=sum+wkqq2(i,k)*ztvinv(k,occ(j),s)
 50            continue
               wkqnm(i,occ(j),s)=sum
 80         continue
 100     continue
         do 150 i=1,q
            sum=dble(0.)
            do 120 j=st,fin
               sum=sum+wkqnm(i,occ(j),s)*delta(j)
 120           continue
            b(i,s)=sum
 150     continue
         do 200 i=1,q
            do 180 j=i,q
               xi(i,j)=xi(i,j)+b(i,s)*b(j,s)
 180        continue
 200     continue
 299     continue
 300  continue
      if(mstar.eq.0) then
         err=1
         msg=3
         goto 999
      endif
      do 350 i=1,q
         do 330 j=i,q
            xi(i,j)=xi(i,j)/(dfloat(mstar)*sigma2)
            if(i.ne.j) xi(j,i)=xi(i,j)
 330     continue
 350  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine gls(ntot,m,ist,ifin,occ,nmax,pcol,pred,iflag,sigma2,
     /     p,xcol,beta,y,delta,w,xtw,xtwx,xtwy,xtwxinv,err)
C Gets GLS estimates for beta and sigma2, using weights in w. Also
C calculates deltai = yi - Xi%*%beta, i=1,...,m.
C If iflag=1 then weights in w are ignored and OLS is used.
      integer ntot,m,ist(m),ifin(m),occ(ntot),nmax,pcol,
     /     iflag,p,xcol(p),st,fin,s,i,j,err
      double precision pred(ntot,pcol),sigma2,beta(p),y(ntot),
     /     delta(ntot),w(nmax,nmax,m),xtw(p,nmax),xtwx(p,p),xtwy(p),
     /     xtwxinv(p,p),sum
      err=0
      do 10 i=1,p
         xtwy(i)=dble(0.)
         do 5 j=i,p
            xtwx(i,j)=dble(0.)
 5       continue
 10   continue
      do 100 s=1,m
         st=ist(s)
         fin=ifin(s)
         call mkxtw(ntot,pcol,pred,p,xcol,occ,st,fin,nmax,w,xtw,s,
     /        m,iflag)
         call mkxtwx(ntot,pcol,pred,p,xcol,occ,st,fin,nmax,xtw,xtwx)
         call mkxtwy(ntot,p,occ,st,fin,nmax,xtw,y,xtwy)
 100  continue
      call chfce(p,p,xtwx,err)
      if(err.eq.1) goto 999
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
      call mkdel(ntot,pcol,pred,p,xcol,y,beta,delta)
      sigma2=dble(0.)
      do 400 s=1,m
         st=ist(s)
         fin=ifin(s)
         if(iflag.ne.1) then
            do 300 i=st,fin
               sum=dble(0.)
               do 150 j=st,i
                  sum=sum+delta(j)*w(occ(j),occ(i),s)
 150           continue
               do 200 j=i+1,fin
                  sum=sum+delta(j)*w(occ(i),occ(j),s)
 200           continue
               sigma2=sigma2+sum*delta(i)
 300        continue
         else
            do 350 i=st,fin
               sigma2=sigma2+delta(i)**2
 350        continue
         endif
 400  continue
      sigma2=sigma2/dfloat(ntot)
 999  continue
      return
      end
C***********************************************************************
      subroutine mkdel(ntot,pcol,pred,p,xcol,y,beta,delta)
      integer ntot,p,xcol(p),pcol,i,j
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
      subroutine mkxtwy(ntot,p,occ,st,fin,nmax,xtw,y,xtwy)
C increments xtwy for subject s
      integer ntot,p,occ(ntot),st,fin,nmax,i,j
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
      integer ntot,pcol,p,xcol(p),occ(ntot),st,fin,nmax,i,j,k
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
      subroutine mkxtw(ntot,pcol,pred,p,xcol,occ,st,fin,nmax,w,xtw,s,
     /     m,iflag)
C Calculates t(x)%*%w for subject s.
C If iflag=1 then w is ignored.
      integer ntot,pcol,p,xcol(p),occ(ntot),st,fin,nmax,s,m,iflag,
     /     i,j,k
      double precision pred(ntot,pcol),w(nmax,nmax,m),xtw(p,nmax),sum
      if(iflag.ne.1) then
         do 500 i=1,p
            do 400 j=st,fin
               sum=dble(0.)
               do 200 k=st,j
                  sum=sum+pred(k,xcol(i))*w(occ(k),occ(j),s)
 200           continue
               do 300 k=j+1,fin
                  sum=sum+pred(k,xcol(i))*w(occ(j),occ(k),s)
 300           continue
               xtw(i,occ(j))=sum
 400        continue
 500     continue
      else
         do 700 i=1,p
            do 600 j=st,fin
               xtw(i,occ(j))=pred(j,xcol(i))
 600        continue
 700     continue
      endif
      return
      end
C***********************************************************************
      subroutine preecme2(ntot,subj,m,ist,ifin,occ,nmax,vmax,wknnm,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,p,xcol,ztvinvx,
     /     iflag,ldv,err)
C Preliminary manipulations for ECME-RML. After execution, 
C     ist  = vector of starting positions for subject i, i=1,...,m
C     ifin = vector of finishing positions for subject i, i=1,...,m
C     vinv   = inverse of V_i, i=1,...,m
C     ztvinv  = t(Z_i)%*% inverse of V_i, i=1,...,m
C     ztvinvz = t(Z_i)%*%inv(V_i)%*%(Z_i), i=1,...,m
C     ztvinvx = t(Z_i)%*%inv(V_i)%*%(X_i), i=1,...,m
C     ldv = .5*sum of log det(V_i)
C Needs a workspace wknnm.
C Note: if V_i's are all identity (indicated by iflag=1) then 
C vmax, wknnm, vinv, are ignored and occ is filled in with 1,...,ni.
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,pcol,
     /     q,zcol(q),p,xcol(p),iflag,err
      double precision vmax(nmax,nmax),wknnm(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),ztvinv(q,nmax,m),
     /     ztvinvz(q,q,m),ztvinvx(q,p,m),ldv
      call istfin(ntot,subj,m,ist,ifin)
      if(iflag.ne.1) then
         call mkv(m,nmax,vmax,ntot,occ,ist,ifin,wknnm)
         call chv(nmax,m,wknnm,ntot,occ,ist,ifin,ldv,err)
         if(err.eq.1) goto 999
         call bkv(nmax,m,wknnm,ntot,occ,ist,ifin)
         call mmulv(nmax,m,wknnm,vinv,ntot,occ,ist,ifin)
      else
         call mkocc(ntot,occ,m,ist,ifin)
         ldv=dble(0.)
      endif
      call mmu(ntot,pcol,pred,q,zcol,nmax,m,wknnm,occ,ist,ifin,ztvinv,
     /     iflag)
      call mmtm(q,nmax,m,ztvinv,ntot,occ,ist,ifin,ztvinvz)
      if(iflag.ne.1) then
         call mml(ntot,q,nmax,m,wknnm,occ,ist,ifin,ztvinv)
      endif
      call mkztvix(ntot,q,nmax,m,occ,ist,ifin,ztvinv,pcol,pred,
     /     p,xcol,ztvinvx)
 999  continue
      return
      end
C***********************************************************************
      subroutine mkv(m,nmax,vmax,ntot,occ,ist,ifin,v)
C Stores V_i, the submatrices of vmax, in v
      integer m,nmax,ntot,occ(ntot),ist(m),ifin(m),s,i,j
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
      subroutine chv(nmax,m,v,ntot,occ,ist,ifin,ldv,err)
C Overwrites v = (V_i, i=1,...,m) with its cholesky factors.
C Also calculates ldv = .5*sum of log-determinants of V_i. 
C If any of the V_i's is not positive definite, err is set to 1.
      integer nmax,m,ntot,occ(ntot),ist(m),ifin(m),s,err,i,j,k
      double precision v(nmax,nmax,m),sum,ldv
      err=0
      ldv=dble(0.)
      do 100 s=1,m
         do 50 i=ist(s),ifin(s)
            sum=dble(0.)
            do 10 k=ist(s),i-1
               sum=sum+v(occ(k),occ(i),s)**2
 10         continue
            if(sum.ge.v(occ(i),occ(i),s)) then
               err=1
               goto 199
            endif
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
         do 60 i=ist(s),ifin(s)
            ldv=ldv+dlog(v(occ(i),occ(i),s))
 60      continue
 100  continue
 199  continue
      return
      end
C***********************************************************************
      subroutine bkv(nmax,m,v,ntot,occ,ist,ifin)
C inverts upper-triangular portions of v by backsolve
      integer nmax,m,ntot,occ(ntot),ist(m),ifin(m),s,st,fin,i,j,k
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
      subroutine mmulv(nmax,m,vinvh,vinv,ntot,occ,ist,ifin)
C calculates upper tri part of vinv = vinvh%*%t(vinvh), where 
C vinvh is upper tri
C(works in layers s=1,...,m)
      integer nmax,m,ntot,occ(ntot),ist(m),ifin(m),s,st,fin,i,j,k
      double precision vinvh(nmax,nmax,m),vinv(nmax,nmax,m),sum
      do 300 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 200 i=st,fin
            do 100 j=i,fin
               sum=dble(0.)
               do 50 k=max(i,j),fin
                  sum=sum+vinvh(occ(i),occ(k),s)*vinvh(occ(j),occ(k),s)
 50            continue
               vinv(occ(i),occ(j),s)=sum
 100        continue
 200        continue
 300  continue
      return
      end
C***********************************************************************
      subroutine mmu(ntot,pcol,pred,q,zcol,nmax,m,v,occ,ist,ifin,
     /     ztvinv,iflag)
C multiplies z^t by upper-triangular v, stores result in ztvinv
      integer ntot,pcol,q,zcol(q),nmax,m,occ(ntot),ist(m),ifin(m),s,
     /     st,fin,iflag,i,j,k
      double precision pred(ntot,pcol),v(nmax,nmax,m),ztvinv(q,nmax,m),
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
                  ztvinv(i,occ(j),s)=sum
 100           continue
 200        continue
 300     continue
      else
         do 600 s=1,m
            st=ist(s)
            fin=ifin(s)
            do 500 i=1,q
               do 400 j=st,fin
                  ztvinv(i,occ(j),s)=pred(j,zcol(i))
 400           continue
 500        continue
 600     continue
      endif
      return
      end
C***********************************************************************
      subroutine mmtm(q,nmax,m,ztvinv,ntot,occ,ist,ifin,ztvinvz)
C multiply ztvinv by its transpose, store result (upper-tri part)
C in ztvinvz
      integer q,nmax,m,ntot,occ(ntot),ist(m),ifin(m),s,i,j,k
      double precision ztvinv(q,nmax,m),ztvinvz(q,q,m),sum
      do 500 s=1,m
         do 400 i=1,q
            do 300 j=i,q
               sum=dble(0.)
               do 200 k=ist(s),ifin(s)
                  sum=sum+ztvinv(i,occ(k),s)*ztvinv(j,occ(k),s)
 200           continue
               ztvinvz(i,j,s)=sum
 300        continue
 400     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine mkztvix(ntot,q,nmax,m,occ,ist,ifin,ztvinv,pcol,pred,
     /     p,xcol,ztvinvx)
C multiplies ztvinv by x, stores result in ztvinvx
      integer ntot,q,nmax,m,occ(ntot),ist(m),ifin(m),s,st,fin,i,j,k,
     /     pcol,p,xcol(p)
      double precision ztvinv(q,nmax,m),pred(ntot,pcol),
     /     ztvinvx(q,p,m),sum
      do 300 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 200 i=1,q
            do 100 j=1,p
               sum=dble(0.)
               do 50 k=st,fin
                  sum=sum+ztvinv(i,occ(k),s)*pred(k,xcol(j))
 50            continue
               ztvinvx(i,j,s)=sum
 100        continue
 200     continue
 300  continue
      return
      end
C***********************************************************************
      subroutine mml(ntot,q,nmax,m,v,occ,ist,ifin,ztvinv)
C multiplies ztvinv by lower-triangular v, stores result in ztvinv
      integer ntot,q,nmax,m,occ(ntot),ist(m),ifin(m),s,st,fin,i,j,k
      double precision v(nmax,nmax,m),ztvinv(q,nmax,m),sum
      do 300 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 200 i=1,q
            do 100 j=st,fin
               sum=dble(0.)
               do 50 k=j,fin
                  sum=sum+ztvinv(i,occ(k),s)*v(occ(j),occ(k),s)
 50            continue
               ztvinv(i,occ(j),s)=sum
 100        continue
 200     continue
 300  continue
      return
      end
C***********************************************************************
      subroutine chfce(p,pw,s,err)
C Overwrites s (upper tri) with cholesky factor
c If s is not positive definite, err is set to one.
      integer p,pw,err,i,j,k
      double precision s(p,p),sum
      err=0
      do 50 i=1,pw
         sum=dble(0.)
         do 10 k=1,i-1
            sum=sum+s(k,i)**2
 10      continue
         if(sum.ge.s(i,i)) then
            err=1
            goto 999
         endif
         s(i,i)=dsqrt(s(i,i)-sum)
         do 40 j=i+1,pw
            sum=dble(0.)
            do 30 k=1,i-1
               sum=sum+s(k,i)*s(k,j)
 30         continue
            s(i,j)=(s(i,j)-sum)/s(i,i)
 40      continue
 50   continue
 999  continue
      return
      end
C***********************************************************************
      subroutine chle(p,pw,m,s,l,err)
C overwrites lth layer of s (upper tri) with cholesky factor
C If it fails, err is set to one.
      integer p,pw,m,l,err,i,j,k
      double precision s(p,p,m),sum
      err=0
      do 50 i=1,pw
         sum=dble(0.)
         do 10 k=1,i-1
            sum=sum+s(k,i,l)**2
 10      continue
         if(sum.ge.s(i,i,l)) then
            err=1
            goto 999
         endif
         s(i,i,l)=dsqrt(s(i,i,l)-sum)
         do 40 j=i+1,pw
            sum=dble(0.)
            do 30 k=1,i-1
               sum=sum+s(k,i,l)*s(k,j,l)
 30         continue
            s(i,j,l)=(s(i,j,l)-sum)/s(i,i,l)
 40      continue
 50   continue
 999  continue
      return
      end
C***********************************************************************
      subroutine bdiag(q,m,sig)
C fills in elements of each layer of sig below the diagonal
      integer q,m,s,i,j
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
C Algorithm FAST-MODE: modification of FAST-RML to find posterior
C mode of variance components.
      subroutine fastmode(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,vinv,
     /     pcol,pred,q,zcol,ztvinv,ztvinvz,isflags,err,msg,u,iter,
     /     sigma2,p,xcol,beta,y,delta,xtw,xtwx,xtwy,xtwxinv,
     /     wkqq1,wkqq2,xi,wkqnm,b,cvgd,oxi,maxits,llvec,eps,
     /     xiecme,g,reject,ztvinvx,a,wkqp,wkg,wkgg,abc,dinv)
C set sflag=1 if starting values supplied.
C Returned error codes in msg:
C     0 = no errors
C     1 = supplied V_i matrix is not positive definite
C     2 = GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank
C     3 = inadequate information to obtain starting value of xi
C     4 = value of xi or inv(Ui) became non-pos.def. during iterations
C     5 = t(X)%*%W%*%X became non-pos.def. during iterations
C    10 = non-positive definite xi at input to scoring step
C    11 = non-positive definite wkgg in scoring step
C
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,
     /     pcol,q,zcol(q),iflag,err,msg,iter,sflag,p,xcol(p),cvgd,
     /     maxits,c2,c3,i,j,g,reject(maxits),notconcave,isflags(2)
      double precision vmax(nmax,nmax),w(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),
     /     ztvinv(q,nmax,m),ztvinvz(q,q,m),ldv,u(q,q,m),sigma2,
     /     beta(p),y(ntot),delta(ntot),xtw(p,nmax),xtwx(p,p),xtwy(p),
     /     xtwxinv(p,p),wkqq1(q,q),wkqq2(q,q),xi(q,q),
     /     wkqnm(q,nmax,m),b(q,m),oxi(q,q),osigma2,ldu,ldxi,
     /     ll,llvec(maxits),eps,xiecme(q,q),
     /     ztvinvx(q,p,m),a(q,q,m),wkqp(q,p),ldxtwx,
     /     sig2ecme,osig2ecme,wkg(0:g),wkgg(0:g,0:g),dinv(q,q),
     /     abc(3),trdx
      iflag=isflags(1)
      sflag=isflags(2)
      msg=0
      iter=0
      notconcave=0
      call prefstrm(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,p,xcol,ztvinvx,
     /     iflag,ldv,err)
      if(err.eq.1) then
         msg=1
         goto 999
      endif
      if(sflag.ne.1) then
         call stval1(ntot,m,ist,ifin,occ,nmax,vinv,pcol,pred,
     /        q,ztvinv,ztvinvz,iflag,err,msg,sigma2,p,xcol,
     /        beta,y,delta,xtw,xtwx,xtwy,xtwxinv,wkqq1,wkqq2,xi,
     /        wkqnm,b)
      endif
      cvgd=0
C********* start of main iteration *****************
 1    continue
         iter=iter+1
         reject(iter)=0
 3       continue
         osigma2=sigma2
C        *** calculate U_i, U_i%*%t(Z_i)%*%inv(V_i), and W_i ****
         call mku(q,xi,m,ztvinvz,u,wkqq1,wkqq2,ldxi,ldu,err)
         if(err.eq.1) then
            msg=4
            goto 999
         endif
         call trdixi(q,trdx,wkqq1,dinv)
         call mkwkqnm(q,m,u,nmax,ztvinv,wkqnm,ntot,occ,ist,ifin)
         call mkw(q,nmax,m,ist,ifin,wkqnm,ztvinv,vinv,w,ntot,occ,
     /        iflag)
C        *** save current value of sigma2ecme ******************
         osig2ecme=sig2ecme
C        *** calculate beta, delta_i, and sigma2ecme ***********
         call gls(ntot,m,ist,ifin,occ,nmax,pcol,pred,0,sig2ecme,
     /        p,xcol,beta,y,delta,w,xtw,xtwx,xtwy,xtwxinv,err)
         if(err.eq.1) then
            msg=5
            goto 999
         endif
         sig2ecme= (abc(1) + sig2ecme*dfloat(ntot) + trdx)/
     /        (abc(2)+dfloat(ntot-p-2)+abc(3)*dfloat(q))
C        *** evaluate log-posterior ******************************
C        *** note that xtwx contains the square root of xtwxinv ***
C        ***** first we get its log-determinant *******************
         ldxtwx=dble(0.)
         do 15 i=1,p
            ldxtwx=ldxtwx+dlog(xtwx(i,i))
 15      continue
         ll = -.5*(dfloat(ntot-p-2)+abc(2)+abc(3)*dfloat(q))
     /        *dlog(osigma2) + (dfloat(m-q-1)+abc(3))*ldxi + ldu 
     /        + ldxtwx -.5*(dfloat(ntot-p-2)+abc(2)+abc(3)*dfloat(q))
     /        *sig2ecme/osigma2
         llvec(iter)=ll
         if(iter.gt.1) then
C           *** if log-poterior has gone down, replace the scoring ****
C           ******  estimate of xi by the ecme estimate ****************
            if(reject(iter-1).eq.0) then
               if(llvec(iter).lt.llvec(iter-1)) then
                  sigma2=osig2ecme
                  do 17 i=1,q
                     do 16 j=1,q
                        xi(i,j)=xiecme(i,j)
 16                  continue
 17               continue
                  reject(iter-1)=1
                  goto 3
               endif
            endif
         endif
C        *** check for convergence **********************************
         if(iter.gt.1) then
            c2=0
            do 20 i=1,q
               do 19 j=i,q
                  if(dabs(xi(i,j)-oxi(i,j)).gt.(eps*dabs(oxi(i,j))))
     /                 c2=1
 19            continue
 20         continue
            c3=0
            if(dabs(sigma2-osigma2).gt.(eps*dabs(osigma2))) c3=1
            if((c2.eq.0).and.(c3.eq.0)) then
               cvgd=1
               goto 50
            endif
         endif
C        *** calculate b_i ********************************
         call mkb(q,nmax,m,wkqnm,ntot,delta,b,occ,ist,ifin)
C        *** save current parameters **********************
         osigma2=sigma2
         do 27 i=1,q
            do 26 j=1,q
               oxi(i,j)=xi(i,j)
 26         continue
 27      continue
C        *** perform scoring step, putting result in xi and sigma2,**
C        ****** and put ECME estimate of xi in xiecme ***************
         call fscovm2(m,q,b,u,a,xi,oxi,xiecme,wkqq1,wkqq2,
     /        p,xtwxinv,xtwx,wkqp,g,wkg,wkgg,sigma2,
     /        msg,osigma2,ntot,sig2ecme,ztvinvx,dinv,abc)
         if(msg.eq.10) goto 999
         if(msg.eq.11) then
C        *** wkgg not pos.def., substitute ECME value instead ****
            notconcave=1
            do 29 i=1,q
               do 28 j=1,q
                  xi(i,j)=xiecme(i,j)
 28            continue
 29         continue
            sigma2=sig2ecme
            reject(iter)=1
         endif
         if(iter.lt.maxits) goto 1
C********* end of main iteration *****************
 50   continue
      if(notconcave.eq.1) msg=11
      call bdiag(q,m,u)
      do 70 i=1,p
         do 60 j=i+1,p
            xtwxinv(j,i)=xtwxinv(i,j)
 60      continue
 70   continue
 999  continue
      return
      end
C***********************************************************************
      subroutine fscovm2(m,q,b,u,a,xi,oxi,xiecme,wkqq1,wkqq2,
     /     p,xtwxinv,wkpp,wkqp,g,wkg,wkgg,sigma2,msg,
     /     osigma2,ntot,sig2ecme,ztvinvx,dinv,abc)
C Modification of fscovr2 for posterior mode.
C Like fscovm, but it performs scoring on eta^* rather than eta.
C       oxi = old value of xi
C        xi = new value
C    xiecme = ecme version of xi
C Error messages returned through msg:
C     0 = no error
C    10 = non-positive definite xi
C    11 = non-positive definite wkgg
      integer m,q,p,g,msg,ntot
      double precision b(q,m),u(q,q,m),a(q,q,m),
     /     xi(q,q),oxi(q,q),xiecme(q,q),wkqq1(q,q),wkqq2(q,q),
     /     xtwxinv(p,p),wkpp(p,p),wkqp(q,p),wkg(0:g),
     /     wkgg(0:g,0:g),sigma2,sum,deflate,
     /     osigma2,sig2ecme,ztvinvx(q,p,m),oldtau,tau,
     /     dinv(q,q),abc(3)
      integer s,i,j,k,ii,jj,jjmin,gi,gj,err
      double precision trahah,trahaj,trajaj
      msg=0
C     *** initialize xiecme *****
      do 2 i=1,q
         do 1 j=i,q
            xiecme(i,j)=dinv(i,j)/osigma2
 1       continue
 2    continue
C     *** store cholesky of xtwxinv in wkpp ******
      do 5 i=1,p
         do 4 j=i,p
            wkpp(i,j)=xtwxinv(i,j)
 4       continue
 5    continue
      call chfce(p,p,wkpp,err)
C     *** initialize the workspaces wkg and wkgg *******
      do 7 i=0,g
         wkg(i)=dble(0.)
 7    continue
      gi=0
      do 20 i=1,q
         do 19 j=i,q
            gi=gi+1
            if(i.eq.j) then
               wkgg(0,gi)=dinv(i,i)/osigma2
            else
               wkgg(0,gi)=dble(2.)*dinv(i,j)/osigma2
            endif
            gj=gi-1
            do 18 ii=i,q
               if(ii.eq.i) then
                  jjmin=j
               else
                  jjmin=ii
               endif
               do 17 jj=jjmin,q
                  gj=gj+1
                  if(i.eq.j) then
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=trahah(q,oxi,i,ii)*
     /                       (abc(3)-dfloat(q+1))
                     else
                        wkgg(gi,gj)=trahaj(q,oxi,i,ii,jj)*
     /                       (abc(3)-dfloat(q+1))
                     endif
                  else
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=trahaj(q,oxi,ii,i,j)*
     /                       (abc(3)-dfloat(q+1))
                     else
                        wkgg(gi,gj)=trajaj(q,oxi,i,j,ii,jj)*
     /                       (abc(3)-dfloat(q+1))
                     endif
                  endif
 17            continue
 18         continue
 19      continue
 20   continue
C     *** main loop to accumulate wkgg and xiecme ******
      do 300 s=1,m
C        *** symmetrize Ui ******************************************
         do 25 i=1,q
            do 24 j=i+1,q
               u(j,i,s)=u(i,j,s)
 24         continue
 25      continue
C        *** store U_i%*%t(Z_i)%*%inv(V_i)%*%X_i in wkqp ***********
         do 30 i=1,q
            do 29 j=1,p
               sum=dble(0.)
               do 27 k=1,q
                  sum=sum+u(i,k,s)*ztvinvx(k,j,s)
 27            continue
               wkqp(i,j)=sum
 29         continue
 30      continue
C        *** multiply wkqp by lower cholesky of xtwxinv *****
         do 60 i=1,q
            do 55 j=1,p
               sum=dble(0.)
               do 52 k=j,p
                  sum=sum+wkqp(i,k)*wkpp(j,k)
 52            continue
               wkqp(i,j)=sum
 55         continue
 60      continue
C        *** multiply wkqp by its transpose, using symmetry *****
C        ***** to get A_i ***************************************
         do 80 i=1,q
            do 70 j=i,q
               sum=dble(0.)
               do 68 k=1,p
                  sum=sum+wkqp(i,k)*wkqp(j,k)
 68            continue
               a(i,j,s)=sum
               if(i.ne.j) a(j,i,s)=sum
 70         continue
 80      continue
C        *** store (xi-U_i) in wkqq1 and ******************
C        ***** accumulate xiecme **********************
         do 85 i=1,q
            do 84 j=i,q
               sum=oxi(i,j)-u(i,j,s)
               wkqq1(i,j)=sum
               if(i.ne.j) wkqq1(j,i)=wkqq1(i,j)
               xiecme(i,j)=xiecme(i,j)+a(i,j,s)+u(i,j,s)
     /              +b(i,s)*b(j,s)/osigma2
 84         continue
 85      continue
C        *** accumulate zeroeth row of wkgg ************
         gi=0
         do 100 i=1,q
            do 90 j=i,q
               gi=gi+1
               if(i.eq.j) then
                  wkgg(0,gi)=wkgg(0,gi)+wkqq1(i,i)
               else
                  wkgg(0,gi)=wkgg(0,gi)+dble(2.)*wkqq1(i,j)
               endif
 90         continue
 100     continue
C        *** accumulate rest of wkgg ********************
         gi=0
         do 200 i=1,q
            do 190 j=i,q
               gi=gi+1
               gj=gi-1
               do 180 ii=i,q
                  if(ii.eq.i) then
                     jjmin=j
                  else
                     jjmin=ii
                  endif
                  do 170 jj=jjmin,q
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahah(q,wkqq1,i,ii)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,i,ii,jj)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,ii,i,j)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trajaj(q,wkqq1,i,j,ii,jj)
                        endif
                     endif
 170              continue
 180           continue
 190        continue
 200     continue
 300  continue
C     *** finish computing xiecme **************
      do 302 i=1,q
         do 301 j=i,q
            xiecme(i,j)=xiecme(i,j)/(dfloat(m-q-1)+abc(3))
            if(i.ne.j) xiecme(j,i)=xiecme(i,j)
 301     continue
 302  continue
C     *** put the inverse of oxi into wkqq1 ****
      do 308 i=1,q
         do 307 j=i,q
            wkqq2(i,j)=oxi(i,j)
 307     continue
 308  continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) then
         msg=10
         goto 999
      endif
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,wkqq1)
C     ************ finish off wkgg ***********************
      wkgg(0,0)=(abc(2)+dfloat(ntot-p-2)+abc(3)*dfloat(q))/dble(2.)
      gi=0
      do 325 i=1,q
         do 324 j=i,q
            gi=gi+1
            wkgg(0,gi)=wkgg(0,gi)/dble(2.)
            if(i.eq.j) wkgg(0,gi)=wkgg(0,gi)*wkqq1(i,i)
            wkgg(gi,0)=wkgg(0,gi)
            gj=gi-1
            do 320 ii=i,q
               if(ii.eq.i) then
                  jjmin=j
               else
                  jjmin=ii
               endif
               do 318 jj=jjmin,q
                  gj=gj+1
                  wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                  if(i.eq.j) then
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(i,i)*
     /                       wkqq1(ii,ii)
                     else
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(i,i)
                     endif
                  else
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(ii,ii)
                     endif
                  endif
                  wkgg(gj,gi)=wkgg(gi,gj)
 318           continue
 320        continue
 324     continue
 325  continue
C     **** compute wkg *******************************
      sum=dble(0.)
      gi=0
      do 344 i=1,q
         do 343 j=i,q
            gi=gi+1
            if(i.eq.j) then
               sum=sum+wkgg(0,gi)*dlog(wkqq1(i,i))
            else
               sum=sum+wkgg(0,gi)*wkqq1(i,j)
            endif
 343     continue
 344  continue
      wkg(0)=(dfloat(ntot-p-2)+abc(2)+abc(3)*dfloat(q))
     /     *dble(.5)*(dble(1.)-sig2ecme/osigma2)
     /     - wkgg(0,0)*dlog(osigma2) + sum
      gi=0
      do 352 i=1,q
         do 351 j=i,q
            gi=gi+1
            if(i.eq.j) then
               wkg(gi)=(abc(3)+dfloat(m-q-1))*dble(.5)*(oxi(i,i)
     /              -xiecme(i,i))*wkqq1(i,i)-wkgg(0,gi)*dlog(osigma2)
            else
               wkg(gi)=(abc(3)+dfloat(m-q-1))*(oxi(i,j)-xiecme(i,j))
     /              - wkgg(0,gi)*dlog(osigma2)
            endif
            gj=0
            sum=dble(0.)
            do 350 ii=1,q
               do 349 jj=ii,q
                  gj=gj+1
                  if(ii.eq.jj) then
                     sum=sum+wkgg(gi,gj)*dlog(wkqq1(ii,ii))
                  else
                     sum=sum+wkgg(gi,gj)*wkqq1(ii,jj)
                  endif
 349           continue
 350        continue
            wkg(gi)=wkg(gi)+sum
 351     continue
 352  continue
C     *** solve the linear system, storing the result in wkg ******
C     **** wkgg is converted to its inverse square root ***********
      call chfce(g+1,g+1,wkgg,err)
      if(err.eq.1) then
         msg=11
         goto 999
      endif
      call bkslv(g+1,g+1,wkgg)
      do 360 i=g,0,-1
         sum=dble(0.)
         do 359 j=0,i
            sum=sum+wkgg(j,i)*wkg(j)
 359     continue
         wkg(i)=sum
 360  continue
      do 370 i=0,g
         sum=dble(0.)
         do 369 j=i,g
            sum=sum+wkgg(i,j)*wkg(j)
 369     continue
         wkg(i)=sum
 370  continue
C     *** invert wkg, putting the result into xi **************
C     ****** and get new sigma2 *******************************
C     *** step-halving is used here if wkg is not pos.def. ****
      deflate=dble(1.)
      oldtau=-dlog(osigma2)
      do 420 i=1,q
         wkqq1(i,i)=dlog(wkqq1(i,i))
 420  continue
 425  continue
      tau=oldtau+deflate*(wkg(0)-oldtau)
      gi=0
      do 430 i=1,q
         do 429 j=i,q
            gi=gi+1
            wkqq2(i,j)=wkqq1(i,j)+deflate*(wkg(gi)-wkqq1(i,j))
 429     continue
         wkqq2(i,i)=dexp(wkqq2(i,i))
 430  continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) then
         deflate=deflate/dfloat(2)
         goto 425
      endif
      sigma2=dexp(-tau)
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,xi)
      do 440 i=1,q
         do 439 j=i+1,q
            xi(j,i)=xi(i,j)
 439     continue
 440  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine fscovm(m,q,b,u,a,xi,oxi,xiecme,wkqq1,wkqq2,
     /     p,xtwxinv,wkpp,wkqp,g,wkg,wkgg,sigma2,msg,
     /     osigma2,ntot,sig2ecme,ztvinvx,dinv,abc)
C Modification of fscovr for posterior mode 
C       oxi = old value of xi
C        xi = new value
C    xiecme = ecme version of xi
C Error messages returned through msg:
C     0 = no error
C    10 = non-positive definite xi
C    11 = non-positive definite wkgg
      integer m,q,p,g,msg,ntot
      double precision b(q,m),u(q,q,m),a(q,q,m),
     /     xi(q,q),oxi(q,q),xiecme(q,q),wkqq1(q,q),wkqq2(q,q),
     /     xtwxinv(p,p),wkpp(p,p),wkqp(q,p),wkg(0:g),
     /     wkgg(0:g,0:g),sigma2,sum,deflate,
     /     osigma2,sig2ecme,ztvinvx(q,p,m),oldtau,tau,
     /     dinv(q,q),abc(3)
      integer s,i,j,k,ii,jj,jjmin,gi,gj,err
      double precision trahah,trahaj,trajaj
      msg=0
C     *** initialize xiecme *****
      do 2 i=1,q
         do 1 j=i,q
            xiecme(i,j)=dinv(i,j)/osigma2
 1       continue
 2    continue
C     *** store cholesky of xtwxinv in wkpp ******
      do 5 i=1,p
         do 4 j=i,p
            wkpp(i,j)=xtwxinv(i,j)
 4       continue
 5    continue
      call chfce(p,p,wkpp,err)
C     *** initialize the workspaces wkg and wkgg *******
      do 7 i=0,g
         wkg(i)=dble(0.)
 7    continue
      gi=0
      do 20 i=1,q
         do 19 j=i,q
            gi=gi+1
            if(i.eq.j) then
               wkgg(0,gi)=dinv(i,i)/osigma2
            else
               wkgg(0,gi)=dble(2.)*dinv(i,j)/osigma2
            endif
            gj=gi-1
            do 18 ii=i,q
               if(ii.eq.i) then
                  jjmin=j
               else
                  jjmin=ii
               endif
               do 17 jj=jjmin,q
                  gj=gj+1
                  if(i.eq.j) then
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=trahah(q,oxi,i,ii)*
     /                       (abc(3)-float(q+1))
                     else
                        wkgg(gi,gj)=trahaj(q,oxi,i,ii,jj)*
     /                       (abc(3)-float(q+1))
                     endif
                  else
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=trahaj(q,oxi,ii,i,j)*
     /                       (abc(3)-float(q+1))
                     else
                        wkgg(gi,gj)=trajaj(q,oxi,i,j,ii,jj)*
     /                       (abc(3)-float(q+1))
                     endif
                  endif
 17            continue
 18         continue
 19      continue
 20   continue
C     *** main loop to accumulate wkgg and xiecme ******
      do 300 s=1,m
C        *** symmetrize Ui ******************************************
         do 25 i=1,q
            do 24 j=i+1,q
               u(j,i,s)=u(i,j,s)
 24         continue
 25      continue
C        *** store U_i%*%t(Z_i)%*%inv(V_i)%*%X_i in wkqp ***********
         do 30 i=1,q
            do 29 j=1,p
               sum=dble(0.)
               do 27 k=1,q
                  sum=sum+u(i,k,s)*ztvinvx(k,j,s)
 27            continue
               wkqp(i,j)=sum
 29         continue
 30      continue
C        *** multiply wkqp by lower cholesky of xtwxinv *****
         do 60 i=1,q
            do 55 j=1,p
               sum=dble(0.)
               do 52 k=j,p
                  sum=sum+wkqp(i,k)*wkpp(j,k)
 52            continue
               wkqp(i,j)=sum
 55         continue
 60      continue
C        *** multiply wkqp by its transpose, using symmetry *****
C        ***** to get A_i ***************************************
         do 80 i=1,q
            do 70 j=i,q
               sum=dble(0.)
               do 68 k=1,p
                  sum=sum+wkqp(i,k)*wkqp(j,k)
 68            continue
               a(i,j,s)=sum
               if(i.ne.j) a(j,i,s)=sum
 70         continue
 80      continue
C        *** store (xi-U_i) in wkqq1 and ******************
C        ***** accumulate xiecme **********************
         do 85 i=1,q
            do 84 j=i,q
               sum=oxi(i,j)-u(i,j,s)
               wkqq1(i,j)=sum
               if(i.ne.j) wkqq1(j,i)=wkqq1(i,j)
               xiecme(i,j)=xiecme(i,j)+a(i,j,s)+u(i,j,s)
     /              +b(i,s)*b(j,s)/osigma2
 84         continue
 85      continue
C        *** accumulate zeroeth row of wkgg ************
         gi=0
         do 100 i=1,q
            do 90 j=i,q
               gi=gi+1
               if(i.eq.j) then
                  wkgg(0,gi)=wkgg(0,gi)+wkqq1(i,i)
               else
                  wkgg(0,gi)=wkgg(0,gi)+dble(2.)*wkqq1(i,j)
               endif
 90         continue
 100     continue
C        *** accumulate rest of wkgg ********************
         gi=0
         do 200 i=1,q
            do 190 j=i,q
               gi=gi+1
               gj=gi-1
               do 180 ii=i,q
                  if(ii.eq.i) then
                     jjmin=j
                  else
                     jjmin=ii
                  endif
                  do 170 jj=jjmin,q
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahah(q,wkqq1,i,ii)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,i,ii,jj)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,ii,i,j)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trajaj(q,wkqq1,i,j,ii,jj)
                        endif
                     endif
 170              continue
 180           continue
 190        continue
 200     continue
 300  continue
C     *** finish computing xiecme **************
      do 302 i=1,q
         do 301 j=i,q
            xiecme(i,j)=xiecme(i,j)/(dfloat(m)+abc(3)-dfloat(q+1))
            if(i.ne.j) xiecme(j,i)=xiecme(i,j)
 301     continue
 302  continue
C     ************ finish off wkgg ***********************
      wkgg(0,0)=(dfloat(ntot-p-2)+abc(2)+abc(3)*dfloat(q))
     /     *(osigma2**2)/dble(2.)
      gi=0
      do 325 i=1,q
         do 324 j=i,q
            gi=gi+1
            wkgg(0,gi)=wkgg(0,gi)*osigma2/dble(2.)
            wkgg(gi,0)=wkgg(0,gi)
            do 322 gj=gi,g
                wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                if(gi.ne.gj) wkgg(gj,gi)=wkgg(gi,gj)
 322        continue
 324      continue
 325   continue
C     *** put the inverse of oxi into wkqq1 ****
      do 335 i=1,q
         do 334 j=i,q
            wkqq2(i,j)=oxi(i,j)
 334     continue
 335  continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) then
         msg=10
         goto 999
      endif
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,wkqq1)
C     **** compute wkg *******************************
      sum=dble(0.)
      gi=0
      do 344 i=1,q
         do 343 j=i,q
            gi=gi+1
            sum=sum+wkgg(0,gi)*wkqq1(i,j)
 343     continue
 344  continue
      wkg(0)=(dfloat(ntot-p-2)+abc(2)+abc(3)*dfloat(q))
     /     *(osigma2-sig2ecme/dble(2.))+sum
      gi=0
      do 352 i=1,q
         do 351 j=i,q
            gi=gi+1
            if(i.eq.j) then
               wkg(gi)=(abc(3)+dfloat(m-q-1))*(oxi(i,i)-xiecme(i,i))
     /              /dfloat(2) + wkgg(0,gi)/osigma2
            else
               wkg(gi)=(abc(3)+dfloat(m-q-1))*(oxi(i,j)-xiecme(i,j))
     /              + wkgg(0,gi)/osigma2
            endif
            gj=0
            sum=dble(0.)
            do 350 ii=1,q
               do 349 jj=ii,q
                  gj=gj+1
                  sum=sum+wkgg(gi,gj)*wkqq1(ii,jj)
 349           continue
 350        continue
            wkg(gi)=wkg(gi)+sum
 351     continue
 352  continue
C     *** solve the linear system, storing the result in wkg ******
C     **** wkgg is converted to its inverse square root ***********
      call chfce(g+1,g+1,wkgg,err)
      if(err.eq.1) then
         msg=11
         goto 999
      endif
      call bkslv(g+1,g+1,wkgg)
      do 360 i=g,0,-1
         sum=dble(0.)
         do 359 j=0,i
            sum=sum+wkgg(j,i)*wkg(j)
 359     continue
         wkg(i)=sum
 360  continue
      do 370 i=0,g
         sum=dble(0.)
         do 369 j=i,g
            sum=sum+wkgg(i,j)*wkg(j)
 369     continue
         wkg(i)=sum
 370  continue
C     *** invert wkg, putting the result into xi **************
C     ****** and get new sigma2 *******************************
C     *** step-halving is used here if wkg is not pos.def. ****
C     ****** or if sigma2 is negative *************************
      deflate=dble(1.)
      oldtau=dble(1.)/osigma2
 425  continue
      tau=oldtau+deflate*(wkg(0)-oldtau)
      gi=0
      do 430 i=1,q
         do 429 j=i,q
            gi=gi+1
            wkqq2(i,j)=wkqq1(i,j)+deflate*(wkg(gi)-wkqq1(i,j))
 429     continue
 430  continue
      call chfce(q,q,wkqq2,err)
      if((err.eq.1).or.(tau.le.dble(0))) then
         deflate=deflate/dfloat(2)
         goto 425
      endif
      sigma2=dble(1.)/tau
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,xi)
      do 440 i=1,q
         do 439 j=i+1,q
            xi(j,i)=xi(i,j)
 439     continue
 440  continue
 999  continue
      return
      end
C***********************************************************************
C The following subroutines are used by fastml().
C***********************************************************************
C Algorithm FAST-ML: Finds ML estimates when V_i are known.
      subroutine fastml(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,vinv,
     /     pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,err,msg,u,iter,
     /     sflag,sigma2,p,xcol,beta,y,delta,xtw,xtwx,xtwy,
     /     xtwxinv,wkqq1,wkqq2,xi,wkqnm,b,cvgd,oxi,maxits,llvec,
     /     eps,xiecme,g,reject,wkg,wkgg)
C set sflag=1 if starting values supplied.
C Returned error codes in msg:
C     0 = no errors
C     1 = supplied V_i matrix is not positive definite
C     2 = GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank
C     3 = inadequate information to obtain starting value of xi
C     4 = value of xi or inv(Ui) became non-pos.def. during iterations
C     5 = t(X)%*%W%*%X became non-pos.def. during iterations
C    10 = non-positive definite xi at input to scoring step
C    11 = non-positive definite wkgg in scoring step
C
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,
     /     pcol,q,zcol(q),iflag,err,msg,iter,sflag,p,xcol(p),cvgd,
     /     maxits,c2,c3,i,j,g,reject(maxits),notconcave
      double precision vmax(nmax,nmax),w(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),
     /     ztvinv(q,nmax,m),ztvinvz(q,q,m),ldv,u(q,q,m),sigma2,
     /     beta(p),y(ntot),delta(ntot),xtw(p,nmax),xtwx(p,p),
     /     xtwy(p),xtwxinv(p,p),wkqq1(q,q),wkqq2(q,q),xi(q,q),
     /     wkqnm(q,nmax,m),b(q,m),oxi(q,q),osigma2,ldu,ldxi,
     /     ll,llvec(maxits),eps,xiecme(q,q),sig2ecme,osig2ecme,
     /     wkg(0:g),wkgg(0:g,0:g)
      msg=0
      iter=0
      notconcave=0
      call prefstml(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,ldv,err)
      if(err.eq.1) then
         msg=1
         goto 999
      endif
      if(sflag.ne.1) then
         call stval1(ntot,m,ist,ifin,occ,nmax,vinv,pcol,pred,
     /        q,ztvinv,ztvinvz,iflag,err,msg,sigma2,p,xcol,
     /        beta,y,delta,xtw,xtwx,xtwy,xtwxinv,wkqq1,wkqq2,xi,
     /        wkqnm,b)
      endif
      cvgd=0
C********* start of main iteration *****************
 1    continue
         iter=iter+1
         reject(iter)=0
 3       continue
         osigma2=sigma2
C        *** calculate U_i, U_i%*%t(Z_i)%*%inv(V_i), and W_i ****
         call mku(q,xi,m,ztvinvz,u,wkqq1,wkqq2,ldxi,ldu,err)
         if(err.eq.1) then
            msg=4
            goto 999
         endif
         call mkwkqnm(q,m,u,nmax,ztvinv,wkqnm,ntot,occ,ist,ifin)
         call mkw(q,nmax,m,ist,ifin,wkqnm,ztvinv,vinv,w,ntot,occ,
     /        iflag)
C        *** save current value of sigma2ecme ******************
         osig2ecme=sig2ecme
C        *** calculate beta, delta_i, and sigma2 ********************
         call gls(ntot,m,ist,ifin,occ,nmax,pcol,pred,0,sig2ecme,
     /        p,xcol,beta,y,delta,w,xtw,xtwx,xtwy,xtwxinv,err)
         if(err.eq.1) then
            msg=5
            goto 999
         endif
C        *** find current loglikelihood *****************************
         ll = - dble(.5)*dfloat(ntot)*dlog(osigma2)
     /        + dfloat(m)*ldxi + ldu 
     /        - dble(.5)*dfloat(ntot)*sig2ecme/osigma2
         llvec(iter)=ll
         if(iter.gt.1) then
C        *** if loglikelihood has gone down, replace the scoring ****
C        ******  estimates of xi and sigma2 by the ecme estimates ***
            if(reject(iter-1).eq.0) then
               if(llvec(iter).lt.llvec(iter-1)) then
                  sigma2=osig2ecme
                  do 17 i=1,q
                     do 16 j=1,q
                        xi(i,j)=xiecme(i,j)
 16                  continue
 17               continue
                  reject(iter-1)=1
                  goto 3
               endif
            endif
         endif
C        *** check for convergence **********************************
         if(iter.gt.1) then
            c2=0
            do 20 i=1,q
               do 19 j=i,q
                  if(dabs(xi(i,j)-oxi(i,j)).gt.(eps*dabs(oxi(i,j))))
     /                 c2=1
 19            continue
 20         continue
            c3=0
            if(dabs(sigma2-osigma2).gt.(eps*dabs(osigma2))) c3=1
            if((c2.eq.0).and.(c3.eq.0)) then
               cvgd=1
               goto 50
            endif
         endif
C        *** calculate b_i ******************************************
         call mkb(q,nmax,m,wkqnm,ntot,delta,b,occ,ist,ifin)
C        *** save current parameters **********************
         osigma2=sigma2
         do 27 i=1,q
            do 26 j=1,q
               oxi(i,j)=xi(i,j)
 26         continue
 27      continue
C        *** perform scoring step, putting result in xi and sigma2,**
C         ****** and put ECME estimate of xi in xiecme ***************
         call fscov2(m,q,b,u,xi,oxi,wkqq1,wkqq2,
     /        g,wkg,wkgg,sigma2,msg,xiecme,osigma2,ntot,sig2ecme)
         if(msg.eq.10) goto 999
         if(msg.eq.11) then
C        *** wkgg not pos.def., substitute ECME values instead ****
            notconcave=1
            do 29 i=1,q
               do 28 j=1,q
                  xi(i,j)=xiecme(i,j)
 28            continue
 29         continue
            sigma2=sig2ecme
            reject(iter)=1
         endif
         if(iter.lt.maxits) goto 1
C********* end of main iteration *****************
 50   continue
      call bdiag(q,m,u)
      do 70 i=1,p
         do 60 j=i+1,p
            xtwxinv(j,i)=xtwxinv(i,j)
 60      continue
 70   continue
      if(notconcave.eq.1) msg=11
 999  continue
      end
C***********************************************************************
      subroutine fscov2(m,q,b,u,xi,oxi,wkqq1,wkqq2,
     /     g,wkg,wkgg,sigma2,msg,xiecme,osigma2,ntot,sig2ecme)
C Fisher scoring algorithm to perform ML estimation for xi and sigma2
C Like fscov, except that scoring is performed on etastar rather 
C than eta.
C     oxi = old value of xi
C      xi = new value
C  xiecme = ecme estimate of xi
C Error messages returned through msg:
C     0 = no error
C    10 = non-positive definite xi
C    11 = non-positive definite wkgg
      integer m,q,g,msg,ntot
      double precision b(q,m),u(q,q,m),xi(q,q),
     /     oxi(q,q),wkqq1(q,q),wkqq2(q,q),wkg(0:g),
     /     wkgg(0:g,0:g),sigma2,sum,deflate,xiecme(q,q),
     /     osigma2,sig2ecme,oldtau,tau
      integer s,i,j,ii,jj,jjmin,gi,gj,err
      double precision trahah,trahaj,trajaj
      msg=0
C     *** initialize xiecme ****************************
      do 4 i=1,q
         do 3 j=i,q
            xiecme(i,j)=dble(0.)
 3       continue
 4    continue
C     *** initialize the workspaces wkg and wkgg *******
      do 8 i=0,g
         wkg(i)=dble(0.)
         do 7 j=0,g
            wkgg(i,j)=dble(0.)
 7       continue
 8    continue
C     *** main loop to accumulate wkgg and xiecme *******************
      do 300 s=1,m
C        *** symmetrize Ui ******************************************
         do 25 i=1,q
            do 24 j=i+1,q
               u(j,i,s)=u(i,j,s)
 24         continue
 25      continue
         do 35 i=1,q
            do 34 j=1,q
C              *** this line accumulates xiecme *********************
               xiecme(i,j)=xiecme(i,j)+u(i,j,s)+b(i,s)*b(j,s)/osigma2
C              ******************************************************
 34         continue
 35      continue
C        *** store (xi-U_i) in wkqq1 **********************
         do 45 i=1,q
            do 44 j=i,q
               wkqq1(i,j)=oxi(i,j)-u(i,j,s)
               if(i.ne.j) wkqq1(j,i)=wkqq1(i,j)
 44         continue
 45      continue
C        *** now we're ready to accumulate wkgg ***********
         gi=0
         do 200 i=1,q
            do 190 j=i,q
               gi=gi+1
               if(i.eq.j) then
                  wkgg(0,gi)=wkgg(0,gi)+wkqq1(i,i)
               else
                  wkgg(0,gi)=wkgg(0,gi)+dble(2.)*wkqq1(i,j)
               endif
               gj=gi-1
               do 180 ii=i,q
                  if(ii.eq.i) then
                     jjmin=j
                  else
                     jjmin=ii
                  endif
                  do 170 jj=jjmin,q
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahah(q,wkqq1,i,ii)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,i,ii,jj)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,ii,i,j)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trajaj(q,wkqq1,i,j,ii,jj)
                        endif
                     endif
 170              continue
 180           continue
 190        continue
 200     continue
 300  continue
C     *** finish off xiecme ***************************
      do 302 i=1,q
         do 301 j=i,q
            xiecme(i,j)=xiecme(i,j)/dfloat(m)
            if(i.ne.j) xiecme(j,i)=xiecme(i,j)
 301     continue
 302  continue
C     *** put the inverse of oxi into wkqq1 ****
      do 304 i=1,q
         do 303 j=i,q
            wkqq2(i,j)=oxi(i,j)
 303     continue
 304  continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) then
         msg=10
         goto 999
      endif
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,wkqq1)
C     *** finish off wkgg ****************************
      wkgg(0,0)=dfloat(ntot)/dble(2.)
      gi=0
      do 310 i=1,q
         do 309 j=i,q
            gi=gi+1
            wkgg(0,gi)=wkgg(0,gi)/dble(2.)
            if(i.eq.j) wkgg(0,gi)=wkgg(0,gi)*wkqq1(i,i)
            wkgg(gi,0)=wkgg(0,gi)
            gj=gi-1
            do 308 ii=i,q
               if(ii.eq.i) then
                  jjmin=j
               else
                  jjmin=ii
               endif
               do 307 jj=jjmin,q
                  gj=gj+1
                  wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                  if(i.eq.j) then
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(i,i)*
     /                       wkqq1(ii,ii)
                     else
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(i,i)
                     endif
                  else
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(ii,ii)
                     endif
                  endif
                  wkgg(gj,gi)=wkgg(gi,gj)
 307           continue
 308        continue
 309     continue
 310  continue
C     **** compute wkg *******************************
      sum=dble(0.)
      gi=0
      do 344 i=1,q
         do 343 j=i,q
            gi=gi+1
            if(i.eq.j) then
               sum=sum+wkgg(0,gi)*dlog(wkqq1(i,i))
            else
               sum=sum+wkgg(0,gi)*wkqq1(i,j)
            endif
 343     continue
 344  continue
      wkg(0)=dfloat(ntot)*dble(.5)*(dble(1.)-sig2ecme/osigma2)
     /     - wkgg(0,0)*dlog(osigma2) + sum
      gi=0
      do 352 i=1,q
         do 351 j=i,q
            gi=gi+1
            if(i.eq.j) then
               wkg(gi)=dfloat(m)*dble(.5)*(oxi(i,i)-xiecme(i,i))
     /              *wkqq1(i,i) - wkgg(0,gi)*dlog(osigma2)
            else
               wkg(gi)=dfloat(m)*(oxi(i,j)-xiecme(i,j))
     /              - wkgg(0,gi)*dlog(osigma2)
            endif
            gj=0
            sum=dble(0.)
            do 350 ii=1,q
               do 349 jj=ii,q
                  gj=gj+1
                  if(ii.eq.jj) then
                     sum=sum+wkgg(gi,gj)*dlog(wkqq1(ii,ii))
                  else
                     sum=sum+wkgg(gi,gj)*wkqq1(ii,jj)
                  endif
 349           continue
 350        continue
            wkg(gi)=wkg(gi)+sum
 351     continue
 352  continue
C     *** solve the linear system, storing the result in wkg ******
C     **** wkgg is converted to its inverse square root ***********
      call chfce(g+1,g+1,wkgg,err)
      if(err.eq.1) then
         msg=11
         goto 999
      endif
      call bkslv(g+1,g+1,wkgg)
      do 375 i=g,0,-1
         sum=dble(0.)
         do 374 j=0,i
            sum=sum+wkgg(j,i)*wkg(j)
 374     continue
         wkg(i)=sum
 375  continue
      do 380 i=0,g
         sum=dble(0.)
         do 379 j=i,g
            sum=sum+wkgg(i,j)*wkg(j)
 379     continue
         wkg(i)=sum
 380  continue
C     *** invert wkg, putting the result into xi **************
C     ****** and get new sigma2 *******************************
C     *** step-halving is used here if wkg is not pos.def. ****
C     ****** or if sigma2 is negative *************************
      deflate=dble(1.)
      oldtau=-dlog(osigma2)
      do 420 i=1,q
         wkqq1(i,i)=dlog(wkqq1(i,i))
 420  continue
 425  continue
      tau=oldtau+deflate*(wkg(0)-oldtau)
      gi=0
      do 430 i=1,q
         do 429 j=i,q
            gi=gi+1
            wkqq2(i,j)=wkqq1(i,j)+deflate*(wkg(gi)-wkqq1(i,j))
 429     continue
         wkqq2(i,i)=dexp(wkqq2(i,i))
 430  continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) then
         deflate=deflate/dfloat(2)
         goto 425
      endif
      sigma2=dexp(-tau)
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,xi)
      do 440 i=1,q
         do 439 j=i+1,q
            xi(j,i)=xi(i,j)
 439     continue
 440  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine fscov(m,q,b,u,xi,oxi,wkqq1,wkqq2,
     /     g,wkg,wkgg,sigma2,msg,xiecme,osigma2,ntot,sig2ecme)
C Fisher scoring algorithm to perform ML estimation for xi and sigma2
C     oxi = old value of xi
C      xi = new value
C  xiecme = ecme estimate of xi
C Error messages returned through msg:
C     0 = no error
C    10 = non-positive definite xi
C    11 = non-positive definite wkgg
      integer m,q,g,msg,ntot
      double precision b(q,m),u(q,q,m),xi(q,q),
     /     oxi(q,q),wkqq1(q,q),wkqq2(q,q),wkg(0:g),
     /     wkgg(0:g,0:g),sigma2,sum,deflate,xiecme(q,q),
     /     osigma2,sig2ecme,oldtau,tau
      integer s,i,j,ii,jj,jjmin,gi,gj,err
      double precision trahah,trahaj,trajaj
      msg=0
C     *** initialize xiecme ****************************
      do 4 i=1,q
         do 3 j=i,q
            xiecme(i,j)=dble(0.)
 3       continue
 4    continue
C     *** initialize the workspaces wkg and wkgg *******
      do 8 i=0,g
         wkg(i)=dble(0.)
         do 7 j=0,g
            wkgg(i,j)=dble(0.)
 7       continue
 8    continue
C     *** main loop to accumulate wkgg and xiecme *******************
      do 300 s=1,m
C        *** symmetrize Ui ******************************************
         do 25 i=1,q
            do 24 j=i+1,q
               u(j,i,s)=u(i,j,s)
 24         continue
 25      continue
         do 35 i=1,q
            do 34 j=1,q
C              *** this line accumulates xiecme *********************
               xiecme(i,j)=xiecme(i,j)+u(i,j,s)+b(i,s)*b(j,s)/osigma2
C              ******************************************************
 34         continue
 35      continue
C        *** store (xi-U_i) in wkqq1
         do 45 i=1,q
            do 44 j=i,q
               wkqq1(i,j)=oxi(i,j)-u(i,j,s)
               if(i.ne.j) wkqq1(j,i)=wkqq1(i,j)
 44         continue
 45      continue
C        *** now we're ready to accumulate wkgg ***********
         gi=0
         do 200 i=1,q
            do 190 j=i,q
               gi=gi+1
               if(i.eq.j) then
                  wkgg(0,gi)=wkgg(0,gi)+wkqq1(i,i)
               else
                  wkgg(0,gi)=wkgg(0,gi)+dble(2.)*wkqq1(i,j)
               endif
               gj=gi-1
               do 180 ii=i,q
                  if(ii.eq.i) then
                     jjmin=j
                  else
                     jjmin=ii
                  endif
                  do 170 jj=jjmin,q
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahah(q,wkqq1,i,ii)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,i,ii,jj)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(q,wkqq1,ii,i,j)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trajaj(q,wkqq1,i,j,ii,jj)
                        endif
                     endif
 170              continue
 180           continue
 190        continue
 200     continue
 300  continue
C     *** finish off xiecme ***************************
      do 305 i=1,q
         do 304 j=i,q
            xiecme(i,j)=xiecme(i,j)/dfloat(m)
            if(i.ne.j) xiecme(j,i)=xiecme(i,j)
 304     continue
 305  continue
C     *** finish off wkgg ****************************
      wkgg(0,0)=dfloat(ntot)*(osigma2**2)/dfloat(2)
      do 310 i=1,g
         wkgg(0,i)=wkgg(0,i)*osigma2/dfloat(2)
         wkgg(i,0)=wkgg(0,i)
         do 309 j=i,g
            wkgg(i,j)=wkgg(i,j)/dfloat(2)
            if(j.ne.i) wkgg(j,i)=wkgg(i,j)
 309     continue
 310  continue
C     *** put the inverse of oxi into wkqq1 ****
      do 312 i=1,q
         do 311 j=i,q
            wkqq2(i,j)=oxi(i,j)
 311     continue
 312  continue
      call chfce(q,q,wkqq2,err)
      if(err.eq.1) then
         msg=10
         goto 999
      endif
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,wkqq1)
C     **** compute wkg *******************************
      sum=dble(0.)
      gi=0
      do 314 i=1,q
         do 313 j=i,q
            gi=gi+1
            sum=sum+wkgg(0,gi)*wkqq1(i,j)
 313     continue
 314  continue
      wkg(0)=dfloat(ntot)*(osigma2-sig2ecme/dfloat(2))+sum
      gi=0
      do 322 i=1,q
         do 321 j=i,q
            gi=gi+1
            if(i.eq.j) then
               wkg(gi)=dfloat(m)*(oxi(i,i)-xiecme(i,i))/dfloat(2)
     /              + wkgg(0,gi)/osigma2
            else
               wkg(gi)=dfloat(m)*(oxi(i,j)-xiecme(i,j))
     /              + wkgg(0,gi)/osigma2
            endif
            gj=0
            sum=dble(0.)
            do 318 ii=1,q
               do 317 jj=ii,q
                  gj=gj+1
                  sum=sum+wkgg(gi,gj)*wkqq1(ii,jj)
 317           continue
 318        continue
            wkg(gi)=wkg(gi)+sum
 321     continue
 322  continue
C     *** solve the linear system, storing the result in wkg ******
C     **** wkgg is converted to its inverse square root ***********
      call chfce(g+1,g+1,wkgg,err)
      if(err.eq.1) then
         msg=11
         goto 999
      endif
      call bkslv(g+1,g+1,wkgg)
      do 325 i=g,0,-1
         sum=dble(0.)
         do 324 j=0,i
            sum=sum+wkgg(j,i)*wkg(j)
 324     continue
         wkg(i)=sum
 325  continue
      do 330 i=0,g
         sum=dble(0.)
         do 329 j=i,g
            sum=sum+wkgg(i,j)*wkg(j)
 329     continue
         wkg(i)=sum
 330  continue
C     *** invert wkg, putting the result into xi **************
C     ****** and get new sigma2 *******************************
C     *** step-halving is used here if wkg is not pos.def. ****
C     ****** or if sigma2 is negative *************************
      deflate=dble(1.)
      oldtau=dble(1.)/osigma2
 425  continue
      tau=oldtau+deflate*(wkg(0)-oldtau)
      gi=0
      do 430 i=1,q
         do 429 j=i,q
            gi=gi+1
            wkqq2(i,j)=wkqq1(i,j)+deflate*(wkg(gi)-wkqq1(i,j))
 429     continue
 430  continue
      call chfce(q,q,wkqq2,err)
      if((err.eq.1).or.(tau.le.dble(0))) then
         deflate=deflate/dfloat(2)
         goto 425
      endif
      sigma2=dble(1.)/tau
      call bkslv(q,q,wkqq2)
      call mm(q,q,wkqq2,xi)
      do 440 i=1,q
         do 439 j=i+1,q
            xi(j,i)=xi(i,j)
 439     continue
 440  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine prefstml(ntot,subj,m,ist,ifin,occ,nmax,vmax,wknnm,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,ldv,err)
C Preliminary manipulations for FASTML. After execution, 
C     ist  = vector of starting positions for subject i, i=1,...,m
C     ifin = vector of finishing positions for subject i, i=1,...,m
C     vinv   = inverse of V_i, i=1,...,m
C     ztvinv  = t(z_i)%*% inverse of V_i, i=1,...,m
C     ztvinvz = t(z_i)%*%inv(V_i)%*%(z_i), i=1,...,m
C     ldv = .5*sum of log det(V_i)
C Needs a workspace wknnm.
C Note: if V_i's are all identity (indicated by iflag=1) then 
C vmax, wknnm, vinv, are ignored and occ is filled in with 1,...,ni.
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,pcol,
     /     q,zcol(q),iflag,err
      double precision vmax(nmax,nmax),wknnm(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),ztvinv(q,nmax,m),
     /     ztvinvz(q,q,m),ldv
      call istfin(ntot,subj,m,ist,ifin)
      if(iflag.ne.1) then
         call mkv(m,nmax,vmax,ntot,occ,ist,ifin,wknnm)
         call chv(nmax,m,wknnm,ntot,occ,ist,ifin,ldv,err)
         if(err.eq.1) goto 999
         call bkv(nmax,m,wknnm,ntot,occ,ist,ifin)
         call mmulv(nmax,m,wknnm,vinv,ntot,occ,ist,ifin)
      else
         call mkocc(ntot,occ,m,ist,ifin)
         ldv=dble(0.)
      endif
      call mmu(ntot,pcol,pred,q,zcol,nmax,m,wknnm,occ,ist,ifin,ztvinv,
     /     iflag)
      call mmtm(q,nmax,m,ztvinv,ntot,occ,ist,ifin,ztvinvz)
      if(iflag.ne.1) then
         call mml(ntot,q,nmax,m,wknnm,occ,ist,ifin,ztvinv)
      endif
C     *** fill in elements of ztvinvz below the diagonal ***
      call bdiag(q,m,ztvinvz)
 999  continue
      return
      end
C***********************************************************************
C The following subroutines are used by ecmeml().
C***********************************************************************
C Algorithm ECME-ML: Finds ML estimates when V_i are known.
      subroutine ecmeml(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,err,msg,
     /     u,iter,sflag,sigma2,p,xcol,beta,y,delta,xtw,xtwx,xtwy,
     /     xtwxinv,wkqq1,wkqq2,xi,wkqnm,b,cvgd,obeta,oxi,maxits,llvec,
     /     eps)
C set sflag=1 if starting values supplied.
C Returned error codes in msg:
C     0 = no errors
C     1 = supplied V_i matrix is not positive definite
C     2 = GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank
C     3 = inadequate information to obtain starting value of xi
C     4 = value of xi or inv(Ui) became non-pos.def. during iterations
C     5 = t(X)%*%W%*%X became non-pos.def. during iterations
      integer ntot,subj(ntot),m,ist(m),ifin(m),occ(ntot),nmax,
     /     pcol,q,zcol(q),iflag,err,msg,iter,sflag,p,xcol(p),cvgd,
     /     maxits,c1,c2,c3,i,j
      double precision vmax(nmax,nmax),w(nmax,nmax,m),
     /     vinv(nmax,nmax,m),pred(ntot,pcol),
     /     ztvinv(q,nmax,m),ztvinvz(q,q,m),ldv,u(q,q,m),sigma2,
     /     beta(p),y(ntot),delta(ntot),xtw(p,nmax),xtwx(p,p),xtwy(p),
     /     xtwxinv(p,p),wkqq1(q,q),wkqq2(q,q),xi(q,q),
     /     wkqnm(q,nmax,m),b(q,m),obeta(p),oxi(q,q),osigma2,ldu,ldxi,
     /     ll,llvec(maxits),eps
      msg=0
      iter=0
      call preecme1(ntot,subj,m,ist,ifin,occ,nmax,vmax,w,
     /     vinv,pcol,pred,q,zcol,ztvinv,ztvinvz,iflag,ldv,err)
      if(err.eq.1) then
         msg=1
         goto 999
      endif
      if(sflag.ne.1) then
         call stval1(ntot,m,ist,ifin,occ,nmax,vinv,pcol,pred,
     /        q,ztvinv,ztvinvz,iflag,err,msg,sigma2,p,xcol,
     /        beta,y,delta,xtw,xtwx,xtwy,xtwxinv,wkqq1,wkqq2,xi,
     /        wkqnm,b)
      endif
      cvgd=0
C********* start of main iteration *****************
 1    continue
         iter=iter+1
         do 2 i=1,p
            obeta(i)=beta(i)
 2       continue
         do 10 i=1,q
            do 5 j=i,q
               oxi(i,j)=xi(i,j)
 5          continue
 10      continue
         osigma2=sigma2
         call mku(q,xi,m,ztvinvz,u,wkqq1,wkqq2,ldxi,ldu,err)
         if(err.eq.1) then
            msg=4
            goto 999
         endif
         call mkwkqnm(q,m,u,nmax,ztvinv,wkqnm,ntot,occ,ist,ifin)
         call mkw(q,nmax,m,ist,ifin,wkqnm,ztvinv,vinv,w,ntot,occ,
     /        iflag)
         call gls(ntot,m,ist,ifin,occ,nmax,pcol,pred,0,sigma2,
     /        p,xcol,beta,y,delta,w,xtw,xtwx,xtwy,xtwxinv,err)
         if(err.eq.1) then
            msg=5
            goto 999
         endif
         call mkb(q,nmax,m,wkqnm,ntot,delta,b,occ,ist,ifin)
         call mkxi(q,m,b,u,xi,osigma2)
         ll = - dble(.5)*dfloat(ntot)*dlog(osigma2)
     /        + dfloat(m)*ldxi + ldu -.5*dfloat(ntot)*sigma2/osigma2
         llvec(iter)=ll
         c1=0
         do 30 i=1,p
            if(dabs(beta(i)-obeta(i)).gt.(eps*dabs(obeta(i)))) c1=1
 30      continue
         c2=0
         do 40 i=1,q
            do 35 j=i,q
               if(dabs(xi(i,j)-oxi(i,j)).gt.(eps*dabs(oxi(i,j))))
     /              c2=1
 35         continue
 40      continue
         c3=0
         if(dabs(sigma2-osigma2).gt.(eps*dabs(osigma2))) c3=1
         if((c1.eq.0).and.(c2.eq.0).and.(c3.eq.0)) cvgd=1
         if((cvgd.eq.0).and.(iter.lt.maxits)) goto 1
C********* end of main iteration *****************
      call bdiag(q,m,u)
      do 70 i=1,p
         do 60 j=i+1,p
            xtwxinv(j,i)=xtwxinv(i,j)
 60      continue
 70   continue
 999  continue
      return
      end
C***********************************************************************
      subroutine mkxi(q,m,b,u,xi,sigma2)
C For ECME-ML: calculates new estimate of xi from b, u, and sigma2
      integer q,m,s,i,j
      double precision b(q,m),u(q,q,m),xi(q,q),sigma2
      do 10 i=1,q
         do 5 j=i,q
            xi(i,j)=dble(0.)
 5       continue
 10   continue
      do 100 s=1,m
         do 90 i=1,q
            do 80 j=i,q
               xi(i,j)=xi(i,j)+sigma2*u(i,j,s)+b(i,s)*b(j,s)
 80         continue
 90      continue
 100  continue
      do 110 i=1,q
         do 105 j=i,q
            xi(i,j)=xi(i,j)/(dfloat(m)*sigma2)
            if(i.ne.j) xi(j,i)=xi(i,j)
 105     continue
 110  continue
      return
      end
C***********************************************************************
