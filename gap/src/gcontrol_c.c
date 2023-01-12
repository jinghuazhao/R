#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
/*#include <malloc/malloc.h> necessary for Mac OS X 10.3 (Panther)*/
#include <string.h>

#define FALSE 0
#define TRUE 1
#ifndef MAX_INPUT
#define MAX_INPUT 80
#endif
#define NUM_VRTS 2
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
/*#define PI 3.14159265*/

int getnum(FILE *, double *);
double poz(double);
double median(double *,int *);
double rnorm(long *),runif(long int *),dnorm(double),rchisq(int,long *);
double chi(int *,int *,int,int);
void sort(double *,int *,double *);
void prog(double *,int,double *,long *,double *,int *);
void datnull(int,double,long *,int *);
void datmult(int,double,long *,int *);
void datdis(int,double,long *,int *);
void datagen(int,int,int,double,double,double,int,double *,double *,long *);
void rgamma(int *,double *,double *,double *,long *);

double chi(int *o, int *e, int R, int N)
{
  int n1,n2,r1,r2;
  double num,den;

  n1=o[1]+e[1];
  n2=o[2]+e[2];
  r1=o[1];
  r2=o[2];
  num=N*pow(N*(r1+2.*r2)-R*(n1+2.*n2),2.);
  den=R*(N-R)*(N*(n1+4.*n2)-pow(n1+2.*n2,2.));
  return(sqrt(num/den));
}

void datagen(int n,int R,int S,double gamma,double F1,double F2,int oo,double *x,double *xt,long *idum)
{
  int k,o[3],e[3],N;
  N=R+S;
  for (k=0;k<(n-oo);k++)
  {
    datnull(R,F1,idum,o);
    datnull(S,F2,idum,e);
    xt[k]=chi(o,e,R,N);
  }
  for (k=(n-oo);k<n;k++)
  {
    datdis(R,gamma,idum,o);
    datnull(S,F2,idum,e);
    xt[k]=chi(o,e,R,N);
  }
  k=n-oo;
  sort(xt,&k,x);
  sort(xt+k,&oo,x+k);
  return;
}

void datdis(int m, double gamma, long *idum, int *rslt)
{
  double p1;

  p1=gamma/(1.+gamma);
  datmult(m,p1,idum,rslt);
  return;
}

void datmult(int m, double p, long *idum, int *rslt)
{
  int i,j,gen;

  for (i=0;i<3;i++) rslt[i]=0;
  for (i=0;i<m;i++)
  {
    gen=0;
    for (j=0;j<2;j++) gen+=(p<runif(idum));
    rslt[gen]++;
  }
  return;
}

void datnull(int m, double theta, long *idum, int *rslt)
{
  int num=NUM_VRTS,i;

  double alphadot,alpha[2*NUM_VRTS],c_1[2*NUM_VRTS]={1.,1.,1.,1.},
  p,
  z[2*NUM_VRTS];
  alphadot=(1.-theta)/theta;
  for (i=0;i<2*num;i++) alpha[i]=1./2.*alphadot;
  rgamma(&num,alpha,c_1,z,idum);
  p=z[0]/(z[0]+z[1]);
  datmult(m,p,idum,rslt);
  return;
}

int getnum(FILE *f, double *x)
{
  static char buffer[MAX_INPUT+1];
  if (fgets(buffer,MAX_INPUT,f)==NULL) return(0);
  return(sscanf(buffer,"%*s %lf %lf %lf %lf %lf %lf",x,x+1,x+2,x+3,x+4,x+5));
}

double median(double *x,int *n)
{
  int i_1;
  static double hold,xmed,*y;
  static int nmid;
  static int i, iflag;
  static int nmidp1;

  --x;
  if (*n < 1) goto L50;
  if (*n == 1) goto L55;
  hold = x[1];
  i_1 = *n;
  for (i = 2; i <= i_1; ++i) if (x[i] != hold) goto L90;
  xmed = x[1];
  goto L101;
L50:
  REprintf("Invalid vector length in median routine");
  error("%d",1);
L55:
  xmed = x[1];
  goto L101;
L90:
  if((y=(double *)malloc((*n)*sizeof(double)))==NULL)
  {
    REprintf("I can't allocate memory: median routine");
    error("%d",1);
  }
  sort(&x[1], n, y);
  iflag = *n - (*n / 2 << 1);
  nmid = *n / 2;
  nmidp1 = nmid + 1;
  if (iflag == 0) xmed = (y[nmid - 1] + y[nmidp1 - 1]) / 2.;
  if (iflag == 1) xmed = y[nmidp1 - 1];
L101:
  return xmed;
}

void sort(double *x,int *n,double *y)
{
  static double amed;
  static int i, j, k, l, m, il[36], iu[36],i_1;
  static double tt;
  static int ip1, nm1, mid, jmi, lmi, jmk;

  --x;
  --y;
  i_1 = *n;
  for (i = 1; i <= i_1; ++i) y[i] = x[i];

  nm1 = *n - 1;
  i_1 = nm1;
  for (i = 1; i <= i_1; ++i)
  {
    ip1 = i + 1;
    if (y[i] <= y[ip1]) goto L200;
    goto L250;
    L200:
    ;
  }
  return;
L250:
  m = 1;
  i = 1;
  j = *n;
L305:
  if (i >= j) goto L370;
L310:
  k = i;
  mid = (i + j) / 2;
  amed = y[mid];
  if (y[i] <= amed) goto L320;
  y[mid] = y[i];
  y[i] = amed;
  amed = y[mid];
L320:
  l = j;
  if (y[j] >= amed) goto L340;
  y[mid] = y[j];
  y[j] = amed;
  amed = y[mid];
  if (y[i] <= amed) goto L340;
  y[mid] = y[i];
  y[i] = amed;
  amed = y[mid];
  goto L340;
L330:
  y[l] = y[k];
  y[k] = tt;
L340:
  --l;
  if (y[l] > amed) goto L340;
  tt = y[l];
L350:
  ++k;
  if (y[k] < amed) goto L350;
  if (k <= l) goto L330;
  lmi = l - i;
  jmk = j - k;
  if (lmi <= jmk) goto L360;
  il[m - 1] = i;
  iu[m - 1] = l;
  i = k;
  ++m;
  goto L380;
L360:
  il[m - 1] = k;
  iu[m - 1] = j;
  j = l;
  ++m;
  goto L380;
L370:
  --m;
  if (m == 0) return;
  i = il[m - 1];
  j = iu[m - 1];
L380:
  jmi = j - i;
  if (jmi >= 11) goto L310;
  if (i == 1) goto L305;
  --i;
L390:
  ++i;
  if (i == j) goto L370;
  amed = y[i + 1];
  if (y[i] <= amed) goto L390;
  k = i;
L395:
  y[k + 1] = y[k];
  --k;
  if (amed < y[k]) goto L395;
  y[k + 1] = amed;
  goto L390;
}

void prog(double *x,int n,double *deltot,long *idum,double *A,int *delta)
{
  double xsumsq=0.,lambda,p,num,r=0,fun1,fun2,cc,dd,rsig,zeta,kappa,tau2,epsilon;
  int i,j,nu=0,df=0,ngib,burn;

  zeta=1000.;kappa=4.;tau2=1.;epsilon=.01;ngib=500;burn=50;
  rsig=median(x,&n)/.675;
  kappa=kappa*rsig;
  tau2=tau2*rsig*rsig;
  for (i=0;i<n;i++) xsumsq+=x[i]*x[i];
  lambda=xsumsq/n;
  dd=1./(1./tau2+1./lambda);
  for (i=0;i<n;i++)
  {
    cc=dd*(x[i]/lambda+kappa/tau2);
    A[i]=dd*rnorm(idum)+cc;
    fun1=(x[i]-A[i])/sqrt(lambda);
    fun2=x[i]/sqrt(lambda);
    p=dnorm(fun1)*epsilon/(dnorm(fun1)*epsilon+dnorm(fun2)*(1.-epsilon));
    delta[i]=(p>runif(idum));
  }
  for (i=0;i<burn+ngib;i++)
  {
    num=(double)nu*zeta;
    df=nu;
    for (j=0;j<n;j++)
    {
      num+=pow(x[j]-delta[j]*A[j],2);
      df+=(1-delta[j]);
    }
    lambda=num/rchisq((long)df,idum);
    dd=1/(1/tau2+1/lambda);
    r=0;
    for (j=0;j<n;j++)
    {
      fun1=(x[j]-A[j])/sqrt(lambda);
      fun2=x[j]/sqrt(lambda);
      p=dnorm(fun1)*epsilon/(dnorm(fun1)*epsilon+dnorm(fun2)*(1.-epsilon));
      delta[j]=(p>runif(idum));
      if (i>=burn) deltot[j]+=p;
      cc=dd*(x[j]/lambda+kappa/tau2);
      A[j]=delta[j]*(sqrt(dd)*rnorm(idum)+cc)+
        (1-delta[j])*(sqrt(tau2)*rnorm(idum)+kappa);
      r+=delta[j];
    }
  }
  for (j=0;j<n;j++) deltot[j]=deltot[j]/ngib;
  REprintf("%14.9f\n",r);
  return;
}

double dnorm(double x)
{
  return(exp(-.5*x*x)/sqrt(2.*PI));
}

double rchisq(int ia, long *idum)
{
  int j;
  double am,e,s,v1,v2,x,y;
  if (ia < 6)
  {
    x=1.0;
    for (j=1;j<=ia;j++) x *= runif(idum);
    x = -log(x);
  }
  else
  {
    do
    {
      do
      {
        do
        {
          v1=runif(idum);
          v2=2.0*runif(idum)-1.0;
        } while (v1*v1+v2*v2 > 1.0);
        y=v2/v1;
        am=ia-1;
        s=sqrt(2.0*am+1.0);
        x=s*y+am;
      } while (x <= 0.0);
      e=(1.0+y*y)*exp(am*log(x/am)-s*y);
    } while (runif(idum) > e);
  }
  return x;
}

void rgamma(int *m,double *na,double *nb,double *x,long *idum)
{
  int na_dim1, na_offset, nb_dim1, nb_offset, x_dim1, x_offset, i_1;
/*int s_stop();*/
  double work[6];
  extern double rgkm3_(double *,double *,int *,long *);
  int i, j, k;
  extern double rgs_(double *,long *);

  na_dim1 = *m;
  na_offset = na_dim1 + 1;
  na -= na_offset;
  nb_dim1 = *m;
  nb_offset = nb_dim1 + 1;
  nb -= nb_offset;
  x_dim1 = *m;
  x_offset = x_dim1 + 1;
  x -= x_offset;
  work[0] = -1.;
  i_1 = *m;
  for (i = 1; i <= i_1; ++i)
  {
    for (j = 1; j <= 2; ++j)
    {
      if (na[i + j * na_dim1] > 0. && nb[i + j * nb_dim1] > 0.)
      {
        L10:
        if (na[i + j * na_dim1] <= 1.)
          x[i + j * x_dim1] = rgs_(&na[i + j * na_dim1], idum);
        else
          x[i+j*x_dim1] = rgkm3_(&na[i+j*na_dim1], work, &k, idum);
        if (x[i + j * x_dim1] <= DBL_EPSILON) goto L10;
      }
      else
      {
        REprintf(" Error in rgamma routine\n");
        error("%d",1);
      }
      x[i + j * x_dim1] *= nb[i + j * nb_dim1];
    }
  }
  return;
}

double rgs_(double *alp,long *idum)
{
  double ret_val;
  double log(double), exp(double);
  static double b, p, x, u1, u2;
  extern double runif(long int *);
L100:
  u1 = runif(idum);
  b = (*alp + 2.718281828) / 2.718281828;
  p = b * u1;
  u2 = runif(idum);
  if (p > 1.) goto L300;
  x = exp(log(p) / *alp);
  if (u2 > exp(-x)) goto L100;
  ret_val = x;
  return ret_val;
L300:
  x = -log((b - p) / *alp);
  if (log(u2) > (*alp - 1.) * log(x)) goto L100;
  ret_val = x;
  return ret_val;
}

double rgkm3_(double *alp,double *work,int *k,long *idum)
{
  double ret_val;
  double sqrt(double), log(double);
  static double u, w, u1, u2;
  extern double runif(long int *);

  --work;
  if (work[1] == *alp) goto L100;
  work[1] = *alp;
  *k = 1;
  if (*alp > 2.5) *k = 2;
  work[2] = *alp - 1.;
  work[3] = (*alp - 1. / (*alp * 6.)) / work[2];
  work[4] = 2. / work[2];
  work[5] = work[4] + 2.;
  switch (*k)
  {
    case 1:  goto L1;
    case 2:  goto L11;
  }
L1:
  u1 = runif(idum);
  u2 = runif(idum);
  goto L20;
L11:
  work[6] = sqrt(*alp);
L15:
  u1 = runif(idum);
  u = runif(idum);
  u2 = u1 + (1. - u * 1.86) / work[6];
  if (u2 <= 0. || u2 >= 1.) goto L15;
L20:
  w = work[3] * u1 / u2;
  if (work[4] * u2 - work[5] + w + 1. / w <= 0.) goto L200;
  if (work[4] * log(u2) - log(w) + w - 1. < 0.) goto L200;
L100:
  switch (*k)
  {
    case 1:  goto L1;
    case 2:  goto L15;
  }
L200:
  ret_val = work[2] * w;
  return ret_val;
}

double rnorm(long *idum)
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  if  (iset == 0)
  {
    do
    {
      v1=2.0*runif(idum)-1.0;
      v2=2.0*runif(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  }
  else
  {
    iset=0;
    return gset;
  }
}

double runif(long *idum)
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  if (*idum < 0 || iff == 0)
  {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++)
    {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++)
    {
      ma[i] -= ma[1+(i+30) % 55];
      if (ma[i] < MZ) ma[i] += MBIG;
    }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

#ifdef gendata

#define r 1
int main()
{
  double *x,*xt,*pp,level[r]={.50},
  gamma, F1, F2, *dwork, out1[r+1], out2[r+1];
  int i,j,iter,*iwork,n,R,S,oo,reps;
  long idum;

  R=100;S=100;gamma=2.5;F1=.001;F2=.0001;oo=10;reps=1;n=400;

  if ((dwork=(double *)malloc(n*sizeof(double)))==NULL)
  {
    error("I can't allocate memory. (perform.c)\n");
  }
  if ((iwork=(int *)malloc(n*sizeof(int)))==NULL)
  {
    error("I can't allocate memory. (perform.c)\n");
  }
  if ((x=(double *)malloc(n*sizeof(double)))==NULL)
  {
    error("I can't allocate memory. (perform.c)\n");
  }
  if ((xt=(double *)malloc(n*sizeof(double)))==NULL)
  {
    error("I can't allocate memory. (perform.c)\n");
  }
  if ((pp=(double *)malloc(n*sizeof(double)))==NULL)
  {
    error("I can't allocate memory. (perform.c)\n");
  }
  for (i=0;i<r+1;i++) out1[i]=out2[i]=0.;
  for (iter=0;iter<reps;iter++)
  {
    datagen(n,R,S,gamma,F1,F2,oo,x,xt,&idum);
    for (j=0;j<n;j++) pp[j]=0.;
    prog(x,n,pp,&idum,dwork,iwork);
    for (i=0;i<r;i++)
    {
      for (j=0;j<(n-oo);j++) out2[i]+=(pp[j]>level[i]);
      for (j=(n-oo);j<n;j++) out1[i]+=(pp[j]>level[i]);
    }
    for (j=0;j<(n-oo);j++) out2[r]+=(x[j]>3.66);
    for (j=(n-oo);j<n;j++) out1[r]+=(x[j]>3.66);
  }
  for (i=0;i<r;i++) REprintf("%5.2f %16.8f %16.8f\n",level[i],out1[i]/(double)reps,out2[i]/(double)reps);
  REprintf("%5.2f %16.8f %16.8f\n",0.,out1[r]/(double)reps,out2[r]/(double)reps);

  return 0;
}
#endif
#ifdef executable
int main()
{
  double  jnk[6];
  double  xsumsq=0.,lambda,p,num,r,fun1,fun2,cc,dd,rsig,
          *x,*A,*out,*deltot,zeta,kappa,tau2,epsilon;
  FILE *pinfile, *rslttxt;
  int i,j,R,N,lcnt=0,n=0,nu=0,*delta,df=0,itmp=0,o[3],e[3],ngib,burn;
  long idum=2348;

  zeta=1000.;kappa=4.;tau2=1.;epsilon=.01;ngib=500;burn=50;
  rslttxt=fopen("rslt.txt","w+");
  if ((pinfile=fopen("kdata.txt","rt"))==NULL)
  {
    error("I can't open file.");
  }
  n=getnum(pinfile,jnk);
  while (!feof(pinfile)) n+=getnum(pinfile,jnk);
  rewind(pinfile);

  if ((x=(double *)malloc(n*sizeof(double)))==NULL)
  {
    error("I can't allocate memory.");
  }
  if ((A=(double *)malloc(n*sizeof(double)))==NULL)
  {
    error("I can't allocate memory.\n");
  }
  if ((deltot=(double *)malloc(n*sizeof(double)))==NULL)
  {
    error("I can't allocate memory.\n");
  }
  if ((delta=(int *)malloc(n*sizeof(int)))==NULL)
  {
    error("I can't allocate memory.\n");
  }
  if ((out=(double *)malloc((burn+ngib)*sizeof(double)))==NULL)
  {
    error("I can't allocate memory.\n");
  }
  for (i=0;i<n;i++)
  {
    x[i]=0.; A[i]=0.; deltot[i]=0.; delta[i]=0;;
  }
  rewind(pinfile);
  n=0;
  while (!feof(pinfile))
  {
    itmp=getnum(pinfile,jnk);
    if (itmp!=6) REprintf("Bad data (n=%d) at line %d\n",itmp,n+1);
    else
    {
      R=0;N=0;
      for(i=0;i<3;i++)
      {
        o[i]=(int)jnk[i];
        e[i]=(int)jnk[i+3];
        R+=o[i];
        N+=o[i]+e[i];
      }
      x[n++]=chi(o,e,R,N);
    }
  }
  fclose(pinfile);
  REprintf("\n%d lines read.\n\n",n);
  rsig=median(x,&n)/.675;
  kappa=kappa*rsig;
  tau2=tau2*rsig*rsig;
  for (i=0;i<n;i++) xsumsq+=x[i]*x[i];
  lambda=xsumsq/n;
  dd=1./(1./tau2+1./lambda);
  for (i=0;i<n;i++)
  {
    cc=dd*(x[i]/lambda+kappa/tau2);
    A[i]=dd*rnorm(&idum)+cc;
    fun1=(x[i]-A[i])/sqrt(lambda);
    fun2=x[i]/sqrt(lambda);
    p=dnorm(fun1)*epsilon/(dnorm(fun1)*epsilon+dnorm(fun2)*(1.-epsilon));
    delta[i]=(p>runif(&idum));
  }
  for (i=0;i<burn+ngib;i++)
  {
    num=(double)nu*zeta;
    df=nu;
    for(j=0;j<n;j++)
    {
      num+=pow(x[j]-delta[j]*A[j],2);
      df+=(1-delta[j]);
    }
    lambda=num/rchisq((long)df,&idum);
    dd=1/(1/tau2+1/lambda);
    r=0;
    for(j=0;j<n;j++)
    {
      fun1=(x[j]-A[j])/sqrt(lambda);
      fun2=x[j]/sqrt(lambda);
      p=dnorm(fun1)*epsilon/(dnorm(fun1)*epsilon+dnorm(fun2)*(1.-epsilon));
      delta[j]=(p>runif(&idum));
      if (i>=burn) deltot[j]+=p;
      if(i>=burn && j==0) /*printf("%14.9f \n",p)*/;
      cc=dd*(x[j]/lambda+kappa/tau2);
      A[j]=delta[j]*(sqrt(dd)*rnorm(&idum)+cc)+
        (1-delta[j])*(sqrt(tau2)*rnorm(&idum)+kappa);
      r+=delta[j];
    }
  }
  for(j=0;j<n;j++)
  {
    if(j==0) REprintf("%d %14.9f \n",ngib,deltot[j]);
    deltot[j]=deltot[j]/ngib;
    fprintf(rslttxt,"%14.9f\n",deltot[j]);
  }
  fclose(rslttxt);
  return 0;
}
#endif

void gcontrol_c(double *kdata, int *nkdata, double *zeta, double *kappa, 
              double *tau2, double *epsilon, int *ngib, int *burn, 
              int *idumR, double *deltot, double *x, double *A)
{
  double jnk[6], xsumsq=0.,lambda,p,num,r=0,fun1,fun2,cc,dd,rsig,*out;
  int    i,j,R,N,n=0,nu=0,*delta,df=0,o[3],e[3];
  long   idum=*idumR;

  n=*nkdata;
  if ((delta=(int *)malloc(n*sizeof(int)))==NULL)
  {
    error("I can't allocate memory.\n");
  }
  if ((out=(double *)malloc((*burn+(*ngib))*sizeof(double)))==NULL)
  {
    error("I can't allocate memory.\n");
  }
  for (i=0;i<n;i++) delta[i]=0;;
  n=0;
  while (n<*nkdata/6)
  {
    for(i=0;i<6;i++) jnk[i]=kdata[n*6+i];
    R=0;N=0;
    for(i=0;i<3;i++)
    {
      o[i]=(int)jnk[i];
      e[i]=(int)jnk[i+3];
      R+=o[i];
      N+=o[i]+e[i];
    }
    x[n++]=chi(o,e,R,N);
  }
  rsig=median(x,&n)/.675;
  *kappa=*kappa*rsig;
  *tau2=*tau2*rsig*rsig;
  for (i=0;i<n;i++) xsumsq+=x[i]*x[i];
  lambda=xsumsq/n;
  dd=1./(1./(*tau2)+1./lambda);
  for (i=0;i<n;i++)
  {
    cc=dd*(x[i]/lambda+(*kappa)/(*tau2));
    A[i]=dd*rnorm(&idum)+cc;
    fun1=(x[i]-A[i])/sqrt(lambda);
    fun2=x[i]/sqrt(lambda);
    p=dnorm(fun1)*(*epsilon)/(dnorm(fun1)*(*epsilon)+dnorm(fun2)*(1.-*epsilon));
    delta[i]=(p>runif(&idum));
  }
  for (i=0;i<(*burn)+(*ngib);i++)
  {
    num=(double)nu*(*zeta);
    df=nu;
    for(j=0;j<n;j++)
    {
      num+=pow(x[j]-delta[j]*A[j],2);
      df+=(1-delta[j]);
    }
    lambda=num/rchisq((long)df,&idum);
    dd=1/(1/(*tau2)+1/lambda);
    r=0;
    for(j=0;j<n;j++)
    {
      fun1=(x[j]-A[j])/sqrt(lambda);
      fun2=x[j]/sqrt(lambda);
      p=dnorm(fun1)*(*epsilon)/(dnorm(fun1)*(*epsilon)+dnorm(fun2)*(1.-*epsilon));
      delta[j]=(p>runif(&idum));
      if (i>=*burn) deltot[j]+=p;
      /*if(i>=*burn && j==0) printf("%14.9f \n",p)*/;
      cc=dd*(x[j]/lambda+(*kappa)/(*tau2));
      A[j]=delta[j]*(sqrt(dd)*rnorm(&idum)+cc)+
        (1-delta[j])*(sqrt(*tau2)*rnorm(&idum)+(*kappa));
      r+=delta[j];
    }
  }
  for(j=0;j<n;j++)
  {
    deltot[j]=deltot[j]/(*ngib);
    REprintf("%14.9f %14.9f\n",x[j],deltot[j]);
  }
  REprintf("%14.9f\n",r);
}
