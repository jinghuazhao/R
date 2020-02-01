/*                                                           */
/* GENECOUNTING - haplotype analysis with missing genotypes  */
/* (C) Copyright JH Zhao 2001-2005 University College London */
/*                                                           */
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ignore.h"

#define version 2.2
#define mxloc 60
#define mxalleles 50

typedef struct
{
  long gid;
  short sex;
  double count,prob;
  union {
  short h[mxloc];
  short d[mxloc][2];
  } hd;
} pat;

typedef struct hnode_type
{
  long int id;
  int n;
  struct hnode_type *left, *right;
  short l[mxloc];
} hnode;

enum {MALE=1,FEMALE};
static double pl,eps,tol;
static int handlemissing,maxit,convll;
static int nloci,nloci2,nalleles,loci[mxloc],obscom,nlocim,locim[mxloc+1];
static double pp[mxloc][mxalleles];
static double d_sample,d_msample;
static int hapall;
static double *c,*h,*hs,*h0,*hm;
static int idm[mxloc];
   
void *xmalloc(long);

hnode *hrt=0;
hnode *hitree(hnode*,long int,short*);
void hrtree(hnode*);

/*ordinary data*/

void hilo(short*,short*);
void getp(pat *),geth(pat *),ch(pat *,char **),counting(pat *, long int);
void genp(pat *,long int);
void digit2(int,short*,int),digitm(short*,short*,int);
int linenum(int*,short*),linenums(int*,short*);
double phasep(short*,short*);
double ll(pat *);

long hapallm;

/*X chromosome data*/
static int nhaploid,ndiploid,h_sample,h_msample,xdata;
static double pm[mxloc][mxalleles];

void gcx(int *verbose,
int *Rhandlemissing,
int *Rconvll,
double *Reps,
double *Rtol,
int *Rmaxit,
double *Rpl,
double *precis,
int *gid,
int *Rnloci,
int *Rloci,
int *Robscom,
int *Rhapall,
int *genotype,
int *count,
int *Rxdata,
int *sex,
int *hapid,
double *prob,
double *Rh0,
double *Rh1,
double *lnl0,
double *lnl1,
int *npusr, 
int *npdat, 
int *iter,  
int *converge,char **assignment  
)
{
time_t now;
short la,ua;
pat *table,tt;
int it,i,j,k;
long int io;
short loci1[mxloc];
double s,lnls,x2mm;

handlemissing=*Rhandlemissing;
eps=*Reps;
tol=*Rtol;
maxit=*Rmaxit;
xdata=*Rxdata;
pl=*Rpl;
convll=*Rconvll;
nloci=*Rnloci;
for(i=0;i<nloci;i++) loci[i]=Rloci[i];
hapall=*Rhapall;
h0=Rh0;
h=Rh1;
obscom=*Robscom;
nloci2=2*nloci;
nalleles=1; 
for(i=0;i<nloci;i++)
{
  if(loci[i]>nalleles) nalleles=loci[i];
}
hs=(double*)xmalloc(hapall*sizeof(double));
c=(double*)xmalloc(hapall*sizeof(double));
if(handlemissing) hm=(double*)xmalloc(hapall*sizeof(double));
nhaploid=ndiploid=h_sample=h_msample=0;
table=(pat*)xmalloc(obscom*sizeof(pat));
for(i=0;i<obscom;i++)
{
  tt.gid=gid[i];
  tt.count=count[i];
  tt.sex=sex[i];
  k=0;
  for(j=0;j<nloci;j++)
  {
    if(!xdata||(xdata&&tt.sex==FEMALE))
    {
      la=genotype[i*nloci2+2*j];
      ua=genotype[i*nloci2+2*j+1];
      if ((la>loci[j]||ua>loci[j])||(la<=0||ua<=0)) ++k;
      if (la>ua) hilo(&la,&ua);
      tt.hd.d[j][0]=la;
      tt.hd.d[j][1]=ua;
    } else
    if(xdata&&tt.sex==MALE)
    {
      la=genotype[i*nloci2+2*j+1];
      if(la>loci[j]||la<=0) ++k;
      tt.hd.h[j]=la;
    }
  }
  table[i]=tt;
  if(xdata)
  {
    if(tt.sex==MALE) nhaploid++;
    else ndiploid++;
  }
}
getp(table);
for(i=0;i<2;i++) npusr[i]=npdat[i]=0;
k=1;
for(i=0;i<nloci;i++)
{
  npusr[0]+=loci[i]-1;
  k*=loci[i];
}
npusr[1]=k-1;
k=1;
for(i=0;i<nloci;i++)
{
  loci1[i]=loci[i];
  for(j=0;j<loci[i];j++) if(pp[i][j]<*precis) --loci1[i];
  npdat[0]+=loci1[i]-1;
  k*=loci1[i];
}
npdat[1]=k-1;
for(io=0;io<obscom;io++) genp(table,io);
*lnl0=ll(table);
lnls=*lnl0;
it=1;
do
{
  for(i=0;i<hapall;i++) hs[i]=h[i];
  geth(table);
  if(!convll)
  {
    s=0;
    for(i=0;i<hapall;i++) s+=fabs(hs[i]-h[i]);
  }
  else
  {
    for(io=0;io<obscom;io++) genp(table,io);
    *lnl1=ll(table);
    s=*lnl1-lnls;
    lnls=*lnl1;
  }
} while(s>eps&&it++<maxit);
if(s>eps) *converge=0; else *converge=1;
*iter=it;
if(!convll)
{
  for(io=0;io<obscom;io++) genp(table,io);
  *lnl1=ll(table);
}
for(io=0;io<obscom;io++) prob[io]=table[io].prob;
if(*verbose)
{
  Rprintf("\nGENECOUNTING, %s (%s) version %.2f \n\t\tJH Zhao 03/01--07/05\n\n",
  (!handlemissing)?"ordinary":"missing-value handling",
  (!xdata)?"autosome data":"X chromosome data",version);
  time(&now);
  Rprintf("%s\n\n",ctime(&now));
  Rprintf("max. loci=%d, max. alleles=%d\n\n",mxloc,mxalleles);
  Rprintf("\nlog-likelihood assuming linkage equilibrium %.2f\n",*lnl0);
  Rprintf("# parameters %d (%d, by # of specified alleles)\n\n",npdat[0],npusr[0]);
  Rprintf("\nlog-likelihood assuming linkage disequilibrium =%.2f\n",*lnl1);
  Rprintf("# parameters %d (%d, by # of specified alleles)\n",npdat[1],npusr[1]);
  x2mm=2*(*lnl1-*lnl0);
  Rprintf("\nTest of marker-marker disequilibrium, chi-square=%.2f\n",x2mm);
}
free(hs);
free(c);
if(handlemissing) free(hm);
ch(table,assignment);
/*
hrtree(hrt);
*/
free(table);
}

void getp(pat *table)
/*allele frequencies*/
{
int i,j,k,l,u;
short d[mxloc+1],loci1[mxloc],d1[mxloc+1];
double s[mxloc],n[mxloc],ss;

for(i=0;i<nloci;i++)
{
  s[i]=n[i]=0;
  for(j=0;j<nalleles;j++) pp[i][j]=pm[i][j]=0;
}
ss=d_msample=0;
for(i=0;i<obscom;i++)
{
  if(xdata&&table[i].sex==MALE) goto x;
  k=0;
  for(j=0;j<nloci;j++)
  {
    l=table[i].hd.d[j][0]-1;
    u=table[i].hd.d[j][1]-1;
    if((l>=loci[j]||u>=loci[j])||(l<0||u<0))
    {
      ++k;
      continue;
    }
    pp[j][l]+=table[i].count;
    pp[j][u]+=table[i].count;
    s[j]+=table[i].count;
  }
  if(k==0) ss+=table[i].count;
  else d_msample+=table[i].count;
  continue;
x:
  k=0;
  for(j=0;j<nloci;j++)
  {
    l=table[i].hd.h[j]-1;
    if(l>=0&&l<=loci[j])
    {
      ++k;
      pm[j][l]+=table[i].count;
      n[j]+=table[i].count;
    }
  }
  if(k==nloci) h_sample+=table[i].count;
  else h_msample+=table[i].count;
}
d_sample=ss;
for(i=0;i<nloci;i++)
{
  for(j=0;j<loci[i];j++)
  {
    if(!xdata) pp[i][j]/=2*s[i];
    else pp[i][j]=(pp[i][j]+pm[i][j])/(2*s[i]+n[i]);
  }
  loci1[i]=loci[i]-1;
}
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<hapall;i++)
{
  ss=1;
  for(j=0;j<nloci;j++) ss*=pp[j][d[j]];
  for(j=0;j<=nloci;j++) d1[j]=d[j]+1;
  l=linenum(loci,d1)-1;
  h[l]=h0[l]=ss;
  digitm(loci1,d,0);
}
}

void counting(pat *table, long int io)
{
long int j,k,l,k1,k2,nhet2;
short la[mxloc],ua[mxloc],hetid[mxloc],nhet,d[mxloc+1];
double s,ej;

l=0;
for(j=0;j<nloci;j++) hetid[j]=0;
for(j=0;j<nloci;j++)
{
  if(table[io].hd.d[j][0]!=table[io].hd.d[j][1])
  {
    hetid[l]=j;
    ++l;
  }
}
nhet=l;
if(nhet>0)
{
  nhet2=(int)pow(2,nhet-1);
  s=0;
  for(j=0;j<=nloci;j++) d[j]=0;
  for(j=0;j<nhet2;j++)
  {
    for(k=nloci-1;k>=0;k--)
    {
      la[k]=table[io].hd.d[k][0];
      ua[k]=table[io].hd.d[k][1];
    }
    for(l=0;l<nhet;l++) if(d[l]==1)
    {
      k=la[hetid[l]];
      la[hetid[l]]=ua[hetid[l]];
      ua[hetid[l]]=k;
    }
    k1=linenum(loci,la)-1;
    k2=linenum(loci,ua)-1;
    s+=2*h[k1]*h[k2];
    digit2(1,d,0);
  }
  for(j=0;j<=nloci;j++) d[j]=0;
  for(j=0;j<nhet2;j++)
  {
    for(k=0;k<nloci;k++)
    {
      la[k]=table[io].hd.d[k][0];
      ua[k]=table[io].hd.d[k][1];
    }
    for(l=0;l<nhet;l++) if(d[l]==1)
    {
      k=la[hetid[l]];
      la[hetid[l]]=ua[hetid[l]];
      ua[hetid[l]]=k;
    }
    k1=linenum(loci,la)-1;
    k2=linenum(loci,ua)-1;
    ej=2*h[k1]*h[k2]/s;
    c[k1]+=ej*table[io].count;
    c[k2]+=ej*table[io].count;
    digit2(1,d,0);
  }
}
else
{
  for(j=0;j<nloci;j++) la[j]=table[io].hd.d[j][0];
  k=linenum(loci,la)-1;
  c[k]+=2*table[io].count;
}
}

void geth(pat *table)
/*haplotype frequencies with missing data (conditioned)*/
{
long int io,i,j,k,l,k1,k2,cycle,ncycle,nhet2;
short loci1[mxloc],l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
short la[mxloc],ua[mxloc],hetid[mxloc],nhet,d[mxloc+1];
double s,tc,ej;

if(!handlemissing)
{
  for(i=0;i<hapall;i++) c[i]=0;
  for(i=0;i<obscom;i++) 
  {
     if((!xdata)||(xdata&&table[i].sex==FEMALE)) counting(table,i);
     else
     {
       k=linenum(loci,table[i].hd.h)-1;
       c[k]+=table[i].count;
     }
  }
  for(i=0;i<hapall;i++) h[i]=c[i]/(2*d_sample+h_sample);
}
else
{
  for(i=0;i<hapall;i++) c[i]=hm[i]=0;
  for(io=0;io<obscom;io++)
  {
    if(xdata&&table[io].sex==MALE)
    {
      j=0;
      cycle=1;
      for(i=0;i<nloci;i++)
      {
         idm[i]=0;
         locim[i]=0;
         k=table[io].hd.h[i];
         if(k<1||k>loci[i])
         {
           idm[i]=1;
           locim[j]=loci[i];
           cycle*=loci[i];
           j++;
         }
      }
      ncycle=cycle;
      nlocim=j;
      if(nlocim==0||!handlemissing)
      {
        k=linenum(loci,table[io].hd.h)-1;
        c[k]+=table[io].count;
        continue;
      }
      tc=0;
      for(i=0;i<=nloci;i++) d[i]=0;
      for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;
      for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;
      do
      {
        for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;
        j=0;
        for(i=0;i<nloci;i++)
        {
          lq[i]=table[io].hd.h[i];
          if(idm[i]==1)
          {
            lq[i]=lk[j];
            j++;
          }
        }
        k=linenum(loci,lq)-1;
        tc+=h[k];
        digitm(loci1,d,0);
      } while(--cycle>0);
      for(i=0;i<=nloci;i++) d[i]=0;
      for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;
      for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;
      do
      {
        for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;
        j=0;
        for(i=0;i<nloci;i++)
        {
          lq[i]=table[io].hd.h[i];
          if(idm[i]==1)
          {
            lq[i]=lk[j];
            j++;
          }
        }
        if(tc==0) s=0;
        else s=table[io].count/tc;
        k=linenum(loci,lq)-1;
        hm[k]+=h[k]*s;
        digitm(loci1,d,0);
      } while(--ncycle>0);
      continue;
    }
    l=0;
    cycle=1;
    for(j=0;j<nloci;j++)
    {
      idm[j]=0;
      locim[j]=0;
      lk[j]=lq[j]=0;
      k1=table[io].hd.d[j][0];
      k2=table[io].hd.d[j][1];
      if((k1<1||k1>loci[j])||(k2<1||k2>loci[j]))
      {
        cycle*=loci[j]*(loci[j]+1)/2;
        idm[j]=1;
        locim[l]=loci[j];
        lk[l]=lq[l]=1;
        ++l;
      }
    }
    ncycle=cycle;
    nlocim=l;
    if(nlocim==0)
    {
      counting(table,io);
      continue;
    }
    /*obtain total weight*/
    tc=0;
    do
    {
      l=0;
      for(j=0;j<nloci;j++)
      {
        l0[j]=table[io].hd.d[j][0];
        l1[j]=table[io].hd.d[j][1];
        if(idm[j]==1)
        {
          l0[j]=lk[l];
          l1[j]=lq[l];
          ++l;
        }
      }
      tc+=phasep(l0,l1);
      if(nlocim==1)
      {
        ++lk[0];
        if(lk[0]>lq[0])
        {
          lk[0]=1;
          ++lq[0];
        }
      }
      else if(nlocim>1)
      {
        ++lk[nlocim-1];
        for(i=nlocim-1;i>0;i--) if(lk[i]>lq[i])
        {
          ++lq[i];
          lk[i]=1;
          if(lq[i]>locim[i])
          {
            lq[i]=1;
            ++lk[i-1];
            if((i==1)&&(lk[i-1]>lq[i-1]))
            {
              ++lq[i-1];
              lk[i-1]=1;
            }
          }
        }
      }
    } while(--cycle>0);
    l=0;
    for(j=0;j<nloci;j++)
    {
      lk[j]=lq[j]=0;
      if(idm[j]==1)
      {
        lk[l]=lq[l]=1;
        ++l;
      }
    }
    /*update haplotype counts*/
    do
    {
      k1=k2=0;
      for(j=0;j<nloci;j++)
      {
        l0[j]=table[io].hd.d[j][0];
        l1[j]=table[io].hd.d[j][1];
        if(idm[j]==1)
        {
          l0[j]=lk[k1];
          l1[j]=lq[k1];
          ++k1;
        }
        hetid[j]=0;
        if(l0[j]!=l1[j])
        {
          hetid[k2]=j;
          ++k2;
        }
      }
      nhet=k2;
      if(tc==0) s=0;
      else s=table[io].count/tc;
      if(nhet>0)
      {
        nhet2=(int)pow(2,nhet-1);
        for(j=0;j<=nloci;j++) d[j]=0;
        for(j=0;j<nhet2;j++)
        {
          for(k=0;k<nloci;k++)
          {
            la[k]=l0[k];
            ua[k]=l1[k];
          }
          for(l=0;l<nhet;l++) if(d[l]==1)
          {
            k=la[hetid[l]];
            la[hetid[l]]=ua[hetid[l]];
            ua[hetid[l]]=k;
          }
          k1=linenum(loci,la)-1;
          k2=linenum(loci,ua)-1;
          ej=2*h[k1]*h[k2];
          hm[k1]+=ej*s;
          hm[k2]+=ej*s;
          digit2(1,d,0);
        }
      }
      else
      {
        k=linenum(loci,l0)-1;
        hm[k]+=2*h[k]*h[k]*s;
      }
      if(nlocim==1)
      {
        ++lk[0];
        if(lk[0]>lq[0])
        {
          lk[0]=1;
          ++lq[0];
        }
      }
      else if(nlocim>1)
      {
        ++lk[nlocim-1];
        for(i=nlocim-1;i>0;i--) if(lk[i]>lq[i])
        {
          ++lq[i];
          lk[i]=1;
          if(lq[i]>locim[i])
          {
            lq[i]=1;
            ++lk[i-1];
            if((i==1)&&(lk[i-1]>lq[i-1]))
            {
              ++lq[i-1];
              lk[i-1]=1;
            }
          }
        }
      }
    } while(--ncycle>0);
  }
  tc=0;
  for(i=0;i<hapall;i++) tc+=hm[i];
  for(i=0;i<hapall;i++) h[i]=(c[i]+hm[i])/(d_sample*2+tc+h_sample);
}
/* if calling geth() for observed data first
tc*=0.5;
for(i=0;i<hapall;i++) h[i]=(h[i]*d_sample+hm[i]/2)/(d_sample+tc);
*/
}

double phasep(short *l0,short *l1)
/*phase probability*/
{
long int k,k1,k2,nhet2;
short i,j,l,nhet,hetid[mxloc],d[mxloc+1],la[mxloc],ua[mxloc];
double s;

for(j=0;j<nloci;j++) hetid[j]=0;
l=0;
for(i=0;i<nloci;i++) if(l0[i]!=l1[i])
{
  hetid[l]=i;
  ++l;
}
nhet=l;
if(nhet>0)
{
  nhet2=(int)pow(2,nhet-1);
  s=0;
  for(j=0;j<=nloci;j++) d[j]=0;
  for(j=0;j<nhet2;j++)
  {
    for(k=0;k<nloci;k++)
    {
      la[k]=l0[k];
      ua[k]=l1[k];
    }
    for(l=0;l<nhet;l++) if(d[l]==1)
    {
      k=la[hetid[l]];
      la[hetid[l]]=ua[hetid[l]];
      ua[hetid[l]]=k;
    }
    k1=linenum(loci,la)-1;
    k2=linenum(loci,ua)-1;
    s+=2*h[k1]*h[k2];
    digit2(1,d,0);
  }
}
else
{
  k=linenum(loci,l0)-1;
  s=h[k]*h[k];
}

return s;
}

void genp(pat *table,long int io)
/*get genotype probability */
{
long int i,j,k,l,k1,k2,nhet2;
long int cycle;
short loci1[mxloc],l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
short la[mxloc],ua[mxloc],hetid[mxloc],nhet,d[mxloc+1];
double tc,s;

if(xdata&&table[io].sex==MALE)
{
  j=0;
  cycle=1;
  for(i=0;i<nloci;i++)
  {
     idm[i]=0;
     locim[i]=0;
     k=table[io].hd.h[i];
     if(k<1||k>loci[i])
     {
       idm[i]=1;
       locim[j]=loci[i];
       cycle*=loci[i];
       j++;
     }
  }
  nlocim=j;
  if(nlocim==0)
  {
    k=linenum(loci,table[io].hd.h)-1;
    tc=h[k];
  }
  else
  {
    tc=0;
    for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;
    for(i=0;i<=nloci;i++) d[i]=0;
    for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;
    do
    {
      for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;
      j=0;
      for(i=0;i<nloci;i++)
      {
        lq[i]=table[io].hd.h[i];
        if(idm[i]==1)
        {
          lq[i]=lk[j];
          j++;
        }
      }
      k=linenum(loci,lq)-1;
      tc+=h[k];
      digitm(loci1,d,0);
    } while(--cycle>0);
  }
  table[io].prob=tc;
  return;
}
l=0;
cycle=1;
for(j=0;j<nloci;j++)
{
  idm[j]=0;
  locim[j]=0;
  lk[j]=lq[j]=0;
  k1=table[io].hd.d[j][0];
  k2=table[io].hd.d[j][1];
  if((k1<1||k1>loci[j])||(k2<1||k2>loci[j]))
  {
    cycle*=loci[j]*(loci[j]+1)/2;
    idm[j]=1;
    locim[l]=loci[j];
    lk[l]=lq[l]=1;
    ++l;
  }
}
nlocim=l;
if(nlocim>0)
{
  s=0;
  do
  {
    l=0;
    for(j=0;j<nloci;j++)
    {
      l0[j]=table[io].hd.d[j][0];
      l1[j]=table[io].hd.d[j][1];
      if(idm[j]==1)
      {
        l0[j]=lk[l];
        l1[j]=lq[l];
        ++l;
      }
    }
    s+=phasep(l0,l1);
    if(nlocim==1)
    {
      ++lk[0];
      if(lk[0]>lq[0])
      {
        lk[0]=1;
        ++lq[0];
      }
    }
    else if(nlocim>1)
    {
      ++lk[nlocim-1];
      for(i=nlocim-1;i>0;i--) if(lk[i]>lq[i])
      {
        ++lq[i];
        lk[i]=1;
        if(lq[i]>locim[i])
        {
          lq[i]=1;
          ++lk[i-1];
          if((i==1)&&(lk[i-1]>lq[i-1]))
          {
            ++lq[i-1];
            lk[i-1]=1;
          }
        }
      }
    }
  } while(--cycle>0);
}
else
{
  l=0;
  for(j=0;j<nloci;j++) hetid[j]=0;
  for(j=0;j<nloci;j++)
  {
    if(table[io].hd.d[j][0]!=table[io].hd.d[j][1])
    {
      hetid[l]=j;
      ++l;
    }
  }
  nhet=l;
  if(nhet>0)
  {
    nhet2=(int)pow(2,nhet-1);
    s=0;
    for(j=0;j<=nloci;j++) d[j]=0;
    for(j=0;j<nhet2;j++)
    {
      for(k=nloci-1;k>=0;k--)
      {
        la[k]=table[io].hd.d[k][0];
        ua[k]=table[io].hd.d[k][1];
      }
      for(l=0;l<nhet;l++) if(d[l]==1)
      {
        k=la[hetid[l]];
        la[hetid[l]]=ua[hetid[l]];
        ua[hetid[l]]=k;
      }
      k1=linenum(loci,la)-1;
      k2=linenum(loci,ua)-1;
      s+=2*h[k1]*h[k2];
      digit2(1,d,0);
    }
  }
  else
  {
    for(j=0;j<nloci;j++) la[j]=table[io].hd.d[j][0];
    k=linenum(loci,la)-1;
    s=h[k]*h[k];
  }
}
table[io].prob=s;
}

double ll(pat *table)
/*log-likelihood*/
{
int i,j,l;
long int k,k1,k2;
double t,lnl,xlnl;

lnl=xlnl=0;
for(i=0;i<obscom;i++)
{
  if(xdata&&table[i].sex==MALE)
  {
    l=0;
    for(j=0;j<nloci;j++)
    {
      k=table[i].hd.h[j];
      if(k<1||k>loci[j]) ++l;
    }
    if(l>0&&!handlemissing) continue;
    t=table[i].count;
    if(t!=0&&table[i].prob>0) xlnl+=t*log(table[i].prob);
    continue;
  }
  l=0;
  for(j=0;j<nloci;j++)
  {
    k1=table[i].hd.d[j][0];
    k2=table[i].hd.d[j][1];
    if((k1<1||k1>loci[j])||(k2<1||k2>loci[j])) ++l;
  }
  if(l>0&&!handlemissing) continue;
  t=table[i].count;
  if(t!=0&&table[i].prob>0) lnl+=t*log(table[i].prob);
}
return (lnl+xlnl);
}

void ch(pat *table,char **assignment)
/*trace observed haplotypes*/
{
FILE *fo;
long int i,j,l,a1,a2;
long int k,k1,k2,nhet2;
short k0,la[mxloc],ua[mxloc];
short hetid[mxloc],nhet;
short d[mxloc+1],loci1[mxloc];
long int cycle, ncycle;
short l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
double tc,ej;

fo=fopen(*assignment,"w");
if(!fo) {
  REprintf("I cannot open file assign.dat for posterior probabilities\n");
  return;
}
for(i=0;i<obscom;i++)
{
  if(xdata&&table[i].sex==MALE)
  {
    l=0;
    cycle=1;
    for(j=0;j<nloci;j++)
    {
       idm[j]=0;
       locim[j]=0;
       k=table[i].hd.h[j];
       if(k<1||k>loci[j])
       {
         idm[j]=1;
         locim[l]=loci[j];
         cycle*=loci[j];
         l++;
       }
    }
    nlocim=l;
    ncycle=cycle;
    if(nlocim==0)
    {
      k=linenum(loci,table[i].hd.h)-1;
      fprintf(fo,"%5ld [1]",table[i].gid);
      for(j=0;j<nloci;j++) fprintf(fo," %2d",table[i].hd.h[j]);
      fprintf(fo," %f %ld\n",1.0,k+1);
/*
      if(!hrt) hrt=hitree(hrt,k+1,table[i].hd.h);
      htr_table[i][k]=1.0;
*/
    }
    else
    {
      tc=0;
      for(j=0;j<nloci;j++) lk[j]=loci1[j]=0;
      for(j=0;j<=nloci;j++) d[j]=0;
      for(j=0;j<nlocim;j++) loci1[j]=locim[nlocim-j-1]-1;
      do
      {
        for(j=nlocim-1;j>=0;j--) lk[nlocim-j-1]=d[j]+1;
        l=0;
        for(j=0;j<nloci;j++)
        {
          lq[j]=table[i].hd.h[j];
          if(idm[j]==1)
          {
            lq[j]=lk[l];
            l++;
          }
        }
        k=linenum(loci,lq)-1;
        tc+=h[k];
        digitm(loci1,d,0);
      } while(--cycle>0);
      for(j=0;j<=nloci;j++) d[j]=0;
      for(j=0;j<nloci;j++) lk[j]=loci1[j]=0;
      for(j=0;j<nlocim;j++) loci1[j]=locim[nlocim-j-1]-1;
      do
      {
        for(j=nlocim-1;j>=0;j--) lk[nlocim-j-1]=d[j]+1;
        l=0;
        for(j=0;j<nloci;j++)
        {
          lq[j]=table[i].hd.h[j];
          if(idm[j]==1)
          {
            lq[j]=lk[l];
            l++;
          }
        }
        k=linenum(loci,lq)-1;
        if(tc==0) ej=0;
        else ej=h[k]/tc;
        if(ej>pl)
        {
          fprintf(fo,"%5ld [1]",table[i].gid);
          for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",lq[a1]);
          fprintf(fo," %f %ld\n",ej,k+1);
/*
          htr_table[i][k+1]+=ej;
*/
        }
        digitm(loci1,d,0);
      } while(--ncycle>0);
    }
    continue;
  }
  l=0;
  cycle=1;
  for(j=0;j<nloci;j++)
  {
    idm[j]=0;
    locim[j]=0;
    lk[j]=lq[j]=0;
    k1=table[i].hd.d[j][0];
    k2=table[i].hd.d[j][1];
    if((k1<1||k1>loci[j])||(k2<1||k2>loci[j]))
    {
      cycle*=loci[j]*(loci[j]+1)/2;
      idm[j]=1;
      locim[l]=loci[j];
      lk[l]=lq[l]=1;
      ++l;
    }
  }
  nlocim=l;
  tc=table[i].prob;
  if(nlocim>0) do
  {
    k1=k2=0;
    for(j=0;j<nloci;j++)
    {
      l0[j]=table[i].hd.d[j][0];
      l1[j]=table[i].hd.d[j][1];
      if(idm[j]==1)
      {
        l0[j]=lk[k1];
        l1[j]=lq[k1];
        ++k1;
      }
      hetid[j]=0;
      if(l0[j]!=l1[j])
      {
        hetid[k2]=j;
        ++k2;
      }
    }
    nhet=k2;
    if(nhet>0)
    {
      nhet2=(int)pow(2,nhet-1);
      for(j=0;j<=nloci;j++) d[j]=0;
      for(j=0;j<nhet2;j++)
      {
        for(k=0;k<nloci;k++)
        {
          la[k]=l0[k];
          ua[k]=l1[k];
        }
        for(l=0;l<nhet;l++) if(d[l]==1)
        {
          k=la[hetid[l]];
          la[hetid[l]]=ua[hetid[l]];
          ua[hetid[l]]=k;
        }
        k1=linenum(loci,la)-1;
        k2=linenum(loci,ua)-1;
        if(tc==0) ej=0;
        else ej=2*h[k1]*h[k2]/tc;
        if(ej>pl)
        {
          fprintf(fo,"%5ld [1]",table[i].gid);
          for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",la[a1]);
          fprintf(fo," %f %ld\n",ej,k1+1);
          fprintf(fo,"%5ld [2]",table[i].gid);
          for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",ua[a2]);
          fprintf(fo," %f %ld\n",ej,k2+1);
/*
          htr_table[i][k1+1]+=ej;
          htr_table[i][k2+1]+=ej;
*/
        }
        digit2(1,d,0);
      }
    }
    else
    {
      k=linenum(loci,l0)-1;
      if(tc==0) ej=0;
      else ej=h[k]*h[k]/tc;
      if(ej>pl)
      {
        fprintf(fo,"%5ld [1]",table[i].gid);
        for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",l0[a1]);
        fprintf(fo," %f %ld\n",ej,k+1);
        fprintf(fo,"%5ld [2]",table[i].gid);
        for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",l0[a2]);
        fprintf(fo," %f %ld\n",ej,k+1);
/*
        htr_table[i][k+1]+=2.0*ej;
*/
      }
    }
    if(nlocim==1)
    {
      ++lk[0];
      if(lk[0]>lq[0])
      {
        lk[0]=1;
        ++lq[0];
      }
    }
    else if(nlocim>1)
    {
      ++lk[nlocim-1];
      for(k=nlocim-1;k>0;k--) if(lk[k]>lq[k])
      {
        ++lq[k];
        lk[k]=1;
        if(lq[k]>locim[k])
        {
          lq[k]=1;
          ++lk[k-1];
          if((k==1)&&(lk[k-1]>lq[k-1]))
          {
            ++lq[k-1];
            lk[k-1]=1;
          }
        }
      }
    }
  } while(--cycle>0);
  else
  {
    l=0;
    for(j=0;j<nloci;j++)
    {
      if(table[i].hd.d[j][0]!=table[i].hd.d[j][1])
      {
        hetid[l]=j;
        ++l;
      }
    }
    nhet=l;
    if(nhet>0)
    {
      nhet2=(int)pow(2,nhet-1);
      for(j=0;j<=nloci;j++) d[j]=0;
      for(j=0;j<nhet2;j++)
      {
        for(k=nloci-1;k>=0;k--)
        {
          la[k]=table[i].hd.d[k][0];
          ua[k]=table[i].hd.d[k][1];
        }
        for(l=0;l<nhet;l++) if(d[l]==1)
        {
          k0=la[hetid[l]];
          la[hetid[l]]=ua[hetid[l]];
          ua[hetid[l]]=k0;
        }
        k1=linenums(loci,la)-1;
        k2=linenums(loci,ua)-1;
        if(tc==0) ej=0;
        else ej=2.0*h[k1]*h[k2]/tc;
        if(ej>pl)
        {
          fprintf(fo,"%5ld [1]",table[i].gid);
          for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",la[a1]);
          fprintf(fo," %f %ld\n",ej,k1+1);
          fprintf(fo,"%5ld [2]",table[i].gid);
          for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",ua[a2]);
          fprintf(fo," %f %ld\n",ej,k2+1);
/*
          htr_table[i][k1+1]+=ej;
          htr_table[i][k2+1]+=ej;
*/
        }
/*
        if(!hrt) hrt=hitree(hrt,k1+1,la);
        else hitree(hrt,k1+1,la);
        if(!hrt) hrt=hitree(hrt,k2+1,ua);
        else hitree(hrt,k2+1,ua);
*/
        digit2(1,d,0);
      }
    } else {
      for(j=0;j<nloci;j++) la[j]=table[i].hd.d[j][0];
      k=linenums(loci,la)-1;
      fprintf(fo,"%5ld [1]",table[i].gid);
      for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",la[a1]);
      fprintf(fo," %f %ld\n",1.0,k+1);
      fprintf(fo,"%5ld [2]",table[i].gid);
      for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",la[a2]);
      fprintf(fo," %f %ld\n",1.0,k+1);
/*
      htr_table[i][k+1]=2.0;
      if(!hrt) hrt=hitree(hrt,k+1,la);
      else hitree(hrt,k+1,la);
      if(!hrt) hrt=hitree(hrt,k+1,ua);
      else hitree(hrt,k+1,ua);
*/
    }
  }
}
fclose(fo);
}

hnode *hitree(hnode *r,long int id,short l[mxloc])
/*insert*/
{
int i;
if (r==NULL)
{
  r=malloc(sizeof(hnode));
  r->left=r->right=NULL;
  r->id=id;
  r->n=0;
  for(i=0;i<nloci;i++) r->l[i]=l[i];
}
else
if (id<r->id) r->left=hitree(r->left,id,l);
else if (id>r->id) r->right=hitree(r->right,id,l);
else r->n++;
return r;
}

void hptree(FILE *fo, hnode *r,long int *l)
/*inorder*/
{
int i;

if (!r) return;
hptree(fo,r->left,l);
++*l;
fprintf(fo," %.6f [%.12f]",h0[r->id-1],h0[r->id-1]);
fprintf(fo," %.6f [%.12f]",h[r->id-1],h[r->id-1]);
for(i=0;i<nloci;i++) fprintf(fo," %2hd",r->l[i]);
fprintf(fo," %ld\n",r->id);
hptree(fo,r->right,l);
}

void hrtree(hnode *t)
/*postorder*/
{
if (!t) return;
hrtree(t->left);
hrtree(t->right);
free(t);
}

void digit2(int radix, short d[], int i)
/*binary routine*/
{
if(d[i]<radix)
{
  ++d[i];
  return;
}
else
{
  d[i]=0;
  ++d[i+1];
  if(d[i+1]<=radix) return;
}
digit2(radix,d,i+1);
}

void digitm(short radix[], short d[], int i)
/*mixed-radix routine*/
{
if(d[i]<radix[i])
{
  ++d[i];
  return;
}
else
{
  d[i]=0;
  ++d[i+1];
  if(d[i+1]<=radix[i+1]) return;
}
digitm(radix,d,i+1);

}

int linenum(int *loci, short *ai)
/*haplotype array index*/
{
int loc,j;
loc=0;
for(j=1;j<=nloci;j++)
  if(j==nloci) loc+=ai[j-1];
  else loc=(loc+ai[j-1]-1)*loci[j];
return(loc);
}

int linenums(int *loci, short *ai)
/*array index of existing haplotypes*/
{
int loc,j;
loc=0;
for(j=1;j<=nloci;j++)
  if(j==nloci) loc+=ai[j-1];
  else loc=(loc+ai[j-1]-1)*loci[j];
return(loc);
}

void hilo(short *a, short *b)
{
short temp;
temp = *a;*a = *b;*b = temp;
}

/*code for dynamic allocation*/

/*
 * pass various tests after s/t initialized
 * JH Zhao 16-17/5/2000, 22/03/2001
 */

void *xmalloc(long len)
{
void *mem;

mem=malloc(len);
if(mem) return mem;
REprintf("Sorry, but I cannot allocate memory\n");
error("%d",-1);

}

typedef struct {
  int id,count;
  short l[mxloc];
  double p;
} haploid;
typedef struct {
  int id,count;
  short l[mxloc][2];
} diploid;

int xgeth(haploid *);
void xgenp(haploid *);
double xll(haploid *);

#define checkid() {\
    Rprintf("[%5d]",men[io].id); \
    for(i=0;i<nlocim;i++) Rprintf(" %2d",lk[i]);Rprintf("\n");}
#define lli() {\
  l=0;\
  for(j=0;j<nloci;j++)\
  {\
    k=htable[i].l[j];\
    if(k<1||k>loci[j]) ++l;\
  }\
  if(l>0&&!handlemissing) continue;\
  t=htable[i].count;\
  if(t!=0&&htable[i].p>tol) xlnl+=t*log(htable[i].p);}
#define xgethi() {\
for(io=0;io<nhaploid;io++)\
{\
  j=0;\
  cycle=1;\
  for(i=0;i<nloci;i++)\
  {\
     idm[i]=0;\
     locim[i]=0;\
     k=htable[io].l[i];\
     if(k<1||k>loci[i])\
     {\
       idm[i]=1;\
       locim[j]=loci[i];\
       cycle*=loci[i];\
       j++;\
     }\
  }\
  ncycle=cycle;\
  nlocim=j;\
  if(nlocim==0||!handlemissing)\
  {\
    k=linenum(loci,htable[io].l)-1;\
    c[k]+=htable[io].count;\
    continue;\
  }\
  tc=0;\
  for(i=0;i<=nloci;i++) d[i]=0;\
  for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;\
  for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;\
  do\
  {\
    for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;\
    j=0;\
    for(i=0;i<nloci;i++)\
    {\
      l[i]=htable[io].l[i];\
      if(idm[i]==1)\
      {\
        l[i]=lk[j];\
        j++;\
      }\
    }\
    k=linenum(loci,l)-1;\
    tc+=h[k];\
    digitm(loci1,d,0);\
  } while(--cycle>0);\
  for(i=0;i<=nloci;i++) d[i]=0;\
  for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;\
  for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;\
  do\
  {\
    for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;\
    j=0;\
    for(i=0;i<nloci;i++)\
    {\
      l[i]=htable[io].l[i];\
      if(idm[i]==1)\
      {\
        l[i]=lk[j];\
        j++;\
      }\
    }\
    if(tc==0) s=0;\
    else s=htable[io].count/tc;\
    k=linenum(loci,l)-1;\
    hm[k]+=h[k]*s;\
    digitm(loci1,d,0);\
  } while(--ncycle>0);\
}\
tc=0;\
for(i=0;i<hapall;i++) tc+=hm[i];\
for(i=0;i<hapall;i++) h[i]=(c[i]+hm[i])/(h_sample+tc);}

#define xgenpi() {\
  j=0;\
  cycle=1;\
  for(i=0;i<nloci;i++)\
  {\
     idm[i]=0;\
     locim[i]=0;\
     k=htable[io].l[i];\
     if(k<1||k>loci[i])\
     {\
       idm[i]=1;\
       locim[j]=loci[i];\
       cycle*=loci[i];\
       j++;\
     }\
  }\
  nlocim=j;\
  if(nlocim==0)\
  {\
    k=linenum(loci,htable[io].l)-1;\
    tc=h[k];\
  }\
  else\
  {\
    tc=0;\
    for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;\
    for(i=0;i<=nloci;i++) d[i]=0;\
    for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;\
    do\
    {\
      for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;\
      j=0;\
      for(i=0;i<nloci;i++)\
      {\
        l[i]=htable[io].l[i];\
        if(idm[i]==1)\
        {\
          l[i]=lk[j];\
          j++;\
        }\
      }\
      k=linenum(loci,l)-1;\
      tc+=h[k];\
      digitm(loci1,d,0);\
    } while(--cycle>0);\
  }\
  htable[io].p=tc;}

int xehm(char *xgcdat)
{
FILE *fi;
char record[301],eol[301];
int idn,count,sex,i,j,k,a1,a2,iter;
double s,ss,ft,lnl0,lnl1,lnls,n[mxloc];
short loci1[mxloc],d[mxloc+1],d1[mxloc+1];
haploid *htable;
diploid *dtable;

lnl0=lnl1=lnls=0;
fi=fopen(xgcdat,"r");
if(!fi)
{
  perror("I cannot open the file");
  return 1;
}
ignore_char(fgets(record,300,fi));
nloci=0;
while(sscanf(record,"%d %[^\n]",&loci[nloci++],eol)>1) strcpy(record,eol);
nhaploid=ndiploid=0;
while(fgets(record,300,fi)
      &&sscanf(record,"%d %d %d %[^\n]",&idn,&count,&sex,eol)>3)
if(sex==0) nhaploid++;
else ndiploid++;
rewind(fi);
ignore_char(fgets(record,300,fi));
htable=(haploid*)malloc(nhaploid*sizeof(haploid));
dtable=(diploid*)malloc(ndiploid*sizeof(diploid));
if(!htable||!dtable)
{
  perror("can not allocate memory");
  return 1;
}
for(i=0;i<nloci;i++) n[i]=0;
for(i=0;i<mxloc;i++) for(j=0;j<loci[i];j++) pm[i][j]=0;
i=j=0;
h_sample=h_msample=0;
while(fgets(record,300,fi)
      &&sscanf(record,"%d %d %d %[^\n]",&idn,&count,&sex,eol)>3)
{
  if(sex==0)
  {
    htable[i].id=idn;
    htable[i].count=1;
  }
  else
  {
    dtable[j].id=idn;
  }
  strcpy(record,eol);
  iter=0;
  for(k=0;k<nloci;k++)
  {
    if(sex==0)
    {
      sscanf(record,"%d %[^\n]",&a1,eol);
      htable[i].l[k]=a1;
      if(a1>=1&&a1<=loci[k])
      {
        ++iter;
        pm[k][a1-1]+=htable[i].count;
        n[k]+=htable[i].count;
      }
    }
    else
    {
      sscanf(record,"%d/%d %[^\n]",&a1,&a2,eol);
      dtable[j].l[k][0]=a1;
      dtable[j].l[k][1]=a2;
    }
    strcpy(record,eol);
  }
  if(iter==nloci) h_sample+=htable[i].count;
  else h_msample+=htable[i].count;
  if(sex==0) i++;
  else j++;
}
fclose(fi);
Rprintf("\nNo of nonmissing individuals at each locus\n\n");
for(i=0;i<nloci;i++) Rprintf("%3d %.0f\n",i+1,n[i]);
Rprintf("\n%d/%d individuals with F/P information\n",h_sample,h_msample);
for(i=0;i<nloci;i++) for(j=0;j<loci[i];j++) pm[i][j]/=n[i];
Rprintf("\nAllele frequencies\n\nLocus/Frequecies\n\n");
for(i=0;i<nloci;i++)
{
  Rprintf("%3d",i+1);
  for(j=0;j<loci[i];j++) Rprintf(" %.4f",pm[i][j]);
  Rprintf("\n");
}
hapall=1;
for(i=0;i<nloci;i++) hapall*=loci[i];
hs=(double*)xmalloc(hapall*sizeof(double));
c=(double*)xmalloc(hapall*sizeof(double));
if(handlemissing) hm=(double*)xmalloc(hapall*sizeof(double));
for(i=0;i<hapall;i++) h[i]=h0[i]=0;
for(i=0;i<nloci;++i) loci1[i]=loci[i]-1;
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<hapall;i++)
{
  ss=1;
  for(j=0;j<nloci;j++) ss*=pm[j][d[j]];
  for(j=0;j<=nloci;j++) d1[j]=d[j]+1;
  j=linenum(loci,d1)-1;
  h[j]=h0[j]=ss;
  digitm(loci1,d,0);
}
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<nloci;i++) loci1[i]=loci[nloci-i-1]-1;
k=0;
lnl0=xll(htable);
Rprintf("\nlog-likelihood assuming linkage equilibrium = %.2f\n\n",lnl0);
k=1;
iter=1;
do
{
  for(i=0;i<hapall;i++) hs[i]=h[i];
  xgeth(htable);
  if(!convll)
  {
    s=0;
    for(i=0;i<hapall;i++) s+=fabs(hs[i]-h[i]);
  }
  else
  {
    Rprintf("Iteration %3d, ",k++);
    lnl1=xll(htable);
    Rprintf("log-likelihood=%.2f\n",lnl1);
    s=lnl1-lnls;
    lnls=lnl1;
  }
} while(s>eps&&iter++<maxit);
if(!convll) lnl1=xll(htable);
Rprintf("\nHaplotype frequencies under linkage disequilibrium\n\n");
s=h_sample+h_msample;
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<nloci;i++) loci1[i]=loci[nloci-i-1]-1;
k=0;
for(i=0;i<hapall;i++)
{
  if(h[i]>tol)
  {
    ft=sqrt(h[i]*s)+sqrt(h[i]*s+1)-sqrt(4*h0[i]*s+1);
    Rprintf("%6d %.6f %.12f %.6f %6.2f",k+1,h[i],h[i],h0[i],ft);
    for(j=nloci-1;j>=0;j--) Rprintf("%3d",d[j]+1);
    Rprintf(" %d\n",i+1);
    k++;
  }
  digitm(loci1,d,0);
}
Rprintf("\n%d haplotypes >%f\n\n",k,tol);
Rprintf("log-likelihood assuming linkage disequilibrium = %.2f\n",lnl1);

return 0;
}

int xgeth(haploid *htable)
/*handle missing genotype*/
{
int idm[mxloc],locim[mxloc],lk[mxloc];
int io,i,j,k,cycle,ncycle;
double tc,s;
short loci1[mxloc],l[mxloc],d[mxloc+1];

for(i=0;i<hapall;i++) c[i]=hm[i]=0;
xgethi();

return 0;
}

void xgenp(haploid *htable)
/*obtain genotype probability*/
{
int idm[mxloc],locim[mxloc],lk[mxloc];
int io,i,j,k,cycle;
double tc;
short loci1[mxloc],l[mxloc],d[mxloc+1];

for(io=0;io<nhaploid;io++) xgenpi();

}

double xll(haploid *htable)
/*log-likelihood including missing genotype*/
{
double xlnl=0,t;
int i,j,k,l;

xgenp(htable);
for(i=0;i<nhaploid;i++) lli();

return xlnl;
}
