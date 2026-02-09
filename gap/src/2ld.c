#include <R.h>
#include <stdio.h>
#include <math.h>
void tbyt(double *h, double *haplotypes, double *D, double *VarD,
          double *Dmax, double *VarDmax, double *Dprime, double *VarDprime,
          double *x2, double *lor, double *vlor)
{
double p,q,u,v,t=0;
double a,b,c,d,or,xi;
/*
double ED;
double EDmax;
*/
p=h[0]+h[1];
q=h[2]+h[3];
u=h[0]+h[2];
v=h[1]+h[3];

*D=h[0]-p*u;
/*
ED=(*D)*(*haplotypes-1)/(*haplotypes);
*/
*VarD=(p*q*u*v+(*D)*(q-p)*(v-u)-(*D)*(*D))/(*haplotypes);
if(*D<0) {
  if(p*u<q*v) {
    *Dmax=p*u;
    xi=h[0];
  } else {
    *Dmax=q*v;
    xi=h[3];
  }
} else {
  if(p*v<q*u) {
    *Dmax=p*v;
    xi=h[1];
  } else {
    *Dmax=q*u;
    xi=h[2];
  }
}
/*
EDmax=*Dmax*(*haplotypes-1)/(*haplotypes);
*/
*Dprime=*D/(*Dmax);
if(*Dprime<0) {
   a=v;
   b=u;
} else {
   a=u;
   b=v;
}
*VarDmax=*Dmax*(p*a+q*b-2*fabs(*D))/(*haplotypes);
if(fabs(*Dprime)==1) *VarDprime=0;
else {
  t=*haplotypes**VarD-fabs(*Dprime)**Dmax*(p*a+q*b-2*fabs(*D));
 *VarDprime=((1-fabs(*Dprime))*t+fabs(*Dprime)*xi*(1-xi))/(*haplotypes)/(*Dmax)/(*Dmax);
}
*x2=*haplotypes*(*D)*(*D)/p/q/u/v;
a=*haplotypes*h[0]+0.5;
b=*haplotypes*h[1]+0.5;
c=*haplotypes*h[2]+0.5;
d=*haplotypes*h[3]+0.5;
or=a*d/b/c;
*lor=log(or);
*vlor=1/a+1/b+1/c+1/d;
return;

}

#define maxalleles 1000
#define maxgenotypes maxalleles*(maxalleles+1)/2
#define minval(a,b) ((a<b)?a:b)

enum {PEARSON,TSCHUPROW,CRAMER};

static int alleles1,alleles2;
static double p[maxalleles],q[maxalleles];
static double sample_size,z1,z2;
static int Dmaxij[maxalleles*maxalleles];

void kbylem(double *,int *,int *,double *,double *,double *);
void abp(int,int,double*,double*,double*,double*);

void kbyl(int *nalleles1, int *nalleles2, double *h, double *haplotypes,
          double *Dp, double *VarDp, double *Dijtable, double *VarDijtable,
          double *X2table,
          double *Dmaxtable, double *Dijptable, double *VarDijptable,
          double *x2, double *seX2, double *rho, double *seR, int *optrho,
          double *klinfo, int *verbose)
{
double Dij,VarDij,Dmax,Xij,Dijp,VarDijp,a,b,Eijtable[maxalleles*maxalleles];
double ai,aip,bj,bjp,ak,akp,bl,blp;
double AI,AIP,BJ,BJP;
double W=0,t=0,tt[4]={0,0,0,0},sgn=0;
double Pij,Pil,Pkj,Dkj,Dkl,Dil,Eij,Eil,Ekj,Ekl;
int i,j,k,l,ij,ik,jl,il,jk;
double po,pe;

double phi2=0,p2[maxalleles],q2[maxalleles];
double s0,s1,s2,s3,s4;

if(*verbose==1) Rprintf("Maximum number of alleles = %d\n",maxalleles);
alleles1=(*nalleles1);
alleles2=(*nalleles2);
for(i=0;i<alleles1;i++) {
  t=0;
  for(j=0;j<alleles2;j++) {
    t+=h[i*alleles2+j];
  }
  p[i]=t;
}
for(j=0;j<alleles2;j++) {
  t=0;
  for(i=0;i<alleles1;i++) {
    t+=h[i*alleles2+j];
  }
  q[j]=t;
}

*Dp=0;
for(i=0;i<alleles1;i++) {
  for(j=0;j<alleles2;j++) {
     ij=i*alleles2+j;
     Dijtable[i*alleles2+j]=Dij=h[ij]-p[i]*q[j];
     t=p[i]*(1-p[i])*q[j]*(1-q[j]); /*+Dij*(1-2*p[i])*(1-2*q[j])-Dij*Dij;*/
     VarDijtable[i*alleles2+j]=VarDij=t/(*haplotypes);
     tt[0]=p[i]*q[j];
     tt[1]=(1-p[i])*(1-q[j]);
     tt[2]=p[i]*(1-q[j]);
     tt[3]=q[j]*(1-p[i]);
     if(Dij<0) {
       a=1-q[j];
       b=q[j];
       if(tt[0]<tt[1]) k=0;
       else k=1;
     } else {
       a=q[j];
       b=1-q[j];
       if(tt[2]<tt[3]) k=2;
       else k=3;
     }
     Dmax=tt[k];
     Dmaxij[i*alleles2+j]=k;
     t=h[ij];
     switch(k) {
       case 0:Xij=t;break;
       case 1:Xij=1-p[i]-q[j]+t;break;
       case 2:Xij=p[i]-t;break;
       case 3:Xij=q[j]-t;break;
       default:Xij=0;
     }
     abp(i,j,&ai,&bj,&aip,&bjp);
     Dijp=VarDijp=0;
     if(Dmax!=0)
     {
       Dijp=Dij/Dmax;
       if(h[ij]==0||1-p[i]-q[j]+h[ij]==0) Dijp=-1;
       if(p[i]-h[ij]==0||q[j]-h[ij]==0) Dijp=1;
     }
     if(fabs(Dijp)==1) VarDijp=0;
     else if(Dmax!=0) {
       t=(*haplotypes)*VarDij-fabs(Dijp)*Dmax*(a*p[i]+b*(1-p[i])-2*fabs(Dij));
       VarDijp=((1-fabs(Dijp))*t+fabs(Dijp)*Xij*(1-Xij))/(*haplotypes)/Dmax/Dmax;
     }
     Dijptable[i*alleles2+j]=Dijp;
     VarDijptable[i*alleles2+j]=VarDijp;
     *Dp+=p[i]*q[j]*fabs(Dijp);
     Eijtable[i*alleles2+j]=(p[i]*aip*bj+q[j]*ai*bjp)*Dij+ai*bj*(Dij-p[i]*q[j]);
   }
}
for(j=0;j<alleles2;j++) 
for(i=0;i<alleles1;i++) {
  for(j=0;j<alleles2;j++) {
     tt[0]=p[i]*q[j];
     tt[1]=(1-p[i])*(1-q[j]);
     tt[2]=p[i]*(1-q[j]);
     tt[3]=q[j]*(1-p[i]);
     Dmaxtable[i*alleles2+j]=tt[Dmaxij[i*alleles2+j]];
   }
}

*VarDp=0;
for(i=0;i<alleles1;i++)
for(j=0;j<alleles2;j++)
for(k=0;k<alleles1;k++)
for(l=0;l<alleles2;l++) {
  Pij=h[i*alleles2+j];
  Pil=h[i*alleles2+l];
  Pkj=h[k*alleles2+j];
  Dij=Dijtable[i*alleles2+j];
  Dil=Dijtable[i*alleles2+l];
  Dkj=Dijtable[k*alleles2+j];
  Dkl=Dijtable[k*alleles2+l];
  Eij=Eijtable[i*alleles2+j];
  Eil=Eijtable[i*alleles2+l];
  Ekj=Eijtable[k*alleles2+j];
  Ekl=Eijtable[k*alleles2+l];
  abp(i,j,&ai,&bj,&aip,&bjp);
  if(i==k&&j==l) {
    t=Pij*pow((aip*bj+ai*bjp)*Dij+ai*bj*(1-p[i]-q[j]),2)
     +(p[i]-Pij)*pow(aip*bj*Dij-ai*bj*q[j],2)
     +(q[j]-Pij)*pow(ai*bjp*Dij-ai*bj*p[i],2)-Eij*Eij;
     sgn=1;
  }
  else if(i==k&&j!=l) {
    abp(i,l,&AI,&bl,&AIP,&blp);
    t=Pij*((aip*bj+ai*bjp)*Dij+ai*bj*(1-p[i]-q[j]))*(AIP*bl*Dil-AI*bl*q[l])
     +Pil*(aip*bj*Dij-ai*bj*q[j])*((AIP*bl+AI*blp)*Dil+AI*bl*(1-p[i]-q[l]))
     +(p[i]-Pij-Pil)*(aip*bj*Dij-ai*bj*q[j])*(AIP*bl*Dil-AI*bl*q[l])-Eij*Eil;
     if((Dij>0&&Dil>0)||(Dij<0&&Dil<0)) sgn=1;
     else sgn=-1;
  }
  else if(i!=k&&j==l) {
    abp(k,j,&ak,&BJ,&akp,&BJP);
    t=Pij*((aip*bj+ai*bjp)*Dij+ai*bj*(1-p[i]-q[j]))*(ak*BJP*Dkj-ak*BJ*p[k])
     +Pkj*(ai*bjp*Dij-ai*bj*p[i])*((akp*BJ+ak*BJP)*Dkj+ak*BJ*(1-p[k]-q[j]))
     +(q[j]-Pij-Pkj)*(ai*bjp*Dij-ai*bj*p[i])*(ak*BJP*Dkj-ak*BJ*p[k])-Eij*Ekj;
     if((Dij>0&&Dkj>0)||(Dij<0&&Dkj<0)) sgn=1;
     else sgn=-1;
  }
  else if(i!=k&&j!=l) {
    abp(k,l,&ak,&bl,&akp,&blp);
    t=Pkj*(akp*bl*Dkl-ak*bl*q[l])*(ai*bjp*Dij-ai*bj*p[i])
     +Pil*(aip*bj*Dij-ai*bj*q[j])*(ak*blp*Dkl-ak*bl*p[k])-Eij*Ekl;
     if((Dij>0&&Dkl>0)||(Dij<0&&Dkl<0)) sgn=1;
     else sgn=-1;
  }
  *VarDp+=sgn*t;
}
*VarDp/=(*haplotypes);
W=0;
for(i=0;i<alleles1;i++) {
  for(j=0;j<alleles2;j++) {
     if(p[i]==0||q[j]==0) continue;
     else {
       Dij=Dijtable[i*alleles2+j];
       X2table[i*alleles2+j]=Dij*Dij/VarDijtable[i*alleles2+j];
       W+=Dij*Dij/p[i]/q[j];
     }
  }
}

*x2=(*haplotypes)*W;
W=sqrt(W);
z1=0;
for(i=0;i<alleles1;i++) if(p[i]!=0) z1++;
z2=0;
for(j=0;j<alleles2;j++) if(p[j]!=0) z2++;
t=(z1<z2)?z1:z2;

for(i=0;i<alleles1;i++) p2[i]=0;
for(j=0;j<alleles2;j++) q2[j]=0;

s1=s2=s3=s4=0;

phi2=0;
for(i=0;i<alleles1;i++)
  for(j=0;j<alleles2;j++) {
  if(p[i]==0||q[j]==0) s0=0;
  else s0=pow(h[i*alleles2+j],2)/p[i]/q[j];
  p2[i]+=s0;
  q2[j]+=s0;
  phi2+=s0;
  if(p[i]==0||q[j]==0) s0=0;
  else s0=pow(h[i*alleles2+j],3)/pow(p[i]*q[j],2);
  s1+=s0;
}
phi2-=1;

for(i=0;i<alleles1;i++) if(p[i]!=0) s2+=p2[i]*p2[i]/p[i];
for(j=0;j<alleles2;j++) if(q[j]!=0) s3+=q2[j]*q2[j]/q[j];

for(i=0;i<alleles1;i++)
  for(j=0;j<alleles2;j++) {
    if(p[i]!=0&&q[j]!=0) s4+=h[i*alleles2+j]/p[i]/q[j]*p2[i]*q2[j];
}
t=sqrt((4*s1-3*s2-3*s3+2*s4)/(*haplotypes));
*seR=t;
*seX2=t;
switch(*optrho)
{
case PEARSON:
     *rho=sqrt(phi2/(phi2+1));
     *seR*=0.5/sqrt(phi2)/pow(1+phi2,1.5);
     break;
case TSCHUPROW:
     s0=(z1-1)*(z2-1);
     *rho=sqrt(phi2/sqrt(s0));
     *seR*=0.5/s0/(*rho);
     break;
case CRAMER:
     s0=(z1<z2)?(z1-1):(z2-1);
     *rho=sqrt(phi2/s0);
     *seR*=0.5/sqrt(s0)/(*rho);
     break;
default:;
}

*klinfo=0;
for(i=0;i<alleles1;i++) {
  for(j=0;j<=i;j++) {
    for(k=0;k<alleles2;k++) {
      for(l=0;l<=k;l++) {
        ik=i*alleles2+k;
        jl=j*alleles2+l;
        if((i!=j)&&(k!=l)) {
          il = i * alleles2 + l;
          jk = j * alleles2 + k;
          po = 2.0 * (h[ik] * h[jl] + h[il] * h[jk]);
          pe=2.0*(p[i]*q[k]*p[j]*q[l]+p[i]*q[l]*p[j]*q[k]);
        } else {
          if((i==j)&&(k==l)) {
            po = h[ik]*h[ik];
            pe = p[i]*q[k]*p[i]*q[k];
          }
          else {
            po = 2.0*h[ik]*h[jl];
            pe = 2.0*p[i]*q[k]*p[j]*q[l];
          }
        }
        if(po!=0&&pe!=0) *klinfo+=po*log(po/pe);
      }
    }
  }
}

return;
}

void abp(int i, int j, double *a, double *b, double *ap, double *bp)
{
double pi;
double qj;

pi=p[i];
qj=q[j];
switch (Dmaxij[i*alleles2+j]) {
case 0:
     *a=*b=1;
     *ap=*bp=0;
     break;
case 1:
     *a=pi/(1-pi);*ap=1/(1-pi)/(1-pi);
     *b=qj/(1-qj);*bp=1/(1-qj)/(1-qj);
     break;
case 2:
     *a=1;*ap=0;
     *b=qj/(1-qj);*bp=1/(1-qj)/(1-qj);
     break;
case 3:
     *a=pi/(1-pi);*ap=1/(1-pi)/(1-pi);
     *b=1;*bp=0;
     break;
default:;
}
}

/*
 obtain haplotype frequencies and log-likelihood for two loci
 changed from ASSOCIAT.PAS but trade efficiency for clarification

 19/02/2001 fix bug on k1,k2
 */

void kbylem(double *obs,int *nalleles1,int *nalleles2,double *Rh,double *l0,double *l1)
{
int g2,iter,i,j,k,l,j1,j2,k1,k2,ik,jl,il,jk;
double o12,pobs,e1,e2,r1,r2;
double *h,hc[maxalleles*maxalleles];

alleles1=(*nalleles1);
alleles2=(*nalleles2);
g2=alleles2*(alleles2+1)/2;
sample_size=0;
j1=0;
for(i=0;i<alleles1;i++) {
   for(j=0;j<=i;j++) {
      j2=0;
      for(k=0;k<alleles2;k++) {
         for(l=0;l<=k;l++) {
            o12=obs[j1*g2+j2];
            sample_size+=o12;
            p[i]+=o12;
            p[j]+=o12;
            q[k]+=o12;
            q[l]+=o12;
            j2++;
         }
      }
      j1++;
   }
}
for(i=0;i<alleles1;i++) p[i]/=sample_size*2;
for(j=0;j<alleles2;j++) q[j]/=sample_size*2;
h=Rh;
k=0;
for(i=0;i<alleles1;i++) for(j=0;j<alleles2;j++) {
  h[i*alleles2+j]=p[i]*q[j];
  hc[k++]=0;
}
iter=0;
do {
  *l1=0;
  k1=0;
  for(i=0;i<alleles1;i++) {
    for(j=0;j<=i;j++) {
      k2 = 0;
      for(k=0;k<alleles2;k++) {
        for(l=0;l<=k;l++) {
          ik = i * alleles2 + k;
          jl = j * alleles2 + l;
          o12 = obs[k1*g2+k2];
          if((i!=j)&&(k!=l)) {
            il = i * alleles2 + l;
            jk = j * alleles2 + k;
            r1 = 2.0 * h[ik] * h[jl];
            r2 = 2.0 * h[il] * h[jk];
            pobs=r1+r2;
            if(obs[k1*g2+k2]>0) {
               e1=r1/pobs;
               e2=r2/pobs;
               hc[ik] += e1*o12;
               hc[il] += e2*o12;
               hc[jk] += e2*o12;
               hc[jl] += e1*o12;
            }
          } else {
            if((i==j)&&(k==l)) {
              pobs = h[ik]*h[ik];
              hc[ik] += 2*o12;
            } else {
              pobs = 2.0*h[ik]*h[jl];
              hc[ik] += o12;
              hc[jl] += o12;
            }
          }
          if (o12) *l1+= o12 * log(pobs);
          ++k2;
        }
      }
      ++k1;
    }
  }
  for(i=0;i<alleles1*alleles2;i++) {
    h[i]=hc[i]/sample_size/2;
    hc[i]=0;
  }
  if(iter==0) *l0=*l1;
} while(iter++<15);
}

/*
enum {EHOUTPUT,RAWDATA,CONTINGTABLE};
static int dfobs;
static double x2obs,x2lrt;
static double obs[maxgenotypes][maxgenotypes];

int getobs(char *obsfile)
{
FILE *fp;
char line[301],id[20];
int a1,a2,b1,b2;
int i,j,k,l,l1,l2,u1,u2,nmiss;
double t,rt[maxgenotypes],ct[maxgenotypes];

for(i=0;i<maxalleles;i++) p[i]=0;
for(j=0;j<maxalleles;j++) q[j]=0;
for(i=0;i<maxgenotypes;i++) rt[i]=ct[i]=0;
fp=fopen(obsfile,"r");
if(!fp) {
  Rprintf("Sorry, I cannot open file %s\n",obsfile);
  error("%d",1);
}
if(fgets(line,301,fp)&&sscanf(line,"%s %d %d %d %d",id,&a1,&a2,&b1,&b2)>4)
{
  filetype=RAWDATA;
  rewind(fp);
  alleles1=alleles2=0;
  sample_size=0;
  nmiss=0;
  for(i=0;i<maxgenotypes;i++) for(j=0;j<maxgenotypes;j++) obs[i][j]=0;
  goto rawcodes;
}
rewind(fp);
if(fgets(line,301,fp)&&sscanf(line,"%d%d%f",&alleles1,&alleles2,&sample_size)>2)
filetype=EHOUTPUT;
else {
  rewind(fp);
  if(fgets(line,301,fp)&&sscanf(line,"%d%d",&alleles1,&alleles2)) {
    filetype=CONTINGTABLE;
    sample_size=0;
    l1=0;
    for(i=0;i<alleles1;i++) {
      for(j=0;j<=i;j++) {
        l2=0;
        for(k=0;k<alleles2;k++) {
          for(l=0;l<=k;l++) {
            fscanf(fp,"%f",&obs[l1][l2]);
            sample_size+=obs[l1][l2];
            p[i]+=obs[l1][l2];
            p[j]+=obs[l1][l2];
            q[k]+=obs[l1][l2];
            q[l]+=obs[l1][l2];
            rt[l1]+=obs[l1][l2];
            ct[l2]+=obs[l1][l2];
            l2++;
          }
        }
        l1++;
      }
    }
    goto ok;
  }
  Rprintf("Sorry, but what type of file is this ?\n");
  error("%d",1);
}
if(filetype==EHOUTPUT)
{
  i=1;
  while(fgets(line,301,fp)&&sscanf(line,"%f",&h[i-1])&&(i<alleles1*alleles2)) i++;
  if(i<alleles1*alleles2) Rprintf("I need %d haplotype frequencies\n",alleles1*alleles2);
  Rprintf("%d haplotype frequencies based on a sample of %.0f diploid individuals\n",i,sample_size);
  for(i=0;i<alleles1;i++) {
    t=0;
    for(j=0;j<alleles2;j++) {
      t+=h[i*alleles2+j];
    }
    p[i]=t;
  }
  for(j=0;j<alleles2;j++) {
    t=0;
    for(i=0;i<alleles1;i++) {
      t+=h[i*alleles2+j];
    }
    q[j]=t;
  }
  return 0;
}
rawcodes:
while(fgets(line,301,fp)&&sscanf(line,"%s %d %d %d %d",id,&a1,&a2,&b1,&b2)>3) {
  if((a1!=0)&&(a2!=0)&&(b1!=0)&&(b2!=0)) {
    l1=(a1<a2)?a1:a2;
    u1=(a1>=a2)?a1:a2;
    p[a1-1]++;
    p[a2-1]++;
    if(u1>alleles1) alleles1=u1;
    l2=(b1<b2)?b1:b2;
    u2=(b1>=b2)?b1:b2;
    q[b1-1]++;
    q[b2-1]++;
    if(u2>alleles2) alleles2=u2;
    obs[l1+u1*(u1-1)/2-1][l2+u2*(u2-1)/2-1]++;
    sample_size++;
  } else nmiss++;
}
fclose(fp);
Rprintf("%.0f individuals with full genotypic information and %d not.\n",sample_size,nmiss);
l1=0;
for(i=0;i<alleles1;i++) {
  for(j=0;j<=i;j++) {
    l2=0;
    for(k=0;k<alleles2;k++) {
      for(l=0;l<=k;l++) {
        rt[l1]+=obs[l1][l2];
        ct[l2]+=obs[l1][l2];
        l2++;
      }
    }
    l1++;
  }
}
ok:
x2obs=x2lrt=0;
l1=0;
for(i=0;i<alleles1;i++) {
  for(j=0;j<=i;j++) {
    l2=0;
    for(k=0;k<alleles2;k++) {
      for(l=0;l<=k;l++) {
        if(rt[l1]>0&&ct[l2]>0) {
          t=rt[l1]*ct[l2]/sample_size;
          x2obs+=pow(obs[l1][l2]-t,2)/t;
          if(obs[l1][l2]>0)
          x2lrt+=obs[l1][l2]*log(obs[l1][l2]/rt[l1]/ct[l2]*sample_size);
        }
        l2++;
      }
    }
    l1++;
  }
}
x2lrt*=2;
z1=0;
z2=0;
for(i=0;i<alleles1;i++) {
  p[i]/=sample_size*2;
  if(p[i]!=0) z1++;
}
for(j=0;j<alleles2;j++) {
  q[j]/=sample_size*2;
  if(q[j]!=0) z2++;
}
dfobs=(z1*(z1+1)/2-1)*(z2*(z2+1)/2-1);

return 0;
}

  information from [ij][kl] based on genotype frequencies

double cellinfo(int nhet,int i,int j,int k,int l)
{
double po,pe,ci;
int ik,il,jk,jl;

ik=i*alleles2+k;
il=i*alleles2+l;
jk=j*alleles2+k;
jl=j*alleles2+l;

po=pe=ci=0;
switch(nhet) {
case 0:
      po=h[ik]*h[ik];
      pe=p[i]*q[k]*p[i]*q[k];
      break;
case 1:
      po=2*h[ik]*h[jl];
      pe=2*p[i]*q[k]*p[j]*q[l];
      break;
case 2:
      po=2*(h[ik]*h[jl]+h[il]*h[jk]);
      pe=2*(p[i]*q[k]*p[j]*q[l]+p[i]*q[l]*p[j]*q[k]);
      break;
default:;
}

if(po!=0&&pe!=0) ci=po*log(po/pe);

return ci;
}


*/
