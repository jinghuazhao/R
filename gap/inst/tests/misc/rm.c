#include <stdio.h>
#include <math.h>
#define ERROR /*Science 273: 1516-17 13SEP1996*/
#undef ERROR  /*Science 275: 1327-30 28FEB1997*/

typedef struct normal{float m,v;} norm;
int fbsize(float,float),pbsize();
float n(norm,int);
float y,h,h1,h2,p_A,lambda_o,lambda_s,n1,n2,n3;
float models[12][2]={
    4.0, 0.01,
    4.0, 0.10,
    4.0, 0.50,
    4.0, 0.80,
    2.0, 0.01,
    2.0, 0.10,
    2.0, 0.50,
    2.0, 0.80,
    1.5, 0.01,
    1.5, 0.10,
    1.5, 0.50,
    1.5, 0.80};

main()
{
float gamma,p;
int i;

/*Risch & Merikangas 1996*/
printf("\nThe family-based result: \n");
printf("\ngamma   p      Y     N_asp    P_A    Het   N_tdt   Het N_asp/tdt  L_o   L_s\n\n");
for(i=0;i<12;i++) {
  gamma=models[i][0];
  p=models[i][1];
  fbsize(gamma,p);
  if((i+1)%4==0) printf("\n");
}
/*APOE-4, Scott WK, Pericak-Vance, MA & Haines JL*/
/*Genetic analysis of complex diseases 1327*/
gamma=4.5;
p=0.15;
printf("Alzheimer's:\n\n");
fbsize(gamma,p);
/*Long AD, Grote MN & Langley CH 1328*/
pbsize();
}

int fbsize(float gamma,float p)
/*Jing Hua Zhao 30-12-98*/
{
norm nl, al, na, aa;
float k,va,vd,q,w;

q=1-p;
k=pow(p*gamma+q,2);
va=2*p*q*pow((gamma-1)*(p*gamma+q),2);
vd=pow(p*q*pow(gamma-1,2),2);

#ifdef DEBUG
printf("K=%f VA=%f VD=%f\n",k,va,vd);
#endif

w=p*q*pow(gamma-1,2)/pow(p*gamma+q,2);
y=(1+w)/(2+w);
lambda_s=pow(1+0.5*w,2);
lambda_o=1+w;
h=h1=p*q*(gamma+1)/(p*gamma+q);
p_A=gamma/(gamma+1);

/*ASP*/
nl.m=0;
nl.v=1;
al.m=2*y-1;
#ifdef ERROR
al.v=0;
#else
al.v=4*y*(1-y);
#endif

n1=n(al,1);

/*TDT*/
na.m=0;
na.v=1;
aa.m=sqrt(h)*(gamma-1)/(gamma+1);
aa.v=1-h*pow((gamma-1)/(gamma+1),2);

n2=n(aa,2);

/*ASP-TDT*/
h=h2=p*q*pow(gamma+1,2)/(2*pow(p*gamma+q,2)+p*q*pow(gamma-1,2));
aa.m=sqrt(h)*(gamma-1)/(gamma+1);
aa.v=1-h*pow((gamma-1)/(gamma+1),2);

n3=n(aa,3);

printf("%4.1f%5.2f|%6.3f%10.0f|%6.3f|%6.3f%8.0f%6.3f%8.0f|%6.2f%6.2f\n",
       gamma,p,y,n1,p_A,h1,n2,h2,n3,lambda_o,lambda_s);
return 0;
}

float n(norm all,int op)
/*m=0,v=1 under the null hypotheses*/
{
float z_alpha,z1_beta;
float s,m,v;

m=all.m;
v=all.v;

z1_beta=-0.84;                       /*1-beta=0.8,-.84162123*/
switch(op) {
case 1:z_alpha=3.72;break;           /*alpha=1E-4,3.7190165*/
default:z_alpha=5.33;break;          /*alpha=5E-8,5.3267239*/
}
s=pow((z_alpha-sqrt(v)*z1_beta)/m,2)/2; /*shared/transmitted for each parent*/
if(op==3) s/=2;                         /*the sample size is halved*/

return (s);
}

#define x2_alpha 29.72 /*Q=29.7168*/
#define z_alpha 5.45   /*5.4513104*/

int pbsize()
{
int i,j;
float z1_beta,gamma,lambda,pi,p,q,n,k,s,kp[3]={0.01,0.05,0.10};
/*
\alpha=5e-8,\beta=0.8, \alpha would give 5\% genome-wide significance level
\lambda is the NCP from the marginal table
\pi is the pr(Affected|aa)
*/
gamma=4.5; /*again for Alzheimer's*/
p=0.15;
q=1-p;
pi=0.065;  /*0.07 generates 163, equivalent to ASP*/
z1_beta=-0.84;

k=pi*pow(gamma*p+q,2);
s=(1-pi*pow(gamma,2))*pow(p,2)+(1-pi*gamma)*2*p*q+(1-pi)*pow(q,2);

/*LGL formula*/
lambda=pi*(pow(gamma,2)*p+q-pow(gamma*p+q,2))/(1-pi*pow(gamma*p+q,2));
/*my own*/
lambda=pi*p*q*pow(gamma-1,2)/(1-pi*pow(gamma*p+q,2));

/*not sure about +/-!*/
n=pow(z1_beta+z_alpha,2)/lambda;

/*may be used to correct for population prevalence*/
printf("\nThe population-based result: Kp=%f Kq=%f n=%10.0f\n",k,s,n);

printf("\nRandom ascertainment with disease prevalence\n");
printf("\n         1%%         5%%        10%%\n\n");
for(i=0;i<12;i++) {
  gamma=models[i][0];
  p=models[i][1];
  q=1-p;
  for(j=0;j<3;j++) {
    k=kp[j];
    pi=k/pow(gamma*p+q,2);
    lambda=pi*p*q*pow(gamma-1,2)/(1-pi*pow(gamma*p+q,2));
    n=pow(z1_beta-z_alpha,2)/lambda;
    printf(" %10.f",n);
  } printf("\n");
  if((i+1)%4==0) printf("\n");
}
printf("This is only an approximation, a more accurate result\n");
printf("can be obtained by Fisher's exact test\n");

return 0;
}

ms_model()
/*Ebers G et al (1996) on multiple sclerosis Nat Genet 13:472*/
/*for modest \lambda_s, two models are the same*/
/*\lambda_s=sibling risk ratio associated with the locus*/
{
float z2,z1,z0,y;

/*model with moderate amount of dominance variance*/
y=1-0.5/sqrt(lambda_s);
z2=y*y;
z1=2*y*(1-y);
z0=pow(1-y,2);

/*additive model, dominance variance eq 0*/
z2=0.5-0.25/lambda_s;
z1=0.5;
z0=0.25/lambda_s;
}
