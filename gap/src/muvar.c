/*
  Additive, Dominance and epistatic variances for single and two-locus QTL model.
  Part of the results was reported in Sham P(1998) "Statistics in Human Genetics", Edward Arnold.

  Adapted from the anova.cpp program by JH Zhao on Dec 17 1996

  3/1/2004
 */

#include <R.h>
#include <math.h>

void onelocus (float *y1, float *p1)
{
  float pa,pA;
  float yAA,yAa,yaa;
  float muA,mua;
  float vD,vA;

  pa=*p1;
  pA=1-pa;
  yAA=y1[0];
  yAa=y1[1];
  yaa=y1[2];

  mua=pA*yAa+pa*yaa;
  muA=pA*yAA+pa*yAa;
  vD=pow(pA*pa*(yAA-2*yAa+yaa),2);
  vA=2*pA*pa*pow((pA*(yAA-yAa)+pa*(yAa-yaa)),2);
  Rprintf("Single locus model:\n");
  Rprintf("===================\n");
  Rprintf("mean mua=%f muA=%f\n",mua,muA);
  Rprintf("Dominant/additive variance = %f %f \n\n",vD,vA);
}

void twolocus (float *y12, float *p1, float *p2)
{
  float pa,pA,pb,pB;
  float yAABB,yAaBB,yaaBB;
  float yAABb,yAaBb,yaaBb;
  float yAAbb,yAabb,yaabb;
  float ydddd;
  float yAddd,yddBd,ydadd,ydddb;
  float yAAdd,yAadd,yaadd;
  float yddBB,yddBb,yddbb;
  float SS1,SS2,SS3;
  float vA,vD,vI,vAE,vDE,vAA,vAD,vDD;
  float yAABd,yAAdb,yAaBd,yAadb,yaaBd,yaadb,
        yAdBB,ydaBB,yAdBb,ydaBb,yAdbb,ydabb;
  float yAdBd,yAddb,ydaBd,ydadb;
  float mu,mAA,mAD,mDD;
  float muA,mua,muAA,muAa,muaa;
  float muAB,muab,muAb,muaB;
  float muAAB,muAAb,muAaB,muAab,muaaB,muaab,muABB,muaBB,muABb,muaBb,muAbb,muabb;
  float muB,mub,muBB,muBb,mubb;
  float muAABB,muAABb,muAAbb,muAaBB,muAaBb,muAabb,muaaBB,muaaBb,muaabb;

  pa=*p1;
  pA=1-pa;
  pb=*p2;
  pB=1-pb;
  yAABB=y12[0];yAaBB=y12[3];yaaBB=y12[6];
  yAABb=y12[1];yAaBb=y12[4];yaaBb=y12[7];
  yAAbb=y12[2];yAabb=y12[5];yaabb=y12[8];

  Rprintf("Two locus model: \n");
  Rprintf("================ \n");
  Rprintf("The model parameters \n\n");
  Rprintf("pA=%f pa=%f \npB=%f pb=%f\n yAABB=%f, yAaBB=%f, yaaBB=%f \n yAABb=%f, yAaBb=%f, yaaBb=%f\n yAAbb=%f, yAabb=%f, yaabb=%f\n",
  pA, pa, pB, pb, yAABB, yAaBB, yaaBB, yAABb, yAaBb, yaaBb, yAAbb, yAabb, yaabb);

  ydddd=pow(pA*pB,2)       * yAABB
       +2*pow(pA,2)*pB*pb  * yAABb
       +pow(pA*pb,2)       * yAAbb
       +2*pA*pa*pow(pB,2)  * yAaBB
       +4*pA*pa*pB*pb      * yAaBb
       +2*pA*pa*pow(pb,2)  * yAabb
       +pow(pa*pB,2)       * yaaBB
       +2*pow(pa,2)*pB*pb  * yaaBb
       +pow(pa*pb,2)       * yaabb;

  yAABB=yAABB-ydddd;yAaBB=yAaBB-ydddd;yaaBB=yaaBB-ydddd;
  yAABb=yAABb-ydddd;yAaBb=yAaBb-ydddd;yaaBb=yaaBb-ydddd;
  yAAbb=yAAbb-ydddd;yAabb=yAabb-ydddd;yaabb=yaabb-ydddd;

  Rprintf("The overall mean = %f \n", ydddd);
  Rprintf("The centered means \n yAABB=%f, yAaBB=%f, yaaBB=%f\n yAABb=%f, yAaBb=%f, yaaBb=%f\n yAAbb=%f, yAabb=%f, yaabb=%f\n",
          yAABB, yAaBB, yaaBB, yAABb, yAaBb, yaaBb, yAAbb, yAabb, yaabb);

  yAAdd=pow(pB,2)*yAABB+2*pB*pb*yAABb+pow(pb,2)*yAAbb;
  yAadd=pow(pB,2)*yAaBB+2*pB*pb*yAaBb+pow(pb,2)*yAabb;
  yaadd=pow(pB,2)*yaaBB+2*pB*pb*yaaBb+pow(pb,2)*yaabb;

  yddBB=pow(pA,2)*yAABB+2*pA*pa*yAaBB+pow(pa,2)*yaaBB;
  yddBb=pow(pA,2)*yAABb+2*pA*pa*yAaBb+pow(pa,2)*yaaBb;
  yddbb=pow(pA,2)*yAAbb+2*pA*pa*yAabb+pow(pa,2)*yaabb;

  yAddd=pA*yAAdd+pa*yAadd;
  ydadd=pA*yAadd+pa*yaadd;
  yddBd=pB*yddBB+pb*yddBb;
  ydddb=pB*yddBb+pb*yddbb;

  SS1=pow(pA*pB,2)       * pow(yAABB,2)
     +2*pow(pA,2)*pB*pb  * pow(yAABb,2)
     +pow(pA*pb,2)       * pow(yAAbb,2)
     +2*pA*pa*pow(pB,2)  * pow(yAaBB,2)
     +4*pA*pa*pB*pb      * pow(yAaBb,2)
     +2*pA*pa*pow(pb,2)  * pow(yAabb,2)
     +pow(pa*pB,2)       * pow(yaaBB,2)
     +2*pow(pa,2)*pB*pb  * pow(yaaBb,2)
     +pow(pa*pb,2)       * pow(yaabb,2);

  SS2=pow(pA*pB,2)       * pow(yAAdd+yddBB,2)
     +2*pow(pA,2)*pB*pb  * pow(yAAdd+yddBb,2)
     +pow(pA*pb,2)       * pow(yAAdd+yddbb,2)
     +2*pA*pa*pow(pB,2)  * pow(yAadd+yddBB,2)
     +4*pA*pa*pB*pb      * pow(yAadd+yddBb,2)
     +2*pA*pa*pow(pb,2)  * pow(yAadd+yddbb,2)
     +pow(pa*pB,2)       * pow(yaadd+yddBB,2)
     +2*pow(pa,2)*pB*pb  * pow(yaadd+yddBb,2)
     +pow(pa*pb,2)       * pow(yaadd+yddbb,2);

  SS3=pow(pA*pB,2)       * pow(2*yAddd+2*yddBd,2)
     +2*pow(pA,2)*pB*pb  * pow(2*yAddd+yddBd+ydddb,2)
     +pow(pA*pb,2)       * pow(2*yAddd+2*ydddb,2)
     +2*pA*pa*pow(pB,2)  * pow(yAddd+ydadd+2*yddBd,2)
     +4*pA*pa*pB*pb      * pow(yAddd+ydadd+yddBd+ydddb,2)
     +2*pA*pa*pow(pb,2)  * pow(yAddd+ydadd+2*ydddb,2)
     +pow(pa*pB,2)       * pow(2*ydadd+2*yddBd,2)
     +2*pow(pa,2)*pB*pb  * pow(2*ydadd+yddBd+ydddb,2)
     +pow(pa*pb,2)       * pow(2*ydadd+2*ydddb,2);

  Rprintf("SS1, SS2, SS3 = %f %f %f \n\n", SS1, SS2, SS3);

  vA=SS3;vD=SS2-SS3;vI=SS1-SS2;

  yAABd=pB*yAABB+pb*yAABb;
  yAAdb=pB*yAABb+pb*yAAbb;
  yAaBd=pB*yAaBB+pb*yAaBb;
  yAadb=pB*yAaBb+pb*yAabb;
  yaaBd=pB*yaaBB+pb*yaaBb;
  yaadb=pB*yaaBb+pb*yaabb;
  yAdBB=pA*yAABB+pa*yAaBB;
  ydaBB=pA*yAaBB+pa*yaaBB;
  yAdBb=pA*yAABb+pa*yAaBb;
  ydaBb=pA*yAaBb+pa*yaaBb;
  yAdbb=pA*yAAbb+pa*yAabb;
  ydabb=pA*yAabb+pa*yaabb;

 /*
  A.B.=AAB.+AaB.=AABB+AABb  + AaBB+AaBb
  A.b.=AAb.+Aab.=AAbB+AAbb  + AabB+Aabb
  a.B.=aAB.+aaB.=aABB+aABb  + aaBB+aaBb
  a.b.=aAb.+aab.=aAbB+aAbb  + aabB+aabb
  */

  yAdBd=pA*pB*yAABB+pA*pb*yAABb+pa*pB*yAaBB+pa*pb*yAaBb;
  yAddb=pA*pB*yAABb+pA*pb*yAAbb+pa*pB*yAaBb+pa*pb*yAabb;
  ydaBd=pA*pB*yAaBB+pA*pb*yAaBb+pa*pB*yaaBB+pa*pb*yaaBb;
  ydadb=pA*pB*yAaBb+pA*pb*yAabb+pa*pB*yaaBb+pa*pb*yaabb;

  mu=yAddd*pA+ydadd*pa;
  muA=yAddd;muB=yddBd;
  mua=ydadd;mub=ydddb;

  muAA=yAAdd-2*muA;
  muAa=yAadd-muA-mua;
  muaa=yaadd-2*mua;
  muBB=yddBB-2*muB;
  muBb=yddBb-muB-mub;
  mubb=yddbb-2*mub;
  muAB=yAdBd-muA-muB;
  muAb=yAddb-muA-mub;
  muaB=ydaBd-mua-muB;
  muab=ydadb-mua-mub;

  muAAB=yAABd-muAA-2*muAB-2*muA-muB;
  muAAb=yAAdb-muAA-2*muAb-2*muA-mub;
  muAaB=yAaBd-muAB-muAa-muaB-muA-mua-muB;
  muAab=yAadb-muAa-muAb-muab-muA-mua-mub;
  muaaB=yaaBd-muaa-2*muaB-2*mua-muB;
  muaab=yaadb-muaa-2*muab-2*mua-mub;
  muABB=yAdBB-muBB-2*muAB-2*muB-muA;
  muaBB=ydaBB-muBB-2*muaB-2*muB-mua;
  muABb=yAdBb-muBb-muAB-muAb-muB-mub-muA;
  muaBb=ydaBb-muBb-muaB-muab-muB-mub-mua;
  muAbb=yAdbb-mubb-2*muAb-2*mub-muA;
  muabb=ydabb-mubb-2*muab-2*mub-mua;

  muAABB=yAABB-2*muAAB-2*muABB-4*muAB-muAA-muBB-2*muA-2*muB;
  muAABb=yAABb-muAAB-muAAb-2*muABb-2*muAB-2*muAb-muAA-muBb-2*muA-muB-mub;
  muAAbb=yAAbb-2*muAAb-2*muAbb-4*muAb-muAA-mubb-2*muA-2*mub;
  muAaBB=yAaBB-2*muAaB-muABB-muaBB-2*muAB-2*muaB-muAa-muBB-muA-mua-2*muB;
  muAaBb=yAaBb-muAaB-muAab-muABb-muaBb-muAB-muAb-muaB-muab-muAa-muBb-muA-mua-muB-mub;
  muAabb=yAabb-2*muAab-muAbb-muabb-2*muAb-2*muab-muAa-mubb-muA-mua-2*mub;
  muaaBB=yaaBB-2*muaaB-2*muaBB-4*muaB-muaa-muBB-2*mua-2*muB;
  muaaBb=yaaBb-muaaB-muaab-2*muaBb-2*muaB-2*muab-muaa-muBb-2*mua-muB-mub;
  muaabb=yaabb-2*muaab-2*muabb-4*muab-muaa-mubb-2*mua-2*mub;

  Rprintf("Additive, Dominant, and Epistatic variances: %f %f %f\n\n",vA,vD,vI);
  vAE=pow(pA*pB,2)       * pow(2*muA+2*muB,2)
     +2*pow(pA,2)*pB*pb  * pow(2*muA+muB+mub,2)
     +pow(pA*pb,2)       * pow(2*muA+2*mub,2)
     +2*pA*pa*pow(pB,2)  * pow(muA+mua+2*muB,2)
     +4*pA*pa*pB*pb      * pow(muA+mua+muB+mub,2)
     +2*pA*pa*pow(pb,2)  * pow(muA+mua+2*mub,2)
     +pow(pa*pB,2)       * pow(2*mua+2*muB,2)
     +2*pow(pa,2)*pB*pb  * pow(2*mua+muB+mub,2)
     +pow(pa*pb,2)       * pow(2*mua+2*mub,2);
  vDE=pow(pA*pB,2)       * pow(muAA+muBB,2)
     +2*pow(pA,2)*pB*pb  * pow(muAA+muBb,2)
     +pow(pA*pb,2)       * pow(muAA+mubb,2)
     +2*pA*pa*pow(pB,2)  * pow(muAa+muBB,2)
     +4*pA*pa*pB*pb      * pow(muAa+muBb,2)
     +2*pA*pa*pow(pb,2)  * pow(muAa+mubb,2)
     +pow(pa*pB,2)       * pow(muaa+muBB,2)
     +2*pow(pa,2)*pB*pb  * pow(muaa+muBb,2)
     +pow(pa*pb,2)       * pow(muaa+mubb,2);
  vAA=pow(pA*pB,2)       * pow(4*muAB,2)
     +2*pow(pA,2)*pB*pb  * pow(2*muAB+2*muAb,2)
     +pow(pA*pb,2)       * pow(4*muAb,2)
     +2*pA*pa*pow(pB,2)  * pow(2*muAB+2*muaB,2)
     +4*pA*pa*pB*pb      * pow(muAB+muAb+muaB+muab,2)
     +2*pA*pa*pow(pb,2)  * pow(2*muAb+2*muab,2)
     +pow(pa*pB,2)       * pow(4*muaB,2)
     +2*pow(pa,2)*pB*pb  * pow(2*muaB+2*muab,2)
     +pow(pa*pb,2)       * pow(4*muab,2);
  vAD=pow(pA*pB,2)       * pow(2*muAAB+2*muABB,2)
     +2*pow(pA,2)*pB*pb  * pow(2*muABb+muAAB+muAAb,2)
     +pow(pA*pb,2)       * pow(2*muAAb+2*muAbb,2)
     +2*pA*pa*pow(pB,2)  * pow(2*muAaB+muABB+muaBB,2)
     +4*pA*pa*pB*pb      * pow(muAaB+muAab+muABb+muaBb,2)
     +2*pA*pa*pow(pb,2)  * pow(2*muAab+muAbb+muabb,2)
     +pow(pa*pB,2)       * pow(2*muaaB+2*muaBB,2)
     +2*pow(pa,2)*pB*pb  * pow(muaaB+muaab+2*muaBb,2)
     +pow(pa*pb,2)       * pow(2*muaab+2*muabb,2);
  vDD=pow(pA*pB,2)       * pow(muAABB,2)
     +2*pow(pA,2)*pB*pb  * pow(muAABb,2)
     +pow(pA*pb,2)       * pow(muAAbb,2)
     +2*pA*pa*pow(pB,2)  * pow(muAaBB,2)
     +4*pA*pa*pB*pb      * pow(muAaBb,2)
     +2*pA*pa*pow(pb,2)  * pow(muAabb,2)
     +pow(pa*pB,2)       * pow(muaaBB,2)
     +2*pow(pa,2)*pB*pb  * pow(muaaBb,2)
     +pow(pa*pb,2)       * pow(muaabb,2);

  mAA=pow(pA*pB,2)       * (4*muAB)
     +2*pow(pA,2)*pB*pb  * (2*muAB+2*muAb)
     +pow(pA*pb,2)       * (4*muAb)
     +2*pA*pa*pow(pB,2)  * (2*muAB+2*muaB)
     +4*pA*pa*pB*pb      * (muAB+muAb+muaB+muab)
     +2*pA*pa*pow(pb,2)  * (2*muAb+2*muab)
     +pow(pa*pB,2)       * (4*muaB)
     +2*pow(pa,2)*pB*pb  * (2*muaB+2*muab)
     +pow(pa*pb,2)       * (4*muab);
  mAD=pow(pA*pB,2)       * (2*muAAB+2*muABB)
     +2*pow(pA,2)*pB*pb  * (2*muABb+muAAB+muAAb)
     +pow(pA*pb,2)       * (2*muAAb+2*muAbb)
     +2*pA*pa*pow(pB,2)  * (2*muAaB+muABB+muaBB)
     +4*pA*pa*pB*pb      * (muAaB+muAab+muABb+muaBb)
     +2*pA*pa*pow(pb,2)  * (2*muAab+muAbb+muabb)
     +pow(pa*pB,2)       * (2*muaaB+2*muaBB)
     +2*pow(pa,2)*pB*pb  * (muaaB+muaab+2*muaBb)
     +pow(pa*pb,2)       * (2*muaab+2*muabb);
  mDD=pow(pA*pB,2)       * (muAABB)
     +2*pow(pA,2)*pB*pb  * (muAABb)
     +pow(pA*pb,2)       * (muAAbb)
     +2*pA*pa*pow(pB,2)  * (muAaBB)
     +4*pA*pa*pB*pb      * (muAaBb)
     +2*pA*pa*pow(pb,2)  * (muAabb)
     +pow(pa*pB,2)       * (muaaBB)
     +2*pow(pa,2)*pB*pb  * (muaaBb)
     +pow(pa*pb,2)       * (muaabb);

  Rprintf("VA, VD, VI = %f %f %f\n",vAE,vDE, vAA+vAD+vDD);
  Rprintf("VI = VAA(%f) + VAD(%f) + VDD(%f) \n\n",vAA,vAD,vDD);
  Rprintf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
  Rprintf("Mean (AA, AD, DD) check: %f %f %f\n\n",mAA,mAD,mDD);

  Rprintf("Additive effects %f %f %f %f\n",muA,mua,muB,mub);
  vAE=2*pA*pow(muA,2)+2*pa*pow(mua,2)+2*pB*pow(muB,2)+2*pb*pow(mub,2);
  Rprintf("sum = %f\n\n",vAE);
  Rprintf("Dominance \n%f %f %f %f %f %f\n",muAA,muAa,muaa,muBB,muBb,mubb);
  vDE=pow(pA*muAA,2)+2*pA*pa*pow(muAa,2)+pow(pa*muaa,2)
     +pow(pB*muBB,2)+2*pB*pb*pow(muBb,2)+pow(pb*mubb,2);
  Rprintf("sum = %f\n\n",vDE);
  Rprintf("AA interactions \n%f %f %f %f\n",muAB,muAb,muaB,muab);
  vAA=4*pA*pB*pow(muAB,2)+4*pA*pb*pow(muAb,2)+4*pa*pB*pow(muaB,2)+4*pa*pb*pow(muab,2);
  Rprintf("sum = %f \n\n", vAA);
  Rprintf("AD interactions \n");
  Rprintf("%f %f %f",  muAAB,muAAb,muAaB);
  Rprintf(" %f %f %f",  muAab,muaaB,muaab);
  Rprintf(" %f %f %f",  muABB,muABb,muAbb);
  Rprintf(" %f %f %f\n",muABB,muaBb,muabb);
  vAD=2*pow(pA,2)*pB*pow(muAAB,2)+2*pow(pA,2)*pb*pow(muAAb,2)
     +4*pA*pa*pB*pow(muAaB,2)+4*pA*pa*pb*pow(muAab,2)
     +2*pow(pa,2)*pB*pow(muaaB,2)+2*pow(pa,2)*pb*pow(muaab,2)
     +2*pA*pow(pB,2)*pow(muABB,2)+4*pA*pB*pb*pow(muABb,2)
     +2*pA*pow(pb,2)*pow(muAbb,2)+2*pa*pow(pB,2)*pow(muaBB,2)
     +4*pa*pB*pb*pow(muaBb,2)+2*pa*pow(pb,2)*pow(muabb,2);
  Rprintf("sum = %f \n\n", vAD);
  Rprintf("DD interactions\n");
  Rprintf("%f %f %f",  muAABB,muAABb,muAAbb);
  Rprintf(" %f %f %f",  muAaBB,muAaBb,muAabb);
  Rprintf(" %f %f %f\n",muaaBB,muaaBb,muaabb);
  vDD=pow(pA*pB*muAABB,2)+2*pow(pA,2)*pB*pb*pow(muAABb,2)+pow(pA*pb*muAAbb,2)
     +2*pA*pa*pow(pB*muAaBB,2)+4*pA*pa*pB*pb*pow(muAaBb,2)+2*pA*pa*pow(pb*muAabb,2)
     +pow(pa*pB*muaaBB,2)+2*pow(pa,2)*pB*pb*pow(muaaBb,2)+pow(pa*pb*muaabb,2);
  Rprintf("sum = %f\n\n",vDD);
  Rprintf("VA, VD, VI = %f %f %f\n",vAE,vDE, vAA+vAD+vDD);
  Rprintf("VI = VAA(%f) + VAD(%f) + VDD(%f) \n\n",vAA,vAD,vDD);
}
