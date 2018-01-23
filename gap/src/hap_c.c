#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include "hap.h"

static int n_subjects, n_loci=0, n_phase, *alleles, all_snps=0;
long n_warn = 0;
double **af,**pq;
static FILE *outfile;
char of1name[MAX_FILENAME_LEN],of2name[MAX_FILENAME_LEN],
     tempname[MAX_FILENAME_LEN];

void hap_c(
int *nobs,
char **idstr,
char **gdata,
int *nloci,
int *loci,
char **names,
double *mb,
double *pr,
double *po,
double *to,
double *th,
int *maxitt,
int *n,
int *sst,
int *rst,
int *rp,
int *ro,
int *rv,
double *sd,
int *mm,
int *mi,
int *mc,
double *ds,
double *de,
int *q,
double *l1,
int *niter,
int *converged,char **hapfile, char **assignfile
)
{
  CODE *coding;
  HAP **so_list, **ho_list, **unique, *h1, *h2, **h;
  char cc[3];
  double memmax, logl, logl_s, lastl=0, best_logl=0, df, df_step=0, *p_unique;
  int n_subject, *i_subject, mpl, mps, mph, *order = (int *) 0;
  int  res, j, k, iter, carry_on, nlen=0;
  long n_hap, hap_n, max_haps, i;
  short a1, a2, max_alleles;
  time_t now;

  double memory = *mb; /* Maximum memory to be consumed by task, in Mbyte*/
  double min_posterior = *po; /* Minimum posterior probability */
  double min_prior = *pr;  /* Minimum prior probability */
  double tol = *to;        /* Log likelihood tolerance */
  int maxit = *maxitt;     /* Maximum EM iterations */
  int num = *n;            /* Force numeric output */
  int ss = *sst;           /* Tab delimited ouput for spreadsheet */
  int rs = *rst;           /* Random starting point */
  long seed = *sd;         /* Random number seed */
  int rano = *ro;          /* Add loci in random order */
  int revo = *rv;          /* Add loci in reverse order */
  int mmle = *mm;          /* Multiple maximizations at last step */
  int mimp = *mi;          /* Multiple random imputations */
  int mcmc = *mc;          /* Number of MCMC steps between samples */
  double dfs = *ds;        /* Starting df parameter for Dirichlet (*N) */
  double dfe = *de;        /* Finishing df parameter for Dirichlet (*N) */
  int ranp = 0;            /* Restart from random prior probabilities */
  int quiet = *q;          /* Suppress screen progress report */
  double of_max = 1.0;     /* Posterior threshold for listing */
  int maxs = 10000;        /* Maximum subjects (doesn't matter yet) */
  double head_room = 0.5;  /* Dynamic memory headroom in Mb */

  GetRNGstate();
  time(&now);
  Rprintf("\nHAP, Version %s", VERSION);
  if (now != -1) {
    Rprintf(" run on  %s", asctime(localtime(& now)));
    Rprintf("hap warnings: %s", asctime(localtime(& now)));
  }
  else Rprintf(": time not available\n");
  memmax = memavail(1024)/1.0e6;
  if (memory<head_room || memory> memmax) {
    Rprintf("Dynamic memory available %.1f Mbyte\n", memmax);
    memory = memmax;
  }
  strcpy(of1name,*hapfile);
  strcpy(of2name,*assignfile);
  n_loci=*nloci;
  n_subjects=*nobs;
  all_snps=1;
  alleles = calloc(n_loci, sizeof(int));
  max_alleles=-1;
  for (j=0; j<n_loci; j++) {
      alleles[j]=loci[j];
      if(alleles[j]>2) all_snps=0;
      if(alleles[j]>max_alleles) max_alleles=alleles[j];
  }
  i_subject = calloc(n_loci, sizeof(int));
  if (!i_subject) goto no_room;
  for (i=0;i<n_loci;i++) i_subject[i]=0;
  af=allocateU(alleles);
  pq=allocateU(alleles);
  mpl = sizeof(CODE);
  if (names) mpl += sizeof(char *) + nlen;
  if (revo || rano) mpl += sizeof(int);
  mps = 0;
  mph = 3*sizeof(HAP *) + sizeof(HAP);
  if (mmle > 0 || mimp > 0) mph += sizeof(double);
  max_haps = ((memory-head_room)*1.e6 - n_loci*mpl - maxs*mps)/mph;
  Rprintf("Maximum memory usage            : %10.1f Mb\n", memory);
  Rprintf("Maximum haplotype instances     : %10ld\n", max_haps);
  if (rano || rs || mimp>0) {
    if (!seed) seed = (long) now;
    Rprintf("Seed for random numbers         : %10ld\n", seed);
/*
    SEED (seed);
*/
  }
  if (rano) {
    order = (int *) calloc(n_loci, sizeof(int));
    if (order) {
      ranord(n_loci, order);
      Rprintf("Loci added in random order\n");
    }
    else goto no_room;
  }
  if (revo) {
    order = (int *) calloc(n_loci, sizeof(int));
    if (order) for (i=0; i<n_loci; i++) order[i] = n_loci-i-1;
    else goto no_room;
  }
  if (rs) rs = 1;
  so_list = calloc(max_haps, sizeof(HAP*));
  if (!so_list) goto no_room;
  ho_list = calloc(max_haps, sizeof(HAP*));
  if (!ho_list) goto no_room;
  coding = (CODE *) calloc(n_loci, sizeof(CODE));
  if (!coding) goto no_room;
  else {
    for (j=0; j<n_loci; j++) {
      coding[j].anum = 0;
      strcpy(coding[j].one,"");
      strcpy(coding[j].two,"");
    }
  }
  n_subject = n_hap = 0;
  h = so_list;
  for(j=0;j<n_subjects;j++) {
    res=gt_read(j,idstr,gdata,order,coding,&h1,&h2);
    switch (res) {
    case 1:
      n_subject++;
      for(i=0;i<n_loci;i++) {
        a1=(*h1).loci[i];
        a2=(*h2).loci[i];
        if(a1>0) {
          af[i][a1-1]++;
          i_subject[i]++;
        }
        if(a2>0) {
          af[i][a2-1]++;
          i_subject[i]++;
        }
      }
      *h++ = h1;
      *h++ = h2;
      n_hap += 2;
      break;
    case 2:
      REprintf("Skipping this subject(%s)\n", (*h1).id);
      break;
    case 3:
      REprintf("Exiting while reading subject %d\n", n_subject);
      PutRNGstate();
      return;
    }
  }
  Rprintf("Thresholds for trimming\n");
  Rprintf("\tPrior probability       : %10.2f\n", min_prior);
  Rprintf("\tPosterior probability   : %10.2f\n", min_posterior);
  Rprintf("EM convergence criteria\n");
  Rprintf("\tMaximum iterations      : %10d\n", maxit);
  Rprintf("\tTolerated LLH change    : %10.2f\n", tol);
  if (mmle > 0) {
    Rprintf("Repetions of final EM iteration : %10d\n", mmle);
  }
  if (mimp > 0) {
    Rprintf("Multiple imputed datasets       : %10d\n", mimp);
    Rprintf("\tMCMC steps              : %10d\n", mcmc);
    dfs *= (2*n_subject);
    dfe *= (2*n_subject);
    df_step = (dfs - dfe)/(double) mcmc;
    Rprintf("\tStarting Dirichlet df   : %10.2f\n", dfs);
    Rprintf("\tFinishing Dirichlet df  : %10.2f\n", dfe);
  }
  if (of1name[0]) {
    Rprintf("Haplotype frequencies output to : %10s", of1name);
    if (mimp>0) Rprintf("(.*)\n");
    else Rprintf("\n");
  }
  if (of2name[0]) {
    Rprintf("Haplotype assignments output to : %10s", of2name);
    if (mimp>0) Rprintf("(.*)\n");
    else Rprintf("\n");
  }
  Rprintf("\nAllele frequencies:\n");
  Rprintf("\nMarker\\Allele number\n\n");
  Rprintf("    ");
  for(i=0;i<max_alleles;i++) Rprintf("%9ld",i+1);
  Rprintf("\n\n");
  for(i=0;i<n_loci;i++) {
    Rprintf("%3ld ",i+1);
    for(j=0;j<alleles[i];j++) {
      af[i][j]/=(double)i_subject[i];
      Rprintf(" %f",af[i][j]);
    }
    Rprintf("\n");
  }
  Rprintf("\nSubjects (%% of %d) with full genotype at each locus:\n\n",n_subject);
  for(i=0;i<n_loci;i++)
    Rprintf("%3ld %5d %6.2f\n",i+1,i_subject[i]/2,50.0*i_subject[i]/n_subject);
  if (!quiet) {
    Rprintf("\nProgress of maximum likelihood estimation:\n");
    Rprintf("\nLoci Phased\tht instances\tIteration\tLog Likelihood");
    Rprintf("\n-----------\t------------\t---------\t--------------\n");
  }
  n_phase = 0;
  while (n_phase < n_loci) {
    n_hap = hap_expand(n_hap, max_haps, so_list, rs);
    if (!n_hap) goto no_room;
    if (!quiet) {
      Rprintf("\r%11d\t%12ld\t%9s", n_phase, n_hap, "(sorting)");
      R_FlushConsole();
    }
    for (i=0; i<n_hap; i++) ho_list[i] = so_list[i];
    qsort(ho_list, n_hap, sizeof(HAP *), cmp_hap);
#if 0
    /* debugging information JH Zhao */
    hap_list(stdout,n_hap,coding,so_list);
    Rprintf("\nafter sorting\n");
    hap_list(stdout,n_hap,coding,ho_list);
#endif
    carry_on =  1;
    iter = 0;
    *converged=0;
    while (carry_on) {
      iter++;
      hap_prior(n_hap, ho_list, min_prior);
      hap_n = hap_posterior(n_hap, so_list, min_posterior, &logl, 0);
      carry_on = (iter==1) || ((logl - lastl)>tol);
      if (carry_on && (iter==maxit)) {
        carry_on = 0;
        REprintf("No convergence in %d iterations", iter);
        REprintf("\t(at %d-locus step)\n",n_phase);
        n_warn++;
      }
      else *converged=1;
      *niter=iter;
      lastl = logl;
      if (!quiet) {
        Rprintf("\r%11d\t%12ld\tEM%7d\t%14.3f", n_phase, n_hap, iter, logl);
        R_FlushConsole();
      }
    }
    n_hap = hap_posterior(n_hap, so_list, min_posterior, &logl, 1);
  }
  logl_s=lastl;
  if (!quiet) {
    Rprintf("\r\t\t%12ld", n_hap);
    R_FlushConsole();
  }
  for (i=0; i<n_hap; i++) ho_list[i] = so_list[i];
  qsort(ho_list, n_hap, sizeof(HAP *), cmp_hap);
  hap_prior(n_hap, ho_list, min_prior);
  hap_n = n_unique_haps(n_hap, ho_list);
  unique = (HAP **) calloc(hap_n, sizeof(HAP *));
  if (!unique) goto no_room;
  unique_haps(n_hap, ho_list, unique);
  p_unique=(double*)0;
  if (mmle>0 || mimp>0) {
    p_unique = (double *) calloc(hap_n, sizeof(double));
    if (!p_unique) goto no_room;
    for (i=0; i<hap_n; i++) p_unique[i] = unique[i]->prior;
  }
  if (mmle > 0) {
    best_logl = logl;
    for (j=0; j<mmle; j++) {
      carry_on =  1;
      iter = 0;
      *converged=0;
      if (ranp) {
        hap_prior_restart(n_hap, ho_list);
        hap_posterior(n_hap, so_list, min_posterior, &logl, 0);
      }
      else hap_posterior_restart(n_hap, so_list);
      while (carry_on) {
        iter++;
        hap_prior(n_hap, ho_list, min_prior);
        hap_posterior(n_hap, so_list, min_posterior, &logl, 0);
        carry_on = (iter==1) || ((logl - lastl)>tol);
        if (carry_on && (iter==maxit)) {
          carry_on = 0;
          REprintf("No convergence in %d iterations", iter);
          REprintf("\t(at %d-locus step)\n",n_phase);
          n_warn++;
        }
        else *converged=1;
        *niter=iter;
        lastl = logl;
        if (!quiet) {
          Rprintf("\r%6d/%04d\t%12ld\tEM%7d\t%14.3f", n_loci, j+1,
               n_hap, iter, logl);
          R_FlushConsole();
        }
      }
      if (logl > best_logl) {
        best_logl = logl;
        for (i=0; i<hap_n; i++) p_unique[i] = unique[i]->prior;
        if (!quiet) {
          Rprintf("\r\t\t\t\t\t\t\t\t(%.3f)", best_logl);
          R_FlushConsole();
        }
      }
    }
    hap_prior_restore(n_hap, ho_list, p_unique);
    if (quiet) Rprintf("\n\nBest solution has log likelihood %.3f", best_logl);
  }
  qsort(unique, hap_n, sizeof(HAP *), more_probable);
  if (!quiet) Rprintf("\n");
  Rprintf("\n");
  if (of1name[0]) {
    outfile = fopen(of1name, "w");
    if (!outfile) goto open_error;
    if (!ss) {
       fprintf(outfile, "hap listing: %s", asctime(localtime(& now)));
       fprintf(outfile,"\n");
       fprintf(outfile,"Log likelihood = %.3f\n",logl_s);
       if(mmle>0) fprintf(outfile,"The best log likelihood = %.3f\n",best_logl);
       fprintf(outfile,"\n");
    }
    j = hap_write(outfile, n_loci, names, coding, order, hap_n, unique, 0, 0.0, num, ss);
    fclose(outfile);
    Rprintf("%d haplotypes written to output file %s\n", j, of1name);
  }
  *l1=logl_s;
  if(mmle>0) *l1=best_logl;
  if (of2name[0]) {
    outfile = fopen(of2name, "w");
    if (!outfile) goto open_error;
    j = hap_write(outfile, n_loci, names, coding, order, n_hap, so_list, 1, of_max, num, ss);
    fclose(outfile);
    Rprintf("%d possible assignments to %d subjects written to output file %s\n",
           j, n_subject, of2name);
  }
  if (mimp > 0 && !quiet) {
    Rprintf("\nProgress of multiple imputation:\n");
    Rprintf("\nImputation\tIP step\t    Prior df");
    Rprintf("\n----------\t-------\t    --------\n");
  }
  for (j=1; j<=mimp; j++) {
    hap_prior_restore(n_hap, ho_list, p_unique);
    df = dfs;
    for (iter=0; iter<mcmc; iter++) {
      df -= df_step;
      if (!quiet) {
        Rprintf("\r%10d\t%7d\t%12.2f", j, iter+1, df);
        R_FlushConsole();
      }
      sample_posterior(n_hap, so_list); /* I step */
      sample_prior(n_hap, ho_list, df); /* P step */
    }
    if (of1name[0]) {
      sprintf(tempname,"%s.%03d", of1name, j);
      outfile = fopen(tempname, "w");
      hap_posterior(n_hap, so_list, min_posterior, &logl, 0);
      if(!ss) {
        fprintf(outfile,"\n");
        fprintf(outfile,"Log likelihood = %.3f\n",logl);
        fprintf(outfile,"\n");
      }
      hap_write(outfile, n_loci, names, coding, order, hap_n, unique, 0, 0.0, num, ss);
      fclose(outfile);
    }
    if (of2name[0]) {
      sprintf(tempname,"%s.%03d", of2name, j);
      outfile = fopen(tempname, "w");
      hap_write(outfile, n_loci, names, coding, order, n_hap, so_list, 1, 0.0, num, ss);
      fclose(outfile);
    }
  }
  if (!quiet) Rprintf("\n");
  if (mimp>0 && of1name[0]) {
    Rprintf("\nSamples from posterior distribution of haplotype frequencies ");
    Rprintf("written to \n\tfiles %s.001 ... %s.%03d\n", of1name, of1name,
           mimp);
  }
  if (mimp>0 && of2name[0]) {
    Rprintf("Multiply imputed datasets written to files %s.001 ... %s.%03d\n",
           of2name, of2name, mimp);
  }

  if (num) {
    Rprintf("\nNumerical recoding of alleles was forced:\n");
    for (j=k=0; j<n_loci; j++) {
      if (coding[j].anum == 2) {
        k++;
        Rprintf("%s: \t", names[j]);
        strcpy(cc, coding[j].one);
        Rprintf("1=%s", cc);
        strcpy(cc,coding[j].two);
        Rprintf(", 2=%s", cc);
        Rprintf("\n");
      }
    }
    if (!k) Rprintf("\t... but no recoding was necessary\n");
  }
  if (n_warn) Rprintf("\nNote: %ld messages\n", n_warn);
  else Rprintf("\n");
  free(i_subject);
  free(alleles);
  freeU(af);
  freeU(pq);
  PutRNGstate();
  return;

 no_room:
  REprintf("Insufficient memory\n");
  PutRNGstate();
  return;

 open_error:
  REprintf("Error opening output file\n");
  PutRNGstate();
  return;

}

HAP* new_hap(char *id, double prior, double posterior, char *loci){
  HAP *result;
  short *loc;
  int i;
  result = (HAP *) malloc(sizeof(HAP));
  if (result) {
    result->id = id;
    result->prior = prior;
    result->posterior = posterior;
    loc = malloc(n_loci*sizeof(short));
    if (loc) {
      result->loci = loc;
      for (i=0; i<n_loci; i++)
        loc[i] = loci? loci[i] : 0;
    }
    else {
      free(result);
      result = 0;
    }
  }
  return result;
}

HAP* cpy_hap(HAP *old) {
  HAP *result;
  short *loc;
  int i;
  result = (HAP *) malloc(sizeof(HAP));
  if (result) {
    result->id = old->id;
    result->prior = old->prior;
    result->posterior = old->posterior;
    loc = malloc(n_loci*sizeof(short));
    if (loc) {
      result->loci = loc;
      for (i=0; i<n_loci; i++)
        loc[i] = old->loci[i];
    }
    else {
      free(result);
      result = 0;
    }
  }
  return result;
}

void kill_hap(HAP* old){
  free(old->loci);
  free(old);
}

int cmp_hap(const void *to_one, const void *to_two){
  int i;
  short *loc1, *loc2;
  short a1, a2;
  HAP *one, *two;
  one = * (HAP **) to_one;
  two = * (HAP **) to_two;
  loc1 = one->loci;
  loc2 = two->loci;
  for (i=0; i<n_phase; i++) {
    a1 = loc1[i];
    a2 = loc2[i];
    if (a1<a2) return -1;
    if (a1>a2) return +1;
  }
  return 0;
}

int more_probable(const void *to_one, const void *to_two) {
  HAP *one, *two;
  double p1, p2;
  one = * (HAP **) to_one;
  two = * (HAP **) to_two;
  p1 = one->prior;
  p2 = two->prior;
  if (p1<p2) return 1;
  if (p1>p2) return -1;
  return 0;
}

/* List haplotypes */

int allele_code(int i, CODE code) {
#if 0
  if (code.anum==1)
    return i; /* '0'+i, Numeric coding */
  else if (code.anum==2)
    return i? (i==1? code.one : code.two) : ' ';
  else
    return '?';
#endif
  return i;
}

void hap_list(FILE *out, long n_hap, CODE *code, HAP **list) {
  long i;
  int j, k;
  HAP *h;

  for (i=0; i<n_hap; i++) {
    h = list[i];
    fprintf(out, "%12s %12.5f %12.5f  ",
            h->id, h->prior, h->posterior);
    for (j=0; j<n_loci; j++) {
      if(code[j].anum==1)
        fprintf(out, "%2d", allele_code(h->loci[j], code[j]));
      else {
        k=h->loci[j];
        fprintf(out,  "%2s", k? (k==1? code[j].one : code[j].two) : " ");
      }
    }
    fprintf(out, "\n");
  }
}

/* Expand unphased data by possible haplotype */

long hap_expand(long n_hap, long max_haps, HAP** so_list, int random_start) {
  double p1, p2, pq11, pq12, pq21, pq22;
  double f, den, ppqq[99][99];
  double u1, u2, ut;
  long i, new_hap;
  short a1, a2, het;
  HAP **h, **hnew;
  HAP *h1, *h2, *nh1, *nh2;
  long j, k;

  /*  First pass */

  f = den = 0.0;
  for (h=so_list, new_hap=i=0; i<n_hap; i+=2) {
    h1 = *h++;
    h2 = *h++;
    a1 = (int) h1->loci[n_phase];
    a2 = (int) h2->loci[n_phase];
    het  = (cmp_hap(&h1, &h2)!=0);
    if (!a1)
      j = (het? 2*alleles[n_phase]:alleles[n_phase]+1)*alleles[n_phase];
    else {
      den += 2;
      if(a1>0) pq[n_phase][alleles[n_phase]-a1]++;
      if(a2>0) pq[n_phase][alleles[n_phase]-a2]++;
/*
      if(a1>0) pq[n_phase][a1-1]++;
      if(a2>0) pq[n_phase][a2-1]++;
*/
      het = (het && (a1!=a2));
      j = (het? 4: 2);
    }
    new_hap += j;
  }
  for(i=0;i<alleles[n_phase];i++) if(den>0.0) pq[n_phase][i]/=den;
    else pq[n_phase][i]=1.0/alleles[n_phase];

  /* Second pass */

  if (new_hap > max_haps)
    return 0; /* Exceed memory use */

  for (h = so_list+n_hap-1, hnew= so_list+new_hap-1, i=n_hap-1; i>=0; i-= 2) {
    h2 = *h--;
    h1 = *h--;
    a1 = (int) h1->loci[n_phase];
    a2 = (int) h2->loci[n_phase];
    het  = (cmp_hap(&h1, &h2)!=0);
    if (!a1) {
      /* new += (het? 4: 3); */
      p2 = h2->posterior;
      p1 = h1->posterior;
      if (random_start) {
        ut=0.0;
        for(j=0;j<alleles[n_phase];j++) {
           for(k=0;k<alleles[n_phase];k++) {
             ppqq[j][k] = URAN();
             ut += ppqq[j][k];
           }
        }
        for(j=0;j<alleles[n_phase];j++)
           for(k=0;k<alleles[n_phase];k++)
             ppqq[j][k] /= ut;
      }
      h1->loci[n_phase]=h2->loci[n_phase]=1;
      if(random_start) f=ppqq[0][0];
      else f=pq[n_phase][0]*pq[n_phase][0];
      h2->posterior=p2*f;
      h1->posterior=p1*f;
      *hnew-- = h2;
      *hnew-- = h1;
      for(j=2;j<=alleles[n_phase];j++) {
        if(random_start) f=ppqq[j-1][j-1];
        else f=pq[n_phase][j-1]*pq[n_phase][j-1];
        nh2 = cpy_hap(h2);
        nh1 = cpy_hap(h1);
        if (!nh1 || !nh2)
          return 0;
        nh1->loci[n_phase]=nh2->loci[n_phase]=j;
        nh2->posterior=p2*f;
        nh1->posterior=p1*f;
        *hnew-- = nh2;
        *hnew-- = nh1;
      }
      for(j=2;j<=alleles[n_phase];j++) {
        for(k=1;k<j;k++) {
           if(random_start) {
             pq11=ppqq[j-1][j-1];
             pq12=ppqq[j-1][k-1];
             pq21=ppqq[k-1][j-1];
             pq22=ppqq[k-1][k-1];
           }
           else {
             pq11=pq[n_phase][j-1]*pq[n_phase][j-1];
             pq12=pq[n_phase][j-1]*pq[n_phase][k-1];
             pq21=pq[n_phase][k-1]*pq[n_phase][j-1];
             pq22=pq[n_phase][k-1]*pq[n_phase][k-1];
           }
           nh2 = cpy_hap(h2);
           nh1 = cpy_hap(h1);
           if (!nh1 || !nh2)
             return 0;
           nh2->loci[n_phase]=k;
           nh1->loci[n_phase]=j;
           if (het) {
             nh2->posterior = p2*pq12;
             nh1->posterior = p1*pq11;
             *hnew-- = nh2;
             *hnew-- = nh1;
             nh2 = cpy_hap(nh2);
             nh1 = cpy_hap(nh1);
             if (!nh1 || !nh2)
               return 0;
             nh2->loci[n_phase]=j;
             nh1->loci[n_phase]=k;
             nh2->posterior = p2*pq22;
             nh1->posterior = p1*pq21;
             *hnew-- = nh2;
             *hnew-- = nh1;
           }
           else {
             nh2->posterior = p2*(pq12 + pq22);
             nh1->posterior = p1*(pq11 + pq21);
             *hnew-- = nh2;
             *hnew-- = nh1;
           }
        }
      }
    }
    else {
      /* new += (het? 2: 1); */
      het = (het && (a1!=a2));
      if (het) {
        nh2 = cpy_hap(h2);
        nh1 = cpy_hap(h1);
        if (!nh1 || !nh2)
          return 0;
        if (random_start) {
          u1 = URAN ();
          u2 = URAN ();
          ut = u1+u2;
          h1->posterior *= u1/ut;
          nh1->posterior *= u2/ut;
          u1 = URAN ();
          u2 = URAN ();
          ut = u1+u2;
          h2->posterior *= u1/ut;
          nh2->posterior *= u2/ut;
        }
        else {
          h2->posterior /= 2.0;
          h1->posterior /= 2.0;
          nh2->posterior /= 2.0;
          nh1->posterior /= 2.0;
        }
        nh2->loci[n_phase] = h1->loci[n_phase];
        nh1->loci[n_phase] = h2->loci[n_phase];
        *hnew-- = nh2;
        *hnew-- = nh1;
      }
      *hnew-- = h2;
      *hnew-- = h1;
    }
  }
  n_phase++;
  return new_hap;
}

/* Functions to work on list in haplotype order */

void hap_prior(long n_hap, HAP** ho_list, double min_prior) {
  double total, subtotal;
  HAP **hs, **he, **hn, **h;
  /* Calculate total posterior probability */
  hs = ho_list;
  he = ho_list + n_hap;
  for (total=0.0, h=hs; h<he; h++)
    total += (*h)->posterior;
  /* Sum posterior  probability for each haplotype and hence calculate prior */
  while (hs < he) {
    h = hs;
    subtotal = 0.0;
    do {
      subtotal += (*h)->posterior;
      h++;
    } while ((h<he) && (cmp_hap(hs, h)==0));
    hn = h;
    subtotal /=  total;
    if (subtotal < min_prior)
      subtotal = 0.0;
    for (h=hs; h<hn; h++)
      (*h)->prior = subtotal;
    hs = hn;
  }
}

void hap_prior_restart(long n_hap, HAP** ho_list) {
  double total, u;
  HAP **hs, **he, **h;
  hs = ho_list;
  he = ho_list + n_hap;
  total = 0.0;
  total = u = URAN ();
  for (h=hs; h<he; h++) {
    (*h)->prior = u;
    if ((h<he) && (cmp_hap(h, h+1)!=0)) {
      u = URAN ();
      total += u;
    }
  }
  for (h=hs; h<he; h++) {
    (*h)->prior /= total;
  }
}

void sample_prior(long n_hap, HAP** ho_list, double prior_df) {
  double total, gv, post_df;
  int count;
  HAP **hs, **he, **hn, **h;
  hs = ho_list;
  he = ho_list + n_hap;
  /* Count occurrences of each haplotype */
  total = 0.0;
  while (hs < he) {
    h = hs;
    count = 0;
    do {
      if ((*h)->posterior)
        count ++;
      h++;
    } while ((h<he) && (cmp_hap(hs, h)==0));
    hn = h;
    post_df = prior_df + (double) count;
    gv = post_df > 0 ? rangamma(post_df):  0.0;
    total += gv;
    for (h=hs; h<hn; h++)
      (*h)->prior = gv;
    hs = hn;
  }
  for (h=ho_list; h<he; h++) {
    (*h)->prior /= total;
  }
}

long n_unique_haps(long n_hap, HAP **ho_list) {
  HAP **hs, **he, **h;
  long res = 0;
  hs = ho_list;
  he = ho_list + n_hap;
  while (hs < he) {
    h = hs;
    do {
      h++;
    } while ((h<he) && (cmp_hap(hs, h)==0));
    res++;
    hs = h;
  }
  return res;
}


void unique_haps(long n_hap, HAP **ho_list, HAP **unique) {
  HAP **hs, **he, **h;
  hs = ho_list;
  he = ho_list + n_hap;
  while (hs < he) {
    h = hs;
    do {
      h++;
    } while ((h<he) && (cmp_hap(hs, h)==0));
    *unique++ = *hs;
    hs = h;
  }
  return;
}

void hap_prior_restore(long n_hap, HAP **ho_list, double *p_unique) {
  HAP **hs, **he, **h;
  long i;
  hs = ho_list;
  he = ho_list + n_hap;
  i = 0;
  while (hs < he) {
    h = hs;
    do {
      (*(h++))->prior = p_unique[i];
    } while ((h<he) && (cmp_hap(hs, h)==0));
    i ++;
    hs = h;
  }
  return;
}

/* Functions to operate on list in subject order */

long hap_posterior(long n_hap, HAP **so_list, double min_posterior,
                   double *llh, int if_pack) {
  HAP **hs, **he, **hn, **h, **h2,  **n;
  long new_hap, small;
  char *id;
  double subtotal, gp, wllh;
  int any;
  hs = so_list;
  he = so_list + n_hap;
  small = 0;
  wllh = 0.0;
  while (hs < he) {
    h = hs;
    subtotal = 0.0;
    any = 0;
    do {
      id = (*h)->id;
      h2 = h+1;
      gp = (*h)->prior * (*h2)->prior;
      if (cmp_hap(h, h2)!=0)
        gp *= 2;
      subtotal += gp;
      (*h)->posterior = (*h2)->posterior = gp;
      h = h2+1;
    } while ((h<he) && ((*h)->id)==id);
    hn = h;
    /* Scale posterior, detecting if any are too small */
    if (subtotal>0.0) {
      for (h=hs; h<hn; h++) {
        (*h)->posterior /= subtotal;
        if ((*h)->posterior < min_posterior)
          small=1;
        else
          any = 1;
      }
      if (!if_pack || any)
        wllh += log(subtotal);
    }
    else {
      small = 1;
    }
    if (if_pack && !any)  {
      REprintf("Subject %s dropped from data ", id);
      REprintf("\t(at %d-locus step)\n", n_phase);
      n_warn++;
    }
    hs = hn;
  }
  if (if_pack && small) {
    for (new_hap=0, n=h=so_list; h<he; h++) {
      if ((*h)->posterior < min_posterior) {
        kill_hap(*h);
      }
      else {
        *(n++) = (*h);
        new_hap++;
      }
    }
    return hap_posterior(new_hap, so_list, min_posterior, llh, if_pack);
  }
  else {
    *llh = wllh;
    return n_hap;
  }
}

void hap_posterior_restart(long n_hap, HAP** so_list) {
  HAP **h, **h2, **hs, **he, **hn;
  char *id;
  double u, subtotal;
  hs = so_list;
  he = so_list + n_hap;
  while (hs < he) {
    h = hs;
    subtotal = 0.0;
    do {
      id = (*h)->id;
      h2 = h+1;
      u = URAN ();
      subtotal += u;
      (*h)->posterior = (*h2)->posterior = u;
      h = h2+1;
    } while ((h<he) && ((*h)->id)==id);
    hn = h;
    /* Scale posterior */
    for (h=hs; h<hn; h++)
      (*h)->posterior /= subtotal;
    hs = hn;
  }
  return;
}

void sample_posterior(long n_hap, HAP **so_list) {
  HAP **hs, **he, **hn, **h, **h2;
  long small;
  char *id;
  double subtotal, gp;
  int any;
  hs = so_list;
  he = so_list + n_hap;
  small = 0;
  while (hs < he) {
    h = hs;
    subtotal = 0.0;
    any = 0;
    do {
      id = (*h)->id;
      h2 = h+1;
      gp = (*h)->prior * (*h2)->prior;
      if (cmp_hap(h, h2)!=0)
        gp *= 2;
      subtotal += gp;
      (*h)->posterior = (*h2)->posterior = subtotal;
      h = h2+1;
    } while ((h<he) && ((*h)->id)==id);
    hn = h;
    /* Sample posterior */
    subtotal *= URAN ();
    for (h=hs; (*h)->posterior<subtotal; h++)
      (*h)->posterior = 0.0;
    (*(h++))->posterior = 1.0;
    (*(h++))->posterior = 1.0;
    while(h<hn)
      (*(h++))->posterior = 0.0;
    hs = hn;
  }
}

long check_hap(long n_hap, HAP **list) {
  long i, res;
  int j, a;
  res = 0;
  for (i=0; i<n_hap; i++) {
    for (j=0; j<n_phase; j++) {
      a = (int) list[i]->loci[j];
      if (a <0 || a>99) {
        res ++;
        break;
      }
    }
  }
  return res;
}

int encode(char a[3], CODE *code) {
  int an;
  char c1[3], c2[3];
  an = (*code).anum;
  if (!an) {
    if (strcmp(a,"0")<0 || strcmp(a,"9")>0)
      (*code).anum = an = 2; /* alpha coding */
    else
      (*code).anum = an = 1; /* numeric coding */
  }

  /* Do the encoding */

  if (an==1) {
    if (strcmp(a,"0")<0 || strcmp(a,"9")>0)
      return -1; /* Error */
    else
      return atoi(a);
  }
  else if (an==2) {
    if (strcmp(a,"A") && strcmp(a,"C") && strcmp(a,"G") && strcmp(a,"T")) {
      return 0;
    }
    strcpy(c1,(*code).one);
    if (strlen(c1)>0) {
      if (!strcmp(a,c1))
        return 1;
      else {
        strcpy(c2, (*code).two);
        if (strlen(c2)>0) {
          if (!strcmp(a,c2))
            return 2;
          else
            return -1;
        }
        else {
          strcpy((*code).two, a);
          return 2;
        }
      }
    }
    else {
      strcpy((*code).one, a);
      return 1;
    }
  }
  else
    return -1;
}

int gt_read(int i, char **idstr, char **gdata, int *order, CODE *code, HAP **one, HAP **two)
{
  char *id;
  char a1[3], a2[3];
  int j, err, i1, i2;

  id = malloc(1+strlen(idstr[i]));
  if (!id) goto no_room;
  strcpy(id, idstr[i]);
  if(!((*one) = new_hap(id, 0.0, 1.0, 0))) goto no_room;
  if(!((*two) = new_hap(id, 0.0, 1.0, 0))) goto no_room;
  for (err=j=0; j<n_loci; j++) {
      strcpy(a1,gdata[i*n_loci*2+2*j]);
      strcpy(a2,gdata[i*n_loci*2+2*j+1]);
/*
      i1=atoi(a1);
      i2=atoi(a2);
*/
      i1 = encode(a1, code+j);
      i2 = encode(a2, code+j);
      if (i1<0 || i2<0 || (i1&&!i2) || (i2&&!i1)) {
        REprintf("Data error on locus %d: %2s/%2s\n", j+1, a1, a2);
        err = 1;
      }
      if (order) {
        (*one)->loci[order[j]] = i1;
        (*two)->loci[order[j]] = i2;
      }
      else {
        (*one)->loci[j] = i1;
        (*two)->loci[j] = i2;
      }
    }
  if (err) return 2;
  else return 1;
no_room:
  return 3;
}

/* Create vector of 1..n in random order */

void ranord(int n, int *order) {
  int i, j, w;
  double u;
  for (i=0; i<n; i++) {
    u =  URAN ();
    w = (int) (i*u + 0.5); /* where to put next number */
    for (j=i; j>w; j--) {
      order[j] = order[j-1];
    }
    order[w] = i;
  }
}

/* Determine how much dynamic memory is available */

int talloc(long request) {
  char *c;
  c = malloc(request);
  if (c) {
    free(c);
    return 1;
  }
  else {
    return 0;
  }
}

long memavail(int atom) {
  long step, res;
  res = 0;
  step = atom;
  while (talloc(res+step)) {
    res += step;
    step *= 2;
  }
  while (step>atom) {
    step /= 2;
    if (talloc(res+step)) {
      res += step;
    }
  }
  return res;
}


int hap_write(FILE *outfile, int n_loci, char **names, CODE *coding,
               int *order, long n_hap, HAP **list, int so, double of_max,
               int numeric, int tabdel) {
  int j, k, l, nlen, len, chr, nwr;
  long i, out;
  char  c[7], rep;
  char *id;
  double cum, max_post=0;
  HAP **h, **hs, **he, **hn;

  nwr = 0;
  if (of_max > 1.0)
    of_max = 1.0;

  /* Header line(s) */

  if (tabdel) {
    if (so) {
      fprintf(outfile, "id\tchr\t");
    }
    for (j=0; j<n_loci; j++)
      if (names)
        fprintf(outfile, "%s\t", names[j]);
      else
        fprintf(outfile, "locus_%d\t", j+1);
    fprintf(outfile, "Probability\n");
  }
  else {
    fprintf(outfile, "\n");
    if (names) {
      nlen =0;
      for (j=0; j<n_loci; j++) {
        len = strlen(names[j]);
        if (len>nlen)
          nlen = len;
      }
      for (k=1; k<=nlen; k++) {
        if (so) {
          if (k<nlen)
            fprintf(outfile,"             ");
          else
            fprintf(outfile,"Subject id   ");
        }
        for (j=0; j<n_loci; j++) {
          if (k<=strlen(names[j])) {
            rep=*(names[j]+k-1);
            if(rep=='_') rep='|';
            if(all_snps)
              fprintf(outfile, "%c", rep);
              else fprintf(outfile, " %c", rep);
            }
          else {
            if(all_snps) fprintf(outfile, " ");
            else fprintf(outfile, "  ");
            }
        }
        if (k==nlen) {
          fprintf(outfile, "\tProbability");
          if (!so)
            fprintf(outfile, "\tCumulative");
        }
        fprintf(outfile, "\n");
      }
    }
    else {
      if (so) {
/*      if (k<nlen)
          fprintf(outfile,"             ");
        else
*/        fprintf(outfile,"Subject id   ");
      }
      for (j=0; j<n_loci; j++) {
        if (!((j+1)% 10)) {
          if(all_snps) fprintf(outfile, "+");
          else fprintf(outfile, " +");
        }
        else {
          if(all_snps) fprintf(outfile, ".");
          else fprintf(outfile, " .");
        }
      }
      fprintf(outfile, "\tProbability");
      if (!so)
        fprintf(outfile, "\tCumulative");
      fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
  }
  if (!so) { /* List haplotypes and their population frequencies */
    cum = 0.0;
    for (i=0; i<n_hap; i++) {
      for (j=0; j<n_loci; j++){
        k = order? order[j]: j;
        if(numeric)
          sprintf(c,"%d",list[i]->loci[k]);
        else {
/*
          sprintf(c, "%d", allele_code(list[i]->loci[k], coding[j]));
*/
          l=list[i]->loci[k];
          strcpy(c, l? (l==1? coding[j].one : coding[j].two) : " ");
          if(coding[j].anum==1) sprintf(c,"%d",list[i]->loci[k]);
        }
      if (tabdel) {
        if(all_snps) fprintf(outfile,"%s\t",c);
        else fprintf(outfile, "%2s\t", c);
        }
      else {
        if(all_snps) fprintf(outfile, "%s", c);
        else fprintf(outfile, "%2s", c);
        }
      }
      cum += list[i]->prior;
      if (tabdel)
        fprintf(outfile, "%11.6f\n", list[i]->prior);
      else
        fprintf(outfile, "\t%11.6f\t%10.6f\n", list[i]->prior, cum);
      nwr ++;
    }
  }
  else {  /* List assignments to subjects and their posterior probabilities */
    hs = list;
    he = list + n_hap;
    while (hs < he) {
      if (of_max>0.0) {
        h = hs;
        max_post = 0.0;
        do {
          id = (*h)->id;
          if ((*h)->posterior > max_post)
            max_post = (*h)->posterior;
          h += 2;
        } while ((h<he) && ((*h)->id==id));
        max_post *= of_max;
        hn = h;
      }
      else {
        hn = he;
      }
      for (h=hs; h<hn;) {
        out = of_max>0.0? ((*h)->posterior >= max_post):
          ((*h)->posterior > 0.0);
        if (out) {
          for (chr=1; chr<3; chr++, h++) {
            if (tabdel)
              fprintf(outfile,"%s\t%1d\t", (*h)->id, chr);
            else
              fprintf(outfile,"%10s %1d ", (*h)->id, chr);
            for (j=0; j<n_loci; j++) {
              k = order? order[j]: j;
              if(numeric)
                sprintf(c,"%d",(*h)->loci[k]);
              else {
/*
                sprintf(c,"%d",allele_code((*h)->loci[k], coding[j]));
*/
                l=(*h)->loci[k];
                strcpy(c,l? (l==1? coding[j].one : coding[j].two) : " ");
                if(coding[j].anum==1) sprintf(c,"%d",(*h)->loci[k]);
              }
              if (tabdel) {
                if(all_snps) fprintf(outfile, "%s\t", c);
                else fprintf(outfile, "%2s\t", c);
                }
              else {
                if(all_snps) fprintf(outfile, "%s", c);
                else fprintf(outfile, "%2s", c);
                }
            }
            if (tabdel)
              fprintf(outfile, "%10.7f\n", (*h)->posterior);
            else
              fprintf(outfile, "\t%10.7f\n", (*h)->posterior);
          }
          nwr ++;
        }
        else {
          h += 2;
        }
      }
      hs = hn;
    }
  }
  return nwr;
}

/* Random number generation */

static int snd_call=0;
static double snd_save;
static double aprev=0.0, c1, c2, c3, c4, c5;
#define E 2.71828182

/* Pair of N(0,1) deviates (Alg 3.6 of Ripley) */

void norm2(double *g1, double *g2) {
  double u1, u2, w, c;
  do  {
    u1 = 2.0* URAN ()-1.0;
    u2 = 2.0* URAN ()-1.0;
    w = u1*u1 + u2*u2;
  }
  while (w>=1.0);
  c = sqrt(-2.0*log(w)/w);
  *g1 = u1*c;
  *g2 = u2*c;
  return;
}

/* Generate standard normal deviate */

double snd() {
   double this_one;
   snd_call = !snd_call;
   if (snd_call) {
     norm2(&this_one, &snd_save);
     return(this_one);
   }
   else
     return(snd_save);
}

/* Random gamma variates */

double rangamma(double alpha) {
  double i, t, u, ustar, b, p, x, aa, u1, u2, w, a9;

  if (alpha==1.0) {
    /* Unit exponential by Alg 3.7 of Ripley */
    for ( i=0.0; ; i++) {
      t = u = URAN ();
      do
        {
          ustar = URAN ();
          if (u <= ustar)
            return i+t;
          else
            u = URAN ();
        }
      while (u<ustar);
    }
  }
  else if (alpha<1.0) {
  /* Alg 3.19 of Ripley */
    for(b=(alpha+E)/E;;) {
      p = b * URAN ();
      if (p<=1.0){
        x = pow(p, 1./alpha);
        if (x<=-log(URAN ()))
          return(x);
      }
      else{
        x = -log((b-p)/alpha);
        if (pow(x,alpha-1.0)>=URAN ())
          return(x);
      }
    }
  }
  else if (alpha<100.0) {
    /* Alg 3.20 of Ripley */
    if (alpha!=aprev) {
      aprev = alpha;
      c1 = alpha-1.0;
      aa = 1.0/c1;
      c2 = aa*(alpha-1.0/(6.0*alpha));
      c3 = 2.0*aa;
      c4 = c3 + 2.0;
      if (alpha>2.5)
        c5 = 1.0/sqrt(alpha);
    }
    for(;;) {
      do {
        u1 = URAN ();
        u2 = URAN ();
        if (alpha>2.5)
          u1 = u2 + c5*(1.0-1.86*u1);
      }
      while ( (u1>=1.0) || (u1<=0.0) );
      w = c2*u2/u1;
      if ( (c3*u1 + w + 1.0/w)<=c4 || (c3*log(u1) - log(w) + w)<1.0 )
        return(c1*w);
    }
  }
  else {
    /* Gamma variate by cube root approximation as given by Abramowitz & Stegun
       section 26.4.14. Perhaps necessary for very large gamma? */
    a9 = alpha*9.0;
    w = (a9 - 1.0 + sqrt(a9)*snd())/a9;
    return alpha*w*w*w;
  }
}

double **allocateU(int allele[])
{
   int size1, size2; /*sizes for the 2 dimensions*/
   int v, i; /*indices for the 2 dimeansions*/
   double ** newU; /*array to be returned*/

   size1 = n_loci;
   newU = (double**) malloc(size1 * sizeof(double*));
   if (NULL == newU) {
     REprintf("\nCould not allocate first dim of U\n");
     error("%d",EXIT_FAILURE);
   }
   for (v = 0; v < size1; v++) {
     size2 = allele[v];
     newU[v] = (double *) malloc(size2 * sizeof(double));
     if (NULL == newU[v]) {
       REprintf("\nCould not allocate second dim of U level v %d\n ", v);
       error("%d",EXIT_FAILURE);
     }
     for(i = 0; i < size2; i++) newU[v][i] = 0.0;
   }
   return(newU);
}

void freeU(double **oldU)
{
   int v;
   int size1;

   size1 = n_loci;
   for (v = 0; v < size1; v++) free(oldU[v]);
   free(oldU);
}
