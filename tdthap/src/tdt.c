#include <R.h>
#include <stdlib.h>
#include <stdio.h>

double asran();
void asran_seed(int, int, int);

/*
  Do computations for quadratic form tdt test

  nhap:   Number of different haplotypes
  ntran:  Number of transmissions
  haps:   Character array of haplotype descriptions
  tr:     The transmitted haplotypes (numbered 1...nhap)
  un:     The untransmitted haplotypes (numbered 1...nhap)
  nsamp:  The number of Monte Carlo samples
  funct:  If non-zero, use the similarity function. Otherwise do the 
          coventional Pearsonian test
  keep:   If non-zero, Monte Carlo values of test statistic are kept in full
  seeds:  Three integer seeds for random number generator
  res:    Result vector. First value is realised statistic and subsequent
          values are Monte Carlo values if keep is T, or the proportion 
          exceeding the observed value if keep is F
*/


int tdt_quad(int *nhap, int *ntran, char **haps, int *tr, int *un, 
	     int *nsamp, int *funct, int *keep, int *seeds, double *res) {
  /*FILE *scratch_file;*/
  int nh, nt, nl, nmc, i, j, ij, nnz=0, inz; 
  int a, b, c, *hmat, *ha, *hb;
  char *ch, **h;
  double *score, s, test, obst=0.0;
  double tdt_similarity(int *, int *);
  /* Structure type to hold non-zero similarities */ 
  typedef struct s_tag {
    int row;
    int col;
    double similarity;
    struct s_tag *next;
  } similar;
  similar *s_list=NULL, *s_last, *s_this;

  nh = *nhap;
  nt = *ntran;
  nmc = *nsamp;
  a = seeds[0];
  b = seeds[1];
  c = seeds[2];
  asran_seed(a, b, c);
  score = (double *) calloc(nh, sizeof(double));
  /* 
     If using a similarity function, calculate all the non-zero ones ... 
  */
  if (*funct) {
    /* Count loci */
    for (ch = *haps, nl = 1; *ch; ch++) if (*ch=='.') nl++;
    /* Convert haplotypes from text strings to  integer array */
    hmat = (int *) calloc((nl+1)*nh, sizeof(int));
    if (!hmat) goto noroom;
    for (ij=i=0, h=haps; i<nh; i++, h++) {
      for (ch = *h, j=0; j<nl; j++, ij++) {
	for (a=0; *ch && *ch!='.'; ch++) a = a*10+(*ch - '0');
	hmat[ij] = a;
	ch ++;
      }
      hmat[ij++] = 0; /* Terminate haplotype with zero */
    }
    /* Identify all non-zero similarities */
    nnz =0;
    s_last = (similar *) 0;
    for (a=0, ha=hmat; a<nh; a++, ha+=(nl+1)) {
      for (b=a, hb=ha; b<nh; b++, hb+=(nl+1)) {
	s = tdt_similarity(ha, hb);
	if (s!=0.0) {
	  nnz ++;
	  s_this = (similar *) malloc(sizeof(similar));
	  if (!s_this) goto noroom;
	  s_this->row = a; 
	  s_this->col = b;
	  s_this->similarity = s;
	  s_this->next = (similar *) 0;
	  if (s_last) 
	    s_last->next = s_this;
	  else 
	    s_list = s_this;
	  s_last = s_this;
	}
      }
    }
    free(hmat);
  }
  /* 
     Otherwise each haplotype is similar to itself only ... 
  */
  else {
    /* Calculate total number of occurrences of each haplotype */
    for (a=0; a<nh; a++) score[a] = 0.0;
    for (i=0; i<nt; i++) {
      score[tr[i]-1] ++;
      score[un[i]-1] ++;
    }
    s_last = (similar *) 0;
    for (nnz=a=0; a<nh; a++) {
      if (score[a]>0.0) {
	nnz++;
	s_this = (similar *) malloc(sizeof(similar));
	if (!s_this) goto  noroom;
	s_this->row = s_this->col = a;
	s_this->similarity = 1.0/score[a];
	s_this->next = (similar *) 0;
	if (s_last) 
	  s_last->next = s_this;
	else 
	  s_list = s_this;
	s_last = s_this;
      }
    }
  }
  /* 
     This loop does the test; the first time for real 
  */
  for (ij=j=0; j<=nmc; j++) {
    /* Calculate score vector  */
    for (a=0; a<nh; a++) score[a] = 0.0;
    for (i=0; i<nt; i++) {
      a = tr[i]-1; /* Remember C indexing! */ 
      b = un[i]-1;
      /* In Monte Carlo passes, swap transmitted and untransmitted haplotypes*/
      if (j && asran()>0.5) {
	score[b]++;
	score[a]--;
      }
      else {
	score[a]++;
	score[b]--;
      }
    }
    for (test=0.0, inz=0, s_this=s_list; inz<nnz; inz++, s_this=s_this->next) {
      a = s_this->row;
      b = s_this->col;
      s = s_this->similarity;
      test += score[a]*score[b]*s;
      if (a!=b) 
	test += score[b]*score[a]*s;
    }
    /* First time through loop is the observed value */
    if (!j) {
      obst = test;
    }
    else {
      if (*keep) res[j] = test;
      else if (test>=obst) ij++;
    }
  }
  res[0] = obst;
  if (!*keep) res[1] = (double) ij / (double) nmc;
  goto      tidyup;
  
  /* Insufficient dynamic memory */

noroom:
  REprintf("*** tdt.c *** Insufficient memory\n");
  res[0] = -1.0;
  goto tidyup;

  /* Return dynamic storage space */
  
tidyup:
  for (inz = 0, s_this=s_list; inz<nnz && s_this; inz++) {
    s_last = s_this;
    s_this = s_this->next;
    free(s_last);
  } 
  if (score) free(score);
  return 0;

}

/* Applied Statistics random number generator */
  
static int ix, iy, iz;

void asran_seed(int i1, int i2, int i3)
   {
   ix = i1;
   iy = i2;
   iz = i3;
   }

double asran()
   {
   double r;

   ix = (171*ix) % 30269;
   iy = (172*iy) % 30307;
   iz = (170*iz) % 30323;
   r  = (double)ix/30269.0 + (double)iy/30307.0 + (double)iz/30323.0;
   return ( r - (int) r );
   }
