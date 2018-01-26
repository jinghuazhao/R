/* $Id: coxfit6c.c,v 1.5 2003/08/09 21:45:31 therneau Exp $  */
/*
**  When iteration is done, this returns all of the other pieces that
**   were saved.
**
**  u   :    first derivative
**  imatb  : sparse portion of the Cholesky
**  imatr  : dense portion of the Cholesky
**  hinvb  : sparse portion of the invers
**  hinvr  : dense portion of the inverse
**  rank   : rank of imat = number of non-zero diagonal elements
**  ny     : on input: 2=Cox model, 3= (start, stop]
*/
#include "coxmeS.h"
#include "coxfit6.h"
#include "kinproto.h"

/* the next line is just so that I can use "c6.n" instead of "coxfit6.n", etc*/
extern struct coxfit6 c6;

void coxfit6c(double *u,      double *imatb,
	      double *imatr, double *hinvb, double *hinvr,
	      Sint  *rank,    Sint *ny) {
    int i,j,k;
    double *dptr;
    int nvar, nvar2, nvar3;
    int nfac;  /* number of factors */
    int    nf, ns;

    nf   = c6.nfrail;
    nvar = c6.nvar;
    ns   = c6.nsparse;
    nfac = c6.nfactor;
    nvar2= nvar + (nf - nfac);  /* number of cols of X */
    nvar3= nvar + nf; /* total number of coefficients */

    k=0;
    for (i=0; i<nvar3; i++) {
	u[i]    = c6.u[i];
	if (c6.imat[i][i] >0) k++;
	}
    *rank =k;

    for (i=0; i<c6.tblock; i++) 
	imatb[i] = c6.imatb[i];

    dptr = imatr;
    for (i=ns; i<nvar3; i++) {
	for (j=0; j<nvar3; j++) *dptr++ = c6.imat[i][j];
	}

    chinv4(&(c6.imat[ns]), nvar3, c6.nblock, 
	     c6.bsize,  c6.imatb, 1);

    for (i=0; i<c6.tblock; i++) 
	hinvb[i] = c6.imatb[i];

    dptr = hinvr;
    for (i=ns; i<nvar3; i++) {
	for (j=0; j<nvar3; j++) *dptr++ = c6.imat[i][j];
	}

    /*
    ** Release all the memory
    **  Somewhere I read that it's better to do this in the reverse order
    **  to which it was allocated
    */
    if (c6.calc2 ==1) {
	Free(c6.dlag2[0]);
	Free(c6.dlag2);
	Free(c6.dlag1);
	Free(c6.tlist);
	}
   if (*ny ==3) {
	Free(c6.start);
	Free(c6.sort1);
	}
    Free(c6.status);
    Free(c6.a);
    if (nvar2 >0) {
	Free(c6.cmat2[0]);
	Free(c6.cmat2);
	Free(c6.cmat[0]);
	Free(c6.cmat);
	}
    Free(c6.imatb);
    Free(c6.imat);
    if (nvar2 >0) {
	Free(c6.x[0]);
	Free(c6.x);
	}
    if (c6.nfx >0) {
	Free(c6.findex);
	Free(c6.fx);
	}
    if (c6.nblock >0) {
	Free(c6.bsize);
	}
    }
