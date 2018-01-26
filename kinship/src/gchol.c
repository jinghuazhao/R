/* $Id: gchol.c,v 1.4 2003/08/22 14:39:42 therneau Exp $ */
/*
** General cholesky decompostion
*/
#include "bdsS.h"
#include "kinproto.h"
#include "dmatrix.h"

void gchol(Sint *n2, double *matrix, double *toler) {
    int i,j;
    int n/*, flag*/;
    double **mat;

    n = *n2;
    mat = dmatrix(matrix, n, n);

    i = cholesky5(mat, n, *toler);
    *toler = i;

    /* zero out the upper triangle */
    for (i=0; i<n; i++) {
	for (j= i+1; j<n; j++) mat[i][j] =0;
	}
    }

void gchol_solve(Sint *n2, double *matrix, double *y, Sint *flag2) {
    int n;
    double **mat;
    int flag;

    n = *n2;
    flag = *flag2;
    mat = dmatrix(matrix, n, n);

    chsolve5(mat, n, y, flag);
    }
    
void gchol_inv(Sint *n2, double *matrix, Sint *flag2) {
    int n;
    double **mat;
    int i,j;
    int flag;

    n = *n2;
    flag = *flag2;
    mat = dmatrix(matrix, n, n);

    chinv5(mat, n, flag);

    /*
    **  the result of chinv5 has the inverse of L and full inverse
    **  all packed together
    */
    if (flag ==1) {
	/* 
	** return L-inverse, by zeroing out the other part
	*/
	for (i=0; i<n; i++) {
	    mat[i][i] = 1;
	    for (j=i+1; j<n; j++) mat[i][j] =0;
	    }
	}
    else {
	/* 
	** replicate the lower part into the upper one, for a symmetric result
	*/
	for (i=0; i<n; i++) {
	    for (j=i+1; j<n; j++) mat[j][i] = mat[i][j];
	    }
	}
    }
