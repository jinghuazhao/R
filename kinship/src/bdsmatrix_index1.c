/* $Id: bdsmatrix_index1.c,v 1.2 2002/12/26 22:54:48 therneau Exp $ */
/*
** Symmetric block diagonal subscripting
**     given a set of row/col numbers 'zed', return
**		    a. the position vector for A[zed,zed] in the compact 
**                       matrix, with zeros for "off diagonal" elements
**		    b. the postion vector of diag(A[zed,zed]) in compact
**		    c. the position vector of the new lower triangle
**                      in this case, bsize is also fixed up.
**                      This is used when, say, the [10:15,10:15] portion of
**                      a bdsmatrix is requested, the result is a bdsmatrix.
**    If flag[0]=0 'a' isn't done, flag[1]=0 -> b isn't done, etc
**    The row/col numbers are assumed to be in sorted order
**
**  The "1+" on all output is to S-subscripts, starting at 1
*/
#include "bdsS.h"
void bdsmatrix_index1(Sint *nblock, Sint *bsize, Sint *flag,
		     Sint *nrow,   Sint *rows,  Sint *indexa,
		     Sint *indexb, Sint *indexc) {
    int i, j, k; 
    int blocksize;
    int na, nb, nc;      /* current pos in indexa, indexb, or indexc vector */
    int block;           /* block currently being processed */
/*  int ib;            *//* row within block */
    int irow;            /* global row counter */
    int jrow;            /* current progress through the "desired" list */
    int pos;             /* current position in the blocks array */
    int firstrow, lastrow;/*first and last rows of a block */
    int newblock;         /*final size of current block */

    irow=0; jrow=0;
    nb=0; na=0; nc=0;
    pos =0;
    for (block=0; block < *nblock; block++) {
  	blocksize = bsize[block];
	firstrow = irow;
	lastrow  = irow + blocksize -1;
	newblock =0;

	for (i=0; i<blocksize; i++) {
	    if (irow == rows[jrow]) {
		newblock++;
		if (flag[0] ==1) {
		    for (j=jrow; j< *nrow && rows[j] <= lastrow; j++) {
			k = pos + rows[j] - irow;
			/* npos tracks the jrow,jrow in indexa*/
			indexa[na + j - jrow]  = 1+ k;
			indexa[na + (j -jrow)*(*nrow)] =1+ k;
			}
		    }

		if (flag[1]==1) indexb[nb++] = 1+ pos;
		    
		if (flag[2]==1) {
		    for (j=jrow; j< *nrow && rows[j] <= lastrow; j++) {
			k = pos + rows[j] - irow;
			indexc[nc++] = 1+k;
			}
		    }
		    
		na += *nrow +1;
		jrow++;
		if (jrow == *nrow) {
		    bsize[block] = newblock;
		    for (j=block+1; j< *nblock; j++) bsize[j] =0;
		    return;
		    }
		}
	    pos += blocksize -i;
	    irow++;
	    }
	bsize[block] = newblock;
	}
    }
