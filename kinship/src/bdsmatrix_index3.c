/* $Id: bdsmatrix_index3.c,v 1.2 2002/12/26 22:54:49 therneau Exp $ */
/*
** The bsdmatrix creation function allows a user to input the full contents
**   of each block diagonal matrix.  This function gives the indices of the
**   lower triangle portion.
**  The "1+" on all output is to S-subscripts, starting at 1
*/
#include "bdsS.h"
#include <stdio.h>
void bdsmatrix_index3(Sint *nblock, Sint *bsize, Sint *index) {
    int i, j; 
    int blocksize;
    int nc;              /* number returned so far, for the index vector */
    int block;           /* block currently being processed */
    int irow;            /* global row counter */
    int pos;             /* current position in the blocks array */
    int lastrow;

    irow=0; 
    nc=0;
    pos =0;
    for (block=0; block < *nblock; block++) {
	blocksize = bsize[block];
	lastrow = irow + blocksize;

	for (i=0; i<blocksize; i++) {
	    for (j=irow; j< lastrow; j++)
		index[nc++] = 1 +pos + j - irow;

	    pos += blocksize +1;
	    irow++;
	    }
	pos -= blocksize;
	}
    }
