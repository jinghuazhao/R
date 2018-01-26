/* $Id: bdsmatrix_prod3.c,v 1.1 2003/01/07 22:09:46 therneau Exp $ */
/*
**  Product of a gchol.bdsmatrix object and a regular one
**
** nb		number of blocks       for the bdsmatrix
** bsize	the block sizes             ""
** bmat		the vector of blocks	    ""
** rmat         right hand matrix           ""
** ydim         dimensions of the matrix on the left-hand side
** y	        the left hand matrix, which will be overwritten
** temp		a scratch vector of doubles
**
** The S code has already verified that ydim[0] equals the dimension
** of the bdsmatrix, thus the dim of rmat can be found by subtraction
*/
#include "bdsS.h"
double sqrt(double);

void bdsmatrix_prod3(Sint *nb,	  Sint *bsize,     Sint *ydim,
		    double *bmat, double *rmat,    
		    double *temp, double *y) {

    int nblock;
    int nrow, ncol;
    int brow, rrow;
    int i,j, k, col;
    int itemp;
    int nk;
    int blocksize, offset, irow, n, block;
    double x;
    
    nblock = *nb;
    nrow   = ydim[0];
    ncol   = ydim[1];
    
    brow =0;	/* number of rows in the block diagonal portion */
    for (i=0; i<nblock; i++) brow += bsize[i];
    rrow = nrow -brow;   /* this will be 0 if there is no rmat */

    /*
    ** Do one colum of the matrix at a time
    */
    for (col=0; col<ncol; col++) {
	offset = col*nrow;
	/*
	** The indexing below is a little opaque, but fast when done right
	** irow = the current row of bmat being processed
	** n    = the index of the [irow, irow] element
	** k    = irow of the start of this block
	** nk   = the index of the [k,k] element of bmat
	** itemp = the indices of the elements representing bmat[irow,]
	**    Say the first block is 4 by 4
	**     	irow=0, itemp= 0
	**    	irow=1, itemp= 1 4
	**   	irow=2, itemp= 2 5 7
	**	irow=3, itemp= 3 6 8 9
	**
	**  The product of bmat[irow,] and y[,col] is placed in temp[irow]
	*/
	irow =0;
	n=0;
	for (block=0; block < nblock; block++) {
	    blocksize = bsize[block];
	    k = irow;
	    nk =n;
	    for (i=0; i<blocksize; i++) {  /* march down the rows */
		x = 0;
		y[offset+irow] *= sqrt(bmat[n]); /* y = sqrt(D) y */ 
		x = y[offset+irow];              /* L[i,i] * y[i] */
		itemp = nk +i;
		for (j=0; j<i; j++) {
		    x += bmat[itemp] * y[offset+j + k];
		    itemp += blocksize - (j+1);
		    }
		temp[irow] = x;
		irow++;
		n += blocksize -i;
		}
	    }
	
	/* 
	** Add in the rmat part, if present 
	**  n is now the index of the first element of the row of rmat
	*/
	n=0;
	for (i=0; i<rrow; i++) {
	    y[offset + irow] *= sqrt(rmat[n+irow]);
	    x = y[offset + irow];
	    for (j=0; j< irow; j++) 
		x += rmat[n+j] * y[offset+j];
	    temp[irow] = x;
	    irow++;
	    n += nrow;
	    }
	/*
	**  Copy the temp vector back over the top of y
	*/
	for (i=0; i<nrow; i++) y[offset+i] = temp[i];
	offset += nrow;
	}
    }
