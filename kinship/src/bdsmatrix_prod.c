/* $Id: bdsmatrix_prod.c,v 1.2 2002/12/26 22:54:50 therneau Exp $ */
/*
**  Product of a block-diagonal matrix and a regular one
**
** nb		number of blocks       for the bdsmatrix
** bsize	the block sizes             ""
** bmat		the vector of blocks	    ""
** rmat         right hand matrix           ""
** offdiag      off diagonal element        ""
** ydim         dimensions of the matrix on the left-hand side
** y	        the left hand matrix, which will be overwritten
** temp		a scratch vector of doubles
** itemp     	a scratch vector of integers
**
** The S code has already verified that ydim[0] equals the dimension
** of the bdsmatrix, thus the dim of rmat can be found by subtraction
*/
#include "bdsS.h"
void bdsmatrix_prod(Sint *nb,	  Sint *bsize,     Sint *ydim,
		    double *bmat, double *rmat,    double *offdiag,
		    double *temp, Sint *itemp,     double *y) {

    int nblock;
    int nrow, ncol;
    int brow, rrow;
    int i,j, k, col;
    int blocksize, offset, irow, n, block;
    double x;
    double offsum;
    
    nblock = *nb;
    nrow   = ydim[0];
    ncol   = ydim[1];
    
    brow =0;	/* number of rows in the block diagonal portion */
    for (i=0; i<nblock; i++) brow += bsize[i];
    rrow = nrow -brow;   /* this will be 0 if there is no rmat */

    /*
    ** First, deal with special case of a non-zero off diagonal element.
    ** If B is the main portion, we have B = (B-off) + off,
    **  so B%*%A = (B-off)%*% A + (off * 1'1) %*% A, the latter is a constant
    **  per column.  We don't adjust rmat.
    */
    if (*offdiag != 0) {
	n =0;
	for (block=0; block < nblock; block++) {
	    blocksize = bsize[block];
	    j = (blocksize * (blocksize+1))/2;
	    for (i=0; i<j; i++ ) 
		bmat[n++] -= *offdiag;
	    }
	}

    /*
    ** Now do one colum of the matrix at a time
    **   offset = the index if the first element of y, for this column of y
    */
    for (col=0; col<ncol; col++) {
	offset = col*nrow;
	offsum = 0;
	if (*offdiag !=0) {
	    /* Get the sum of this colum of y, multiplied by offdiag */
	    for (i=0; i<brow; i++) offsum += y[offset+i];
	    offsum *= *offdiag;
	    }
	/*
	** The indexing below is a little opaque, but fast when done right
	** irow = the current row of bmat being processed
	** n    = the index of the [irow,irow] element of bmat
	** itemp = the indices of the elements representing bmat[irow,]
	**    Say the first block is 4 by 4
	**     	irow=0, itemp= 0 1 2 3
	**    	irow=1, itemp= 1 4 5 6
	**   	irow=2, itemp= 2 5 7 8
	**	irow=3, itemp= 3 6 8 9
	**
	**  The product of bmat[irow,] and y[,col] is placed in temp[irow]
	*/
	irow =0;
	n=0;
	for (block=0; block < nblock; block++) {
	    blocksize = bsize[block];
	    k = irow;                  /* row number of upper corner of block*/
	    for (i=0; i<blocksize; i++) itemp[i] = n+i;  /* setup */
	    for (i=0; i<blocksize; i++) {  /* march down the rows */
		x = 0;
		for (j=0; j<blocksize; j++) {
		    x += bmat[itemp[j]] * y[offset+j + k];
		    if (j>i) itemp[j] += blocksize - (i+1);
		    else     itemp[j] += 1;
		    }
		temp[irow] = x;
		
		irow++;
		n += blocksize -i;
		}
	    }
	
	/* Add in the rmat part, if present */
	if (rrow >0) {
	    /* First, the pieces on the rhs of the block-diagonal part */
	    for (irow=0; irow<brow; irow++) {
		x =0;
		k =0;
		for (j=0; j<rrow; j++) {
		    x += rmat[irow+k] * y[offset+j+brow];
		    k += nrow;
		    }
		temp[irow] += x;
		}
	    
	    /* Now, the product of rmat and y */
	    k =0;
	    for (i=0; i<rrow; i++) {
		x =0;
		for (j=0; j<nrow; j++) {
		    x += rmat[j+k] * y[offset + j];
		    }
		k += nrow;
		temp[i+brow] = x;
		}
	    }
	
	/*
	**  Copy the temp vector back over the top of y
	*/
	for (i=0; i<brow; i++)    y[offset+i] = temp[i] + offsum;
	for (i=brow; i<nrow; i++) y[offset+i] = temp[i];
	offset += nrow;
	}
    }
