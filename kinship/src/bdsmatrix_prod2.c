/* $Id: bdsmatrix_prod2.c,v 1.1 2002/12/30 21:59:22 therneau Exp $ */
/*
**  Product of a block-diagonal matrix and a vector
**   This is a trivial variant of bdsmatrix_prod, to be called from
**  within C code instead of from S.
**   It assumes that y is a vector, of the appropriate length, and that
**  offdiag=0.
**
** nb		number of blocks       for the bdsmatrix
** bsize	the block sizes             ""
** nrow         total number of rows/cols in the bdsmatrix
** bmat		the vector of blocks	    ""
** rmat         right hand matrix           ""
** y	        the left hand vector
** result	the result of the calcluation -- a vector like "y"
** itemp     	a scratch vector of integers, of length >= max(bsize)
*/
void bdsmatrix_prod2(int nblock,     int *bsize,     int nrow,
		     double *bmat,   double *rmat,  
		     double *y,      double *result, int *itemp) {

    int brow, rrow;
    int i,j, k;
    int blocksize, /*offset, */irow, n, block;
    double x;
/*  double offsum;*/
    
    brow =0;	/* number of rows in the block diagonal portion */
    for (i=0; i<nblock; i++) brow += bsize[i];
    rrow = nrow -brow;   /* this will be 0 if there is no rmat */


    /*
    ** The indexing below is a little opaque, but fast when done right
    ** irow = the current row of bmat being processed
    ** n    = the index of the [irow,irow] element of bmat
    ** itemp = the indices of the elements representing bmat[irow,]
    **    Say the first block is 4 by 4
    **     	irow=0, itemp= 0 1 2 3
    **    	irow=1, itemp= 1 4 5 6
    **   	irow=2, itemp= 2 5 7 8
    **	        irow=3, itemp= 3 6 8 9
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
		x += bmat[itemp[j]] * y[j + k];
		if (j>i) itemp[j] += blocksize - (i+1);
		else     itemp[j] += 1;
    	        }
	    result[irow] = x;
    	
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
		x += rmat[irow+k] * y[j+brow];
		k += nrow;
    	        }
	    result[irow] += x;
    	    }
        
        /* Now, the product of rmat and y */
        k =0;
        for (i=0; i<rrow; i++) {
	    x =0;
	    for (j=0; j<nrow; j++) {
		x += rmat[j+k] * y[j];
    	        }
	    k += nrow;
	    result[i+brow] = x;
    	    }
	}
    }
    

