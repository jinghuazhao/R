/* $Id: cholesky4.c,v 1.3 2003/02/05 19:19:25 therneau Exp $ */
/*
** subroutine to do generalized Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
** This is the second specialized form, for the correlated frailty problem.
**   The matrix C has C[1:m, 1:m] block diagonal and the remainder dense.
**
** arguments are:
**     n         the size of the matrix to be factored
**     nblock    the number of blocks
**     bsize     vector of block sizes (m above is sum(bsize))
**     bd        compressed contents of the block diagonal matrices,
**	            their lower triangles strung together
**                  length(bd) = sum[ bsize*(bsize+1)/2 ]
**     **matrix  a ragged array containing the dense portion
**     toler     the threshold value for detecting "singularity"
**
**  The factorization is returned in the lower triangle, D occupies the
**    diagonal and the upper triangle is unused.
**
**  Return value:  the rank of the matrix 
**
**  If a column is deemed to be redundant, then that diagonal is set to zero.
**
**  For documentation below, A[] is the uncompressed matrix.  We think of
**   ourselves as operating on it, with the compression just an indexing
**   nuisance (a real pain-in-the-ass nuisance, actually).
**
**   Terry Therneau
*/
#include <math.h> 
int cholesky4(double **matrix, int n, int nblock, int *bsize,
	      double *bd, double toler) {

    double temp;
    int  i,j,k, m;
    double eps, pivot;
    int rank;
    int n2;
    int ii, ji, jj, kj, blocksize;
    int block;
    

    m=0;
    eps =0;

    /* Find the max diagonal element of the matrix, for scaling */
    ii =0;
    for (block=0; block<nblock; block++) {
	m += bsize[block];
	for (blocksize= bsize[block]; blocksize>0; blocksize--) {
	    if (fabs(bd[ii]) < eps) eps = bd[ii];  
	    ii += blocksize;
	    }
	}
    n2=n-m;
    for (i=0; i<n2; i++) if (fabs(matrix[i][i+m]) > eps)  
	eps = fabs(matrix[i][i+m]);
    if (eps > 0.0) eps *= toler;
    else eps = toler;              /* just in case diagonal ==0 */

    /*  
    ** Do the Cholesky  for the blocks diagonal portion
    */
    ji=0; 
    rank =0;
    ii =0;
    i =0;
    for (block=0; block<nblock; block++) {  
	/* process block i */
	for (blocksize= bsize[block]; blocksize>0; blocksize--) {
	    jj = ii;
	    pivot = bd[ii];
	    if (fabs(pivot) < eps) {
		for (j=0; j<blocksize; j++) bd[ii+j] =0;
		for (j=0; j<n2; j++) matrix[j][i] =0;
		}
	    else {
		rank++;
		/*
		** sweep out the rest of the block
		*/
		for (j=1; j<blocksize; j++) {
		    ji = ii + j;             /* points to A[j,i] */
		    jj += blocksize +1 -j;   /* points to A[j,j] */
		    temp = bd[ji] /pivot;
		    bd[ji] = temp;
		    bd[jj] -= temp*temp*pivot;
		    for (k=(j+1); k<blocksize; k++) {
			kj = jj + (k-j);
			bd[kj] -= temp*bd[ii+k];
			}
		    for (k=0; k<n2; k++)
			matrix[k][i+j] -= temp*matrix[k][i];
		    }

		/* sweep out dense rows below the block */
		for (j=0; j<n2; j++) {
		    temp = matrix[j][i] / pivot;
		    matrix[j][i] = temp;
		    matrix[j][j+m] -= temp*temp*pivot;
		    for (k=(j+1); k<n2; k++) 
			matrix[k][j+m] -= temp*matrix[k][i];
		    }
		}
	    ii+= blocksize;
	    i += 1;
	    }
	}

    /* 
    ** Now the final corner, which isn't below the block diagonal portion
    */
    for (i=0; i<n2; i++) {
        pivot = matrix[i][i+m];
        if (fabs(pivot) < eps) {
            for (j=i; j<n2; j++) matrix[j][i+m] =0;  /* zero the column */
            }
        else  {
            rank++;
            for (j=(i+1); j<n2; j++) {
                temp = matrix[j][i+m]/pivot;
                matrix[j][i+m] = temp;
                matrix[j][j+m] -= temp*temp*pivot;
                for (k=(j+1); k<n2; k++) matrix[k][j+m] -= temp*matrix[k][i+m];
                }
            }
        }
    return(rank);
    }
