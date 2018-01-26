/* $Id: chinv4.c,v 1.2 2002/12/26 22:54:50 therneau Exp $ */
/*
** Invert a cholesky decompostion of L where C = LDL', L lower triangular.  
**    For each column of the inverse:
**	1. Create the column by backsolving Lb = e_i, where e_i is the
**          ith column of an identity matrix.  The algorithm is exactly the
**	    first half of chsolve4.c
**      2. Replace column i of L with the inverse.  Because it's lower 
**	    triangular, the inverse is lower triangular, and if we process
**          columns from left to right the overlaid portion is one that 
** 	    wouldn't be used again (the next right hand side vector has zeros
** 	    there).
**      3. Trickery: if one is careful, no scratch space for b is needed.
**
**   At the same time invert D, which is stored on the diagonal of L.
**
**   The matrix C has special structure,
**            C[1:m, 1:m] block diagonal and  C[(m+1):n, 1:n)] is dense. 
**
** arguments are:
**     **matrix, which contains the chol decomp of the dense portion
**     n         the size of the matrix to be factored
**     nblock    the number of blocks
**     bsize     vector of block sizes (m above is sum(bsize))
**     bd        compressed contents of the block diagonal matrices,
**                  their lower triangles strung together
**                  length(bd) = sum[ bsize*(bsize+1)/2 ]
**     flag      if 0, return L-inverse.  If 1, return C-inverse =
**                     (L-inverse)' D-inverse (L-inverse)
**
**  Terry Therneau
*/

void chinv4(double **matrix, int n, int nblock, int *bsize, double *bd,
            int flag) {
    int i,j, k, n2;
    int ii, /*jk,*/ block/*, blocksize*/;
    double temp;
    int m;
    int i2, ii2, j2, yi;

    m =0;
    for (i=0; i<nblock; i++) m+= bsize[i];
    n2 = n-m;

    /*
    ** i is the column of the inverse currently being worked upon, and
    ** ii is the pointer to L[i,i] in the sparse storage.
    */
    i=0;
    ii =0;
    for (block=0; block<nblock; block++) {
	for (j=bsize[block]; j>0; j--) {
	    /* Invert D */
	    if (bd[ii] ==0) continue;
	    else bd[ii] = 1/bd[ii];

	    /*
	    ** solve Fb =e , using eq 2.2.2 of A. George and A Liu, Computer 
	    **  Solution of Large Sparse Positive Definite Systems, 
	    **  Prentice-Hall, 1981.
	    ** e = vector of 0's with a 1 at position i
	    ** i2 = column of L currently being "solved against"
	    */
	    /* backsolve wrt the "1" on the diagonal */
	    for (k=1; k<j; k++) bd[ii+k] *= -1;
	    for (k=0; k<n2; k++) matrix[k][i] *= -1;
	    
	    /* backsolve wrt the rest of this block */
	    i2 = i+1;
	    ii2= ii+j;
	    yi = ii+1;   /*points to "y[i]", a la chsolve4 */
	    for (j2=j-1; j2>0; j2--) {
		temp = bd[yi];
		for (k=1; k<j2; k++) 
		    bd[yi+k] -= bd[ii2+k]*temp;
		for (k=0; k<n2; k++)
		    matrix[k][i] -= matrix[k][i2] *temp;
		yi++;
		i2++;
		ii2 += j2;
		}

	    /* backsolve wrt to the dense corner */
	    for (j2=0; j2<n2; j2++) {
		temp = matrix[j2][i];
		for (k=j2+1; k<n2; k++) matrix[k][i] -= temp * matrix[k][j2+m];
		}

	    i++;
	    ii += j;
	    }
	}

    /*
    ** Finish up by doing the inverse for the dense corner
    */
    for (i=0; i<n2; i++) {
	if (matrix[i][i+m] > 0) {
	    matrix[i][i+m] = 1/matrix[i][i+m];
	    for (k=i+1; k<n2; k++) matrix[k][i+m] *= -1;
	    for (j=i+1; j<n2; j++) {
		temp = matrix[j][i+m];
		for (k=j+1; k<n2; k++) matrix[k][i+m] -= temp * matrix[k][j+m];
		}
	    }
	}

    if (flag==0) return;  /* I'm done, no need to multiply */

    i=0;
    ii=0;
    for (block=0; block<nblock; block++) {
	for (j=bsize[block]; j>0; j--) {
	    if (bd[ii] ==0) { /* this column of the inverse is all zeros */
		for (k=0; k<j; k++) bd[ii+k] = 0;
		for (k=0; k<n2; k++)matrix[k][i] =0;
		}
	    else {
		/*
		** compute inner product of cols i and i,
		**  and at the same time scale col i for the further work
		*/
		temp = bd[ii];
		ii2= ii+j;
		for (k=1; k<j; k++) {
		    temp += bd[ii+k]* bd[ii+k] * bd[ii2];
		    bd[ii+k] *= bd[ii2];
		    ii2 += (j-k);
		    }
		for (k=0; k<n2;k++) {
		    temp += matrix[k][i]*matrix[k][i] * matrix[k][k+m];
		    matrix[k][i] *= matrix[k][k+m];
		    }
		bd[ii] = temp;

		/* inner product of cols i and i2, for i2 in this block */
		i2 = i+1;
		ii2= ii+j;
		for (j2=1; j2<j; j2++) {
		    temp = bd[ii+j2];   /* the multiplication by '1' on diag*/
		    for (k=1; k<(j-j2); k++) temp += bd[ii+j2+k]* bd[ii2+k];
		    for (k=0; k<n2; k++) temp += matrix[k][i]*matrix[k][i2];
		    bd[ii+j2] = temp;
		    i2++;
		    ii2 += (j-j2);
		    }

		/* 
		** inner product for i and other columns in the block
		**   diagonal portion is not zero (if n2>0), but we don't have
		**   anywhere to store it and so ignore it
		** finish with inner product of this block and the dense  
		*/
		for (j2=0; j2< n2; j2++) {
		    temp = matrix[j2][i];
		    i2 = j2 + m;
		    for (k=j2+1; k<n2; k++) 
			temp += matrix[k][i]*matrix[k][i2];
		    matrix[j2][i] = temp;
		    }
		}
	    ii += j;
	    i++;
	    }
	}
    /*
    ** Now the dense corner
    */
    for (i=0; i<n2; i++) {
	if (matrix[i][i+m] ==0) { /* zero out the row/column */
	    for (k=i; k<n2; k++) {
		matrix[k][i+m] =0;
		matrix[i][k+m] =0;
		}
	    }
	else {
	    /* compute cols i and i */	
	    temp = matrix[i][i+m];
	    for (k=i+1; k<n2; k++) {
		temp += matrix[k][i+m] * matrix[k][i+m] * matrix[k][k+m];
		matrix[k][i+m] *= matrix[k][k+m];
		}
	    matrix[i][i+m] = temp;
	    
	    /* compute columns i and j */
	    for (j=i+1; j<n2; j++) {
		temp = matrix[j][i+m];
		for (k=j+1; k<n2; k++) 
		    temp += matrix[k][i+m] * matrix[k][j+m];
		matrix[j][i+m] = temp;
		matrix[i][j+m] = temp;
		}
	    }
	}
    } 
