/* $Id: chsolve4.c,v 1.2 2002/12/26 22:54:51 therneau Exp $ */
/*
** Solve the equation Ab = y, where the cholesky decomposition of A and y
**   are the inputs.
** This is a specialized form for the frailty problem.  The matrix C in this
**   case has C[1:m, 1:m] block diagonal and  C[(m+1):n, 1:n)] is dense. 
**
** arguments are:
**     n         the size of the matrix to be factored
**     nblock    the number of blocks
**     bsize     vector of block sizes (m above is sum(bsize))
**     bd        compressed contents of the block diagonal matrices,
**                  their lower triangles strung together
**                  length(bd) = sum[ bsize*(bsize+1)/2 ]
**     **matrix, which contains the chol decomp of the dense portion
**     y[n] contains the right hand side
**
**  y is overwriten with b
**
**  If flag=0, solve Ab=y,
**         =1, solve L sqrt(D)b=y, where A=LDL', L lower triangular
**  Terry Therneau
*/
#include <math.h>

void chsolve4(double **matrix, int n, int nblock, int *bsize,
	      double *bd, double *y, int flag) {
    int i,j, k, /*l,*/ n2;
    int ii, /*jk,*/ block, blocksize;
    double temp;
    int m;
     
    m =0;
    for (i=0; i<nblock; i++) m+= bsize[i];
    n2 = n-m;

    /*
    ** solve Fz =y , using eq 2.2.2 of A. George and A Liu, Computer Solution
    **  of Large Sparse Positive Definite Systems, Prentice-Hall, 1981.
    **
    **  first the block diagonal portion
    */
    i =0;
    ii =0;
    for (block=0; block<nblock; block++) {
        /* process block i */
        for (j=bsize[block]; j>0; j--) {
            temp = y[i];
            for (k=1; k<j; k++) 
		y[i+k] -= bd[ii+k]*temp;
            for (k=0; k<n2; k++)
		y[m+k] -= matrix[k][i] *temp;
            i++;
            ii += j;
            }
        }
    /*
    ** At this point i and ii point to the END of the bd array, and need to
    **  be saved for later use
    ** Finish up by doing the dense portion of Fz=y
    */
    for (j=0; j<n2; j++) {
         temp = y[j+m];
         for (k=j+1; k<n2; k++) y[k+m] -= temp * matrix[k][j+m];
         }

    /*
    ** Now the second half of the work
    */
    if (flag==1) {
        /* solve sqrt(D) b =z */
        i =0;
        ii=0;
        for (block=0; block<nblock; block++) {
            for (j=bsize[block]; j>0; j--) {
		/* ii points to A[i,i] */
		if (bd[ii] >0) y[i] /= sqrt(bd[ii]);
		else           y[i] = 0;
		i++;
		ii += j;
	        }
            }
        /* dense portion */
        for (j=0; j<n2; j++) {
            temp = matrix[j][j+i];
            if (temp>0) y[j+i] /= sqrt(temp);
            else y[i+j] = 0;
            }
        }

    else {
        /*
        ** solve DF'b =z, using equation 2.2.1
        */
        /* dense portion */
        for (j=(n2-1); j>=0; j--) {
       	  if (matrix[j][j+m]==0)  y[j+m] =0;
       	  else {
       	      temp = y[j+m]/matrix[j][j+m];
       	      for (k= j+1; k<n2; k++)
       		   temp -= y[k+m]*matrix[k][j+m];
       	      y[j+m] = temp;
       	      }
       	  }
        /* block diag portion, walking backwards through the blocks */
        for (block=nblock-1; block >=0; block--) {
       	 for (blocksize=1; blocksize <=bsize[block]; blocksize++) {
       	     i--;
       	     ii -= blocksize;
       	     if (bd[ii] ==0) y[i] =0;
       	     else {
       		 temp = y[i] / bd[ii];
       		 for (j=1; j<blocksize; j++)
       		     temp -= y[i+j]*bd[ii+j];
       		 for (j=0; j<n2; j++)
       		     temp -= y[j+m]*matrix[j][i];
       		 y[i] = temp;
       		 }
       	     }
            }
        }
    }
