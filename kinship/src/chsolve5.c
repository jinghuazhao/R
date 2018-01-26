/* $Id: chsolve5.c,v 1.3 2003/10/20 22:27:13 therneau Exp $ */
/*
** Solve the equation Ab = y, where the cholesky decomposition of A and y
**   are the inputs.
**
** Input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**        y[n] contains the right hand side
**
**  y is overwriten with b
**
** This differs from chsolve2 only in the flag, which if =1 gives the
**   "half solution" FD^{1/2}b = y.  This is used to transform
**    data to "uncorrelated" form in linear regression.  This operation
**    requires that D>=0, however.
**
**  Terry Therneau
*/
#include <math.h>

void chsolve5(double **matrix, int n, double *y, int flag)
     {
     int i,j;
     double temp;

     /*
     ** solve Fz =y, where z== DF'b
     */
     for (i=0; i<n; i++) {
	  temp = y[i] ;
	  for (j=0; j<i; j++)
	       temp -= y[j] * matrix[i][j] ;
	  y[i] = temp ;
	  }

     if (flag==1) {
	 /*
	 ** solve D^{1/2}b =z
	 */
	 for (i=0; i<n; i++) {
	     if (matrix[i][i]<=0) y[i]=0;
	     else y[i] /= sqrt(matrix[i][i]);
	     }
	 }
     else {
	 /*
	 ** solve DF'b = z
	 */
	 for (i=(n-1); i>=0; i--) {
	     if (matrix[i][i]==0)  y[i] =0;
	     else {
		 temp = y[i]/matrix[i][i];
		 for (j= i+1; j<n; j++)
		     temp -= y[j]*matrix[j][i];
		 y[i] = temp;
		 }
	     }
	  }
     }
