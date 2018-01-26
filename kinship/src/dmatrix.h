#ifndef DMATRIX_H
#define DMATRIX_H

/* from R 1.9.0 .C/.Fortran calls are forced to associate with package name*/
/*
** the next section represents a solution to a bit of a 
**  conundrum.  
**  1. Both the gchol and the survival routines make use of
**     the subroutines chsolve2 and dmatrix.
**  2. If I include these routines, then S can be confused when I
**      load the "survival" library (as we do locally at Mayo).
**  3. If I don't, then most of the rest of the world can't find the
**      symbols, because Splus re-defines the symbols to S_dmatrix and
**      S_chsolve2 to avoid any future conflicts with me. 
**      R doesn't redefine the names, so R users that load survival would
**      find them, but not those that don't load survival.
** How do we satisfy everyone?
**
** Since the routines are short, I include a copy of them here defined as
** "static".  This uses a trivial amount more of memory for the extra
** copy, but leaves them invisible to outside routines and thus
** avoids conflicts.  The downside, of course, is that I'd have to fix
** both the global copy (for survival) and this one should I find a bug.
** After 5+ years of use, a bug is unlikely however.
*/

static double **dmatrix(double *array, int ncol, int nrow)
    {
/*S_EVALUATOR*/
    int i;
    double **pointer;

    pointer = (double **) ALLOC(nrow, sizeof(double *));
    for (i=0; i<nrow; i++) {
	pointer[i] = array;
	array += ncol;
	}
    return(pointer);
    }

#if 0
static void chsolve2(double **matrix, int n, double *y)
     {
     register int i,j;
     register double temp;

     /*
     ** solve Fb =y
     */
     for (i=0; i<n; i++) {
	  temp = y[i] ;
	  for (j=0; j<i; j++)
	       temp -= y[j] * matrix[i][j] ;
	  y[i] = temp ;
	  }
     /*
     ** solve DF'z =b
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
#endif
#endif
