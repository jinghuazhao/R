/*  $Id: coxfit6a.c,v 1.7 2003/10/20 12:40:48 therneau Exp $ */
/* A reentrant version of the Coxfit program, for random effects modeling
**   with reasonable efficiency (I hope).  The important arrays are saved
**   from call to call so as to speed up the process.  The x-matrix itself
**   is the most important of these.
** This version of the routine uses the bdsmatrix library.
**
** coxfit6a: Lots of setup --- save away the necessary data
** coxfit6b: Iterate to convergence for a given theta vector, including the
**               modification terms in the logliklihood.
** coxfit6c: Compute residuals and release the saved memory.
**
**  the input parameters to coxfit6_a are
**
**       nused        :number of people
**       nvar         :number of covariates 
**       ny           : 2 for cox, 3 for (start,stop] data.  Affects sort
**       y[2,n]       :row 1: time of event or censoring for person i
**                    :row 2: status for the ith person    1=dead , 0=censored
**       covar(--,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**                        for the first dimension see comments in the code
**                        section "allocate covar memory"
**       fcol         : number of columns in the fcovar matrix
**       fx(fcol,n)   : compressed memory for factor variables, associated
**                       with sparse matrix portions.
**       findex(fcol,nfrail): for each fcol, a 0/1 vector showing which
**                      coefficients it corresponds to.
**       strata(nstrat):sizes of the strata, cumulative
**       sort 	     :  sort order for the obs, last to first within strata
**                       if ny==3, then there is a second column of data
**       offset(n)    :offset for the linear predictor
**       weights(n)   :case weights
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       tolerch      :tolerance for the Cholesky routines
**       eps          : iteration tolerance
**       method       : Method 0=Breslow, 1=Efron
**       nblock       : nblock, bsize, and rsize give the number of blocks
**       bsize        :  the vector of blocksizes, and the number of columns
**       rsize        :  for the non-sparse (rmat) component of the bdsmatrix
**                       associated with the random effects for the model.
**                       rsize will often be 0, when all of the random effects
**			 are factor variables
**	standard      : which type of calculation: standard or alternate
**
**  returned parameters
**       means        : vector of colum means of X
**       xscale       : vector of values to scale the columns of X (MAD)
** 
**  "nvar3" will be the total number of coefficients, random + fixed.
**   For computational reasons, they are in the order factor/sparse, 
**    then factor/non-sparse, random/not factor/not sparse, then fixed.  
**   The variable nfactor will be the number of factor terms, nsparse
**    the number of sparse terms, and nfrail the total number of
**    penalized (random) ones,  which is the dimension of the penalty matrix.
**    Only factor terms are ever sparse, but not all factor terms need be so.
**   The variable "nfactor" is the number of frailty covariates that are
**    levels of one of the factor variables, and thus encoded in the fx
**    matrix rather than the X matrix, and "nvar" the number of ordinary
**    covariates, which are of course coded in the X matrix.  
**   Factors which are not sparse correspond to levels of a factor that have
**    too large a fraction of the data in that level.  We have found that
**    ignoring the covariate terms for these can make the NR algorithm take
**    bad steps and get lost.
**   The X matrix thus has (nfrail-nfactor) + nvar columns.  Random effects
**    which are not factors will come first.  An example would be random
**    slopes, one slope per group.
**
**  work arrays: these are allocated using Calloc (so they persist
**                until the call of coxfit6_c)
**       mark(n)      : number of tied deaths at a time point.  This is an
**                       integer, but stored as double since it is a frequent
**                       multiplier in expressions with other doubles.
**       wtave(n)     : average weight for those tied deaths
**       score(n)     : vector of Xbeta values, one per subject
**       u(nv)        : first derivative or score vector
**       a(nv)        : work vector for computing u, contains the weighted 
**                         mean at time t
**       a2(nv)       : same as a, but only for the current tied deaths
**       imat         : information matrix, a bsdmatrix 
**         imat.b  = block portion
**                      we engage in some serious trickery here to make 
**                      subscripting easier   -- see comments below where
**                     imat is allocated
**       cmat, cmat2  : correspond to the dense portions of imat, containing
**                      sum(x^2) to compute variance at time t.  
**                      cmat2 = only the tied deaths
**       coef(nv)     : current coefficient estimate
**       oldbeta(nv)  : one-iteration-ago value of coef, needed for 
**                       step-halving
**       temp(nv)     : scratch vector of doubles
**       dlag3, dsum1, dsum2, dsum3: see the Notes file
**       tlist        : integer, scratch vector of length max(mark)
*/
#include <math.h>
#include "coxmeS.h"
#include "coxfit6.h"

/* the next line is just so that I can use "c6.n" instead of "coxfit6.n", etc*/
extern struct coxfit6 c6;

static double **bmatrix(int nblock, int *bsize, int rcol,
			int nfrail, int size);
static double ** cmatrix(int nrow, int ncol);

void coxfit6a(Sint *nused,      Sint *nvar,      Sint *ny,     
	      double *yy,       double *covar2,  double *offset2,
	      double *weights2, Sint  *nstrat,   Sint *strata,    
	      Sint *sorts,      Sint *fcol,      Sint *fx,
	      Sint *findex,
	      Sint *nblock,     Sint *bsize,     Sint *rsize,
	      double *means,    double *xscale,  Sint *method,
	      double *tolerch,  double *eps,     Sint *standard) { 
/*S_EVALUATOR */

    int i,j,k, p1, p2, istrat;
    int nsparse, fsize;
    int nvar2, nvar3;
    int nvar2b;
    int n;
    
    double  temp, temp2;
    double  ndead;
    double  *dptr;
    int *iptr;

    n          = *nused;
    c6.n       = *nused;
    c6.nvar    = *nvar;
    c6.method  = *method;
    c6.tolerch = *tolerch;
    c6.eps     = *eps;
    c6.nblock  = *nblock;
    nsparse =0;  	/* number of random effects coefficients */
    fsize  =0;          /* memory size of the block diagonal matrices */
    if (c6.nblock >0)  { /* no sparse terms is a possibility*/
	for (i=0; i<c6.nblock; i++)
	    nsparse += bsize[i];

	c6.bsize = (int *) Calloc(c6.nblock + 2*nsparse, int);
	c6.bstart= c6.bsize + c6.nblock;
	c6.bstop  = c6.bstart + nsparse;
	k=0; 
	for (i=0; i<c6.nblock; i++) {
	    c6.bsize[i] = bsize[i];
	    fsize   += (bsize[i] * (1+bsize[i])) /2;
	    for (j=0; j<c6.bsize[i]; j++) {
		c6.bstart[k+j] =k;
		c6.bstop[k+j] = k + c6.bsize[i];
		}
	    k += c6.bsize[i];
	    }
	c6.tblock = fsize;
	}
    else c6.tblock  =  0; 
    c6.nsparse = nsparse;
    c6.nfrail =  nsparse + *rsize;

    c6.nfx = *fcol;
    if (c6.nfx >0) {  /* there are factor variables: save fcol matrix */
	c6.fx  = (int *) Calloc(c6.nfx *n, int);
	c6.findex = (double *) Calloc(c6.nfx * c6.nfrail, double);
	j =0;
	for (i=0; i<(n * c6.nfx); i++)  if (fx[i] > j) j = fx[i];
	c6.nfactor = j+1;  /* max number of levels for any factor */
	for (i=0; i<(c6.nfx*c6.nfrail); i++) {
	    c6.findex[i] = findex[i];  /* save findex matrix too */
	    }
	}
    else c6.nfactor = 0;  /* No factor terms! */
    
    nvar3 = *nvar + c6.nfrail;      /* total number of coefficients */
    nvar2 = nvar3 - c6.nfactor;     /* total columns in X matrix */
    nvar2b= nvar3 - nsparse;        /* non-sparse columns in imat */

    /*
    **  Allocate storage for the arrays and vectors
    **  Since they will be used later, sizes are based on what will be
    **    needed with the random effect terms.
    */
    if (nvar2 >0) {  /* save away the X matrix */
	c6.x =  (double **) Calloc(nvar2, double *);
	dptr =  (double *)  Calloc(n*nvar2, double);
	if (*ny==2) {
	    for (j=0; j<nvar2; j++) {  
		c6.x[j] = dptr;
		for (i=0; i<n; i++) {
		    k = sorts[i];    /* store data in sorted order */
		    *dptr++ = covar2[j*n +k];
		    }
		}
	    }
	else {
	    for (j=0; j<nvar2; j++) {
		c6.x[j] = dptr;
		for (i=0; i<n; i++) {
		    *dptr++ = covar2[j*n +i];
		    }
		}
	    }
	}


    c6.imat  = bmatrix(c6.nblock, c6.bsize, nvar2b, nsparse, fsize);
    c6.imatb = c6.imat[0];  
    if (nvar2b > 0) {
	c6.cmat  = cmatrix(nvar2b, nvar3);
	c6.cmat2 = cmatrix(nvar2b, nvar3);
	}

    c6.a = (double *) Calloc(5*nvar3 + 4*n + c6.nfx, double);
    c6.oldbeta = c6.a + nvar3;
    c6.a2      = c6.oldbeta + nvar3;
    c6.wtave   = c6.a2 + nvar3;
    c6.weights = c6.wtave+ n;
    c6.offset  = c6.weights + n;
    c6.u       = c6.offset + n;
    c6.temp    = c6.u     + nvar3;
    c6.stop    = c6.temp  + nvar3 + c6.nfx;
    c6.status  = Calloc(2*n+ *nstrat + nvar3, int);
    c6.mark    = c6.status +n;
    c6.strata  = c6.mark +n;
    c6.itemp   = c6.strata + *nstrat;

    for (i=0; i<*nstrat; i++) c6.strata[i] = strata[i];

    if (*ny==2) {
	/* store the data in sorted order, saving a trivial amount of time*/
	for (i=0; i<n; i++) {
	    j = sorts[i];
	    c6.weights[i] = weights2[j];
	    c6.offset[i]  = offset2[j];
	    c6.status[i]  = yy[n +j];
	    c6.stop[i]    = yy[j];
	    }
	iptr = c6.fx;
	for (k=0; k<*fcol; k++) {
	    j = k*n;
	    for (i=0; i<n; i++) *iptr++ = fx[j + sorts[i]];
	    }
	}
    else {
	c6.sort1 = (int *) Calloc(2*n, int);
	c6.sort2 = c6.sort1 +n;
	c6.start = (double *) Calloc(n, double);
	for (i=0; i<n; i++) {
	    c6.weights[i] = weights2[i];
	    c6.offset[i]  = offset2[i];
	    c6.start[i]   = yy[i];
	    c6.stop[i]    = yy[n+i];
	    c6.status[i]  = yy[2*n +i];
	    c6.sort1[i]   = sorts[i];
	    c6.sort2[i]   = sorts[n+i];
	    }
	for (i=0; i<(*fcol *n); i++) c6.fx[i] = fx[i];
	}

    /*
    **   Mark(i) contains the number of tied deaths at this point,
    **    for the last person of several tied times. It is zero for
    **    all other points.
    **   Wtave contains the average weight for the deaths, at this
    **    time point.
    */
    temp=0;
    istrat=0;
    ndead =0;
    for (i=0; i<n; i++){
	c6.wtave[i] =0;
	c6.mark[i] =0;
	}
    for (i=0; i<n; i++) {
	if (*ny==2) p1=i; else p1=sorts[i];  /* death order */
	if (c6.status[p1]==1) {
	    ndead++;
	    temp += c6.weights[p1];
	    }
	if ((i+1) == c6.strata[istrat]) {
	    istrat++;
	    c6.mark[p1] = ndead;
	    c6.wtave[p1]= temp/ndead;
	    ndead =0;
	    temp  =0;
	    }
	else {  /* prior statement is always true when i== n-1 */
	    if (*ny==2) p2 = i+1; else p2= sorts[i+1];
	    if (c6.stop[p1] != c6.stop[p2]) {
		c6.mark[p1] = ndead;
		if (ndead >0) c6.wtave[p1]= temp/ndead;  /*avoid a 0/0 */
		ndead =0;
		temp  =0;
	        }
	    }
	}

    /*
    ** Subtract the mean from each x and scale it, as this makes the
    **   regression much more stable.  (And with these models, we sometimes
    **   need all the stability we can get).
    ** Penalized terms, however, must be left in peace.  Changing the size
    **   of their associated coefficients would mess up the penalty.
    */
    k = nvar2 - c6.nvar;   /* number of penalized terms in X matrix */
    for (j=k; j<nvar2; j++) {
	temp=0; temp2=0;
	for (i=0; i<n; i++) temp += c6.x[j][i];
	temp /= n;
/*	temp =0;     debug*/
	means[j-k] = temp;
	for (i=0; i<n; i++) {
	    c6.x[j][i] -= temp;
	    temp2 += fabs(c6.x[j][i]);
	    }
	temp2 /= n;
/*	temp2 =1;   debug */
	xscale[j-k] = temp2;
	for (i=0; i<n; i++) c6.x[j][i] /= temp2;
	}

    /*
    ** Space for the later dsum tricks
    */
    c6.calc2 = *standard;
    if (c6.calc2 ==1) {
	j=0;
	for (i=0; i<n; i++)
	    if (c6.mark[i] >j) j=c6.mark[i];   /* max number of tied deaths */
	c6.tlist = (int *) Calloc(j, int);
	c6.dlag1 = (double *) Calloc(nvar3, double);
	c6.dsum3 = c6.dlag1 + nsparse;

	/* 
	** This is dsum2 of the note (first nfactor rows/cols) and
	**   also dsum3 (next non-factor rows)
	*/
	c6.dlag2  = bmatrix(c6.nblock, c6.bsize, nvar3-nsparse, 
	                     nsparse, fsize);
	}
    }

/*
** Allocate space for a ragged array
*/
static double ** cmatrix(int nrow, int ncol) {
    double **mat;
    double *ptr;
    int i;

    mat = (double **) Calloc(nrow, double *);
    ptr = (double *)  Calloc(nrow*ncol, double);
    for (i=0; i<nrow; i++) {
	mat[i] = ptr;
	ptr += ncol;
	}
    return(mat);
    }

/*
** Allocate space for a bdsmatrix 
**  and cleverly arrange the indexing so that it looks like a single
**  ragged array.  
**  For the sparse portion, x[i[][j] is only legal for j>=i, j within the block
**    (as though we only saved the upper portion of the ragged part)
**  For the non-sparse portion, it looks ordinary.
**
**  nblock, bsize, rcol: parameters of the bdsmatrix object
**            number of blocks, block sizes, #columns of the dense portion
**  nfrail = # sparse coefs = sum(bsize)
**  size   = amount of memory for the sparse portion = sum[bsize*(bsize+1)/2]
** 
**  Yes, the last 2 parameters are redundant, but the calling routine happens
**    to have them at hand, so why not.
*/
static double **bmatrix(int nblock, int *bsize, int rcol,
			int nfrail, int size) {
    int i,j,k;
    double **pointer;
    double *temp;

    pointer = (double **) Calloc(nfrail+rcol, double *);  /* ragged array */

    i = size + rcol*(rcol+nfrail);
    temp = (double *) Calloc(i, double);  /* the data itself */

    /* index the sparse portion */
    k=0;
    for (i=0; i<nblock; i++) {
	for (j=bsize[i]; j>0; j--) {
	    pointer[k] = temp -k;
	    k++;
	    temp += j;
	    }
	}

    /* index the dense portion */
    for (i=nfrail; i< (rcol+nfrail); i++) {
	pointer[i] =temp;
	temp += (rcol + nfrail);
	}

    return(pointer);
    }
