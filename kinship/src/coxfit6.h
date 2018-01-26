/* $Id: coxfit6.h,v 1.4 2003/04/07 18:31:39 therneau Exp $  */
/*
** The common variables for the coxfit6 routines.  These are saved in
**  a single external list, in order to not use up external
**  symbols in the Splus core.  The alternative was to make them all
**  "static external", and combine all the routines into a single .c 
**  source file, but that file was just getting too long to edit
**  easily.  (Not an problem for emacs -- a problem for me!)
** They are needed by all 4 of the routines.  Saving them from call to
**   call is faster than passing everything in from Splus each time.
*/

#ifndef COXFIT6_H
#define COXFIT6_H
struct coxfit6 {
    double eps;         /* convergence criteria */
    double tolerch;     /* tolerance for the cholesky routines */
    double *stop;       /* vector of death times */
    double *start;      /* vector of start times */
    double **x;    	/* the matrix of covariates */
    double *weights;    /* case weights */
    double *offset;     /* offset terms, if present */
    double *a;		/* mean covariate at the current time */
    double *a2;         /* bit of 'a' due to deaths, for Efron approx */
    double *u;          /* first derivative vector, of length p */
    double  *imatb;     /*    block portion of imat */
    double *wtave;      /* average weight at this death time, for Efron */
    double *oldbeta;    /* lagged beta coef, for step-halving */
    double *temp;       /* scratch vector */
    double *findex;     /* which coefficients each sparse term refers to */
                        /* these are integers, but kept as double to be used
			   as an argument to bdsmatrix_prod2 */
    double *dlag1;      /* lagged denominators */
    double *dsum3;      /*  " " */
    double **cmat;	/* sum of squares & cross products, dense part only*/
    double **cmat2;     /* bit of cmat due to deaths, for Efron approx */
    double **imat;      /* information matrix */
    double **dlag2;      /* lagged denominators */
    int    *itemp;      /* scratch vector */
    int    *status;     /* 0/1 status variable */
    int    *mark;       /* number of tied deaths at a time point */
    int    *sort1;      /* only used for (start,stop] data */
    int    *sort2;      /*         ditto */
    int    *strata;     /* contains the cumulative count of number per strata*/
    int    *fx;         /* the factor variables corresponding to sparse terms*/
    int    *bsize;      /* block size structure for imat, cmat, cmat2 */
    int    *tlist;      /* list of factors tied at a given time */
    int    *bstart;     /* index of the start of each block, repeated */
    int    *bstop;      /* index of the end of each block, repeated */
    int    n;           /* number of subjects */
    int    nvar;        /* number of fixed effect covariates */
    int    nfrail;      /* total number of random effects */
    int    nsparse;     /* number of the nfrail which are sparse */
    int    nfactor;     /* number of coefficients represented by fx matrix */
    int    nblock;      /* number of blocks in imat.b */
    int    tblock;      /* total number of elements in the sparse portion */
    int    nfx;         /* number of columns in fx, usually = 1 */
    int    method;      /* 0= Breslow, 1=Efron */
    int    calc2;       /* 0=usual, 1=alternate summation for sparse terms */
    };
#endif
