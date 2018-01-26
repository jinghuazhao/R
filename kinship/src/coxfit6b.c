/* $Id: coxfit6b.c,v 1.8 2003/10/20 12:40:49 therneau Exp $ */
/*
** This call is used for iteration
**   It assumes a fixed penalty matrix, stored as a bdsmatrix of course
**   
**    Input
** iter         : a vector containing (min, max) number of iterations
** beta         : vector of coefficients, random first, then others
**              :    On input contains starting values
** pmatb        : block diagonal portion of the penalty matrix
** pmatr        : dense portion of the penalty matrix
**
**    Output
** maxiter      : number of iterations actually done
** beta         : ending values
** loglik       : loglik for the starting values, and for the final ones
** hdet         : log(determinant) of the sparse portion of H
*/
#include "coxmeS.h"
#include "coxfit6.h"
#include <math.h>
#include <stdio.h>
#include "kinproto.h"

/* the next line is just so that I can use "c6.n" instead of "coxfit6.n", etc*/
extern struct coxfit6 c6;
static void update(int j, int upper);
static double dsum1, dsum2;
static int nvar3;

void coxfit6b(Sint *maxiter,  double *beta,
	      double *loglik, double *pmatb,  double *pmatr,
	      double *hdet) {
    int i,j,k, l, p;
    int ii, istrat;
    int     iter;
    int    nvar, nvar2, nvar2b;
    int    nfrail, ns, nfac;
    int    nfns;    /* number of factors that are not sparse */
    int    halving;

    double  denom=0, zbeta, risk;
    double  temp, temp2;
    double  newlik=0;
    double  d2, efron_wt=0;
    double  ndead;
    double  *dptr;
    double  *psum;
    int     dohalf =0;
    int     ntie, ntie2, dup1, dup2;

    nfrail = c6.nfrail;    /* number of penalized coefficients */
    nvar   = c6.nvar;      /* number of unpenalized coefficients */
    ns     = c6.nsparse;   /* number of factor levels that are sparse */
    nfac   = c6.nfactor;   /* number of factor levels (penalized) */
    nvar2  = nvar + (nfrail - nfac);  /* number of cols of X */
    nvar3  = nvar + nfrail; /* total number of coefficients */
    nvar2b = nvar3 - ns;    /* total number non-sparse terms */
    nfns   = nfac - ns;     /* number of factor levels that are NOT sparse */
    for (i=0; i<nvar3; i++) c6.oldbeta[i] = beta[i];

    /*
    ** Compute the sums of the penalty matrix, used for recentering
    **  the frailty coefficients in a sparse model
    ** Save the result in "psum"
    */
    psum = c6.temp + nvar3;
    for (i=0; i<c6.nfx; i++) {
	/* 
	** get the product of the penalty matrix with the 0/1 indicator
	**  variable that represents the ith sparse factor variable.
	** recentering the sparse factors is important for the NR; non-sparse
	**  terms don't have the problem.
	*/
	bdsmatrix_prod2(c6.nblock, c6.bsize, nfrail, pmatb, pmatr,
			c6.findex + i*nfrail, c6.temp, c6.itemp); 
	temp =0;
	for (j=0; j<nfrail; j++) temp += c6.temp[j];
	psum[i] = temp;
	}

    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (iter=0; iter<= maxiter[1]; iter++) {
	/*
	** Initialize things to the value of the penalty,
	**  using c6.temp as a temporary vector for column sums
	** First the information matrix
        */
	for (i=0; i<c6.tblock; i++) 
	    c6.imatb[i] = pmatb[i];
	dptr = pmatr;
	for (i=ns; i<nfrail; i++) {
	    /* dense rows of penalty */
	    for (j=0; j<nfrail; j++) c6.imat[i][j] = *dptr++;
	    for (j=nfrail; j<nvar3; j++) c6.imat[i][j] =0;
	    }
	for (i=nfrail; i<nvar3; i++) {
	    /* unpenalized part */
	    for (j=0; j<nvar3; j++) c6.imat[i][j] =0;
	    }

	/* form the product of penalty times beta, save in c6.temp */
	bdsmatrix_prod2(c6.nblock, c6.bsize, nfrail, pmatb,
			pmatr, beta, c6.temp, c6.itemp);
	
	/* u and penalized loglik */
	temp =0;
	for (i=0; i<nfrail; i++) {
	    c6.u[i] = -c6.temp[i];
	    temp += c6.temp[i]*beta[i];
	    }
	newlik = -temp/2;  /* -(1/2) b' \sigma^{-1}b */
	for (i=nfrail; i<nvar3; i++) c6.u[i] =0;
	
	/*
	** Now loop through the data, and compute the loglik and the
	**  first and second derivatives (u and imat).
	** One would think that we could be faster computing only
	**  the loglik, then deciding if another iter was needed,
	**  then only computing the first and second derivs if we needed
	**  them.  But the final iter needs the determinant of imat,
	**  so the only savings would come when we are step halving.
	**  This, we expect, happens very rarely.
	*/ 
	istrat=0;
	ntie =0;  ntie2=0;
	for (p=0; p<c6.n; p++) {  /* p = person */
	    if (p==0 || p==c6.strata[istrat]) {
		if (p>0) {
		    istrat++;
		    if (c6.calc2==1) {
			for (j=0; j<ns; j++) update(j, 0);
			}
		    }

		efron_wt =0;
		denom = 0;
		for (i=0; i<nvar3; i++) {
		    c6.a[i] = 0;
		    c6.a2[i]=0 ;
		    for (j=0; j<nvar2b; j++) {
			c6.cmat[j][i] = 0;
			c6.cmat2[j][i]= 0;
                        }
		    }
		if (c6.calc2 ==1) {
		    dsum1 = 0; dsum2 =0;
		    for (i=0; i<nvar2b; i++) c6.dsum3[i] =0;
		    for (i=0; i<ns; i++) c6.dlag1[i] =0;
		    for (i=0; i<c6.tblock; i++) c6.dlag2[0][i] =0; 
		    for (i=ns; i<nvar3; i++) {
			for (j=0; j<=i; j++) c6.dlag2[i][j] =0;
			}
		    }
		}

	    /*
	    ** Form the linear predictor zbeta, and the risk score
	    */
	    zbeta = c6.offset[p];
	    for (i=0; i<c6.nfx; i++) {
		j = c6.fx[p + i*c6.n];  /* level of covariate i */
		zbeta = zbeta + beta[j];
		}
	    for (i=0; i<nvar2; i++)
		zbeta += beta[i+nfac]* c6.x[i][p];
	    risk = exp(zbeta) * c6.weights[p];
	    denom += risk;

	    /*
	    ** Compute the a vector (sums) and 
	    **  the c matrix (sums of squares and cross products)
	    ** There are no cross products between sparse factors, as
	    **  no space was left for them in the matrices.  Within
	    **  a factor, each row of data has a single "1", so cmat is
	    **  not needed.
	    ** It's an oddity of the indexing, due to how bdsmatrix
	    **  objects were first defined, that I appear to be
	    **  using the above diagonal part for the sparse, and
	    **  the below diagonal part for the dense part of imat
	    */
	    for (i=0; i< c6.nfx; i++) {
		j = c6.fx[p + i*c6.n];   /* jth covariate is = to 1 */
		/* first, update u and imat based on the OLD a[j] */
		if (c6.calc2==1 && j<ns) update(j,1);
		    
		/* Now update a and cmat */
		c6.a[j] += risk;
		if (j>=ns) c6.cmat[j-ns][j] += risk;
		for (k=i+1; k<c6.nfx; k++) /* crossed factors */ 
		    c6.cmat[c6.fx[p+ k*c6.n] -ns][j] += risk;
		for (k=0; k<nvar2; k++) /* covariates */ 
		    c6.cmat[k+nfns][j] += risk*c6.x[k][p];
		}

	    for (i=0; i<nvar2; i++) {   /* non-factor variables */
		c6.a[i+nfac] += risk * c6.x[i][p];
		for (j=0; j<=i; j++)
		    c6.cmat[i+nfns][j+nfac] += risk*c6.x[i][p]*c6.x[j][p];
		}

	    /*
	    ** Extra terms for the deaths
	    */
	    if (c6.status[p]==1) {
		newlik += c6.weights[p] *zbeta;
		efron_wt += risk;
		for (i=0; i< c6.nfx; i++) {
		    j = c6.fx[p + i*c6.n];   /* jth covariate is = to 1 */
		    c6.u[j] += c6.weights[p];
		    c6.a2[j] += risk;	
		    if (j>=ns) c6.cmat2[j-ns][j] += risk;
		    for (k=i+1; k<c6.nfx; k++) /* crossed factors */
			c6.cmat2[c6.fx[p+ k*c6.n] -ns][j] += risk;
		    for (k=0; k<nvar2; k++)
			c6.cmat2[k+nfns][j] += risk*c6.x[k][p];

		    /* 
		    ** Add this factor variable to the list of "it changed
		    **  values at this death time"
		    ** The Efron imat calculations will require that we
		    **  update all the rows for this factor variable
		    **  It's useful to keep both all rows, and all unique
		    **  blocks; ntie2 < ntie if two rows from the same block
		    **  occur.  The first ntie are from unique blocks.
		    */
		    if (j<ns && c6.calc2==1 && c6.method==1) {
			dup1=0; dup2=0;
			for (k=0; k<ntie; k++) {
			    if (c6.tlist[k] ==j) {
				dup1 =1;  /* exact duplicate */
				break;
				}
			    if ((c6.bstart[c6.tlist[k]] <= j) &&
				(c6.bstop[c6.tlist[k]]  >  j))  dup2=1;
			    }
			if (dup1==0) {
			    if (dup2==1) c6.tlist[ntie++] =j;
			    else {
				c6.tlist[ntie++] = c6.tlist[ntie2]; 
				c6.tlist[ntie2++] = j;
				}
			    }
			}
		    }

		for (i=0; i<nvar2; i++) {  /* non-factor terms */
		    c6.u[i+nfac] += c6.weights[p] *c6.x[i][p];
		    c6.a2[i+nfac] +=  risk*c6.x[i][p];
		    for (j=0; j<=i; j++)
			c6.cmat2[i+nfns][j+nfac] += risk*c6.x[i][p]*c6.x[j][p];
   		    }
		}

	    if (c6.mark[p] >0) {  /* once per unique death time */
		/* use cmat, cmat2, a, and a2 to update u and imat */
		ndead = c6.mark[p];
		if (c6.method==0 || ndead==1)  {
		    /*
		    ** Breslow approx -- we can ignore a2 and cmat2
		    */
		    temp = c6.wtave[p] * ndead;
		    newlik -= temp *log(denom);

		    if (c6.calc2==1) {
			ii = ns;
			dsum1 += temp/denom;
			dsum2 += temp/(denom * denom);
		
			for (i=ns; i<nvar3; i++) {  /* update u */
			    c6.temp[i] = c6.a[i]/ denom;
			    c6.u[i] -= temp *c6.temp[i];
			    c6.dsum3[i-ns] += c6.temp[i] * temp/denom;
			    }
			}
		    else {
			ii =0;
			for (i=0; i<nvar3; i++) {
			    c6.temp[i] = c6.a[i]/ denom;
			    c6.u[i] -= temp *c6.temp[i];
			    }
			for (i=0; i<ns; i++) {
			    c6.imat[i][i] += temp * c6.temp[i];
			    for (j=i; j<c6.bstop[i]; j++)
				c6.imat[i][j] -= temp*c6.temp[j]*c6.temp[i];
			    }
			}
		    /* non-sparse variables - factor or continuous*/
		    for (i=0; i<nvar2b; i++) {
			k = i+ns;  /* i=row number in cmat, k= row in imat */
			for (j=ii; j<=k; j++) 
			    c6.imat[k][j] +=  temp *(
				c6.cmat[i][j] /denom - c6.temp[k]*c6.temp[j]);
			}
		    }
		
		else {
		    /* 
		    ** Do the Efron approx 
		    ** In this case we update the non-sparse, along with
		    **  those sparse factors which got changed at this death
		    **  time (those with a2 != 0)
		    */
		    for (temp2=0; temp2<ndead; temp2++) {
			temp = temp2/ ndead;
			d2= denom - temp*efron_wt;
			newlik -= c6.wtave[p] *log(d2);

			if (c6.calc2==1) {
			    ii = ns;
			    dsum1 += c6.wtave[p]/d2;
			    dsum2 += c6.wtave[p]/(d2*d2);

			    for (i=ns; i<nvar3; i++) {  /* update u */
				c6.temp[i] = (c6.a[i] - temp*c6.a2[i])/d2;
				c6.u[i] -= c6.wtave[p] *c6.temp[i];
				c6.dsum3[i-ns] += c6.temp[i] * c6.wtave[p]/d2;
			        }

			    for (i=0; i<ntie2; i++) {
				for (j=c6.bstart[c6.tlist[i]]; 
				     j<c6.bstop[c6.tlist[i]]; j++) {
				    c6.temp[j] = (c6.a[j] - temp*c6.a2[j])/d2;
				    }
				}
			    for (i=0; i<ntie; i++) {
				j = c6.tlist[i];
				c6.u[j] -= c6.wtave[p] *c6.temp[j];
				c6.imat[j][j] +=  c6.wtave[p] *c6.temp[j];
				/*
				** Update imat[k,j] for all k, unless k<j 
				**  and k is also on tlist (no double updates!)
				*/
				for (k=c6.bstart[j]; k<j; k++) {
				    dup1=0;
				    for (l=0; l<ntie; l++)
					if (c6.tlist[l]==k) dup1=1;
				    if (dup1==0) 
					c6.imat[k][j] -= c6.temp[j]*c6.temp[k]
					               * c6.wtave[p];
				    }
				for (k=j; k<c6.bstop[j]; k++) 
				    c6.imat[j][k] -= c6.temp[j]*c6.temp[k]
					               * c6.wtave[p];
				for (k=ns; k<nvar3; k++) 
				    c6.imat[k][j] += c6.wtave[p]* (
					     (c6.cmat[k-ns][j] - 
					         temp*c6.cmat2[k-ns][j])/d2 -
                                              c6.temp[k]*c6.temp[j]);
				}
			    }
			else {
			    ii=0;
			    for (i=0; i<nvar3; i++) {
				c6.temp[i] = (c6.a[i] - temp*c6.a2[i])/d2;
				c6.u[i] -= c6.wtave[p] *c6.temp[i];
				}
			    for (i=0; i<ns; i++) {	
				c6.imat[i][i] += c6.wtave[p] *c6.temp[i];
				for (j=i; j< c6.bstop[i]; j++)
				    c6.imat[i][j] -= c6.wtave[p] *
					              c6.temp[i] * c6.temp[j];
				}
			    }

			/*
			** Update the non-sparse part of imat
			*/
			for (i=0; i<nvar2b; i++) {
			    k = i+ns;  
			    for (j=ii; j<=k; j++) {
				c6.imat[k][j] +=  c6.wtave[p]*(
				    (c6.cmat[i][j] - temp*c6.cmat2[i][j]) /d2 -
                                          c6.temp[k]*c6.temp[j]);
			        }
			    }
		        }
		    
		    if (c6.calc2 == 1) { /* update denominators */
			for (i=0; i<ntie; i++) {
			    j = c6.tlist[i];
			    c6.dlag1[j] = dsum1;
			    for (k=c6.bstart[j]; k<j; k++)
				    c6.dlag2[k][j] = dsum2;
			    for (k=j; k<c6.bstop[j]; k++)
				    c6.dlag2[j][k] = dsum2;
			    for (k=ns; k <nvar3; k++)
				    c6.dlag2[k][j] = c6.dsum3[k-ns];
			    }
			}
		    } /* end of Efron loop */
			 
		/* rezero temps */
		efron_wt =0;
		ntie =0; ntie2=0;
		for (i=0; i<nvar3; i++) {
		    c6.a2[i]=0;
		    for (j=0; j<nvar2b; j++)  c6.cmat2[j][i]=0;
		    }
		}   /* matches "if (mark[p] >0)"  */
	    } /* end  of accumulation loop  */

	/* 
	** Finish up any deferred sums for sparse terms
	*/
	if (c6.calc2==1) {
	    for (j=0; j<ns; j++) update(j, 0);
	    }
	    
	/* 
	**   Am I done?
	** Note, when doing "minimum" iterations, don't allow step halving at
	**  the tail end of the iterations.  
	*/
	if (iter==0) loglik[0] = newlik;
	if (iter>0 && newlik < loglik[1] && 
	           fabs(1-(loglik[1]/newlik)) > c6.eps)  {  
	    /*it is not converging ! */
	    halving =1;
	    dohalf = iter;
	    for (i=0; i<nvar3; i++)
		beta[i] = (c6.oldbeta[i] + beta[i]) /2; 
	    continue;
	    }

	halving =0;
	cholesky4(&(c6.imat[ns]), nvar3, c6.nblock, 
				  c6.bsize,  c6.imatb, c6.tolerch);

	if (iter >= maxiter[0] && fabs(1-(loglik[1]/newlik)) <= c6.eps) break;
	loglik[1] = newlik;
	if (iter < maxiter[1]) {
	    chsolve4(&(c6.imat[ns]), nvar3, c6.nblock, 
				  c6.bsize,  c6.imatb, c6.u, 0);
	    for (i=0; i<nvar3; i++) {
		c6.oldbeta[i] = beta[i];
		beta[i] += c6.u[i];
		}
	    
	    /*
	    ** Impose the constraint, mean frailty for any factor term
	    **  is 0.  If the problem is not sparse, this happens
	    **  automatically with the NR iteration.  If it is sparse,
	    **  this helps efficiency of the maximizer.
	    ** c6.a is used as a scratch variable.  Each call to prod2
	    **   is penalty matrix[rows/cols for this factor] %*% beta[this
	    **   factor].  The divisor is penalty[same] %*% rep(1, nrows)
	    */
	    for (i=0; i<c6.nfx; i++) {
		for (j=0; j<nfrail; j++)
		    c6.a[j] = beta[j] * c6.findex[j + i*nfrail]; 
		bdsmatrix_prod2(c6.nblock, c6.bsize, nfrail, pmatb, pmatr,
				c6.a, c6.temp, c6.itemp);
		temp =0;
		for (j=0; j<nfrail; j++) temp += c6.temp[j];
		temp /= psum[i];  /* the mean */
		for (j=0; j<nfrail; j++) {
		    if (c6.findex[j + i*nfrail] ==1) beta[j] -= temp;
		    }
		}
	    }
	}   /* return for another iteration */

    temp =0;
    for (i=0; i<nfrail; i++) temp += log(c6.imat[i][i]);
    *hdet = temp;
    if (maxiter[1] > iter) maxiter[1] = iter;
    loglik[1] = newlik;
    return;
    }

static void update(int j, int upper) {
    double temp;
    int k;

    if (dsum1 == c6.dlag1[j]) return;  /* all the terms below just add a zero*/

    if (c6.a[j] > 0) {  /* for 1 factor/obs, this saves half the evals! */
	temp = c6.a[j] * (dsum1 - c6.dlag1[j]);
	c6.u[j] -= temp;
	c6.imat[j][j] += temp; 

	if (upper==1) {
	    for (k=c6.bstart[j]; k<j; k++) 
		c6.imat[k][j] -= c6.a[j]*c6.a[k] * (dsum2 - c6.dlag2[k][j]);
	    }

	for (k=j; k<c6.bstop[j]; k++) 
	    c6.imat[j][k] -= c6.a[j]*c6.a[k] * (dsum2 - c6.dlag2[j][k]);
	for (k=c6.nsparse; k<nvar3; k++) 
	    c6.imat[k][j] +=  c6.cmat[k-c6.nsparse][j]*(dsum1 - c6.dlag1[j]) -
	                   c6.a[j] *(c6.dsum3[k-c6.nsparse] - c6.dlag2[k][j]);
	}

    c6.dlag1[j] = dsum1;
    if (upper==1) for (k=c6.bstart[j]; k<j; k++) c6.dlag2[k][j] = dsum2;
    for (k=j; k<c6.bstop[j]; k++) c6.dlag2[j][k] = dsum2;
    for (k=c6.nsparse; k<nvar3; k++)   c6.dlag2[k][j]=c6.dsum3[k-c6.nsparse];
    }
