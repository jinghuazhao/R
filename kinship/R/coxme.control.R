# $Id: coxme.control.s,v 1.5 2003/08/21 20:30:53 Therneau Exp $
#
# Gather all of the control parameters for coxme into one spot
#
coxme.control <- function(eps=1e-5, 
                          toler.chol = .Machine$double.eps ^ .75, 
                          toler.ms = 1e-2,
			  inner.iter=4,
			  iter.max=10,
			  simplex =0 ,
			  lower = 0,
			  upper = Inf,
			  sparse.calc=NULL) {
    if (iter.max <0) stop("Invalid value for iterations")
    if (inner.iter<1) stop("Invalid value for inner iterations")
    if (eps <=0) stop ("Invalid convergence criteria")
    if (eps <= toler.chol) 
	    warning("For numerical accuracy, tolerance should be < eps")
    if (toler.ms <=eps)
        warning("For numerical accuracy, eps should be < toler.ms")
    if (!is.null(sparse.calc)) {
        if (sparse.calc !=0 && sparse.calc !=1)
            stop("Invalid value for sparse.calc option")
        }
    if (simplex <0 || simplex != floor(simplex))
	    stop("Number of simplex iterations must be a positive integer")
    list(eps=eps, toler.chol=toler.chol, iter.max=iter.max,
	 inner.iter=inner.iter, sparse.calc=sparse.calc,
	 simplex=simplex, lower=lower, upper=upper, toler.ms=toler.ms)
    }
