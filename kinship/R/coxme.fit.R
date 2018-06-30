# $Id: coxme.fit.s,v 1.15 2003/10/31 19:07:42 Therneau Exp $
#
# Fit with a Gaussian penalty and one or more variance matrices
#
coxme.fit <- function(x, y, strata, offset, init, control,
			weights, ties, rownames, 
			fmat, varlist, ntheta, theta, vinit) {
    time0 <- proc.time()
    eps <- control$eps
    n <-  nrow(y)
    if (is.matrix(x)) nvar <- ncol(x)
    else  if (length(x)==0) nvar <- 0
    else nvar <-1
    
    if (missing(offset) || is.null(offset)) offset <- rep(0.0,n)
    if (missing(weights)|| is.null(weights))weights<- rep(1.0,n)
    else {
	if (any(weights<=0)) stop("Invalid weights, must be >0")
	}

    # Get the list of sort indices, but don't sort the data itself
    if (ncol(y) ==3) {
	if (length(strata) ==0) {
	    sorted <- cbind(order(-y[,2], y[,3]), 
			    order(-y[,1]))
	    newstrat <- n
	    }
	else {
	    sorted <- cbind(order(strata, -y[,2], y[,3]),
			    order(strata, -y[,1]))
	    newstrat  <- cumsum(table(strata))
	    }
	status <- y[,3]
        routines <- 'agfit6b'
        }
    else {
	if (length(strata) ==0) {
	    sorted <- order(-y[,1], y[,2])
	    newstrat <- n
	    }
	else {
	    sorted <- order(strata, -y[,1], y[,2])
	    strata <- (as.numeric(strata))[sorted]
	    newstrat <-  cumsum(table(strata))
	    }
	status <- y[,2]
        routines <- 'coxfit6b'
        }
    
    #
    # The following lines are a failsafe, for the rare case that this
    #   routine is called by a user, and not by coxme()
    # 
    if (!is.list(varlist)) stop("variance matrix list isn't a list!")
    if (is.matrix(fmat) && ncol(fmat >1)) {
	ncluster <- ncol(fmat)
	clnames <- dimnames(fmat)[[2]]
	}
    else ncluster <- 1
    if (ncluster != length(varlist))
        stop("Lengths of variance list and of fmat disagree")
    if (length(ntheta)!= ncluster) stop("Wrong length for ntheta")
#    if (any(theta <0)) stop ("Invalid value of for theta")
    if (length(theta) != sum(ntheta)) stop ("Wrong length for theta")

    #
    # Check out fmat:
    #   each column should be integer
    #   each col must contain numbers from 1 to k for some k
    #   column 1 should contain at least one "1" (the C routine
    if (any(fmat != floor(fmat))) stop("fmat must be integers")
    if (any(fmat <1)) stop("fmat must be >0")
    fmat <- as.matrix(fmat)
    nfrail <- apply(fmat, 2, max)
    temp <- apply(as.matrix(fmat), 2, function(x) length(unique(x)))
    if (any(nfrail != temp)) 
	stop("fmat must be a set of integer indices, with none missing")
    
    #
    # Check the variance array is ok: it must be a list of length
    #  ncluster, and each element must be a list of variance matrices of
    #  dimension equal to fcoef.n
    if (any(sapply(varlist, class) != 'list'))
        stop("varlist must be a list of lists")
    if (any(sapply(varlist, length) != ntheta))
        stop("lengths of varlist disagree with ntheta")
    for (i in 1:ncluster) {
        for (j in varlist[[i]]) { #grab each element of the varlist
            if (!is.matrix(j)) stop("varlist subelements must be matrices")
            if (any(dim(j) != rep(nfrail[i],2))) 
                stop("Invalid dimension in a varlist matrix")
            }
        }
    # End of the redundant checks.  (Although it uses lots of lines of code,
    #  the compute burden is neglible.)

    #
    # For the C code, all of the random effects are glued together into one
    #  large bdsmatrix object "kmat"; 'ikmat', it's inverse, is what is
    #  actually passed to C.
    # kfun() creates kmat; it uses a fixed theta vector of full length.  The
    #  maximization over theta is done by nlmin, and uses a shorter vector.  On
    #  input, any theta=0 is maximized, and all others are fixed.
    #    'itheta' is the vector of "iterated thetas"
    #    'tindex' maps from itheta into theta
    # The total coefficient vector is c(random effects, fixed effects)
    #
    # If fmat has multiple columns, e.g., crossed random effects, then
    #  we only allow one sparse matrix in the final kmat.  For crossed
    #  effects the off-diagonal terms can often be key to reasonable
    #  convergence speed.  And as well -- will someone really have multiple
    #  factors, crossed, where both have a very large number of levels?
    # Whichever factor is first in the list, and only that one, is allowed
    #  to be sparse.  Thus, computation can change just by reordering the
    #  random effects formula.
    #
    # 'findex' is needed by the C routines to impose the 0-sum constraint.
    #   to do this C needs to know which set of columns belongs to each term
    #   Note: this constraint only applies to factors (fmat), not to other 
    #   random effects, and is in fact redundant if the full Hessian matrix
    #   is used for iteration -- the constraint takes care of itself in that
    #   case.
    tindex <- which(theta==0)
	
    # Make sure the first matrix is sparse, and turn all others into
    #  ordinary matrices.  
    # nfrail and tsize are vectors, one element per random term
    nfrail <- unlist(lapply(varlist, function(x) nrow(x[[1]])))

    if (!inherits(varlist[[1]][[1]], 'bdsmatrix')) {
	temp <- varlist[[1]]
	temp <- lapply(temp, function(x) 
		       bdsmatrix(blocks=x, blocksize=nrow(x),
				 dimnames=dimnames(x)))
	varlist[[1]] <- temp
	}
    
    if (ncluster >1) {
	j <- nfrail[1]
	for (i in 2:ncluster) {
	    varlist[[i]] <- lapply(varlist[[i]], as.matrix)
	    fmat[,i] <- fmat[,i] + j
	    j <- j + nfrail[i]
	    }
	}

    nsparse <- rep(0,ncluster)  #number of sparse terms in each matrix
    nsparse[1] <- sum(varlist[[1]][[1]]@blocksize)
    findex <- matrix(0, sum(nfrail), ncluster)
    for (i in 1:ncluster) {
	findex[unique(fmat[,i]),i] <- 1
	}

    kfun <- function(theta, varlist, bsize, rcol) {
        # function to create kmat
        n <- sum((bsize *(bsize+1))/2)
        if (rcol==0)
            kmat <- bdsmatrix(blocksize=bsize, blocks=rep(0., n))
        else kmat<- bdsmatrix(blocksize=bsize, blocks=rep(0., n),
                              rmat=matrix(0., ncol=rcol, nrow=rcol+sum(bsize)))
        
	# Process the first crossed term, which may be sparse
	#  and will be a bdsmatrix
	tlist <- varlist[[1]]
	ntheta <- length(tlist)
	xmat <- tlist[[1]]
	if (length(xmat@rmat) ==0) {
	    krow <-nrow(xmat)
	    kcol <-0
	    }
	else {
	    krow <- nrow(xmat@rmat)
	    kcol <- ncol(xmat@rmat)
	    }
	if (ntheta == 1) {
	    kmat@blocks <- theta[1]*xmat@blocks
	    if (length(xmat@rmat) >0) {
		kmat@rmat[1:krow, 1:kcol] <- theta[1]* xmat@rmat
		}
	    }
	else {
	    #Within a cluster, prior functions have guarranteed that
	    #  all the matrices in the list have exactly the same
	    #  form (both dimensionality and sparse structure)
	    temp <- theta[1]*xmat@blocks
	    rtemp <- theta[1]*xmat@rmat
	    for (j in 2:ntheta) {
		temp <- temp + theta[j]*tlist[[j]]@blocks
		if (length(xmat@rmat)>0) 
		    rtemp <- rtemp + theta[j] * tlist[[j]]@rmat
                
                }
	    kmat@blocks <- temp
	    if (kcol<0) kmat@rmat[1:krow, 1:kcol] <- rtemp
            }
	
	# Now for any other terms (if present)
	#  They will all be ordinary matrices
	ncluster <- length(varlist)
	if (ncluster >1) {
	    ti <- ntheta
	    for (i in 2:ncluster) {
		tlist <- varlist[[i]]
		ntheta <- length(tlist)
		xmat <- tlist[[1]]
		indx <- 1:nrow(xmat)
		if (ntheta == 1) {	 
		    kmat@rmat[krow+indx, kcol+indx] <- theta[ti+1]*xmat
		    }
		else {
		    #Within a cluster, prior functions have guarranteed
		    #  that all the matrices in the list have exactly the 
		    #  same dimensionality 
		    temp <- theta[ti]*xmat@rmat
		    for (j in 2:ntheta) {
			temp <- temp + theta[ti+j-1]*tlist[[j]]
			}
		    kmat@rmat[krow+indx, kcol+indx] <- temp
		    }
		krow <- krow + nrow(xmat)
		kcol <- kcol + nrow(xmat)
		ti <- ti + ntheta
		}
	    }
	kmat
	}
	
    #initialize arguments  -- the way I'm doing this is S-specific.  For
    #  R, use the lines that say "formals"
    #We don't do it for varlist, because I think this makes
    #  a copy of each arg, and that would be a lot of memory.  But doing
    #  these 2 does mean I don't have to do it in the function, nor pass
    #  the args through the nlmin() calling tree
    # kfun[[3]] <- varlist[[1]][[1]]@blocksize 
    # kfun[[4]] <- sum(nfrail-nsparse)
    formals(kfun)[[3]] <- varlist[[1]][[1]]@blocksize
    formals(kfun)[[4]] <- sum(nfrail-nsparse)
    
    # We need to call kfun at least once before the first .C routine
    #  to get the sizes of things.  What should we use as an intial estimate 
    #  of the unknown thetas?
    # Our first use of kmat is to get a starting estimate for the coefs, that
    #  we will use in all the iterations.  From the assumption that
    #  it is easier to shrink random effects than to expand them
    #  (in terms of fast iteration), we first used theta=1 for unknown
    #  random effects, as that is a fairly large value.  However, for some
    #  data sets this is too large a theta, and there are convergence problems!
    #  Now the default for "vinit" in the coxme routine sets the start.
    # Probable long-term solution is addition of Levenberg-Marquardt steps
    #  to the C code, so that it doesn't get lost even for large theta.
    gkmat <- gchol(kfun(ifelse(theta==0,vinit,theta), varlist))
    ikmat <- solve(gkmat)  #inverse of kmat
    if (any(diag(ikmat) <=0)) stop("Random effects variance is not spd")

    #
    # Have C store the data
    # The "-1" in both sorted and kindex is because subscripts start a 0 in C
    # If there are no covariates, make sure that all args have something to
    #   pass.
    # Decide which form of calculation will be done, standard (0) or alternate
    #  (1).  The latter is better with a large number of sparse factors.
    if (nvar==0) {
        x <- 0
        temp.nvar <- 1
        }
    else temp.nvar <- nvar

    rcol <- sum(nfrail-nsparse)
    nevent <- sum(y[,ncol(y)])

    if (is.null(control$sparse.calc)) {
        if ((2*n) > (nevent*(sum(nsparse) - rcol))) control$sparse.calc <- 0
        else control$sparse.calc <- 1 
        }

    ifit <- .C("coxfit6a", 
               as.integer(n),
               as.integer(nvar),
               as.integer(ncol(y)),
               as.double(c(y)),
               as.double(x),
               as.double(offset),
               as.double(weights),
               as.integer(length(newstrat)),
               as.integer(newstrat),
               as.integer(sorted-1),
               as.integer(ncol(fmat)),
               as.integer(fmat-1),
               as.integer(findex),
               as.integer(length(ikmat@blocksize)),
               as.integer(ikmat@blocksize),
               as.integer(rcol),
               means = double(temp.nvar),
               scale = double(temp.nvar),
               as.integer(ties=='efron'),
               as.double(control$toler.chol),
               as.double(control$eps),
               as.integer(control$sparse.calc),
               copy=c(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,T,F,F,F,F), PACKAGE="kinship")
    means   <- ifit$means
    scale   <- ifit$scale

    # Set up initial values for the coefficients
    # The vector of frailty coefficients has those for cluster term 1
    #  followed by those for cluster term 2 followed by ....etc
    # Remember the scaling!
    finit <- rep(0.0, sum(nfrail))
    if (!missing(init) && !is.null(init)) {
	if (length(init) != nvar) {
	    if (length(init) == (nvar+ sum(nfrail))) {
		finit <- init[-(1:nvar)]
		init  <- init[1:nvar]
		}
	    else stop("Wrong length for inital values")
	    }
        init <- c(finit, init*scale)
	}
    else init <- double(nvar+sum(nfrail))

    time1 <- proc.time()
    # The initial fit to the data
    #  If the thetas are fixed, this is the only fit, otherwise it
    #  provides starting estimates for the ms() controlled ones.
    if (ncol(y) ==3)
    fit <- .C("agfit6b",
              iter= as.integer(c(0, control$iter.max)),
              beta = as.double(init),
              loglik = double(2),
              as.double(ikmat@blocks),
              as.double(ikmat@rmat),
              hdet = double(1),
              copy=c(F,T,T,T,F,F,T), PACKAGE="kinship")
    else
    fit <- .C("coxfit6b",
              iter= as.integer(c(0, control$iter.max)),
              beta = as.double(init),
              loglik = double(2),
              as.double(ikmat@blocks),
              as.double(ikmat@rmat),
              hdet = double(1),
              copy=c(F,T,T,T,F,F,T), PACKAGE="kinship")
    log0 <- fit$loglik  #save this for later
    if (fit$iter[2] >1 && fit$iter[2]==control$iter.max)
	    warning("Initial fit did not converge")
    if (fit$iter[2] > control$inner.iter)
           warning(paste("The initial fit took", fit$iter[2], 
                         " steps, which is bigger than the inner.iter ",
                         "paramter.  Consider increasing the latter"))

    if (length(tindex)==0) {
        # all the thetas are fixed -- no more work to do
        ilik <- fit$loglik[2] -
                  .5*(sum(log(diag(gkmat))) + fit$hdet)
        iter <- c(0, fit$iter[2])
        }
    else {
        # Use nlminb to find the optimal parameters
        #    (For R, we need to use optim() instead)
        # For iteration, there are two distinct cases we consider:
        #   fast one -- if there is only one theta, then the
        #     inverse of (theta * kmat) is (1/theta) * kmat^{-1}, and
        #     we only need to invert kmat once for all iterations.
        #     Since shared frailty models with only one random term are
        #     fairly common, it's worthwhile to treat this special.
        #   slow -- otherwise kmat is a sum of terms, and it's inverse
        #     does not necessarily have a simple form
        # At the end of the nlminb calls some of the results of the internal
        #   function are not available to us (i.e., the coefs), so we do
        #   just one more call.
	if (length(vinit)==1) vinit <- rep(vinit, length(tindex))
	else vinit <- vinit[theta==0]

        if (sum(ntheta) ==1) {
            # We are in the special, fast case.
            # Note that if theta is fixed, we never arrive here.
	    # The subtraction of fit0 scales the problem so that the
	    #   relative tolerance criteria of nlminb works better.
            # The "1+" helps when the true solution is near theta=0, beta=0 by
            #   making the min of the function be near -1 rather than 0,
            #   so relative convergence criteria aren't overly fussy
            itemp1 <- ikmat@blocks * vinit  # ikmat for vinit =1
            if (length(ikmat@rmat) >0) itemp2 <- as.vector(ikmat@rmat * vinit)
            else itemp2 <- 0
            logfun1 <- function(x, iblock, rmat, nf, gdet, ofile, fit0, 
                               iter, init) {
                if (x<=0) return(-1)# return a "same as null" fit
                if(ofile=="agfit6b")
                fit <- .C("agfit6b",
                          iter=as.integer(c(iter,iter)),
                          beta = as.double(init),
                          loglik = double(2),
                          as.double(iblock/x),
                          as.double(rmat/x),
                          hdet = double(1),
                          copy=c(F,T,T,T,F,F,T), PACKAGE="kinship")
                else
                fit <- .C("coxfit6b",
                          iter=as.integer(c(iter,iter)),
                          beta = as.double(init),
                          loglik = double(2),
                          as.double(iblock/x),
                          as.double(rmat/x),
                          hdet = double(1),
                          copy=c(F,T,T,T,F,F,T), PACKAGE="kinship")

		kdet <- gdet + nf*log(x)
		ilik <- 1+ fit$loglik[2] -  .5*(kdet + fit$hdet) 
		fit0 -ilik  #ms wants to minimize something, not maximize
		}

            gdet <- sum(log(diag(gkmat))) - nfrail*log(vinit) 
	    mfit <- optim(par=vinit, logfun1, method="L-BFGS-B",
			  lower=control$lower, upper=control$upper,
			  control=list(pgtol=control$toler.ms),    # pgtol as reltol
			  iblock=itemp1, rmat=itemp2, init=fit$beta,
                          nf=nfrail, ofile=routines, gdet=gdet,
                          fit0=fit$loglik[1], iter=control$inner.iter)
            }
        else {
            # slower iteration
            logfun <- function(x, varlist, theta, tindex, init, ofile,
                               fit0, iter, kfun) {
                temp <- theta
                temp[tindex] <- x
                gkmat <- gchol(kfun(temp, varlist))
                ikmat <- solve(gkmat)
		if (any(diag(ikmat) <=0)) { #Not an spd matrix
		    return(0)  # return a "worse than null" fit
		    }
                if (ofile=="agfit6b")
                fit <- .C("agfit6b",
                          iter= as.integer(c(iter,iter)),
                          beta = as.double(init),
                          loglik = double(2),
                          as.double(ikmat@blocks),
                          as.double(ikmat@rmat),
                          hdet = double(1),
                          copy=c(F,T,T,T,F,F,T), PACKAGE="kinship")
                else
                fit <- .C("coxfit6b",
                          iter= as.integer(c(iter,iter)),
                          beta = as.double(init),
                          loglik = double(2),
                          as.double(ikmat@blocks),
                          as.double(ikmat@rmat),
                          hdet = double(1),
                          copy=c(F,T,T,T,F,F,T), PACKAGE="kinship")
		ilik <- 1+ fit$loglik[2] -
                        .5*(sum(log(diag(gkmat))) + fit$hdet)
                        -(ilik - fit0) 
		}
    
	    if (control$simplex <2) {
		mfit <- optim(par=vinit, logfun, method="L-BFGS-B",
			      lower=control$lower, upper=control$upper,
                              control=list(pgtol=control$toler.ms),       # pgtol as reltol
			      varlist=varlist, theta=theta, tindex=tindex,
			      init=fit$beta, ofile=routines, kfun=kfun,
                              fit0=fit$loglik[1], iter=control$inner.iter)
		}
	    else {
		# Use the Nelder-Mead simplex algorithm to constrain the
		#  solution to a reasonable space, before invoking
		#  nlminb
		# Initialize
		nparm <- length(vinit)
		points <- outer(vinit, rep(1,nparm))
		points <- cbind(vinit, vinit + diag(rep(.5, nparm)))
		result <- double(nparm+1)
		for (i in 1:(nparm+1)) 
			result[i] <- logfun(points[,i], varlist, theta,
					    tindex, init=fit$beta)
		ord <- order(result)
		points <- points[,ord]
		result <- result[ord]

		# iterate
		step <- 'reflect'
		for (iter in 1:control$simplex) {
		    if (step=='reflect') {
			# Walk away from the worst (last on list)
			center <- apply(points[,1:nparm], 1,mean) 
			    #center of opposite face
			newx <- 2*center - points[,nparm+1] 
			newy <-  logfun(newx, varlist, theta,
					    tindex, init=fit$beta)
			ii <- iter
			if (newy < result[1]) {
			    step <- 'expand'
			    points <- cbind(newx, points[,1:nparm])
			    result <- c(newy, result[1:nparm])
			    }
			else if (newy > result[nparm]) {
			    step <- 'contract'
			    if (newy > result[nparm+1]) 
				    newx <- points[,nparm+1]
			    }
			else {
			    ord <- sum(newy > result)
			    points <- cbind(points[,1:ord], newx, 
					    points[,(ord+1):nparm])
			    result <- c(result[1:ord], newy, 
					result[(ord+1):nparm])
			    }
			}
		    else if (step=='expand') {
			# Last step worked great: go further
			newx <- 2*newx - center
			newy <-  logfun(newx, varlist, theta,
					    tindex, init=fit$beta)
			
			if (newy < result[1]) { #keep going
			    points[,1] <- newx
			    result[1] <- newy
			    }
			else step <- 'reflect'
			}
	
		    else if (step=='contract') {
			# Last attempt got worse, try to recover
			newx <- .5*(center + newx)
			newy <- logfun(newx, varlist, theta,
					    tindex, init=fit$beta)
			
			if (newy > result[nparm]) {
			    # still the worst, do a shrinkage towards best 
			    #  point so far
			    points[,-1] <- 0.5*(points[,1] + points[,-1])
			    for (i in 1+(1:nparm)) 
				result[i] <- logfun(points[,i], varlist, theta,
					    tindex, init=fit$beta)
			    ord <- order(result)
			    points <- points[,ord]
			    result <- result[ord]
			    }
			else if (newy < result[1]) {
			    points <- cbind(newx, points[,1:nparm])
			    result <- c(newy, result[1:nparm])
			    }
			else {
			    ord <- sum(newy > result)
			    points <- cbind(points[,1:ord], newx, 
					    points[,(ord+1):nparm])
			    result <- c(result[1:ord], newy,
					result[(ord+1):nparm])
			    }
			step <- 'reflect'
			}
                    if (sqrt(var(result)) < control$toler.ms) break
		    }
		vinit <- apply(points,1,mean)  #center of final simplex
		vmin  <- apply(points,1,min)
		vmax  <- apply(points,1,max)
		mfit <- optim(par= vinit, logfun, method="L-BFGS-B",
			      lower=pmax(control$lower, 2*vmin - vinit),
			      upper=pmin(control$upper, 2*vmax - vinit),
			      control=list(pgtol=control$toler.ms),    # pgtol as reltol
			      varlist=varlist, theta=theta, tindex=tindex,
			      init=fit$beta)
		}
		
	    }

        #
        # Now do one more iteration, to get the final coefs
        #
        iter <- c(fit$iter[2], (mfit$counts[1] + mfit$counts[2])*4 + fit$iter[2])
        theta[tindex] <- mfit$par
        gkmat <- gchol(kfun(theta, varlist))
        ikmat <- solve(gkmat,full=T)  #the inverse of K, not of gchol(K)
        if (ncol(y) ==3)
        fit <- .C("agfit6b",
                  iter=as.integer(c(0,control$iter.max)),
                  beta = as.double(fit$beta),
                  loglik = double(2),
                  as.double(ikmat@blocks),
                  as.double(ikmat@rmat),
                  hdet = double(1),
                  copy=c(F,T,T,T,F,F,T), PACKAGE="kinship")
        else
        fit <- .C("coxfit6b",
                  iter=as.integer(c(0,control$iter.max)),
                  beta = as.double(fit$beta),
                  loglik = double(2),
                  as.double(ikmat@blocks),
                  as.double(ikmat@rmat),
                  hdet = double(1),
                  copy=c(F,T,T,T,F,F,T), PACKAGE="kinship")
        ilik <- fit$loglik[2] -
                  .5*(sum(log(diag(gkmat))) + fit$hdet)
        iter[2] <- iter[2] + fit$iter[2]
        }

    #
    # Put together the result vector
    #
    nfrail <- sum(nfrail)
    nsparse<- sum(nsparse)
    nvar2 <- nvar + nfrail -nsparse  #number dense terms
    nvar3 <- nfrail + nvar           #total coefficients
    btot <- length(ikmat@blocks)
    fit3 <- .C("coxfit6c",
               u    = double(nvar3),
               h.b  = double(btot),
               h.r  = double(nvar2*nvar3),
               hi.b = double(btot),
               hi.r = double(nvar2*nvar3),
               hrank= integer(1),
               as.integer(ncol(y)),
               PACKAGE="kinship")
    time2 <- proc.time()
    
    if (nvar2 > 0) {
        rmat1 <- matrix(fit3$h.r, nrow=nvar3)
        rmat2 <- matrix(fit3$hi.r, nrow=nvar3)
        if (nvar >0) {
            #undo the scaling
            # For the ordinary matrix (hinv),  this is essentially
            #   diag(1/scale) %*% hinv %*% diag(1/scale).  But, the scale 
            #   factors are all 1 for the factor variables, so we don't 
            #   have to scale the sparse part of the matrix.
            # For the Cholesky, we have  X=LDL', L lower triangular with 1's
            #   on the diagonal.  The matrix hmat contains L off the diagonal
            #   and D on the diagonal.  We want the cholesky of SXS, where
            #   S=diag(scale).  It turns out that the new diagonal = old
            #   diagonal times S^2, the new L is SL(S-inverse).
            if (nvar==1) { #diag(x) when x is of length 1 causes grief
                rmat1[nvar3,] <- rmat1[nvar3,]/scale
                rmat2[nvar3,] <- rmat2[nvar3,]/scale
                rmat1[,nvar2] <- rmat1[,nvar2]*scale
                rmat2[,nvar2] <- rmat2[,nvar2]/scale
                rmat1[nvar3,nvar2] <- rmat1[nvar3,nvar2]*scale^2
                }
            else {
                temp <- (nfrail+1):nvar3 
                rmat1[temp,] <- (1/scale)*rmat1[temp,] #multiply rows* scale
                rmat2[temp,] <- (1/scale)*rmat2[temp,] 
                temp <- (1+nfrail-nsparse):nvar2        #multiply cols
                rmat1[,temp] <- rmat1[,temp] %*% diag(scale)
                rmat2[,temp] <- rmat2[,temp] %*% diag(1/scale)
                temp <- seq(length=length(scale), to=length(rmat1), by=1+nvar3)
                rmat1[temp] <- rmat1[temp]*(scale^2)    #fix the diagonal
                }
            }
        hmat <- new('gchol.bdsmatrix', .Dim=as.integer(c(nvar3, nvar3)),
                    blocksize=as.integer(ikmat@blocksize), blocks=fit3$h.b,
                    rmat= as.matrix(rmat1), rank=as.integer(fit3$hrank))
        hinv <- bdsmatrix(blocksize=ikmat@blocksize, blocks=fit3$hi.b,
                          rmat=rmat2)
        }
    else {
        hmat <- new('gchol.bdsmatrix', .Dim=as.integer(c(nvar3, nvar3)),
                    blocksize=as.integer(ikmat@blocksize), blocks=fit3$h.b,
                    rmat=as.matrix(numeric(0)), rank=as.integer(fit3$hrank))
        hinv <- bdsmatrix(blocksize=ikmat@blocksize, blocks=fit3$hi.b)
        }

    #
    # Compute the df and REML updates
    #  The latter is returned as an advisory value as a next update to "try"
    #
    tprod <- function(H, pen) {
        # trace of the product of two bdsmatrix objects
        #  The first one, H, might be larger than the penalty one
        # If these were matrices, then the trace = sum(H*t(pen)), using
        #  elementwise multiplication.  
        # Both H and the penalty are symmetric, so no transpose is needed.
        # Since only the lower part of the blocks is kept, we need 2*sum,
        #   but this would count the diagonal elements twice.
        nfrail <- nrow(pen)
        temp1 <- sum(H@blocks * pen@blocks)
        temp2 <- sum(diag(H)[1:nfrail] * diag(pen))
        2*temp1 - temp2
        }

    trace <- tprod(hinv, ikmat)
    df <- (nfrail + nvar) - trace   # formula 5.16 of the book
    fcoef <- fit$beta[1:nfrail]
    b.sigma.b <- sum(fcoef * (as.matrix(ikmat) %*% fcoef))
    idf <- nvar + sum(ntheta)

    time3 <- proc.time()
    timeused <- c(time0[1]+time0[2], time1[1]+time1[2], time2[1]+time2[2],
		  time3[1]+time3[2])
    timeused <- diff(timeused)

    if (nvar > 0) {
	u <- fit3$u
	temp <- seq(to=nvar3, length=length(scale))
	u[temp] <- u[temp] * scale
	list(coefficients=list(fixed=fit$beta[-(1:nfrail)]/scale, 
	                       random=theta),
	     frail=fcoef, penalty=b.sigma.b/2,
	     loglik=c(log0[1], ilik, fit$log[2]), var=hinv,
	     df=c(idf, df), hmat=hmat, iter=iter, control=control,
	     u=u, means=means, scale=scale, timeused=timeused)
	}
    else list(coefficients=list(fixed=NULL, random=theta),
	      frail=fcoef, penalty=b.sigma.b/2,
	      loglik=c(log0[1], ilik, fit$log[2]), var=hinv,
	      df=c(idf, df), hmat=hmat, iter=iter, control=control,
	      u=fit3$u, means=means, scale=scale, timeused=timeused)
    }
