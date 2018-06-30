# $Id: coxme.s,v 1.11 2003/08/22 14:40:37 Therneau Exp $
# This program, CoxMixedEffects model or coxme, is a first version of what
#  I hope will eventually be a more general mixed effects model.
#  At the moment, it supports one more argument: which is a variance
#  structure.

coxme <- function(fixed=formula(data), data=parent.frame(), random, 
	weights, subset, na.action, init, 
	control, ties= c("efron", "breslow", "exact"),
	singular.ok =T, varlist=NULL, variance, vinit=.2, sparse=c(50,.02),
	rescale=T, pdcheck=T, x=F, y=T, shortlabel=T, ...) {
    time0 <- proc.time()
    ties <- match.arg(ties)
    call <- match.call()
    m <- match.call(expand.dots=FALSE)
    temp <- c("", "data", "weights", "subset", "na.action")
    m <- m[ match(temp, names(m), nomatch=0)]

    if (missing(variance)) theta <- NULL
    else  theta <- variance  #We always liked "theta" better as a name, but it
                             # didn't seem as obvious to the user community

    # Be sneaky, and paste the extra variables onto my formula
    # This causes the model frame m to have all the material that it needs
    #  for both fixed and random components, and to deal correctly with missing
    # The first line for temp.fixed gives the character form of the original
    #  formula.  The paste is needed in case it is too long, in which case
    #  deparse returns multiple character strings, one for each line that would
    #  be printed on the terminal.
    gvars <- all.vars(random)
    temp.fixed <- paste(deparse(as.vector(fixed)), collapse='')
    temp.fixed <- paste(temp.fixed, paste(gvars, collapse='+'), sep='+')
    temp.fixed <- as.formula(temp.fixed)
    
    m$formula <- temp.fixed
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    Terms <- terms(fixed)
    if (missing(control)) control <- coxme.control()
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    weights <- model.extract(m, 'weights')
    offset<- attr(Terms, "offset")
    tt <- length(offset)
    n <- nrow(Y)

    # idiot proofing: if more than one offset statement, just
    #   add the offset terms together
    offset <- if(tt == 0)
		    rep(0, nrow(Y))
	      else if(tt == 1)
		      m[[offset]]
	      else {
		    ff <- m[[offset[1]]]
		    for(i in 2:tt)
			    ff <- ff + m[[offset[i]]]
		    ff
		    }

    special <- c("strata", "cluster")
    Terms <- terms(fixed, special)
    attr(Terms,"intercept")<- 1  #Cox model always has \Lambda_0
    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
        # people keep wanting to mix GEE and random effects...
        stop ("A cluster() statement is invalid in coxme")
	}
    if (length(strats)) {
	temp <- untangle.specials(Terms, 'strata', 1)
	dropx <- c(dropx, temp$terms)
	if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
	else strata.keep <- strata(m[,temp$vars], shortlabel=T)
	strats <- as.numeric(strata.keep)
	}

    if (length(dropx)) X <- model.matrix(Terms[-dropx], m)[,-1,drop=FALSE]
    else               X <- model.matrix(Terms, m)[,-1,drop=FALSE]
	
    type <- attr(Y, "type")
    if (type!='right' && type!='counting')
	stop(paste("Cox model doesn't support \"", type,
			  "\" survival data", sep=''))

    if (missing(init)) init <- NULL

    # Check for penalized terms
    pterms <- sapply(m, inherits, 'coxph.penalty')
    if (any(pterms)) {
	stop("You cannot have, as yet, penalized terms in coxme")
	}

    if (ties=='exact') 
	    stop("Exact (or discrete) method for breaking ties is not supported for coxme")
    if( ties!="breslow" && ties !="efron")
	    stop(paste ("Unknown option for ties", ties))
    
    # Check that the "sparse" option was legally set
    if (!missing(sparse)) {
        if (is.logical(sparse)) {
            if (sparse) stop("A value of sparse=T is not valid")
            sparse <- c(Inf, 0) # no sparse terms
            }
        else {
            if (!is.numeric(sparse)) stop("Invalid sparse option")
            if (length(sparse) == 1) {
                if (sparse >1) sparse <- c(sparse, 1)
                else if (sparse >=0)  sparse <- c(Inf, sparse)
                else stop("Invalid value for sparse option")
                }
            else if (length(sparse)==2) {
                if (sparse[1] <2) 
                    stop("Invalid value for first element of sparse option")
                if (sparse[2] >1 || sparse[2] <0)
                    stop("Invalid value for second element of sparse option")
                }
            }
        }

    #
    # Now, process the random formula.
    #  First, check for the common mistake of age~sex instead of ~age|sex
    #  Then get the various crossed terms
    #  Eventually the material following will loop over the crossed
    # terms, finding the grouping and covariates for each
    #
    temp <- terms(random)
    if (attr(temp, 'response') > 0)
        stop("Random formula cannot have a variable before the ~")
    xform <- getCrossedTerms(random)
    if (length(xform) > 1) 
        stop ("Sorry, crossed random effects are still a work in progress")

    # Find out the grouping
    gform <- getGroupsFormula2(xform[[1]])
    gnames <- attr(terms(gform), 'term.labels')
    ncluster <- length(gnames)
    if (ncluster==0) stop("No grouping variables found")
    groups <- getGroups(m, gform)
    # If there is nesting, create nested "names" for the lower
    #   levels of nesting
    if (ncluster >1) {
        temp <- groups
	for (i in 2:ncluster)
	    temp[,i] <- strata2(groups[,1:i], shortlabel=shortlabel,
                                 sep='/')
        groups <- temp
	}

    # Any random slopes?
    cform <- getCovariateFormula2(xform[[1]])
    cterms <- terms(cform)
    if (length(attr(cterms, 'term.labels')) >0) {
        # There are covariates: we need to augment the X matrix
        attr(cterms, 'intercept') <- 0
        tempx <- model.matrix(cterms, m)
        ngroup <- max(kindex)
        x2 <- matrix(0., nrow=length(tempx), ncol=ngroup*ncol(tempx))
        for (i in 1:max(kindex)) {
            j <- i + seq(0, length=ncol(tempx), by=ngroup)
            x2[,j] <- tempx* (kindex==i)
            }
        ntheta[1] <- ntheta[1] + ncol(tempx)
        vmat <- varlist[[1]][[1]]
        stop("Random slopes code not yet finished")
	}
    
    #
    # Do initial checking of the variance matrix list
    #
    temp <- coxme.varcheck(ncluster, varlist, n, gnames, groups, sparse,
                           rescale, pdcheck)
    varlist <- temp$varlist
    kindex  <- temp$kindex
    ntheta <- temp$ntheta
    
    # Set up theta to the right length
    theta.names <- NULL
    for (i in 1:ncluster) {
        if (ntheta[i]==1) theta.names <- c(theta.names, gvars[i])
        else  theta.names <- c(theta.names, 
                               paste(gnames[i], 1:ntheta[i], sep=""))
        }
    if (length(theta)==0) theta <- rep(0., sum(ntheta))
    else if (length(theta) != sum(ntheta)) stop("Wrong length for theta")
    names(theta) <- theta.names

    time1 <- proc.time()
    fit <- coxme.fit(X, Y, strats, offset, init, control, weights=weights,
			    ties=ties, row.names(m),
		            kindex, varlist, ntheta, theta, vinit)
    time2 <- proc.time()

    if (is.character(fit)) {
	fit <- list(fail=fit)
	oldClass(fit) <- 'coxme'
	}
    else {
	fcoef <- fit$coefficients$fixed
        nvar <- length(fcoef)
	if (length(fcoef)>0 && any(is.na(fcoef))) {
	    vars <- (1:length(fcoef))[is.na(fcoef)]
	    msg <-paste("X matrix deemed to be singular; variable",
			   paste(vars, collapse=" "))
	    if (singular.ok) warning(msg)
	    else             stop(msg)
	    }
        if (length(fcoef) >0) {
            names(fcoef) <- dimnames(X)[[2]]
            fit$coefficients <- list(fixed=fcoef, random=fit$coeff$random)
            }

	fit$n <- nrow(Y)
	oldClass(fit) <-  'coxme'
	fit$terms <- Terms
	fit$assign <- attr(X, 'assign')

	#Wald test
	#not for intercept only models, or if test is already done
	if (nvar>0  && is.null(fit$wald.test)) { 
            nfrail <- ncol(fit$var) - nvar
	    nabeta <- !is.na(fcoef)
	    # The init vector might be longer than the betas, for a sparse term
	    if (is.null(init)) temp <- fcoef[nabeta]
	    else temp <- (fcoef - init[1:length(fcoef)])[nabeta]
            wvar <- as.matrix(fit$var[nfrail + 1:nvar, nfrail+1:nvar])
	    fit$wald.test <-  coxph.wtest(wvar[nabeta,nabeta], temp,
					  control$toler.chol)$test
	    }
	na.action <- attr(m, "na.action")
	if (length(na.action)) fit$na.action <- na.action
        if (x)  {
            fit$x <- X
            if (length(strats)) fit$strata <- strata.keep
            }
        if (y)     fit$y <- Y
        }

    if (!is.null(weights) && any(weights!=1)) fit$weights <- weights

    # Label the frailty terms
    #  If there is a single term, leave it as a vector
    #  If there are multiple terms, return a list
    # Also create the linear predictor
    if (ncluster==1) {
	names(fit$frail) <- dimnames(varlist[[1]][[1]])[[1]]
        flinear <- fit$frail[kindex]
	}
    else {
	ftemp <- vector('list', ncluster)
	j <- 0
        flinear <- 0
	for (i in 1:ncluster) {
	    tname <- dimnames(varlist[[i]][[1]])[[1]]
	    nf <- length(tname)
	    temp <- fit$frail[j + 1:nf]
            flinear <- flinear + temp[kindex[,i]]
	    names(temp) <- tname
	    ftemp[[i]]<- temp
	    j <- j+nf
	    }
        names(ftemp) <- gnames
	fit$frail <- ftemp
	}
    if (nvar ==0) fit$linear.predictor <- flinear
    else fit$linear.predictor <- flinear + c(X %*% fit$coef$fixed)

    time3 <- proc.time()
    timeused <- c((time1[1]+ time1[2]) - (time0[1] + time0[2]), fit$timeused,
		  (time3[1]+ time3[2]) - (time2[1] + time2[2]))
    timeused <- c(sum(timeused), timeused)
    names(timeused) <- c("Total", "setup", "fit1", "fit2", "fit3", "finish")
    fit$timeused <- timeused
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    fit$ties <- ties
    fit$kindex <- kindex
    names(fit$loglik) <- c("NULL", "Integrated", "Penalized")
    oldClass(fit) <- 'coxme'
    fit
}
