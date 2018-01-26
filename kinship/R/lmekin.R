# $Id: lmekin.s,v 1.8 2003/04/07 18:32:40 Therneau Exp $
#
# An lme function, specialized to kinship matrices, based on equation 2.14
#  of Pinheiro and Bates (the one they say is awful).
# Fits a model with y ~ X\beta + Zb, where Z is the identity, b has
#  length n (one random effect per subject), and 
#  var(y) = \sigma^2 (I + A), A is \tau^2_1 A_1 + ...
#  with usually only 1 or two components, the kinship matrix K + an ibd.
#
lmekin <- function(fixed, data=parent.frame(), random, 
                   varlist=NULL, variance, sparse=c(20,.05),
                   rescale=T, pdcheck=T,
                   subset, weight, na.action) {
    # start with the standard stuff, stolen from coxme
    # however, we assume a random statement that has just one effect, and
    # match it up with the kinship matrix
    call <- match.call()
    m <- match.call(expand.dots=FALSE)
    temp <- c("", "data", "weights", "subset", "na.action")
    m <- m[ match(temp, names(m), nomatch=0)]

    if (missing(variance)) theta <- NULL
    else  theta <- variance  #We always liked "theta" better as a name, but it
                             # didn't seem as obvious to the user community
    # We are borrowing some tools from lme here
    reSt <- reStruct(random, REML=F, data=NULL)
    gform <- getGroupsFormula(reSt)  #formula for the random effects
    if (is.null(gform)) {
        temp.fixed <- fixed
        gvars <- NULL
        }
    else {
        # Be sneaky, and paste the "random" variables onto my formula.
        # This causes the data frame m to have all necesary variables.
        # If the formula is long, deparse() can be a character
        #  string of length >1 --- would print as more than one line.  Thus
        #  the use of "collapse".  Setting a huge line length would probably
        #  be another way around it.
        gvars <- all.vars(random)
        fvars <- all.vars(formula)
        gvars <- gvars[is.na(match(gvars, fvars))]  #no need for duplicates
        
        temp.fixed <- paste(deparse(as.vector(fixed)), collapse='')
        temp.fixed <- paste(temp.fixed, paste(gvars, collapse='+'),
                                sep='+')
        temp.fixed <- as.formula(temp.fixed) 
        }

    m$formula <- temp.fixed
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    Terms <- terms(fixed)
    X <- model.matrix(Terms, m)
    Y <- model.extract(m, "response")
    n <- length(Y)
    weights <- model.extract(m, 'weights')
    offset<- attr(Terms, "offset")
    tt <- length(offset)

    # idiot proofing: if more than one offset statement, just
    #   add the offset terms together
    offset <- if(tt == 0)
		    rep(0, n)
	      else if(tt == 1)
		      m[[offset]]
	      else {
		    ff <- m[[offset[1]]]
		    for(i in 2:tt)
			    ff <- ff + m[[offset[i]]]
		    ff
		    }

    #
    # Do initial checking of the variance matrix list
    #
    ncluster <- length(gvars)
    if (ncluster==0) stop("No grouping variables found")
    groups <- getGroups(m, gform)
    temp <- coxme.varcheck(ncluster, varlist, n, gvars, groups, sparse,
                           rescale, pdcheck)
    varlist <- temp$varlist
    kindex  <- temp$kindex
    ntheta <- temp$ntheta
    
    # Set up theta to the right length
    theta.names <- NULL
    for (i in 1:ncluster) {
        if (ntheta[i]==1) theta.names <- c(theta.names, gvars[i])
        else  theta.names <- c(theta.names, 
                               paste(gvars[i], 1:ntheta[i], sep=""))
        }
    if (length(theta)==0) theta <- rep(0., sum(ntheta))
    else if (length(theta) != sum(ntheta)) stop("Wrong length for theta")
    names(theta) <- theta.names
    # Set up theta to the right length
    theta.names <- NULL
    for (i in 1:ncluster) {
        if (ntheta[i]==1) theta.names <- c(theta.names, gvars[i])
        else  theta.names <- c(theta.names, 
                               paste(gvars[i], 1:ntheta[i], sep=""))
        }
    if (length(theta)==0) theta <- rep(0., sum(ntheta))
    else if (length(theta) != sum(ntheta)) stop("Wrong length for theta")
    names(theta) <- theta.names
    tindex <- which(theta==0)

    #
    # All the above assumes that I can have multiple random effects, which
    #  coxme can do.  This routine is still stuck with only 1, and it
    #  had better be the right size
    #
    if (ncluster >1) stop("function can have only 1 random effect")
    varlist <- varlist[[1]]  #There is only one element
    kindex <- kindex[,1]     # Ditto
    if (max(kindex) != n)
        stop("The random effect must be 1 per subject")
    ntheta <- ntheta[1]

    # This variable orders the data to match kmat
    kindex2 <- integer(n)
    kindex2[kindex] <- 1:n
    
    logfun <- function(itheta, X, Y, varlist, theta, tindex, center) {
        theta[tindex] <- exp(itheta)
        tkmat <- varlist[[1]]
        tkmat@blocks <- tkmat@blocks * theta[1]
        diag(tkmat) <- diag(tkmat +1)
        if (length(varlist) >1) {
            for (i in 2:length(varlist))
                tkmat@blocks <- varlist[[i]]@blocks * theta[i] +
                                  tkmat@blocks
            }
        #
        # The loglik below is invariant to multiplication of tkmat
        #  by a constant, mathematically.  Keeping the diagonal of tkmat
        #  closer to 1 avoids numerical round-off problems, however.
        # (It multiplies one term and divides the other).
        tkmat@blocks <- tkmat@blocks/ tkmat@blocks[1]
        gk <- gchol(tkmat)
        newx <- solve(gk, X, full=FALSE)
        newy <- solve(gk, Y, full=FALSE)
        resid <- qr.resid(qr(newx), newy)
        n <- length(Y)
        loglik <- (n/2)*(log(mean(resid^2)) - center) + 
                         sum(log(diag(gk)))/2
        loglik
        }
    newX <- X[kindex2,]
    newY <- as.vector(Y[kindex2]) #remove any names
    dimnames(newX) <- NULL        #Not carrying these around saves time
    if (length(tindex) > 0) {
        center <- log(mean((Y-mean(Y))^2))
        nfit <- optim(par= rep(-1, length(tindex)), logfun, method="L-BFGS-B",
                      lower=log(.00001), X=newX, Y=newY, 
                      varlist=varlist, theta=theta, tindex=tindex,
                      center=center)
        iter <- nfit$counts
        theta[tindex] <- exp(nfit$par)
        }
    else iter <- 0

    # Ok, now we have found the ratio of error var to random var.
    # One more iteration of the "logfun" computations above, to solve the
    #  problem.
    # For this iteration we need to scale things properly, since we need both
    #  the se and the loglik for printout.
    tkmat <- varlist[[1]]
    tkmat@blocks <- tkmat@blocks * theta[1]
    diag(tkmat) <- diag(tkmat +1)
    if (length(varlist) >1) {
        for (i in 2:length(varlist))
            tkmat@blocks <- varlist[[i]]@blocks * theta[i] +
                                  tkmat@blocks
        }
    gk <- gchol(tkmat)
#    lfit <- lm.fit(as.matrix(solve(gk, newX, full=F)), 
#                   solve(gk, newY, full=F))
    xok <- as.matrix(solve(gk, newX, full=F))
    yok <-  solve(gk, newY, full=FALSE)
    lfit <- lm(yok~0+xok)
    names(lfit$coefficients) <- dimnames(X)[[2]]
    ls <- summary(lfit)
    resid.var <- mean(lfit$residuals^2)   #differs from ls$sigma, division by N
    theta <- c(theta*resid.var, resid.var)
    names(theta) <- c(theta.names, 'resid')
    
    fitted <- c(X %*% lfit$coef)  #fitted, on the original scale
    residuals <- Y - fitted

    # The next line is wrong and I know it is, leave it as a placeholder.
    frail <- residuals[kindex2] 
    names(frail) <- groups
    
    fcoef <- lfit$coef
    call$fixed <- fixed
    call$random <- random
    fit <- list(coefficients=list(fixed=fcoef, random=frail),
                theta = theta,
                variance= ls$cov.unscaled * ls$sigma^2,
                ctable = ls$coefficients,
                residuals= residuals,
                fitted.values= fitted,
                effects=lfit$effects,
                rank = lfit$rank,
                assign=lfit$assign,
                df.residual = lfit$df.residual - length(theta),
                loglik = (-n/2)*(log(mean(lfit$residuals^2)) + 1+ log(2*pi)) - 
                            sum(log(diag(gk)))/2,
                iter = iter, n=n,
                call = call,
                method='ML')
    # For the life of me, I don't know where the +1 in the loglik comes
    #  from.  But if I put it in, I match lme's result....

    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action

    oldClass(fit) <- c('lmekin')
    fit
}
