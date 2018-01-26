# $Id: print.lmekin.s,v 1.1 2003/01/07 21:39:53 Therneau Exp $
print.lmekin <- function(x, ...) {
    cat("Linear mixed-effects kinship model fit by maximum likelihood\n")
    cat("  Data:", deparse(x$call$data), "\n")
    if(!is.null(x$call$subset)) {
        cat("  Subset:", deparse(x$call$subset), "\n")
        }
    cat("  Log-likelihood =", format(x$loglik), "\n")
    omit <- x$na.action
    if(length(omit))
        cat("  n=", x$n, " (", naprint(omit), ")\n", sep = "")
    else cat("  n=", x$n, "\n\n")

    fixF <- x$call$fixed
    if(inherits(fixF, "formula") || is.call(fixF)) {
        cat("Fixed effects:", deparse(as.vector(x$call$fixed)), "\n")
        }
    else {
        cat("Fixed effects:", deparse(lapply(fixF, function(el)
                                       as.name(deparse(as.vector(el))))), "\n")
        }
    print(x$ctable)
    cat("\n")
    l2 <- dim(x$var)[1]
    if (l2>1) {
        l <- 2:l2
        df <- length(l)
        w.chisq <- x$coef$fixed[l]%*%solve(x$var[l,l])%*%x$coef$fixed[l]
        cat("Wald test of fixed effects = ", w.chisq, "df = ", df, "p = ", 1-pchisq(w.chisq,df))
        cat("\n\n")
    }
    
    randF <- x$call$random
    if(inherits(randF, "formula") || is.call(randF)) {
        cat("Random effects:", deparse(as.vector(randF)), "\n")
        }
    else {
        cat("Random effects:", deparse(lapply(randF, function(el)
                                       as.name(deparse(as.vector(el))))), "\n")
        }
    varlist <- x$call$varlist
    if (!is.null(varlist)) {
        cat(" Variance list:", deparse(varlist), "\n")
        }

    rcoef <- x$theta
    temp <- matrix(c(sqrt(rcoef), rcoef/sum(rcoef)), nrow=2, byrow=T)
    dimnames(temp) <- list(c("Standard Dev:", "% Variance:"), names(rcoef))
    print(temp)
    }
