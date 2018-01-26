# $Id: print.coxme.s,v 1.6 2003/08/09 21:46:42 Therneau Exp $
print.coxme <- function(x, digits=options()$digits, ...) {
    cat("Cox mixed-effects model fit by maximum likelihood\n")
    cat("  Data:", deparse(x$call$data), "\n")
    if(!is.null(x$call$subset)) {
        cat("  Subset:", deparse(x$call$subset), "\n")
	}

    beta <- x$coefficients
    nvar <- length(beta$fixed)
    nfrail<- nrow(x$var) - nvar

    omit <- x$na.action
    if(length(omit))
        cat("  n=", x$n, " (", naprint(omit), ")\n", sep = "")
    else cat("  n=", x$n, "\n")
    temp <- matrix(x$loglik, nrow=1)
    cat("  Iterations=", x$iter, "\n")
    dimnames(temp) <- list("Log-likelihood", 
                           c("NULL", "Integrated", "Penalized"))
    print(temp)
    chi <- 2*diff(x$loglik[c(1,3)]) 
    cat("\n  Penalized loglik: chisq=", format(round(chi,2)), 
        "on", format(round(x$df[2],2)), "degrees of freedom, p=",
        format(signif(1- pchisq(chi,x$df[2]),2)),"\n")
    chi <- 2*diff(x$loglik[1:2]) 
    cat(" Integrated loglik: chisq=", format(round(chi,2)), 
        "on", format(round(x$df[1],2)), "degrees of freedom, p=",
        format(signif(1- pchisq(chi,x$df[1]),2)),"\n\n")

    if (nvar > 0)  { # Not a ~1 model
        fixF <- x$call$fixed
        if(inherits(fixF, "formula") || is.call(fixF)) {
            cat("Fixed effects:", deparse(as.vector(x$call$fixed)), "\n")
            }
        else {
            cat("Fixed effects:", deparse(lapply(fixF, function(el)
                                      as.name(deparse(as.vector(el))))), "\n")
            }

        coef <- beta$fixed
        se <- sqrt(diag(x$var)[nfrail+1:nvar])
        tmp <- cbind(coef, exp(coef), se, round(coef/se,2),
               signif(1 - pchisq((coef/ se)^2, 1), 2))
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
            "se(coef)", "z", "p"))
        prmatrix(tmp)
        }

    cat("\n")
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
    temp <- matrix(beta$random, nrow = 1)
    dimnames(temp) <- list("Variance:", names(beta$random))
    print(temp)
    invisible(x)
    }
