# 
# Defining function for kinship gaussian
#
frailty.kin <- function(x, theta, df, kmat, 
		   method=c("df", "fixed"), ...) {
    nclass <- length(unique(x))
    # Check for consistency of the arguments
    if (missing(method)) {
	if (!missing(theta)) {
	    method <- 'fixed'
	    if (!missing(df)) 
		    stop("Cannot give both a df and theta argument")
	    }
	else if (!missing(df)) {
	    if (df==0) method <- "aic"
	    else       method <- 'df'
	    }
	}
    method <- match.arg(method)
    if (method=='df' && missing(df)) stop("Method = df but no df argument")
    if (method=='fixed' && missing(theta))
	    stop("Method= fixed but no theta argument")
    if (method !='fixed' && !missing(theta)) 
	    stop("Method is not 'fixed', but have a theta argument")

    dmat <- solve(kmat)

    x <- factor(x, levels=dimnames(kmat)[[1]])
    # are there values of x that are not in kmat?
    if (any(is.na(x))) 
	stop("The argument of frailty.kin does not match with kmat")
    # values in kmat that are not in x?
    if (any(table(x) ==0))
	stop("The argument of frailty.kin does not match with kmat")

    oldClass(x) <- "coxph.penalty"
    attr(x,'contrasts') <- function(n,...) contr.treatment(n,F)
    
    if (!missing(theta) & !missing(df)) 
	    stop("Cannot give both a df and theta argument")

    pfun<- function(coef, theta, ndead, dmat){
	if (theta==0) list(recenter=0, penalty=0, flag=T)
	else list(penalty= c(coef %*% dmat %*% coef) /(2*theta),
		  first  = c(dmat %*% coef) /theta ,
		  second = c(dmat /theta),
		  flag=F)
	}

     printfun <- function(coef, var, var2, df, history) {
	if (!is.null(history$history)) 
	     theta <- history$history[nrow(history$history),1]
	else theta <- history$theta
		
	if (is.matrix(var)) test <- coxph.wtest(var, coef)$test
	else 		    test <- sum(coef^2/var)
	df2 <- max(df, .5)      # Stop silly p-values
	list(coef=c(NA, NA, NA, test, df, 1-pchisq(test, df2)),
		 history=paste("Variance of random effect=", format(theta)))
	}

    if (method=='fixed') {
	temp <- list(pfun=pfun,
		     printfun = printfun,
		     pparm = dmat,
		     diag =F,
		     sparse= F,
		     cfun = function(parms, iter, old){
		          list(theta=parms$theta, done=T)},
		     cparm= list(theta=theta, ...))
        }
    else {  #df method
	temp <- list(pfun=pfun,
		     printfun =printfun,
		     pparm=dmat,
		     diag =F,
		     sparse= F,
		     cargs=('df'),
		     cparm=list(df=df, thetas=0, dfs=0,
		                guess=3*df/length(unclass(x)), ...),
                     cfun = frailty.controldf)
	}

    vname <- paste("frail", levels(x), sep=':')
    temp <- c(temp, list(varname=vname))
    
    attributes(x) <- c(attributes(x), temp)
    x
    }

			  
			   
			   
