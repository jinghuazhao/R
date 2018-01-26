# $Id
# Matrix multiplication for symmetric block diagonal (bds) matrices
# renamed from %*%.bdsmatrix.s

# Both S3 and S4 work
#'%*%' <- function(x,y) UseMethod('%*%') 
#'%*%.default' <- base::"%*%" #function(x,y) get("%*%",NULL)(x,y)

setGeneric("%*%")

'%*%.bdsmatrix' <- function(x, y) {
    if (!inherits(x, 'bdsmatrix'))
	stop("First argument must be a bdsmatrix object")
    if (inherits(y, 'bdsmatrix'))
	stop("Product of two bdsmatrices is not yet supported")
    if (!is.numeric(y))
	stop("Matrix multiplication is defined only for numeric objects")
    dy <- dim(y)
    dx <- dim(x)
    ldy <- length(dy)
    if (ldy!=2) dy <- c(length(y), 1)
    if (dx[2] != dy[1]) 
	stop("Number of columns of x should be the same as number of rows of y")

    # Do the multiplication in C code.  Y is replaced by the result
    #  (Since x is a square matrix, the result is the same size as y)
    nblock <- length(x@blocksize)
    temp <- .C("bdsmatrix_prod", 
	       as.integer(nblock),
	       as.integer(x@blocksize),
	       as.integer(dy),
	       as.double(x@blocks),
	       as.double(x@rmat),
	       as.double(x@offdiag),
	       temp = double(dy[1]),
	       itemp= integer(max(x@blocksize)),
	       y =   as.double(y), PACKAGE="kinship")
    z <- matrix(temp$y, nrow=dx[1])

    # Create dimnames for the result, using the dimnames of the input args
    dnx <- dimnames(x)
    dny <- dimnames(y)
    if(!is.null(dnx) || !is.null(dny)) {
	dnz <- list(NULL, NULL)
	if(!is.null(dnx))
	    dnz[1] <- dnx[1]
	if(!is.null(dny))
	    dnz[2] <- dny[2]
	dimnames(z) <- dnz
        }
    z
    }

#
# This allows for multiplication in the other direction
#
setMethod("%*%", signature(x='matrix', y='bdsmatrix'),
    function(x, y) {
	t(y%*% t(x))
	})

setMethod("%*%", signature(x='numeric', y='bdsmatrix'),
    function(x, y) {
	t(y%*% x)
	})

#
# This is essentially the inverse of solve(gchol(x),y, full=F)
#  Needed for residuals in lmekin (at least for now)
#
'%*%.gchol.bdsmatrix' <- function(x, y) {
    if (!inherits(x, 'gchol.bdsmatrix'))
    if (inherits(y, 'bdsmatrix'))
	stop("Product of two bdsmatrices is not yet supported")
    if (!is.numeric(y))
	stop("Matrix multiplication is defined only for numeric objects")
    dy <- dim(y)
    dx <- dim(x)
    ldy <- length(dy)
    if (ldy!=2) dy <- c(length(y), 1)
    if (dx[2] != dy[1]) 
	stop("Number of columns of x should be the same as number of rows of y")

    # Do the multiplication in C code.  Y is replaced by the result
    #  (Since x is a square matrix, the result is the same size as y)
    nblock <- length(x@blocksize)
    temp <- .C("bdsmatrix_prod3", 
	       as.integer(nblock),
	       as.integer(x@blocksize),
	       as.integer(dy),
	       as.double(x@blocks),
	       as.double(x@rmat),
	       temp = double(dy[1]),
	       y =   as.double(y), PACKAGE="kinship")

    if (dy[2]==1) temp$y
    else matrix(temp$y, nrow=dx[1])

    # About the only thing to call this, ever, is lmekin.
    # So don't bother about dimnames
    }

