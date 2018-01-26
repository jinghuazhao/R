\name{as.matrix.bdsmatrix}
\alias{as.matrix.bdsmatrix}
\title{
a function for bdsmatrix
}
\description{
A function for bdsmatrix
}
\usage{
\method{as.matrix}{bdsmatrix}(x,...)
}
\arguments{
\item{x}{a bdsmatrix}
\item{...}{additional arguments to be passed}
}
\examples{
\dontrun{
# The function is currently defined as
function(x)
{
	if(class(x) != "bdsmatrix")
		stop("argument must be a bdsmatrix object")
	dd <- dim(x)
	d3 <- sum(x@blocksize)
	# dim of square portion
	d4 <- sum(x@blocksize^2)
	# size of x@blocks
	newmat <- matrix(0., dd[1], dd[2], dimnames = x@.Dimnames)
	temp <- .C("bdsmatrix_index1",
		as.integer(length(x@blocksize)),
		as.integer(x@blocksize),
		as.integer(c(1, 0, 0)),
		as.integer(d3),
		as.integer(1:d3 - 1),
		indexa = integer(d3 * d3),
		indexb = 0,
		indexc = 0)$indexa
	newmat[x@permute, x@permute] <- c(x@offdiag, x@blocks)[1 + temp]
	if(length(x@rmat) > 0) {
		newmat[,  - (1:d3)] <- x@rmat
		newmat[ - (1:d3),  ] <- t(x@rmat)
	}
	newmat
}
}
}
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 1.21.
