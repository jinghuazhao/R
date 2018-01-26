\name{solve.bdsmatrix}
\alias{solve.bdsmatrix}
\title{
  Solve a matrix equation using the generalized Cholesky decompostion
}
\description{
This function solves the equation Ax=b for x, when
A is a block diagonal sparse matrix 
(an object of class \code{bdsmatrix}).
}
\usage{
\method{solve}{bdsmatrix}(a, b, tolerance=1e-10, full=T, \dots)
}
\arguments{
\item{a}{
a block diagonal sparse matrix object
}
\item{b}{
a numeric vector or matrix, that forms the right-hand side of the equation.
}
\item{tolerance}{
the tolerance for detecting singularity in the a matrix
}
\item{full}{
if true, return the full inverse matrix; if false return only
that portion corresponding to the blocks.  
This argument is ignored if \code{b} is present.
If the bdsmatrix \code{a} has a non-sparse portion,
i.e., if the  \code{rmat} component is present,
then the inverse of \code{a} will not be 
block-diagonal sparse.  
In this case setting full=F returns only a portion of the inverse.
The elements that are returned are those of the full inverse, but the
off-diagonal elements that are not returned would not have been zero.
}
\item{\ldots}{
an argument to achieve compatibility with solve from R base
}
}
\value{
if argument \code{b} is not present, the inverse of
\code{a} is returned, otherwise the solution to 
matrix equation.
The equation is solved using a generalized Cholesky decomposition.
}
\details{
The matrix \code{a} consists of a block diagonal
sparse portion with an optional dense border.
The inverse of \code{a}, which is to be computed if
\code{y} is not provided, will have the same
block diagonal structure as \code{a} only if there
is no dense border, otherwise the resulting matrix will not be sparse.

However, these matrices may often be very large, and a non sparse
version of one of them will require gigabytes of even terabytes of
space.  For one of the
common computations (degrees of freedom in a penalized model) only those
elements of the inverse that correspond to the non-zero part of
\code{a} are required;
the \code{full=F} option returns only that portion
of the (block diagonal portion of) the inverse matrix.
}
\seealso{
\code{\link{bdsmatrix}, \link{gchol}}
}
\examples{
\dontrun{
tmat <- bdsmatrix(c(3,2,2,4), 
              c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
              matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dim(tmat)
solve(tmat, cbind(1:13, rep(1,13)))
}
}
\keyword{algebra}
% docclass is function
% Converted by Sd2Rd version 1.21.
