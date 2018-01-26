\name{solve.gchol}
\alias{solve.gchol}
\title{
  Solve a matrix equation using the generalized Cholesky decompostion
}
\description{
This function solves the equation Ax=b for x, given b and the
generalized Cholesky decompostion of A.
If only the first argument is given, then a G-inverse of A is
returned.
}
\usage{
\method{solve}{gchol}(a, b, full=T, \dots)
}
\arguments{
\item{a}{
a generalized cholesky decompostion of a matrix, as
returned by the \code{gchol} function.
}
\item{b}{
a numeric vector or matrix, that forms the right-hand side of the equation.
}
\item{full}{
solve the problem for the full (orignal) matrix, or for the cholesky
matrix.
}
\item{\ldots}{
an argument to achieve compatibility with solve from R base
}
}
\value{
if argument \code{b} is not present, the inverse of
\code{a} is returned, otherwise the solution to 
matrix equation.
}
\details{
A symmetric matrix A can be decomposed as LDL', where L is a lower
triangular matrix with 1's on the diagonal, L' is the transpose of
L, and D is diagonal.  
This routine solves either the original problem Ay=b 
(\code{full} argument) or the subproblem sqrt(D)L'y=b.
If \code{b} is missing it returns the inverse of
A or L, respectively.
}
\seealso{
\code{\link{gchol}}
}
\examples{
\dontrun{
# Create a matrix that is symmetric, but not positive definite
#   The matrix temp has column 6 redundant with cols 1-5
smat <- matrix(1:64, ncol=8)
smat <- smat + t(smat) + diag(rep(20,8))  #smat is 8 by 8 symmetric
temp <-  smat[c(1:5, 5:8), c(1:5, 5:8)]
ch1  <- gchol(temp)

print(as.matrix(ch1))   # print out L
print(diag(ch1))        # print out D
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
aeq(diag(ch1)[6], 0)    # Check that it has a zero in the proper place

ginv <- solve(ch1)    # see if I get a generalized inverse
aeq(temp \%*\% ginv \%*\% temp, temp)
aeq(ginv \%*\% temp \%*\% ginv, ginv)
}
}
\keyword{algebra}
% docclass is function
% Converted by Sd2Rd version 1.21.
