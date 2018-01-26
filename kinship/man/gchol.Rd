\name{gchol}
\alias{gchol}
\alias{gchol.bdsmatrix}
\title{
  Generalized Cholesky decompostion
}
\description{
Perform the generalized Cholesky decompostion of a real symmetric matrix.
}
\usage{
gchol(x, tolerance=1e-10)
}
\arguments{
\item{x}{
the symmetric matrix to be factored
}
\item{tolerance}{
the numeric tolerance for detection of singular columns in x.
}
}
\value{
an object of class \code{gchol} containing the
generalized Cholesky decompostion.
It has the appearance of a lower triangular matrix.
}
\details{
A symmetric matrix A can be decomposed as LDL', where L is a lower
triangular matrix with 1's on the diagonal, L' is the transpose of
L, and D is diagonal.
The inverse of L is also lower-triangular, with 1's on the diagonal.  
If all elements of D are positive, then A must be symmetric
positive definite (SPD), and the solution can be reduced the usual
Cholesky decomposition U'U where U is upper triangular and
U = sqrt(D) L'.

The main advantage of the generalized form is that it admits
of matrices that are not of full rank: D will contain zeros 
marking the redundant columns, and the rank of A is the
number of non-zero columns.  If all elements of D are zero or
positive, then A is a non-negative definite (NND) matrix.
The generalized form also has the (quite minor) numerical advantage
of not requiring square roots during its calculation.
To extract the components of the decompostion, use the 
\code{diag} and \code{as.matrix}
functions.

The \code{solve} has a method for gchol decompostions,
and there are gchol methods for block diagonal symmetric
(\code{bdsmatrix}) matrices as well.
}
\seealso{
\code{\link{bdsmatrix}, \link{solve.gchol}}
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
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 1.21.
