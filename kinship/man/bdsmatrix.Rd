\name{bdsmatrix}
\alias{bdsmatrix}
\alias{print.bdsmatrix}
\title{
  Create a sparse symmetric block diagonal matrix object
}
\description{
 Sparse block diagonal matrices are used in the the large parameter matrices
that can arise in random-effects coxph and survReg models.  This routine
creates such a matrix.  Methods for these matrices allow them to be 
manipulated much like an ordinary matrix, but the total memory use can be
much smaller.
}
\usage{
bdsmatrix(blocksize, blocks, rmat, dimnames)
}
\arguments{
\item{blocksize}{
 vector of sizes for the matrices on the diagonal
}
\item{blocks}{
 contents of the diagonal blocks, strung out as a vector
}
\item{rmat}{
 the dense portion of the matrix, forming a right and lower border
}
\item{dimnames}{
a list of dimension names for the matrix
}
}
\value{
an object of type bdsmatrix
}
\details{
Consider the following matrix, which has been divided into 4 parts.
\preformatted{
             1  2  0  0  0 | 4  5
	     2  1  0  0  0 | 6  7
	     0  0  3  1  2 | 8  8
	     0  0  1  4  3 | 1  1
	     0  0  2  3  5 | 2  2
	     --------------+-----
	     4  6  8  1  2 | 7  6
	     5  7  8  1  2 | 6  9
}
The upper left is block diagonal, and can be stored in a compressed form
without the zeros.  With a large number of blocks, the zeros can
actually account for over 99\% of a matrix; this commonly happens with
the kinship matrix for a large collection of families (one block/family).
The arguments to this routine would be block sizes of 2 and 3, along with
a 2 by 7 "right hand" matrix.  Since the matrix is symmetrical, the bottom
slice is not needed.
}
\examples{
# The matrix shown above is created by
tmat <- bdsmatrix(c(2,3), c(1,2,1, 3,1,2, 4,3, 5),
                  rmat=matrix(c(4,6,8,1,2,7,6, 5,7,8,1,2,6,9), ncol=2))

# Note that only the lower part of the blocks is needed, however, the
#  entire block set is also allowed, i.e., c(1,2,2,1, 3,1,2,1,4,3,2,3,5)
}
\keyword{manip}
% docclass is function
% Converted by Sd2Rd version 1.21.
