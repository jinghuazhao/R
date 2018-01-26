\name{bdsI}
\alias{bdsI}
\title{
  Sparse identity matrices
}
\description{
This function will create an identitiy matrix, in the sparse
\code{bdsmatrix} format.
}
\usage{
bdsI(id, blocksize)
}
\arguments{
\item{id}{
the identifier list.  
This will become the dimnames of the final matrix, and must be a 
set of unique values.
It's length determines the dimension of the final matrix
}
\item{blocksize}{
the blocksize vector of the final matrix. 
If supplied, the sum of blocksizes must equal the dimension of the matrix. 
By default, the created matrix is as sparse as possible.
}
}
\value{
an identity matrix.
}
\examples{
imat <- bdsI(1:10)
}
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 1.21.
