\name{bdsmatrix.ibd}
\alias{bdsmatrix.ibd}
\title{
  Create a bdsmatrix from a list
}
\description{
Routines that create identity-by-descent (ibd) coefficients
often output their results as a list of values (i, j, x[i,j]),
with unlisted values of the x matrix assumed to be zero.
This routine recasts such a list into 
\code{bdsmatrix} form.
}
\usage{
bdsmatrix.ibd(id1, id2, x, idmap, diagonal)
}
\arguments{
\item{id1}{
row identifier for the value, in the final matrix.
Optionally, \code{id1} can be a 3 column
matrix or data.frame, in which case it is assumed to contain
the first 3 arguments, in order.
}
\item{id2}{
column identifier for the value, in the final matrix.
}
\item{x}{
the value to place in the matrix
}
\item{idmap}{
a two column matrix or data frame.
Sometimes routines create output with integer values for 
\code{id1} and \code{id2},
and then this argument is the mapping from this internal label to
the ``real'' name)
}
\item{diagonal}{
If diagonal elements are not preserved in the list, this value
will be used for the diagonal of the result.
If the argument appears, then the output matrix will contain
an entry for each value in \code{idlist}.
Otherwise only those with an explicit entry appear.
}
}
\value{
a \code{bdsmatrix} object representing a
block-diagonal sparse matrix.
}
\details{
The routine first checks for non-symmetric or otherwise inconsistent input.
It then groups observations together into `families' of
related subjects, which determines the structure of the final matrix.
As with the \code{makekinship} function, 
singletons with no relationships are first in the output matrix, and
then families appear one by one.
}
\seealso{
\code{\link{bdsmatrix}, \link{kinship}, \link{coxme}, \link{lmekin}}
}
\examples{
\dontrun{
ibdmat <- bdsmatrix.ibd(i,j, ibdval, idlist=subject)
}
}
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 1.21.
