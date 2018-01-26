\name{bdsBlock}
\alias{bdsBlock}
\title{
  Block diagonal matrices.
}
\description{
Create a block-diagonal matrix of ones.
}
\usage{
bdsBlock(id, group)
}
\arguments{
\item{id}{
the identifier list.  
This will become the dimnames of the final matrix, and must be a 
set of unique values.
It's length determines the dimension of the final matrix
}
\item{group}{
a vector giving the grouping structure.
All rows/cols belonging to a given group will form a block of 
1's in the final matrix.
}
}
\value{
a block-diagonal matrix of class \code{bdsmatrix}
}
\seealso{
\code{\link{bdsmatrix}, \link{bdsI}}
}
\examples{
\dontrun{
id    <- letters[1:10]
group <- c(1,1,3,2,3,3,2,3,2,4)
bdsBlock(id, group)
    a b d g i c e f h j 
  a 1 1 0 0 0 0 0 0 0 0
  b 1 1 0 0 0 0 0 0 0 0
  d 0 0 1 1 1 0 0 0 0 0
  g 0 0 1 1 1 0 0 0 0 0
  i 0 0 1 1 1 0 0 0 0 0
  c 0 0 0 0 0 1 1 1 1 0
  e 0 0 0 0 0 1 1 1 1 0
  f 0 0 0 0 0 1 1 1 1 0
  h 0 0 0 0 0 1 1 1 1 0
  j 0 0 0 0 0 0 0 0 0 1

# Create the matrices for a sparse nested fit of family within city
group <- paste(mydata$city, mydata$family, sep='/')
mat1 <- bdsI(group)
mat2 <- bdsBlock(group, mydata$city)
fit <- coxme(Surv(time, status) ~ age + sex, data=mydata,
               random= ~1|group, varlist=list(mat1, mat2))
}
}
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 1.21.
