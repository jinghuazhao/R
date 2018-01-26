\name{kinship}
\alias{kinship}
\alias{kindepth}
\title{
  Compute a kinship matrix
}
\description{
Computes the n by n kinship matrix for a set of n related subjects
}
\usage{
kinship(id, father.id, mother.id)
}
\arguments{
\item{id}{
a vector of subject identifiers.  
It may be either numeric or character.
}
\item{father.id }{
for each subject, the identifier of the biological father.
}
\item{mother.id }{
for each subject, the identifier of the biological mother.
}
}
\value{
a matrix of kinship coefficients.
}
\details{
Two genes G1 and G2 are identical by descent (ibd) if they are both physical
copies of the same ancestral gene; two genes are identical by state if
they represent the same allele.  
So the brown eye gene that I inherited from my mother is ibd with hers;
the same gene in an unrelated individual is not.

The kinship coefficient between two subjects is the probability that a
randomly selected allele will be ibd between them.
It is obviously 0 between unrelated individuals. 
If there is no inbreeding in the pedigree, it will be .5 for an individual
with themselves (we could choose the same allele twice), .25 between mother
and child, etc.

The computation is based on a recursive algorithm described in Lange.
It is unfortunately not vectorizable, so the S code is slow.
For studies with multiple disjoint families see the 
\code{makekinship} routine.
}
\section{References}{
K Lange, Mathematical and Statistical Methods for Genetic Analysis,
Springer-Verlag, New York, 1997.
}
\seealso{
\code{\link{makekinship}, \link{makefamid}}
}
\examples{
\dontrun{
test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
round(8*kinship(test1$id, test1$dad, test1$mom))

   1 2 3 4 5 6 7 8 9 10 11 12 13 14 
 1 4 0 0 0 2 2 0 0 1  0  0  0  0  0
 2 0 4 0 0 2 2 0 0 1  2  0  0  0  1
 3 0 0 4 0 0 0 2 2 2  1  0  0  0  0
 4 0 0 0 4 0 0 2 2 0  1  0  0  0  0
 5 2 2 0 0 4 2 0 0 1  1  0  0  0  0
 6 2 2 0 0 2 4 0 0 2  1  0  0  0  0
 7 0 0 2 2 0 0 4 2 1  2  0  0  0  1
 8 0 0 2 2 0 0 2 4 1  1  0  0  0  0
 9 1 1 2 0 1 2 1 1 4  1  0  0  0  0
10 0 2 1 1 1 1 2 1 1  4  0  0  0  2
11 0 0 0 0 0 0 0 0 0  0  4  0  2  1
12 0 0 0 0 0 0 0 0 0  0  0  4  2  1
13 0 0 0 0 0 0 0 0 0  0  2  2  4  2
14 0 1 0 0 0 0 1 0 0  2  1  1  2  4
}
}
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 1.21.
