% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kin.morgan.R
\name{kin.morgan}
\alias{kin.morgan}
\title{kinship matrix for simple pedigree}
\usage{
kin.morgan(ped, verbose = FALSE)
}
\arguments{
\item{ped}{individual's id, father's id and mother's id.}

\item{verbose}{an option to print out the original pedigree.}
}
\value{
The returned value is a list containing:
\itemize{
\item kin the kinship matrix in vector form.
\item kin.matrix the kinship matrix.
}
}
\description{
kinship matrix for simple pedigree
}
\details{
kinship matrix according to Morgan v2.1.
}
\note{
The input data is required to be sorted so that parents preceed their children.
}
\examples{
\dontrun{
# Werner syndrome pedigree
werner<-c(
 1, 0,  0,  1,
 2, 0,  0,  2,
 3, 0,  0,  2,
 4, 1,  2,  1,
 5, 0,  0,  1,
 6, 1,  2,  2,
 7, 1,  2,  2,
 8, 0,  0,  1,
 9, 4,  3,  2,
10, 5,  6,  1,
11, 5,  6,  2,
12, 8,  7,  1,
13,10,  9,  2,
14,12, 11,  1,
15,14, 13,  1)
werner<-t(matrix(werner,nrow=4))
kin.morgan(werner[,1:3])
}

}
\references{
Morgan V2.1 \url{https://sites.stat.washington.edu/thompson/Genepi/MORGAN/Morgan.shtml}
}
\seealso{
\code{\link{gif}}
}
\author{
Morgan development team, Jing Hua Zhao
}
\keyword{datagen}
