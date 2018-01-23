\name{print.hap.score}
\alias{print.hap.score}
\title{
  Print a hap.score object
}
\description{
Method function to print a class of type hap.score
}
\usage{
\method{print}{hap.score}(x, ...)
}
\arguments{
\item{x}{ The object returned from hap.score (which has class hap.score).}
\item{...}{Optional argunents.}
}
\value{
Nothing is returned.
}
\details{
This is a print method function used to print information from
hap.score class, with haplotype-specific information given in a
table. Because hap.score is a class, the generic print function 
can be used, which in turn calls this print.hap.score function.
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM, Poland  GA (2002)
Score tests for association of traits with haplotypes when
linkage phase is ambiguous. Amer J Hum Genet 70:425-34
}
\seealso{
\code{\link[gap]{hap.score}}
}
\examples{
\dontrun{
save <- hap.score(y, geno, trait.type = "gaussian")

# Example illustrating generic print function:
print(save)

# Example illustrating specific method print function:
print.hap.score(save)
}
}
\keyword{print}
