\name{plot.hap.score}
\alias{plot.hap.score}
\title{
  Plot haplotype frequencies versus haplotype score statistics
}
\description{
Method function to plot a class of type hap.score
}
\usage{
\method{plot}{hap.score}(x, ...)
}
\arguments{
\item{x}{The object returned from hap.score (which has class hap.score).}
\item{...}{Optional arguments}
}
\value{
Nothing is returned.
}
\details{
This is a plot method function used to plot haplotype frequencies on
the x-axis and haplotype-specific scores on the y-axis. Because
hap.score is a class, the generic plot function 
can be used, which in turn calls this plot.hap.score function.
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA (2002)
Score tests for association of traits with haplotypes when
linkage phase is ambiguous. Amer J Hum Genet 70:425-34
}
\seealso{
\code{\link[gap]{hap.score}}
}
\examples{
\dontrun{
save <- hap.score(y, geno, trait.type = "gaussian")

# Example illustrating generic plot function:
plot(save)

# Example illustrating specific method plot function:
plot.hap.score(save)
}
}
\keyword{hplot}
