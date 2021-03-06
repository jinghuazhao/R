\name{s2k}
\alias{s2k}
\title{Statistics for 2 by K table}
\usage{s2k(y1,y2)}
\description{
This function calculates one-to-others and maximum accumulated chi-squared
statistics for a 2 by K contingency table.
}
\arguments{
\item{y1}{a vector containing the first row of a 2 by K contingency table}
\item{y2}{a vector containing the second row of a 2 by K contingency table}
}

\value{
The returned value is a list containing:

\item{x2a}{the one-to-other chisquare}
\item{x2b}{the maximum accumulated chisquare}
\item{col1}{the column index for x2a}
\item{col2}{the column index for x2b}
\item{p}{the corresponding p value}
}

\section{References}{
Hirotsu C, Aoki S, Inada T, Kitao Y (2001) An exact test for the association 
between the disease and alleles at highly polymorphic loci with particular interest 
in the haplotype analysis. Biometrics 57:769-778
}

\examples{
\dontrun{
# an example from Mike Neale
# termed 'ugly' contingency table by Patrick Sullivan
y1 <- c(2,15,16,35,132,30,25,7,12,24,10,10,0)
y2 <- c(0, 6,31,49,120,27,15,8,14,25, 3, 9,3)

result <- s2k(y1,y2)
}
}
\author{Chihiro Hirotsu, Jing Hua Zhao}
\note{The lengths of y1 and y2 should be the same}
\keyword{models}
