\name{tdt.quad}
\alias{tdt.quad}
\title{
TDT tests for extended haplotypes
}
\usage{
tdt.quad(hap, nsim=5000, funct=FALSE, keep=TRUE, seeds=c(0, 0, 0))
}
\arguments{
\item{hap}{
A list containing the transmitted and untransmitted haplotypes. This would 
normally be computed using \code{\link{tdt.select}}.
}
\item{nsim}{
The number of Monte Carlo simulations from the null hypothesis.
}
\item{funct}{
If T, a similarity function is used and the test is a Geary-Moran test. 
Otherwise, the Pearsonian test, Sum \eqn{(O-E)^2/E}, is used.
}
\item{keep}{
If TRUE, all simulated values of the test statistic are kept. Otherwise only 
the realised value of the test statistic and the p-value are returned.
}
\item{seeds}{
Three numbers to seed the random number generator. The default is to use 
a different three random numbers each time.
}}
\value{
A list containing, the number of distinct haplotypes (\eqn{n.hap}), the number of 
informative transmissions (\eqn{n.trans}), the test statistic (\eqn{test}), the p-value
(\eqn{p.value}) and, optionally, all the simulated values of the test statistic
under the null hypothesis (\eqn{sim}).
}
\description{
The function calculates the test statistic and then simulates its distribution
under the null hypothesis by randomly transmitting parental haplotypes with 
probability 0.5. The test statistic is recalculated for each simulated dataset.
For Geary-Moran tests in particular this can be quite slow.
}
\references{
Clayton, D. and Jones, H. (1999) Transmission/disequilibrium tests for extended marker 
haplotypes.  Am.J.Hum.Gen., 65:1161-1169.
}
\seealso{
  \code{\link{hap.transmit}}, \code{\link{tdt.select}},
  \code{\link{tdt.rr}}, \code{\link{set.similarity}},
  \code{\link{get.similarity}} 
}
\examples{
\dontrun{
#  Do a Pearsonian test using 10000 simulations and summarise the distribution
#  of the statistic under the null hypothesis


	test <- tdt.quad(hap.use, nsim=10000, keep=T)
	test
	summary(test$sim)
}}
\keyword{htest}
% Converted by Sd2Rd version 0.2-a4.
