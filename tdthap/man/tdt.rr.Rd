\name{tdt.rr}
\alias{tdt.rr}
\title{
Calculate haplotype relative risks in TDT studies
}
\usage{
tdt.rr(hap, prior.df=c(0.5, 0.5), prob=c(0.05, 0.95))
}
\arguments{
\item{hap}{
A list containing the transmitted and untransmitted haplotypes. This would 
normally be computed using \code{\link{tdt.select}}.
}
\item{prior.df}{
a vector of length two containing the degree of freedom parameters for the 
prior distribution of the haplotype relative risk - a beta distribution of the
second kind. 
}
\item{prob}{
The probability levels for Bayesian credibility intervals for the haplotype
relative risks.
}}
\value{
A matrix containing the numbers of transmitted and untransmitted haplotypes, 
the (binomial) p-values, the Bayes estimates of the haplotype relative risks, 
and the lower and upper bounds of the credible interval. The prior estimate 
and credible interval is also shown.  
}
\description{
The p-value is the conventional "exact" test based on the binomial distribution
of transmissions. The estimated relative risks use a Bayesian method, 
recommended because of the multiplicity problem. the prior is a beta 
distribution of the second kind, defined by two "degrees of freedom" 
parameters. Note that the prior mean is prior.df[1]/prior.df[2] and that  
Bayes estimates based on small numbers of transmissions are pulled in towards 
this. A "realistic" choice of these parameters is recommended, and to aid 
this, the function returns credible intervals using the prior alone as well as
the a posteriori interval for each haplotype. 
}
\references{
Spielman R., McGinnis R., and Ewens, W. (1993) Transmission tests for linkage 
disequilibrium. American Journal of Human Genetics, 52, 506-16.
}
\seealso{
\code{\link{hap.transmit}}, \code{\link{tdt.select}}, \code{\link{tdt.quad}}
}
\examples{
\dontrun{
# Select the sub-haplotype made up from the first two markers and 
# print tables of TDT tests and haplotype realtaive risks


	hap.use <- tdt.select(haps, markers=1:2)
	rr <- tdt.rr(hap.use)
	rr
}}
\keyword{htest}
% Converted by Sd2Rd version 0.2-a4.
