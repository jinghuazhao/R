\name{hwe.cc}
\alias{hwe.cc}
\title{A likelihood ratio test of population Hardy-Weinberg equilibrium for case-control studies}
\usage{hwe.cc(model, case, ctrl, k0, initial1, initial2)}

\description{
A likelihood ratio test of population Hardy-Weinberg equilibrium for case-control studies
}
\arguments{
 \item{model}{model specification, dominant, recessive}
 \item{case}{a vector of genotype counts in cases}
 \item{ctrl}{a vector of genotype counts in controls}
 \item{k0}{prevalence of disease in the population}
 \item{initial1}{initial values for beta, gamma, and q}
 \item{initial2}{initial values for logit(p) and log(gamma)}
}

\details{
This is a collection of utility functions. The null hypothesis declares that the proportions of
genotypes are according to Hardy-Weinberg law, while under the alternative hypothesis, the expected
genotype counts are according to the probabilities that particular genotypes are obtained conditional
on the prevalence of disease in the population. In so doing, Hardy-Weinberg equilibrium is considered
using both case and control samples but pending on the disease model such that 2-parameter multiplicative
model is built on baseline genotype \eqn{\alpha}{alpha}, \eqn{\alpha\beta}{alpha beta} and
\eqn{\alpha\gamma}{alpha gamma}.
}

\value{
The returned value is a list with the following components.
\item{Cox}{statistics under a general model}
\item{t2par}{under the null hypothesis}
\item{t3par}{under the alternative hypothesis}
\item{lrt.stat}{the log-likelihood ratio statistic}
\item{pval}{the corresponding p value}
}

\section{References}{
Yu C, Zhang S, Zhou C, Sile S. A likelihood ratio test of population Hardy-Weinberg equilibrium
for case-control studies. Genetic Epidemiology 33:275-280, 2009
}

\seealso{
\code{\link[gap]{hwe}}
}

\examples{
\dontrun{

### Saba Sile, email of Jan 26, 2007, data always in order of GG AG AA, p=Pr(G),
### q=1-p=Pr(A)
case=c(155,27,4)
ctrl=c(408,55,15)
k0=.2
initial1=c(1.0,0.94,0.0904)
initial2=c(logit(1-0.0904),log(0.94))
hwe.cc("recessive",case,ctrl,k0, initial1, initial2)

### John Phillips III, TGFb1 data codon 10: TT CT CC, CC is abnormal and increasing
### TGFb1 activity
case=c(29,78,13)
ctrl=c(17,28,6)
k0 <- 1e-5
initial1 <- c(2.45,2.45,0.34)
initial2 <- c(logit(1-0.34),log(2.45))
hwe.cc("dominant",case,ctrl,k0,initial1,initial2)
}
}
\author{Chang Yu, \url{http://biostat.mc.vanderbilt.edu/wiki/Main/ChangYu}, Li Wang, Jing Hua Zhao}
\keyword{htest}
