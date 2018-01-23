\name{gcontrol}
\alias{gcontrol}
\title{genomic control}
\usage{gcontrol(data,zeta,kappa,tau2,epsilon,ngib,burn,idum)}

\description{
The Bayesian genomic control statistics with the following parameters,

\tabular{ll}{
n \tab number of loci under consideration \cr
lambdahat \tab median(of the n trend statistics)/0.46 \cr
  \tab Prior for noncentrality parameter Ai is \cr
  \tab Normal(sqrt(lambdahat)kappa,lambdahat*tau2) \cr
kappa \tab multiplier in prior above, set at 1.6 * sqrt(log(n)) \cr
tau2  \tab multiplier in prior above \cr
epsilon \tab prior probability a marker is associated, set at 10/n \cr
ngib \tab number of cycles for the Gibbs sampler after burn in \cr
burn \tab number of cycles for the Gibbs sampler to burn in
}

Armitage's trend test along with the posterior probability that each marker is
associated with the disorder is given. The latter is not a p-value but any value
greater than 0.5 (pout) suggests association. 
}

\arguments{
\item{data}{the data matrix}
\item{zeta}{program constant with default value 1000}
\item{kappa}{multiplier in prior for mean with default value 4}
\item{tau2}{multiplier in prior for variance with default value 1}
\item{epsilon}{prior probability of marker association with default value 0.01}
\item{ngib}{number of Gibbs steps, with default value 500}
\item{burn}{number of burn-ins with default value 50}
\item{idum}{seed for pseudorandom number sequence}
}

\source{
\url{http://www.stat.cmu.edu}
}

\value{
The returned value is a list containing:

\item{deltot}{the probability of being an outlier}
\item{x2}{the \eqn{\chi^2}{chi-squared} statistic}
\item{A}{the A vector}
}

\section{References}{
Devlin B, Roeder K (1999) Genomic control for association studies. 
Biometrics 55:997-1004
}

\examples{
\dontrun{
test<-c(1,2,3,4,5,6,  1,2,1,23,1,2, 100,1,2,12,1,1, 
        1,2,3,4,5,61, 1,2,11,23,1,2, 10,11,2,12,1,11)
test<-matrix(test,nrow=6,byrow=T)
gcontrol(test)
}
}

\author{Bobby Jones, Jing Hua Zhao}

\note{Adapted from gcontrol by Bobby Jones and Kathryn Roeder, 
use -Dexecutable for standalone program, function getnum in the original 
code needs \%*s to skip id string}

\keyword{models}
