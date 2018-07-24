\name{hwe.jags}
\alias{hwe.jags}
\title{Hardy-Weinberg equlibrium test for a multiallelic marker using JAGS}
\usage{hwe.jags(k,n,delta=rep(1/k,k),lambda=0,lambdamu=-1,lambdasd=1,
              parms=c("p","f","theta","lambda"),...)}

\description{
Hardy-Weinberg equilibrium test
}
\arguments{
 \item{k}{number of alleles}
 \item{n}{a vector of k(k+1)/2 genotype counts}
 \item{delta}{initial parameter values}
 \item{lambda}{initial parameter values}
 \item{lambdamu}{initial parameter values}
 \item{lambdasd}{initial parameter values}
 \item{parms}{monitored parmeters}
 \item{...}{parameters passed to jags, e.g., n.chains, n.burnin, n.iter}
}

\details{
This function performs Bayesian Hardy-Weinberg equilibrium test, which mirrors \code{\link[gap]{hwe.hardy}},
another implementation for highly polymorphic markers when asymptotic results do not hold.
}

\value{
The returned value is a fitted model from jags().
}

\section{References}{
Wakefield J (2010). Bayesian methods for examining Hardy-Weinberg equilibrium. Biometrics 2010;66:257-65
}

\seealso{
\code{\link[gap]{hwe.hardy}} 
}

\examples{
\dontrun{
ex1 <- hwe.jags(4,c(5,6,1,7,11,2,8,19,26,15))
print(ex1)
ex2 <- hwe.jags(2,c(49,45,6))
print(ex2)
ex3 <- hwe.jags(4,c(0,3,1,5,18,1,3,7,5,2),lambda=0.5,lambdamu=-2.95,lambdasd=1.07)
print(ex3)
ex4 <- hwe.jags(9,c(1236,120,3,18,0,0,982,55,7,249,32,1,0,12,0,2582,132,20,1162,29,1312,6,0,0,4,0,4,0,2,0,0,0,0,0,0,0,115,5,2,53,1,149,0,0,4),
                delta=c(1,1,1,1,1,1,1,1,1),lambdamu=-4.65,lambdasd=0.21)
print(ex4)
ex5 <- hwe.jags(8,n=c(
         3,
         4, 2,
         2, 2, 2,
         3, 3, 2, 1,
         0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 1,
         0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 2, 1, 0, 0, 0))
print(ex5)
}
}
\author{Jing Hua Zhao, Jon Wakefield}
\keyword{htest}
