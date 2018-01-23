\name{hwe}
\alias{hwe}
\title{Hardy-Weinberg equlibrium test for a multiallelic marker}
\usage{hwe(data, data.type="allele", yates.correct=FALSE, miss.val=0)}

\description{
Hardy-Weinberg equilibrium test
}
\arguments{
 \item{data}{A rectangular data containing the genotype, or an array of genotype counts}
 \item{data.type}{An option taking values "allele", "genotype", "count"  if data is alleles, genotype or genotype count}
 \item{yates.correct}{A flag indicating if Yates' correction is used for Pearson \eqn{\chi^2}{chi-squared} statistic}
 \item{miss.val}{A list of missing values}
}

\details{
This function obtains Hardy-Weinberg equilibrium test statistics. It can
handle data coded as allele numbers (default), genotype identifiers (by
setting data.type="genotype") and counts corresponding to individual genotypes
 (by setting data.type="count") which requires that genotype counts for all n(n+1) possible genotypes,
 with n being the number of alleles.

For highly polymorphic markers when asymptotic results do not hold, please resort to hwe.hardy.
}

\value{
The returned value is a list containing:
 \item{allele.freq}{Frequencies of alleles}
 \item{x2}{Pearson \eqn{\chi^2}{chi-square}}
 \item{p.x2}{p value for \eqn{\chi^2}{chi-square}}
 \item{lrt}{Log-likelihood ratio test statistic}
 \item{p.lrt}{p value for lrt}
 \item{df}{Degree(s) of freedom}
 \item{rho}{\eqn{\sqrt{\chi^2/N}}{sqrt{chi-square/N}} the contingency table coefficient}
}

\seealso{
\code{\link[gap]{hwe.hardy}} 
}

\examples{
\dontrun{
a <- c(3,2,2)
a.out <- hwe(a,data.type="genotype")
a.out
a.out <- hwe(a,data.type="count")
a.out
require(haplo.stats)
data(hla)
hla.DQR <- hwe(hla[,3:4])
summary(hla.DQR)
# multiple markers
s <- vector()
for(i in seq(3,8,2))
{
  hwe_i <- hwe(hla[,i:(i+1)])
  s <- rbind(s,hwe_i)
}
s
}
}
\author{Jing Hua Zhao}
\keyword{htest}
