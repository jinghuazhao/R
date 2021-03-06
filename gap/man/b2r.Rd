\name{b2r}
\alias{b2r}
\title{Obtain correlation coefficients and their variance-covariances}
\usage{b2r(b,s,rho,n)}
\description{
This function converts linear regression coefficients of phenotype on
single nucleotide polymorphisms (SNPs) into Pearson correlation coefficients
with their variance-covariance matrix. It is useful as a preliminary step 
for meta-analyze SNP-trait associations at a given region. Between-SNP
correlations (e.g., from HapMap) are required as auxiliary information.
}
\arguments{
  \item{b}{the vector of linear regression coefficients}
  \item{s}{the corresponding vector of standard errors}
  \item{rho}{triangular array of between-SNP correlation}
  \item{n}{the sample size}
}

\value{
The returned value is a list containing:

\item{r}{the vector of correlation coefficients}
\item{V}{the variance-covariance matrix of correlations}
}

\section{References}{
Becker BJ (2004). Multivariate meta-analysis. in Tinsley HEA,
Brown SD (Ed.) Handbook of Applied Multivariate Statistics and
Mathematical Modeling (Chapter 17, pp499-525). Academic Press.

Casella G, Berger RL (2002). Statistical Inference, 2nd Edition, Duxbury.

Elston RC (1975). On the correlation between correlations. Biometrika 62:133-40
}
\seealso{
\code{\link[gap]{mvmeta}}, \code{\link[gap]{LD22}}
}

\examples{
\dontrun{
n <- 10
r <- c(1,0.2,1,0.4,0.5,1)
b <- c(0.1,0.2,0.3)
s <- c(0.4,0.3,0.2)
bs <- b2r(b,s,r,n)
}
}
\author{Jing Hua Zhao}
\keyword{datagen}
