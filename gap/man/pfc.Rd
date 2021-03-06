\name{pfc}
\alias{pfc}
\title{Probability of familial clustering of disease}
\usage{pfc(famdata,enum)}
\description{
To calculate exact probability of familial clustering of disease
}
\arguments{
\item{famdata}{collective information of sib size, number of affected sibs and their frequencies}
\item{enum}{a switch taking value 1 if all possible tables are to be enumerated}
}

\value{
The returned value is a list containing (tailp,sump,nenum are only available if enum=1):

 \item{p}{the probabitly of familial clustering}
 \item{stat}{the deviances, chi-squares based on binomial and hypergeometric distributions, 
the degrees of freedom should take into account the number of marginals used}
 \item{tailp}{the exact statistical significance}
 \item{sump}{sum of the probabilities used for error checking}
 \item{nenum}{the total number of tables enumerated}
}

\section{References}{
Yu C, Zelterman D (2001) Exact inference for family disease clusters. Commun Stat -- Theory
Meth 30:2293-2305

Yu C, Zelterman D (2002) Statistical inference for familial disease clusters. Biometrics
58:481-491
}
\seealso{
\code{\link[gap]{kin.morgan}}
}

\examples{
\dontrun{
# IPF among 203 siblings of 100 COPD patients from Liang KY, SL Zeger,
# Qaquish B. Multivariate regression analyses for categorical data
# (with discussion). J Roy Stat Soc B 1992, 54:3-40

# the degrees of freedom is 15
famtest<-c(
1, 0, 36,
1, 1, 12,
2, 0, 15,
2, 1,  7,
2, 2,  1,
3, 0,  5,
3, 1,  7,
3, 2,  3,
3, 3,  2,
4, 0,  3,
4, 1,  3,
4, 2,  1,
6, 0,  1,
6, 2,  1,
6, 3,  1,
6, 4,  1,
6, 6,  1)
test<-t(matrix(famtest,nrow=3))
famp<-pfc(test)
}
}
\author{Dani Zelterman, Jing Hua Zhao}
\note{Adapted from family.for by Dani Zelterman, 25/7/03}
\keyword{models}
