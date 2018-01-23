\name{pfc.sim}
\alias{pfc.sim}
\title{Probability of familial clustering of disease}
\usage{pfc.sim(famdata,n.sim=1000000,n.loop=1)}
\description{
To calculate probability of familial clustering of disease using Monte Carlo simulation
}
\arguments{
\item{famdata}{collective information of sib size, number of affected sibs and their frequencies}
\item{n.sim}{number of simulations in a single Monte Carlo run}
\item{n.loop}{total number of Monte Carlo runs}
}

\value{
The returned value is a list containing:

\item{n.sim}{a copy of the number of simulations in a single Monte Carlo run}
\item{n.loop}{the total number of Monte Carlo runs}
\item{p}{the observed p value}
\item{tailpl}{accumulated probabilities at the lower tails}
\item{tailpu}{simulated p values}
}

\section{References}{
Yu C and D Zelterman (2001) Exact inference for family disease clusters. Commun Stat -- Theory
Meth 30:2293-2305

}
\seealso{
\code{\link[gap]{pfc}}
}

\examples{
\dontrun{
# Li FP, Fraumeni JF Jr, Mulvihill JJ, Blattner WA, Dreyfus MG, Tucker MA,
# Miller RW. A cancer family syndrome in twenty-four kindreds.
# Cancer Res 1988, 48(18):5358-62. 

# family_size  #_of_affected frequency

famtest<-c(
1, 0, 2,
1, 1, 0,
2, 0, 1,
2, 1, 4,
2, 2, 3,
3, 0, 0,
3, 1, 2,
3, 2, 1,
3, 3, 1,
4, 0, 0,
4, 1, 2,
5, 0, 0,
5, 1, 1,
6, 0, 0,
6, 1, 1,
7, 0, 0,
7, 1, 1,
8, 0, 0,
8, 1, 1,
8, 2, 1,
8, 3, 1,
9, 3, 1)

test<-matrix(famtest,byrow=T,ncol=3)

famp<-pfc.sim(test)
}
}
\author{Chang Yu, Dani Zelterman}
\note{Adapted from runi.for from Change Yu, 5/6/4}
\keyword{models}
