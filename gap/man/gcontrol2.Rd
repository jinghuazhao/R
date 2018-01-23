\name{gcontrol2}
\alias{gcontrol2}
\title{genomic control based on p values}
\usage{gcontrol2(p,col=palette()[4],lcol=palette()[2],...)}

\description{
The function obtains 1-df \eqn{\chi^2}{chi-squared} statistics (observed) according
to a vector of p values, and the inflation factor (lambda) according
to medians of the observed and expected statistics. The latter is based
on the empirical distribution function (EDF) of 1-df \eqn{\chi^2}{chi-squared} statstics.

It would be appropriate for genetic association analysis as of 1-df Armitage trend test
for case-control data; for 1-df additive model with continuous outcome one has to
consider the compatibility with p values based on z-/t- statistics.
}

\arguments{
\item{p}{a vector of observed p values}
\item{col}{colour for points in the Q-Q plot}
\item{lcol}{colour for the diagonal line in the Q-Q plot}
\item{...}{other options for plot}
}

\value{
A list containing:
\item{x}{the expected \eqn{\chi^2}{chi-squared} statistics}
\item{y}{the observed \eqn{\chi^2}{chi-squared} statistics}
\item{lambda}{the inflation factor}
}

\section{References}{
Devlin B, Roeder K (1999) Genomic control for association studies. 
Biometrics 55:997-1004
}

\examples{
\dontrun{
x2 <- rchisq(100,1,.1)
p <- pchisq(x2,1,lower.tail=FALSE)
r <- gcontrol2(p)
print(r$lambda)
}
}

\author{Jing Hua Zhao}

\keyword{models}
