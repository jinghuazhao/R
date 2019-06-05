\name{AE3}
\alias{AE3}
\title{AE model using nuclear family trios}
\usage{AE3(model, random, data, seed=1234, n.sim=50000, verbose=TRUE)}
\description{
This function is adapted from example 7.1 of Rabe-Hesketh et al. (2008). It
also procides heritability estimate and confidence intervals.
}

\arguments{
 \item{model}{a linear mixed model formula, see example below}
 \item{random}{random effect, see exampe below}
 \item{data}{data to be analyzed}
 \item{seed}{random number seed}
 \item{n.sim}{number of simulations}
 \item{verbose}{a flag for printing out results}
}

\value{
The returned value is a list containing:

 \item{lme.result}{the linear mixed model result}
 \item{h2}{the heritability estimate}
 \item{CL}{confidence intervals}

}

\section{References}{
Rabe-Hesketh S, Skrondal A, Gjessing HK. Biometrical modeling of twin and family data using standard mixed model 
software. Biometrics 2008, 64:280-288
}

\examples{
\dontrun{
require(gap.datasets)
AE3(bwt ~ male + first + midage + highage + birthyr,
    list(familyid = pdIdent(~var1 + var2 + var3 -1)), mfblong)
}
}
\author{Jing Hua Zhao}
\note{Adapted from f.mbf.R from the paper}
\keyword{models}
\keyword{htest}
