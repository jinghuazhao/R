\name{mgibbs}
\alias{mgibbs}
\title{
Modified Gibbs sampler for Bayesian inference in linear mixed
models
}

\description{
Simulates posterior draws of parameters in linear mixed models
using a Markov chain Monte Carlo (MCMC) procedure, the
modified Gibbs sampler described by Schafer (1998). This algorithm
may be slow, requiring a large number of cycles to achieve
stationarity. In most cases, "fastmcmc" will perform
better. This algorithm is provided mainly for comparison
against "fastmcmc".

For a description of the model and the prior distribution, see
the "Details" section below.
}
\usage{
mgibbs(y, subj, pred, xcol, zcol, prior, seed, vmax, occ,
   start, iter=1000)
}
\arguments{
\item{y}{
vector of responses. This is simply the individual yi vectors
stacked upon one another. Each element of y represents the 
observed response for a particular subject-occasion, or for a
particular unit within a cluster.
}
\item{subj}{
vector of same length as y, giving the subject (or cluster)
indicators i for the elements of y. For example, suppose 
that y is c(y1,y2,y3,y4) where length(y1)=2, length(y2)=3,
length(y3)=2, and length(y4)=7. Then subj should be
c(1,1,2,2,2,3,3,4,4,4,4,4,4,4).
}
\item{pred}{
matrix of covariates used to predict y. The number of rows
should be length(y). The first column will typically be
constant (one), and the remaining columns correspond to other
variables appearing in Xi and Zi.
}
\item{xcol}{
vector of integers indicating which columns of pred will be
used in Xi. That is, pred[,xcol] is the Xi matrices (stacked
upon one another.
}
\item{zcol}{
vector of integers indicating which columns of pred will be
used in Zi. That is, pred[,zcol] is the Zi matrices (stacked
upon one another).
}
\item{prior}{
A list with four components specifying the hyperparameters of
the prior distribution applied to sigma2 and psi. The
components must be named "a", "b", "c", and "Dinv". All are
scalars except for "Dinv", which is a matrix of dimension
c(length(zcol),length(zcol)).
}
\item{seed}{
Seed for random number generator. This should be a positive
integer. 
}
\item{vmax}{
optional matrix of dimension c(max(occ),max(occ)) from which
the Vi matrices will be extracted. In a longitudinal dataset, 
vmax would represent the Vi matrix for an individual with
responses at all possible occasions 1,2,\dots,nmax=max(occ);
for individuals with responses at only a subset of these
occasions, the Vi will be obtained by extracting the rows
and columns of vmax for those occasions. If no vmax is
specified by the user, an identity matrix is used. In most
applications of this model one will want to have Vi =
identity, so most of the time this argument can be omitted.
}
\item{occ}{
vector of same length as y indicating the "occasions" for the
elements of y. This argument is relevant only if a
non-identity vmax is specified. In a longitudinal dataset
where each individual is measured on at most nmax distinct
occasions, each element of y corresponds to one
subject-occasion, and the elements of occ should be coded as
1,2,\dots,nmax to indicate these occasion labels. (You should
label the occasions as 1,2,\dots,nmax even if they are not
equally spaced in time; the actual times of measurement will
be incorporated into the matrix "pred".) 
}
\item{start}{
optional starting values of the parameters. If this argument
is not given then the function chooses its own starting values.
This argument should be a list of three elements named
"beta", "psi", and "sigma2". Note that "beta" should be a
vector of the same length as "xcol", "psi" should be a
matrix of dimension c(length(zcol),length(zcol)), and
"sigma2" should be a scalar.
}
\item{iter}{
number of cycles of the modified Gibbs sampler to be
performed. 
}}

\value{
a list containing the following components.

\item{beta}{
simulated value of coefficients beta after "iter" cycles of
the modified Gibbs sampler. This is a vector of the same
length as xcol.
}
\item{sigma2}{
simulated value of the residual variance sigma2 after "iter"
cycles of the modified Gibbs sampler.
}
\item{psi}{
simulated value of the between-unit covariance matrix psi
after "iter" cycles of the modified Gibbs sampler.
}
\item{sigma2.series}{
vector of length "iter" containing the entire history of
simulated values of sigma2. That is, sigma2.series[t] contains
the value of sigma2 at cycle t.
}
\item{psi.series}{
array of dimension c(length(zcol),length(zcol),iter)
containing the entire history of simulated values of psi. That
is, psi.series[,,t] contains the value of psi at cycle t.
}}

\details{
The algorithm is described in Section 5 of Schafer (1998).

The model, which is typically applied to longitudinal or
clustered responses, is 

	yi = Xi\%*\%beta + Zi\%*\%bi + ei ,    i=1,\dots,m,

where

	yi    = (ni x 1) response vector for subject
		or cluster i;
	Xi    = (ni x p) matrix of covariates;
	Zi    =	(ni x q) matrix of covariates;
	beta  = (p x 1) vector of coefficients common to the
		population (fixed effects);
	bi    = (q x 1) vector of coefficients specific to
		subject or cluster i (random effects); and
	ei    = (ni x 1) vector of residual errors.

The vector bi is assumed to be normally distributed
with mean zero and unstructured covariance matrix psi, 

        bi \eqn{\sim}{~} N(0,psi) independently for i=1,\dots,m.

The residual vector ei is assumed to be

        ei \eqn{\sim}{~} N(0,sigma2*Vi)

where Vi is a known (ni x ni) matrix. In most applications,
Vi is the identity matrix.

The prior distribution applied to the within-unit residual
variance is scaled inverted-chisquare,

                sigma2 \eqn{\sim}{~} a / chisq(b),

where chisq(b) denotes a chisquare random variable with b
degrees of freedom, and a and b are user-defined
hyperparameters. Values for the hyperparmeters may be chosen
by regarding a/b as a rough prior guess for sigma2, and as
the imaginary degrees of freedom on which this guess is based.

The prior distribution applied to the between-unit covariance
matrix is inverted Wishart,

                   psiinv \eqn{\sim}{~} W(c,D),

where psiinv is the inverse of the between-unit covariance
matrix psi, and W(c,D) denotes a Wishart distribution with
degrees of freedom c and scale matrix D. Values for the
hyperparameters may be chosen by regarding Dinv/c (the inverse
of D divided by c) as a rough prior guess for psi, and c as
the imaginary degrees of freedom on which this guess is based.

An improper uniform prior density function is applied to
the fixed effects beta.
}

\references{
Schafer, J.L. (1998) Some improved procedures for linear mixed
models. Submitted to Journal of the American Statistical
Association.
}

\seealso{
\code{\link[lmm]{ecmeml}}, \code{\link[lmm]{ecmerml}},
\code{\link[lmm]{fastml}}, \code{\link[lmm]{fastrml}},
\code{\link[lmm]{fastmode}}, \code{\link[lmm]{fastmcmc}},
\code{\link[lmm]{example}}
}

\examples{
\dontrun{
For a detailed example, see the file "example.R" distributed
with this library.

}}

\keyword{models}
