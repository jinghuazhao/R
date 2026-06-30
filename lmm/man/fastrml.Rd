\name{fastrml}
\alias{fastrml}
\title{
Rapidly converging algorithm for restricted maximum-likelihood
(RML) estimation in linear mixed models
}

\description{
Computes RML estimates of parameters in linear mixed models
using the rapidly converging procedure described by Schafer
(1998), which combines Fisher scoring with an ECME algorithm.

For a description of the model, see the "Details" section
below. 
}

\usage{
fastrml(y, subj, pred, xcol, zcol, vmax, occ, start,
   maxits=50, eps=0.0001)
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
\item{maxits}{
maximum number of cycles to be performed. The algorithm runs
to convergence or until "maxits" iterations, whichever comes
first.
}
\item{eps}{
convergence criterion. The algorithm is considered to have
converged if the relative differences in all parameters from
one iteration to the next are less than eps--that is, if
all(abs(new-old)<eps*abs(old)).
}}

\value{
a list containing the following components.

\item{beta}{
vector of same length as "xcol" containing estimated fixed
effects.
}
\item{sigma2}{
estimate of residual error variance.
}
\item{psi}{
matrix of dimension c(length(zcol),length(zcol)) containing
estimated variances and covariances of the random effects.
}
\item{converged}{
T if the algorithm converged, F if it did not.
}
\item{iter}{
number of iterations actually performed. Will be equal
to "maxits" if converged=F.
}
\item{reject}{
a logical vector of length iter indicating, for each
iteration, whether the scoring estimates were rejected
and replaced by ECME estimates (T), or whether the scoring
estimates were accepted (F). Scoring estimates are rejected
if they do not increase the loglikelihood.
}
\item{loglik}{
vector of length "iter" reporting the value of the
loglikelihood at each iteration.
}
\item{cov.beta}{
matrix of dimension c(length(xcol),length(xcol)) containing
estimated variances and covariances for elements of "beta".
These are conventional estimates which regard the variance
parameters (sigma2 and psi) as fixed at their RML estimates.
}
\item{b.hat}{
a matrix with length(zcol) rows and m columns, where
b.hat[,i] is an empirical Bayes estimate of bi.
}
\item{cov.b}{
an array of dimension length(zcol) by length(zcol) by m,
where cov.b[,,i] is an empirical Bayes estimate of the
covariance matrix associated with bi. These are
conventional estimates which regard the variance
parameters (sigma2 and psi) as fixed at their RML estimates.
}
\item{cov.beta.new}{
matrix of dimension c(length(xcol),length(xcol)) containing
estimated variances and covariances for elements of "beta".
These are improved estimates which account for uncertainty
in estimating the variance parameters (sigma2 and psi).
}
\item{cov.b.new}{
an array of dimension length(zcol) by length(zcol) by m,
where cov.b.new[,,i] is an estimated covariance matrix 
for bi. These are improved estimates which account for
uncertainty in estimating the variance parameters (sigma2
and psi).
}
\item{cov.b.beta.new}{
an array of dimension length(zcol) by length(xcol) by m,
where cov.b.beta.new[,,i] contains the estimated covariances
bewteen bi and beta. These are improved estimates which account
for uncertainty in estimating the variance parameters (sigma2
and psi). Note that conventional estimates which regard
sigma2 and psi as fixed assume zero covariance between
each bi and beta.
}}

\details{
A full description of the algorithm is given in Section 3 of
Schafer (1998). Scoring is carried out on log(sigma2) and the
nonredundant elements of the inverse of psi/sigma2, taking
logs of the diagonal elements. Improved estimates of variances
and covariances are described in Section 4.

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

        bi \eqn{\sim}{~}  N(0,psi) independently for i=1,\dots,m.

The residual vector ei is assumed to be

        ei \eqn{\sim}{~} N(0,sigma2*Vi)

where Vi is a known (ni x ni) matrix. In most applications,
Vi is the identity matrix.
}

\references{
Schafer, J.L. (1998) Some improved procedures for linear mixed
models. Submitted to Journal of the American Statistical
Association.
}

\seealso{
\code{\link[lmm]{ecmeml}}, \code{\link[lmm]{ecmerml}},
\code{\link[lmm]{fastml}}, \code{\link[lmm]{fastmode}},
\code{\link[lmm]{mgibbs}}, \code{\link[lmm]{fastmcmc}},
\code{\link[lmm]{example}}
}

\examples{
\dontrun{
For a detailed example, see the file "example.R" distributed
with this library.

}}

\keyword{models}
