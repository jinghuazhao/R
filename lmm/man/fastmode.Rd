\name{fastmode}
\alias{fastmode}
\title{
Rapidly converging algorithm for calculating posterior modes
in linear mixed models
}
\description{
Computes the marginal posterior mode of the variance
parameters in linear mixed models using a rapidly converging
procedure described by Schafer (1998), which combines Fisher
scoring with an ECME algorithm. The method is a minor
modification of the restricted maximum-likelihood (RML)
procedure used in "fastrml". The model is identical to
that of "fastrml" with the addition of prior
distributions for the variance parameters.

For a description of the prior distribution, see the "Details"
section below.
}

\usage{
fastmode(y, subj, pred, xcol, zcol, prior, vmax, occ, start,
   maxits=100, eps=0.0001)
}

\arguments{
Identical to those for the function "fastrml", with one
additional required argument:
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
maximum number of cycles to be performed.
The algorithm runs to convergence or until "maxits"
iterations, whichever comes first.
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
effects. This estimate represents the posterior mean for beta,
conditional upon the estimated values of the variance
parameters sigma2 and psi.
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
if they do not increase the log-posterior density.
}
\item{logpost}{
vector of length "iter" reporting the value of the
log-posterior density at each iteration.
}
\item{cov.beta}{
matrix of dimension c(length(xcol),length(xcol)) containing
estimated variances and covariances for elements of "beta".
These are conventional estimates which regard the variance
parameters (sigma2 and psi) as fixed at their estimated
values. 
}}

\details{
The algorithm is described in the appendix of Schafer (1998).
Scoring is carried out on log(sigma2) and the nonredundant
elements of the inverse of psi/sigma2, taking logs of the
diagonal elements. Upon convergence, the estimates  represent
the mode of the joint posterior density of 1/sigma2 and the
inverse of psi, marginalized (i.e. integrated) over beta.

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
\code{\link[lmm]{mgibbs}}, \code{\link[lmm]{fastmcmc}},
\code{\link[lmm]{example}}
}

\examples{
\dontrun{
For a detailed example, see the file "example.R" distributed
with this library.

}}

\keyword{models}
