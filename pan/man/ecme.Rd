\name{ecme}
\alias{ecme}
\title{
ECME algorithm for general linear mixed model
}
\description{
Performs maximum-likelihood estimation for generalized linear
mixed models. The model, which is typically applied to
longitudinal or clustered responses, is

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

\usage{
ecme(y, subj, occ, pred, xcol, zcol=NULL, vmax, start, 
     maxits=1000, eps=0.0001, random.effects=F)

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
that y is in fact c(y1,y2,y3,y4) where length(y1)=2,
length(y2)=3, length(y3)=2, and length(y4)=7. Then subj should
be c(1,1,2,2,2,3,3,4,4,4,4,4,4,4).
}

\item{occ}{
vector of same length as y indicating the "occasions" for the
elements of y. In a longitudinal dataset where each individual
is measured on at most nmax distinct occasions, each element
of y corresponds to one subject-occasion, and the elements of
of occ should be coded as 1,2,\dots,nmax to indicate these
occasion labels. (You should label the occasions as
1,2,\dots,nmax even if they are not equally spaced in time; the
actual times of measurement will be incorporated into the
matrix "pred" below.) In a clustered dataset, the elements 
of occ label the units within each cluster i, using the labels
1,2,\dots,ni.
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
upon one another). If zcol=NULL then the model is assumed
to have no random effects; in that case the parameters are
estimated noniteratively by generalized least squares.
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

\item{start}{
optional starting values of the parameters. If this argument
is not given then ecme() chooses its own starting values.
This argument should be a list of three elements named
"beta", "psi", and "sigma2". Note that "beta" should be a
vector of the same length as "xcol", "psi" should be a
matrix of dimension c(length(zcol),length(zcol)), and
"sigma2" should be a scalar. This argument has no effect if
zcol=NULL.
}

\item{maxits}{
maximum number of cycles of ECME to be performed.
The algorithm runs to convergence or until "maxits"
iterations, whichever comes first.
}

\item{eps}{
convergence criterion. The algorithm is considered to have
converged if the relative differences in all parameters from
one iteration to the next are less than eps--that is, if
all(abs(new-old)<eps*abs(old)).
}

\item{random.effects}{
if TRUE, returns empirical Bayes estimates of all the random
effects bi (i=1,2,\dots,m) and their estimated covariance
matrices.
}}

\value{
a list containing estimates of beta, sigma2, psi, an estimated
covariance matrix for beta, the number of iterations actually
performed, an indicator of whether the algorithm converged,
and a vector of loglikelihood values at each iteration. If
random.effects=T, also returns a matrix of estimated random
effects (bhat) for individuals and an array of corresponding
covariance matrices. 

\item{beta}{
vector of same length as "xcol" containing estimated fixed
effects.
}

\item{sigma2}{
estimate of error variance sigma2.
}

\item{psi}{
matrix of dimension c(length(zcol),length(zcol)) containing
the estimated covariance matrix psi.
}

\item{converged}{
T if the algorithm converged, F if it did not
}

\item{iter}{
number of iterations actually performed. Will be equal
to "maxits" if converged=F.
}

\item{loglik}{
vector of length "iter" reporting the value of the
loglikelihood at each iteration.
}

\item{cov.beta}{
matrix of dimension c(length(xcol),length(xcol)) containing
estimated variances and covariances for elements of "beta".
}

\item{bhat}{
if random.effects=T, a matrix with length(zcol) rows and
m columns, where bhat[,i] is an empirical Bayes estimate
of bi.
}

\item{cov.b}{
if random.effects=T, an array of dimension length(zcol) by
length(zcol) by m, where cov.b[,,i] is an empirical Bayes
estimate of the covariance matrix associated with bi.
}}

\references{
Schafer JL (1997) Imputation of missing covariates under
a multivariate linear mixed model. Technical report 97-04, Dept. of
Statistics, The Pennsylvania State University, 

Schafer JL (2001). Multiple imputation with PAN. Chapter 12, pp357-77.
of New Methods for the Analysis of Change. Edited by Collins LM,
Sayer AG. American Psychological Association, Washington DC.

Schafer JL, Yucel RM (2002). Computational strategies for 
multivariate linear mixed-effects models with missing values. 
Journal of Computational and Graphical Statistics. 11:437-457
}

\examples{
########################################################################
# A simple linear model to these data using ecme(). This will be a
# traditional repeated-measures style additive model with a fixed effect
# for each column (occasion) and a random intercept for each subject.
#
# The data to be used is contained the object marijuana. Since the six
# measurements per subject were not clearly ordered in time, we consider
# a model that has an intercept and five dummy codes to allow the
# population means for the six occasions to be estimated freely together
# with an intercept randomly varied by subject. For a subject i with no
# missing values, the covariate matrices will be
#
#                   1 1 0 0 0 0              1
#                   1 0 1 0 0 0              1
#           Xi =    1 0 0 1 0 0       Zi =   1
#                   1 0 0 0 1 0              1
#                   1 0 0 0 0 1              1
#                   1 0 0 0 0 0              1
#
# When using ecme(), these are combined into a single matrix called
# pred. The pred matrix has length(y) rows. Each column of Xi and Zi
# must be represented in pred. Because Zi is merely the first column
# of Xi, we do not need to enter that column twice. So pred is simply
# the matrices Xi (i=1,...,9), stacked upon each other.
#
data(marijuana)
# we only use the complete data to illustrate
complete <- subset(marijuana,!is.na(y))
attach(complete)
pred <- with(complete,cbind(int,dummy1,dummy2,dummy3,dummy4,dummy5))
xcol <- 1:6
zcol <- 1
# Now we can fit the model.
result <- ecme(y,subj,occ,pred,xcol,zcol)
result

# Now we compare to lmer
if(require(lme4)) {
result <- lmer(y~-1+pred+(1|subj))
result
vcov(result)
detach(complete)
}
########################################################################
}

\keyword{models}
