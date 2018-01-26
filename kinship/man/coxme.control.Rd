\name{coxme.control}
\alias{coxme.control}
\alias{strata2}
\title{
  Control parameters for coxme
}
\description{
Set various control parameters for the \code{coxme}
function.
}
\usage{
coxme.control(eps=0.00001, toler.chol=.Machine$double.eps^0.75, 
toler.ms=.01, inner.iter=4, 
iter.max=10, simplex=0, lower=0, upper=Inf, sparse.calc=NULL )
}
\arguments{

\item{eps}{
convergence criteria for the inner Cox model computations.  
Iteration ceases when the relative change
in the log-likelihood is less than \code{eps}.
}
\item{toler.chol }{
tolerance that is used to detect singularity, i.e., redundant predictor
variables in the model, in the underlying Cholesky decompostion
routines.
}
\item{toler.ms }{
convergence criteria for the minimization of the integrated loglikelihood
over the variance parameters.
Since this ``outer'' iteration uses the Cox iteration as an inner loop,
and the Cox iteration in turn uses the cholesky decomposition as an inner
look, each of these treating the computations below it as if they were
exact, the cholesky tolerance should be tighter than the Cox tolerance,
which in turn should be tighter than that for the variance estimates.
}
\item{inner.iter }{
the number of iterations for the inner iteration loop.
}
\item{iter.max }{
maximum number of iterations for solution of a Cox partial
likelihood, given the values of the random effect variances.
Calls with \code{iter=0} are useful to
evaluate the likelihood for a prespecified parameter vector,
such as in the computation of a profile likelihood.
}
\item{simplex}{
number of iterations for the Nelder-Mead simplex algorithm.  
The simplex method is very good at finding the general neighborhood
of a minimum without getting lost, but can take a very large number
of iterations to narrow in on the final answer; 
opposite strentghs to the standard minimizer 
\code{optim}.
For hard problems, adding 50-100 iterations of the simplex as a
starting estimate for the usual method can be very helpful.
}
\item{lower, upper }{
limits for the variance parameters, used by optim.
}
\item{sparse.calc }{
style of computation for the inner likelihood code.  The results
of the two computations are identical, but can differ
in total compute time.  The optional calculation (calc=1)
uses somewhat more memory, but can be substantially faster when
the total number of random effects is of order n, the total
sample size.  The standard calculation (calc=0) is faster
when the number of random effects is small.
By default, the \code{coxme.fit} function chooses
the method dynamically.  
It may not always do so optimally.
}
}
\value{
a list containing values for each option.
}
\details{
The central computation consists of an outer maximization to
determine the variances of the random effects, performed by
the \code{optim} function.  
Each evaluation for optim, however, itself requires the
solution of a minimization problem; this is the inner loop.
It is important that the inner loop use a fixed number of
iterations, but it is not yet clear what is the minimal
sufficient number for that inner loop.  Making this number
smaller will make the routine faster, but perhaps at the
expense of accuracy.
}
\section{References}{
Therneau and Pankratz
}
\seealso{
\code{\link{coxme}}
}
\keyword{survival}
% docclass is function
% Converted by Sd2Rd version 1.21.
