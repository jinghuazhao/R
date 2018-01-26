\name{coxme}
\alias{coxme}
\alias{coxme.fit}
\alias{coxme.varcheck}
\alias{print.coxme}
\title{
  Fit a mixed-effects Cox model
}
\description{
Returns an object of class \code{coxme} representing the
fitted model.
}
\usage{
coxme(fixed, data=parent.frame(), random, 
   weights, subset, na.action, init, control, 
   ties=c("efron", "breslow", "exact"), singular.ok=T, 
   varlist, variance, vinit=.2, sparse=c(50, .02), rescale=T, pdcheck=T, x=F, y=T, shortlabel=T, \dots)
}
\arguments{
\item{fixed}{
formula decribing the fixed effects part of the model.
}
\item{data}{
a data frame containing the variables.
}
\item{random}{
a one-sided formula describing the random effects part of the model.
}
\item{weights}{
case weights for each observation
}
\item{subset}{
an expression describing the subset of the data that should be used in the fit.
}
\item{na.action }{
a function giving the default action on encountering missing values.  
It is more usual to use the global na.action system option to control this.
}
\item{init}{
initial values for the coefficients for the fixed portion of the model, or 
the frailties followed by the fixed effect coefficients.
}
\item{control}{
the result of a call to \code{coxme.control}
}
\item{ties}{
the approximation to be used for tied death times: either "efron" or
"breslow"
}
\item{singular.ok }{
if TRUE, then redundant coefficients among the fixed effects are set to
NA, if FALSE the program will fail with an error message if there are
redundant variables.
}
\item{varlist}{
variance specifications, often of class \code{bdsmatrix},
decsribing the variance/covariance structure of one or more of the random
effects.  
}
\item{variance}{
fixed values for the variances of selected random effects.
Values of 0 indicate that the final value should be solved for.
}
\item{vinit}{
vector of initial values for variance terms.
It is necessary that the initial variance matrix be symmetric positive
definite.  Normally, a simple sum of the \code{varlist} matrices will
suffice, i.e., a vector of 1s; but not always.
}
\item{sparse}{
determines which levels of random effects factor variables, if any, 
for which the program will use sparse matrix techniques.
If a grouping variable has less than sparse[1] levels, then sparse methods
are not used for that variable.
If it has greater than or equal to sparse[1] unique levels, sparse methods
will be used for those values which represent less than sparse[2] as a 
proportion of the data. 
For instance, if a grouping variable has 4000 levels, but 40\% of the
subjects are in group 1 then 3999 of the levels will be represented sparsely
in the variance matrix. 
A single logical value of F is equivalent to setting sparse[1] to infinity.
}
\item{rescale}{
scale any user supplied variance matrices so as to have a diagonal of 1.0.
}
\item{pdcheck}{
verify that any user-supplied variance matrix is positive definite (SPD).
It has been observed that IBD matrices produced by some software are not
strictly SPD.  Sometimes models with these matrices still work (throughout
the iteration path, the weighted sum of variance matrices was always SPD)
and sometimes they don't.  In the latter case, messages about taking
the log of negative numbers will occur, and the results of the fit are
not necessarily trustworthy.
}
\item{x}{
retain the X matrix in the output.
}
\item{y}{
retain the dependent variable (a Surv object) in the output.
}
\item{shortlabel}{
no comment(s)
}
\item{...}{any other arguments}
}
\value{
an object of class coxme
}
\seealso{
\code{coxph}
}
\examples{
\dontrun{
coxme(Surv(time, status) ~ rx, data=rats, random= ~1|litter)

Cox mixed-effects kinship model fit by maximum likelihood
  Data: rats 
  n= 150 
                    NULL Integrated Penalized 
Log-likelihood -185.6556   -180.849  -173.774

  Penalized loglik: chisq= 23.76 on 13.17 degrees of freedom, p= 0.036 
 Integrated loglik: chisq= 9.61 on 2 degrees of freedom, p= 0.0082 

Fixed effects: Surv(time, status) ~ rx 
        coef exp(coef)  se(coef)        z         p 
rx 0.9132825  2.492491 0.3226856 2.830255 0.0046511

Random effects:  ~ 1 | litter 
             litter 
Variance: 0.4255484
}
}
\keyword{survival}
% docclass is function
% Converted by Sd2Rd version 1.21.
