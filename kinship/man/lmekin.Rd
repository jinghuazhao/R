\name{lmekin}
\alias{lmekin}
\alias{print.lmekin}
\alias{getCovariateFormula2}
\alias{getCrossedTerms}
\alias{getGroupsFormula2}
\title{
  Linear Mixed Effects model using a kinship matrix.
}
\description{
A similar function to lme, but allowing for a complete specification of the
covariance matrix for the random effects.
}
\usage{
lmekin(fixed, data=parent.frame(), random,
    varlist=NULL, variance, sparse=c(20, .05),
    rescale=T, pdcheck=T,
    subset, weight, na.action)
}
\arguments{
\item{fixed}{
model statement for the fixed effects
}
\item{random}{
model statement for the random effects
}
\item{data}{
data frame containing the variables
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
\item{subset}{
selection of a subset of data
}
\item{weight}{
optional case weights
}
\item{na.action }{
the action for missing data values
}
}
\value{
an object of class 'lmekin', sharing similarities with both lm and lme
objects.
}
\details{
The lme function is designed to accept a prototype for the variance matrix
of the random effects, with the same prototype applying to all of the
groups in the data.  For familial genetic random effects, however, each
family has a different covariance pattern, necessitating the input of the
entire set of covariance matrices.  
In return, at present lmekin does not have the prototype abilities of lme.
}
\section{References}{
Pinheiro and Bates, Mixed Effect Models in S and Splus
}
\seealso{
\code{\link{print.lmekin}, lme}
}
\examples{
\dontrun{
#
# Make a kinship matrix for the entire study
# These two functions are NOT fast, the makekinship one in particular
#
cfam <- makefamid(main$gid, main$momid, main$dadid)
kmat <- makekinship(cfam, main$gid, main$momid, main$dadid)

# The kinship matrix for the females only: quite a bit smaller
#
kid <- dimnames(kmat)[[1]]
temp <- main$sex[match(kid, main$gid)] == 'F'
fkmat <- kmat[temp,temp]

# The dimnames on kmat are the gid value, which are necessary to match
#  the appropriate row/col of kmat to the analysis data set
# A look at \%dense tissue on a mammogram, with age at mammogram and
#  weight as covariates, and a familial random effect
#
fit <- lmekin(percdens ~ mammage + weight, data=anal1,
                 random = ~1|gid, varlist=list(fkmat))

Linear mixed-effects kinship model fit by maximum likelihood
  Data: anal1 
  Log-likelihood = -6093.917 
  n= 1535 

Fixed effects: percdens ~ mammage + weight 
 (Intercept)    mammage     weight 
     87.1593 -0.5333198 -0.1948871

Random effects:  ~ 1 | gid 
         Kinship Residual 
StdDev: 7.801603 10.26612
}
}
\keyword{regression}
% docclass is function
% Converted by Sd2Rd version 1.21.
