\name{pan.bd}
\alias{pan.bd}

\title{
Imputation of multivariate panel or cluster data
}

\description{
Implementation of pan() that restricts the covariance matrix
for the random effects to be block-diagonal. This function
is identical to pan() in every way except that psi is now 
characterized by a set of r matrices of dimension q x q.
}

\usage{
pan.bd(y, subj, pred, xcol, zcol, prior, seed, iter=1, start)
}

\arguments{
\item{y}{
See description for pan().
}

\item{subj}{
See description for pan().
}

\item{pred}{
See description for pan().
}

\item{xcol}{
See description for pan().
}

\item{zcol}{
See description for pan().
}

\item{prior}{
Same as for pan() except that the hyperparameters for psi
have new dimensions. The hyperparameter c is now a vector of
length r, where c[j] contains the prior degrees of freedom for
the jth block portion of psi (j=1,\dots,r). The hyperparameter
Dinv is now an array of dimension c(q,q,r), where Dinv[,,j]
contains the prior scale matrix for the jth block portion of
psi (j=1,\dots,r).
}

\item{seed}{
See description for pan().
}

\item{iter}{
See description for pan().
}

\item{start}{
See description for pan().
}}

\value{
A list with the same components as that from pan(), with two
minor differences: the dimension of "psi" is now (q x q x r x
"iter"), and the dimension of "last\$psi" is now (q x q x r).
}

\keyword{models}

