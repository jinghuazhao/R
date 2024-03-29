% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chow.test.R
\name{chow.test}
\alias{chow.test}
\title{Chow's test for heterogeneity in two regressions}
\usage{
chow.test(y1, x1, y2, x2, x = NULL)
}
\arguments{
\item{y1}{a vector of dependent variable.}

\item{x1}{a matrix of independent variables.}

\item{y2}{a vector of dependent variable.}

\item{x2}{a matrix of independent variables.}

\item{x}{a known matrix of independent variables.}
}
\value{
The returned value is a vector containing (please use subscript to access them):
\itemize{
\item F the F statistic.
\item df1 the numerator degree(s) of freedom.
\item df2 the denominator degree(s) of freedom.
\item p the p value for the F test.
}
}
\description{
Chow's test for heterogeneity in two regressions
}
\details{
Chow's test is for differences between two or more regressions.  Assuming that
errors in regressions 1 and 2 are normally distributed with zero mean and
homoscedastic variance, and they are independent of each other, the test of
regressions from sample sizes \eqn{n_1} and \eqn{n_2} is then carried out using
the following steps.  1.  Run a regression on the combined sample with size
\eqn{n=n_1+n_2} and obtain within group sum of squares called \eqn{S_1}.  The
number of degrees of freedom is \eqn{n_1+n_2-k}, with \eqn{k} being the number
of parameters estimated, including the intercept.  2.  Run two regressions on
the two individual samples with sizes \eqn{n_1} and \eqn{n_2}, and obtain their
within group sums of square \eqn{S_2+S_3}, with \eqn{n_1+n_2-2k} degrees of
freedom.  3.  Conduct an \eqn{F_{(k,n_1+n_2-2k)}} test defined by \deqn{F =
\frac{[S_1-(S_2+S_3)]/k}{[(S_2+S_3)/(n_1+n_2-2k)]}} If the \eqn{F} statistic
exceeds the critical \eqn{F}, we reject the null hypothesis that the two
regressions are equal.

In the case of haplotype trend regression, haplotype frequencies from combined
data are known, so can be directly used.
}
\note{
adapted from chow.R.
}
\examples{
\dontrun{
dat1 <- matrix(c(
     1.2, 1.9, 0.9,
     1.6, 2.7, 1.3,
     3.5, 3.7, 2.0,
     4.0, 3.1, 1.8,
     5.6, 3.5, 2.2,
     5.7, 7.5, 3.5,
     6.7, 1.2, 1.9,
     7.5, 3.7, 2.7,
     8.5, 0.6, 2.1,
     9.7, 5.1, 3.6), byrow=TRUE, ncol=3)

dat2 <- matrix(c(
     1.4, 1.3, 0.5,
     1.5, 2.3, 1.3,
     3.1, 3.2, 2.5,
     4.4, 3.6, 1.1,
     5.1, 3.1, 2.8,
     5.2, 7.3, 3.3,
     6.5, 1.5, 1.3,
     7.8, 3.2, 2.2,
     8.1, 0.1, 2.8,
     9.5, 5.6, 3.9), byrow=TRUE, ncol=3)

y1<-dat1[,3]
y2<-dat2[,3]
x1<-dat1[,1:2]
x2<-dat2[,1:2]
chow.test.r<-chow.test(y1,x1,y2,x2)
# from http://aoki2.si.gunma-u.ac.jp/R/
}

}
\references{
\insertRef{chow60}{gap}
}
\seealso{
\code{\link{htr}}
}
\author{
Shigenobu Aoki, Jing Hua Zhao
}
\keyword{htest}
