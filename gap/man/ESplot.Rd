\name{ESplot}
\alias{ESplot}
\title{Effect-size plot}
\usage{ESplot(ESdat,SE=TRUE,logscale=TRUE,alpha=0.05,xlim=c(-2,8),v=1,...)}

\description{
The function accepts parameter estimates and their standard
errors for a range of models.
}

\arguments{
\item{ESdat}{A data frame consisting of model id, parameter estimates and standard errors or confidence limits}
\item{SE}{If TRUE, the third column of ESdata contains the standard error estimates}
\item{logscale}{If TRUE, indicates log-scale as appropriate for odds ratio}
\item{alpha}{Type-I error rate used to construct 100(1-alpha) confidence interval}
\item{xlim}{Lower and upper limits of the horizontal axis, roughly corresponding to confidence limits}
\item{\dots}{Other options for \code{plot}}
\item{v}{Location of the vertical line}
}

\examples{
\dontrun{
# 7-4-2008 MRC-Epid JHZ
options(stringsAsFactors=FALSE)
testdata <- data.frame(
              models=c("Basic model","Adjusted","Moderately adjusted","Heavily adjusted","Other"),
               logOR=log(c(4.5,3.5,2.5,1.5,1)),
             SElogOR=c(0.2,0.1,0.2,0.3,0.2)
            )
ESplot(testdata,v=1)
title("A fictitious plot")
# Outcomes A2, B2, C2 in three columns
# par(mfrow=c(1,3))
# ESplot(snp_effects[c("snpid","A2_b2","A2_se2")], lty=2, xlim=c(0.7,1.4))
# snp_effects["snpid"] <- ""
# ESplot(snp_effects[c("snpid","B2_b2","B2_se2")], lty=2, xlim=c(0.7,1.4))
# ESplot(snp_effects[c("snpid","C2_b2","C2_se2")], lty=2, xlim=c(0.7,1.4))
#
# Quantitative trait, as appropriate for linear regression
# testdata <- data.frame(modelid, beta, se(beta))
# ESplot(testdata, logscale=FALSE)
#
# Other scenarios
# OR with CI
# ESplot(testdata,SE=FALSE)
}
}

\author{Jing Hua Zhao}
\keyword{hplot}
