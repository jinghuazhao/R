\name{ab}
\alias{ab}
\title{Test/Power calculation for mediating effect}
\usage{ab(type,n=25000,a=0.15,sa=0.01,b=log(1.19),sb=0.01,alpha=0.05,fold=1)}
\description{
This function tests for or obtains power of mediating effect based on estimates of
two regression coefficients and their standard errors. Note that for binary outcome
or mediator, one should use log-odds ratio and its standard error.
}
\arguments{
  \item{type}{string option: "test", "power"}
  \item{n}{default sample size to be used for power calculation}
  \item{a}{regression coefficient from indepdendent variable to mediator}
  \item{sa}{SE(a)}
  \item{b}{regression coefficient from mediator variable to outcome}
  \item{sb}{SE(b)}
  \item{alpha}{size of siginficance test for power calculation}
  \item{fold}{fold change for power calculation, as appropriate for a range of sample sizes}
}

\value{
The returned value are z-test and significance level for significant testing or sample size/power for a
given fold change of the default sample size.
}

\section{References}{
Freathy RM, Timpson NJ, Lawlor DA, Pouta A, Ben-Shlomo Y, Ruokonen A, Ebrahim S, Shields B, Zeggini E, Weedon MN, 
Lindgren CM, Lango H, Melzer D, Ferrucci L, Paolisso G, Neville MJ, Karpe F, Palmer CN, Morris AD, Elliott P, Jarvelin 
MR, Smith GD, McCarthy MI, Hattersley AT, Frayling TM. Common variation in the FTO Gene alters diabetes-related 
metabolic traits to the extent expected given its effect on BMI. Diabetes 57:1419-1426, 2008.

Kline RB. Principles and practice of structural equation modeling, Second Edition. The Guilford Press 2005.

MacKinnon DP. Introduction to Statistical Mediation Analysis. Taylor & Francis Group 2008.

Preacher KJ, Leonardelli GJ. Calculation for the Sobel Test-An interactive calculation tool for mediation tests
http://www.people.ku.edu/~preacher/sobel/sobel.htm
}
\seealso{
\code{\link[gap]{ccsize}}
}
\examples{
\dontrun{

ab()
n <- power <- vector()
for (j in 1:10)
{
   z <- ab(fold=j*0.01)
   n[j] <- z[1]
   power[j] <- z[2]
}
plot(n,power,xlab="Sample size",ylab="Power")
title("SNP-BMI-T2D association in EPIC-Norfolk study")
}
}
\author{Jing Hua Zhao}
\keyword{htest}
