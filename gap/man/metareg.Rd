\name{metareg}
\alias{metareg}
\title{Fixed and random effects model for meta-analysis}
\usage{metareg(data, N, verbose="Y", prefixb="b", prefixse="se")}
\description{
Given \eqn{k=n} studies with \eqn{b_1, ..., b_N} being \eqn{\beta}'s and
\eqn{se_1, ..., se_N} standard errors from regression, the fixed effects
model uses inverse variance weighting such that \eqn{w_1=1/se_1^2}, ...,
\eqn{w_N=1/se_N^2} and the combined \eqn{\beta} as the weighted average,
\eqn{\beta_f=(b_1*w_1+...+b_N*w_N)/w}, with \eqn{w=w_1+...+w_N} 
being the total weight, the se for this estimate is \eqn{se_f=\sqrt{1/w}}.
A normal z-statistic is obtained as \eqn{z_f=\beta_f/se_f}, and the
corresponding p value \eqn{p_f=2*pnorm(-abs(z_f))}. For the random effects
model, denote \eqn{q_w=w_1*(b_1-\beta_f)^2+...+w_N*(b_N-\beta_f)^2}
and \eqn{dl=max(0,(q_w-(k-1))/(w-(w_1^2+...+w_N^2)/w))}, corrected
weights are obtained such that \eqn{{w_1}_c=1/(1/w_1+dl)}, ..., 
\eqn{{w_N}_c=1/(1/w_N+dl)}, totaling \eqn{w_c={w_1}_c+...+{w_N}_c}.
The combined \eqn{\beta} and se are then \eqn{\beta_r=(b_1*{w_1}_c+...+b_N*{w_N}_c)/w_c}
and \eqn{se_r=\sqrt(1/w_c)}{se_r=sqrt(1/wc)}, leading to a z-statistic
\eqn{z_r=\beta_r/se_r} and a p-value \eqn{p_r=2*pnorm(-abs(z_r))}. Moreover, a
p-value testing for heterogeneity is \eqn{p_{heter}=pchisq(q_w,k-1,lower.tail=FALSE)}.
}

\arguments{
 \item{data}{Data frame to be used}
 \item{N}{Number of studies}
 \item{verbose}{A control for screen output}
 \item{prefixb}{Prefix of estimate; default value is "b"}
 \item{prefixse}{Prefix of standard error; default value is "se"}
The function accepts a wide format data with estimates as \eqn{b1,...,bN}
and standard errors as \eqn{se1,...,seN}. More generally, they can be specified
by prefixes in the function argument.
}

\value{
The returned value is a data frame with the following variables:
 \item{p_f}{P value (fixed effects model)}
 \item{p_r}{P value (random effects model)}
 \item{beta_f}{regression coefficient}
 \item{beta_r}{regression coefficient}
 \item{se_f}{standard error}
 \item{se_r}{standard error}
 \item{z_f}{z value}
 \item{z_r}{z value}
 \item{p_heter}{heterogeneity test p value}
 \item{i2}{\eqn{I^2}{I^2} statistic}
 \item{k}{No of tests used}
 \item{eps}{smallest double-precision number}
}
\section{References}{
JPT Higgins, SG Thompson, JJ Deeks, DG Altman. Measuring inconsistency in meta-analyses. BMJ 
327:557-60
}
\examples{
\dontrun{
abc <- data.frame(chromosome=1,rsn='abcd',startpos=1234,
                  b1=1,se1=2,p1=0.1,b2=2,se2=6,p2=0,b3=3,se3=8,p3=0.5)
metareg(abc,3)
abc2 <- data.frame(b1=c(1,2),se1=c(2,4),b2=c(2,3),se2=c(4,6),b3=c(3,4),se3=c(6,8))
print(metareg(abc2,3))
}
}
\author{Shengxu Li, Jing Hua Zhao}
\note{Adapted from a SAS macro}
\keyword{models}
