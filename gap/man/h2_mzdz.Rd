\name{h2_mzdz}
\alias{h2_mzdz}
\title{Heritability estimation according to twin correlations}
\usage{h2_mzdz(mzDat=NULL,dzDat=NULL,rmz=NULL,rdz=NULL,nmz=NULL,ndz=NULL,selV=NULL)}

\description{
Heritability and variance estimation according to twin pair correlations.
}
\arguments{
 \item{mzDat}{a data frame for monzygotic twins (MZ)}
 \item{dzDat}{a data frame for dizygotic twins (DZ)}
 \item{rmz}{correlation for MZ twins}
 \item{rdz}{correlation for DZ twins}
 \item{nmz}{sample size for MZ twins}
 \item{ndz}{sample size for DZ twins}
 \item{selV}{names of variables for twin and cotwin}
}

\details{
The example section shows how to obtain bootstrap 95\%CI.
}

\value{
The returned value is a matrix containing heritability and their variance estimations
for "h2","c2","e2","vh","vc","ve".
}

\section{References}{
Keeping ES. Introduction to Statistical Inference, Dover Pulications, Inc. 1995
}

\examples{
\dontrun{

ACE_CI <- function(mzData,dzData,n.sim=5,selV=NULL,verbose=TRUE)
{
ACEr_twinData <- h2(mzDat=mzData,dzDat=dzData,selV=selV)
print(ACEr_twinData)

nmz <- dim(mzData)[1]
ndz <- dim(dzData)[1]
a <- ar <- vector()
set.seed(12345)
for(i in 1:n.sim)
{
  cat("\rRunning # ",i,"/", n.sim,"\r",sep="")
  sampled_mz <- sample(1:nmz, replace=TRUE)
  sampled_dz <- sample(1:ndz, replace=TRUE)
  mzDat <- mzData[sampled_mz,]
  dzDat <- dzData[sampled_dz,]
  ACEr_i <- h2(mzDat=mzDat,dzDat=dzDat,selV=selV)
  if(verbose) print(ACEr_i)
  ar <- rbind(ar,ACEr_i)
}
cat("\n\nheritability according to correlations\n\n")
ar <- as.data.frame(ar)
m <- mean(ar,na.rm=TRUE)
s <- sd(ar,na.rm=TRUE)
allr <- data.frame(mean=m,sd=s,lcl=m-1.96*s,ucl=m+1.96*s)
print(allr)
}

selVars <- c('bmi1','bmi2')

library(mvtnorm)
n.sim <- 500
cat ("\nThe first study\n\n")
mzm <- as.data.frame(rmvnorm(195, c(22.75,22.75),
                     matrix(2.66^2*c(1, 0.67, 0.67, 1), 2)))
dzm <- as.data.frame(rmvnorm(130, c(23.44,23.44),
                     matrix(2.75^2*c(1, 0.32, 0.32, 1), 2)))
mzw <- as.data.frame(rmvnorm(384, c(21.44,21.44),
                     matrix(3.08^2*c(1, 0.72, 0.72, 1), 2)))
dzw <- as.data.frame(rmvnorm(243, c(21.72,21.72),
                     matrix(3.12^2*c(1, 0.33, 0.33, 1), 2)))
names(mzm) <- names(dzm) <- names(mzw) <- names(dzw) <- c("bmi1","bmi2")
ACE_CI(mzm,dzm,n.sim,selV=selVars,verbose=FALSE)
ACE_CI(mzw,dzw,n.sim,selV=selVars,verbose=FALSE)

}
}
\author{Jing Hua Zhao}
\keyword{htest}
