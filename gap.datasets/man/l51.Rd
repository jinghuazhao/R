\name{l51}
\docType{data}
\alias{l51}
\title{An example pedigree data}
\description{
The data contains data on 51 individuals in a pedigree. Below it is used for comparing results from various packages.
}
\usage{data(l51)}
\format{A data frame}
\source{Morgan v3.}
\references{
Morgan v3. \url{http://www.stat.washington.edu/thompson/Genepi/MORGAN/Morgan.shtml}
} 
\examples{\dontrun{
km <- kin.morgan(l51)
k2 <- km$kin.matrix*2

# quantitative trait
library(regress)
r <- regress(qt ~ 1, ~k2, data=l51)
names(r)
r
# qualitative trait
N <- dim(l51)[1]
w <- with(l51,quantile(qt,probs=0.75,na.rm=TRUE))
ped51 <- within(l51, bt <- ifelse(qt<=w,0,1))
d <- regress(bt ~ 1, ~k2, data=ped51)
d
# for other tests not shown here
set.seed(12345)
ped51 <- within(ped51,{r <- rnorm(N); bt[is.na(bt)] <- 0})
library(foreign)
write.dta(ped51,"ped51.dta")
}
}
\keyword{datasets}
