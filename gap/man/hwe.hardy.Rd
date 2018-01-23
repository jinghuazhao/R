\name{hwe.hardy}

\alias{hwe.hardy}

\title{Hardy-Weinberg equilibrium test using MCMC}

\usage{hwe.hardy(a, alleles = 3, seed = 3000, sample = c(1000, 1000, 5000))}

\description{Hardy-Weinberg equilibrium test by MCMC}

\arguments{
%  \item{x}{a genotype object, as genotype.}
  \item{a}{an array containing the genotype counts, as integer.}
  \item{alleles}{number of allele at the locus, greater than or equal to 3, 
                 as integer}
  \item{seed}{pseudo-random number seed, as integer.}
  \item{sample}{optional, parameters for MCMC containing number of chunks, 
                size of a chunk and burn-in steps, as integer.}
%  If \code{x} is given, \code{a} and \code{alleles} are prepaired on the fly.
%  Look in examples.
}

\author{Sun-Wei Guo, Jing Hua Zhao, Gregor Gorjanc}

\source{
  http://www.stat.washington.edu/thompson/Genepi/pangaea.shtml,
}

\value{
  The returned value is a list containing:
  \item{method}{Hardy-Weinberg equilibrium test using MCMC}  
  \item{data.name}{name of used data if \code{x} is given}  
  \item{p.value}{Monte Carlo p value}
  \item{p.value.se}{standard error of Monte Carlo p value}
  \item{switches}{percentage of switches (partial, full and altogether)}
}

\section{References}{
  Guo, S.-W. and E. A. Thompson (1992) Performing the exact test of
  Hardy-Weinberg proportion for multiple alleles. Biometrics. 48:361--372.
}

\seealso{\code{\link[gap]{hwe}}, 
         \code{\link[genetics]{HWE.test}},
         \code{\link[genetics]{genotype}}}

\examples{
  \dontrun{
    # example 2 from hwe.doc:
    a<-c(
    3,
    4, 2,
    2, 2, 2,
    3, 3, 2, 1,
    0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1,
    0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 2, 1, 0, 0, 0)
    ex2 <- hwe.hardy(a=a,alleles=8)

    # example using HLA
    data(hla)
    x <- hla[,3:4]
    y <- pgc(x,handle.miss=0,with.id=1)
    n.alleles <- max(x,na.rm=TRUE)
    z <- vector("numeric",n.alleles*(n.alleles+1)/2)
    z[y$idsave] <- y$wt
    hwe.hardy(a=z,alleles=n.alleles)

    # with use of class 'genotype'
    # this is to be fixed
    library(genetics)
    hlagen <- genotype(a1=x$DQR.a1, a2=x$DQR.a2, 
                       alleles=sort(unique(c(x$DQR.a1, x$DQR.a2))))
    hwe.hardy(hlagen)

    # comparison with hwe
    hwe(z,data.type="count")

    # to create input file for HARDY
    print.tri<-function (xx,n) {
        cat(n,"\n")
        for(i in 1:n) {
            for(j in 1:i) {
                cat(xx[i,j]," ")
            }
        cat("\n")
        }
        cat("100 170 1000\n")
    }
    xx<-matrix(0,n.alleles,n.alleles)
    xxx<-lower.tri(xx,diag=TRUE)
    xx[xxx]<-z
    sink("z.dat")
    print.tri(xx,n.alleles)
    sink()
    # now call as: hwe z.dat z.out
  }
}

\note{Adapted from HARDY, testable with -Dexecutable as standalone program}

\keyword{htest}
