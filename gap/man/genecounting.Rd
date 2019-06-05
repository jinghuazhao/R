\name{genecounting}
\alias{genecounting}
\title{Gene counting for haplotype analysis}
\usage{genecounting(data,weight=NULL,loci=NULL,control=gc.control())}
\description{
Gene counting for haplotype analysis with missing data
}
\arguments{
 \item{data}{genotype table}
 \item{weight}{a column of frequency weights}
 \item{loci}{an array containing number of alleles at each locus}
 \item{control}{is a function with the following arguments: 
      \enumerate{
      \item xdata. a flag indicating if the data involves X chromosome, if so, the
first column of data indicates sex of each subject: 1=male, 2=female. The marker 
data are no different from the autosomal version for females, but for males, two copies
of the single allele present at a given locus.
      \item convll. set convergence criteria according to log-likelihood, if its value set to 1
      \item handle.miss. to handle missing data, if its value set to 1
      \item eps. the actual convergence criteria, with default value 1e-5
      \item tol. tolerance for genotype probabilities with default value 1e-8
      \item maxit. maximum number of iterations, with default value 50
      \item pl. criteria for trimming haplotypes according to posterior probabilities
      \item assignment. filename containing haplotype assignment
      \item verbose. If TRUE, yields print out from the C routine
      }
  }
}

\value{
The returned value is a list containing:

 \item{h}{haplotype frequency estimates under linkage disequilibrium (LD)}
 \item{h0}{haplotype frequency estimates under linkage equilibrium (no LD)}
 \item{prob}{genotype probability estimates}
 \item{l0}{log-likelihood under linkage equilibrium}
 \item{l1}{log-likelihood under linkage disequilibrium}
 \item{hapid}{unique haplotype identifier (defunct, see gc.em)}
 \item{npusr}{number of parameters according user-given alleles}
 \item{npdat}{number of parameters according to observed}
 \item{htrtable}{design matrix for haplotype trend regression (defunct, see gc.em)}
 \item{iter}{number of iterations used in gene counting}
 \item{converge}{a flag indicating convergence status of gene counting}
 \item{di0}{haplotype diversity under no LD, defined as \eqn{1-\sum (h_0^2)}{1-sum (h0^2)}}
 \item{di1}{haplotype diversity under LD, defined as \eqn{1-\sum (h^2))}{1-sum (h^2)}}
 \item{resid}{residuals in terms of frequency weights = o - e}
}

\section{References}{

Zhao, J. H., Lissarrague, S., Essioux, L. and P. C. Sham (2002).
GENECOUNTING: haplotype analysis with missing genotypes.
Bioinformatics 18(12):1694-1695

Zhao, J. H. and P. C. Sham (2003). Generic number systems and haplotype
analysis. Comp Meth Prog Biomed 70: 1-9

Zhao, J. H. (2004). 2LD, GENECOUNTING and HAP: Computer programs for linkage
disequilibrium analysis. Bioinformatics, 20, 1325-1326 

}
\seealso{
\code{\link[gap]{gc.em}}, \code{\link[gap]{LDkl}}
}

\examples{
\dontrun{
require(gap.datasets)
# HLA data
data(hla)
hla.gc <- genecounting(hla[,3:8])
summary(hla.gc)
hla.gc$l0
hla.gc$l1

# ALDH2 data
data(aldh2)
control <- gc.control(handle.miss=1,assignment="ALDH2.out")
aldh2.gc <- genecounting(aldh2[,3:6],control=control)
summary(aldh2.gc)
aldh2.gc$l0
aldh2.gc$l1

# Chromosome X data
# assuming allelic data have been extracted in columns 3-13
# and column 3 is sex
filespec <- system.file("tests/genecounting/mao.dat")
mao2 <- read.table(filespec)
dat <- mao2[,3:13]
loci <- c(12,9,6,5,3)
contr <- gc.control(xdata=TRUE,handle.miss=1)
mao.gc <- genecounting(dat,loci=loci,control=contr)
mao.gc$npusr
mao.gc$npdat
}
}
\author{Jing Hua Zhao}
\note{adapted from GENECOUNTING}
\keyword{models}
