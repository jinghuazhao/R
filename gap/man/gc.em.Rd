\name{gc.em}
\alias{gc.em}
\title{Gene counting for haplotype analysis}
\usage{
gc.em(data, locus.label=NA, converge.eps=1e-06, maxiter=500, 
      handle.miss=0, miss.val=0, control=gc.control())
}

\description{
Gene counting for haplotype analysis with missing data, adapted for hap.score
}
\arguments{
 \item{data}{Matrix of alleles, such that each locus has a  pair  of
   adjacent  columns  of  alleles,  and  the order of columns
   corresponds to the order of  loci  on  a  chromosome.   If
   there  are  K  loci, then ncol(data) = 2*K. Rows represent
   alleles for each subject.}
 \item{locus.label }{Vector of  labels  for  loci,  of  length  K  (see definition of data matrix).}
 \item{converge.eps }{Convergence criterion, based on absolute  change in log likelihood (lnlike).}
 \item{maxiter}{Maximum number of iterations of EM.}
 \item{handle.miss}{a flag for handling missing genotype data, 0=no, 1=yes}
 \item{miss.val}{missing value}
 \item{control}{a function, see \code{\link[gap]{genecounting}}}
}

\value{
List with components:
  \item{converge}{Indicator of convergence of the EM algorithm
  (1=converged, 0 = failed).}
  \item{niter}{Number of iterations completed in the EM alogrithm.}
  \item{locus.info}{A list with  a  component for each locus.  Each
   component is also a list, and  the  items of a locus-
   specific list are the locus name and a vector for the
   unique alleles for the locus.}
  \item{locus.label}{Vector of  labels  for  loci,  of  length  K  (see
    definition of input values).}
  \item{haplotype}{Matrix of unique haplotypes. Each row represents a
   unique  haplotype, and the number of columns is the number of loci.}
  \item{hap.prob}{Vector of mle's of haplotype probabilities.  The ith
   element of hap.prob corresponds to the ith row of haplotype.}
  \item{hap.prob.noLD}{Similar to hap.prob, but assuming no linkage
   disequilibrium.}
  \item{lnlike}{Value of lnlike at last EM iteration (maximum lnlike if converged).}
  \item{lr}{Likelihood ratio statistic to test no linkage disequilibrium among all loci.}
  \item{indx.subj}{Vector for index of subjects, after  expanding  to
   all possible  pairs  of  haplotypes  for  each person. If
   indx=i, then i is the ith row of input matrix data. If the
   ith subject has  n possible  pairs  of haplotypes that
   correspond to their marker phenotype, then i is repeated n times.}
  \item{nreps}{Vector for the count of haplotype pairs that map to
   each subject's marker genotypes.}
  \item{hap1code}{Vector of codes for each subject's first haplotype.
   The values in hap1code are the row numbers of the unique
   haplotypes in the returned matrix haplotype.}
  \item{hap2code}{Similar to hap1code, but for  each  subject's  second haplotype.}
  \item{post}{Vector of posterior probabilities of pairs of
   haplotypes for a person, given thier marker phenotypes.}
  \item{htrtable}{A table which can be used in haplotype trend regression}
}

\section{References}{

Zhao, J. H., Lissarrague, S., Essioux, L. and P. C. Sham (2002).
GENECOUNTING: haplotype analysis with missing genotypes.
Bioinformatics 18(12):1694-1695

Zhao, J. H. and P. C. Sham (2003). Generic number systems and haplotype
analysis. Comp Meth Prog Biomed 70: 1-9

}
\seealso{
\code{\link[gap]{genecounting}}, \code{\link[gap]{LDkl}}
}

\examples{
\dontrun{
data(hla)
gc.em(hla[,3:8],locus.label=c("DQR","DQA","DQB"),control=gc.control(assignment="t"))
}
}

\author{Jing Hua Zhao}
\note{Adapted from GENECOUNTING}
\keyword{models}
