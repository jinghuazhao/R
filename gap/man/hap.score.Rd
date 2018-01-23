\name{hap.score}
\alias{hap.score}
\title{Score statistics for association of traits with haplotypes}
\description{
Compute score statistics to evaluate the association of a trait with haplotypes, when linkage phase is unknown and diploid marker 
phenotypes are observed among unrelated subjects. For now, only autosomal loci are considered. This package haplo.score
which this function is based is greatly acknowledged.
}
\usage{
hap.score(y, geno, trait.type="gaussian", offset=NA, x.adj=NA,
          skip.haplo=0.005, locus.label=NA, miss.val=0, n.sim=0, 
          method="gc", id=NA, handle.miss=0, mloci=NA, sexid=NA)
}

\arguments{
 \item{y}{Vector of trait values. For  trait.type  =  "binomial",  y  must have values of 1 for event, 0 for no event.}
 \item{geno}{Matrix of alleles, such that each locus has a pair of adjacent columns of alleles, and the order of columns corresponds to the order of loci on a chromosome. If there are K loci, then ncol(geno) = 2*K. Rows represent alleles for each subject.}
 \item{trait.type}{Character string  defining  type  of  trait, with values of "gaussian", "binomial", "poisson", "ordinal".}
 \item{offset}{Vector of offset when trait.type = "poisson"}
 \item{x.adj}{Matrix of non-genetic covariates used to adjust the score statistics. Note that intercept should not be included, as it will be added in this function.}
 \item{skip.haplo}{Skip score statistics for haplotypes with frequencies < skip.haplo}
 \item{locus.label}{Vector of labels for loci, of length K (see definition of geno matrix).}
 \item{miss.val}{Vector of codes for missing values of alleles.}
 \item{n.sim}{Number of simulations for empirical p-values.  If n.sim=0, no empirical p-values are computed.}
 \item{method}{method of haplotype frequency estimation, "gc" or "hap"}
 \item{id}{an added option which contains the individual IDs}
 \item{handle.miss}{flag to handle missing genotype data, 0=no, 1=yes}
 \item{mloci}{maximum number of loci/sites with missing data to be allowed in the analysis}
 \item{sexid}{flag to indicator sex for data from X chromosome, i=male, 2=female}
}

\value{
List with the following components:
 \item{score.global}{Global statistic to test association of trait with haplotypes that have frequencies >= skip.haplo.}
 \item{df}{Degrees of freedom for score.global.}
 \item{score.global.p}{P-value of score.global based on chi-square distribution, with degrees of freedom equal to df.}
 \item{score.global.p.sim}{P-value of score.global based on simulations (set equal to NA when n.sim=0).}
 \item{score.haplo}{Vector of score statistics for individual haplotypes that have frequencies >= skip.haplo.}
 \item{score.haplo.p}{Vector of p-values for score.haplo, based on a chi-square distribution with 1 df.}
 \item{score.haplo.p.sim}{Vector of p-values for score.haplo, based on  simulations (set equal to NA when n.sim=0).}
 \item{score.max.p.sim}{P-value  of  maximum  score.haplo, based on simulations (set equal to NA when n.sim=0).}
 \item{haplotype}{Matrix of hapoltypes  analyzed.  The ith row of haplotype corresponds to the ith item of score.haplo, score.haplo.p, and score.haplo.p.sim.}
 \item{hap.prob}{Vector of haplotype probabilies, corresponding to the haplotypes in the matrix haplotype.}
 \item{locus.label}{Vector of labels for loci, of length K (same as input argument).}
 \item{n.sim}{Number of simulations.}
 \item{n.val.global}{Number of valid simulated global statistics.}
 \item{n.val.haplo}{Number of valid simulated score statistics (score.haplo) for individual haplotypes.}
}

\details{This is a version which substitutes haplo.em}

\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA (2002)
Score tests for association of traits with haplotypes when
linkage phase is ambiguous. Amer J Hum Genet 70:425-34
}

\examples{
\dontrun{
data(hla)
y<-hla[,2]
geno<-hla[,3:8]
# complete data
hap.score(y,geno,locus.label=c("DRB","DQA","DQB"))
# incomplete genotype data
hap.score(y,geno,locus.label=c("DRB","DQA","DQB"),handle.miss=1,mloci=1)
unlink("assign.dat")

### note the differences in p values in the following runs
data(aldh2)
# to subset the data since hap doesn't handle one allele missing
deleted<-c(40,239,256)
aldh2[deleted,]
aldh2<-aldh2[-deleted,]
y<-aldh2[,2]
geno<-aldh2[,3:18]
# only one missing locus
hap.score(y,geno,handle.miss=1,mloci=1,method="hap")
# up to seven missing loci and with 10,000 permutations
hap.score(y,geno,handle.miss=1,mloci=7,method="hap",n.sim=10000)

# hap.score takes considerably longer time and does not handle missing data
hap.score(y,geno,n.sim=10000)
}
}
\keyword{models}
\keyword{regression}
