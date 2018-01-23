\name{SNP}
\alias{SNP}
\alias{snp.ES}
\alias{snp.HWE}
\alias{snp.PAR}
\title{Functions for single nucleotide polymorphisms (SNPs)}
\usage{
snp.ES(beta,SE,N)
snp.HWE(g)
snp.PAR(RR,MAF,unit=2)
}
\description{
snp.PAR gives PAR for a particular SNP.

snp.ES gives the effect size estimate based on the linear regression
coefficient and standard error. For logistic regression, we can have
similar idea for log(OR) and log(SE(OR)).

snp.HWE gives an exact Hardy-Weinberg Equilibrium (HWE) test, and -1
in the case of misspecification of genotype counts.

Eventually, this will be a set of functions specifically for single
nucleotide polymorphisms (SNPs), which are biallelic markers. This
is particularly relevant to the genomewide association studies (GWAS)
using GeneChips and in line with the classic generalised single-locus
model. snp.HWE is from Abecasis's website and yet to adapt for 
chromosome X.

Internally, snp.PAR calls for an internal function PARn, whcih
calculates the the population attributable risk (PAR) given a 
set of frequencies and associate relative risks (RR). Other
2x2 table statistics familiar to epidemiologists can be added when
necessary.
}
\arguments{
  \item{MAF}{Minar allele frequency}
  \item{RR}{Relative risk}
  \item{unit}{Unit to exponentiate for homozygote}
  \item{beta}{Regression coefficient}
  \item{SE}{Standard error for beta}
  \item{N}{Sample size}
  \item{g}{Observed genotype vector}
}
\author{Jing Hua Zhao, Shengxu Li}
\keyword{utilities}
