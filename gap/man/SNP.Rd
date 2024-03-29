% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNP.R
\name{snpHWE}
\alias{snpHWE}
\alias{PARn}
\alias{snpPVE}
\alias{snpPAR}
\title{Functions for single nucleotide polymorphisms}
\usage{
snpHWE(g)

PARn(p, RRlist)

snpPVE(beta, se, N)

snpPAR(RR, MAF, unit = 2)
}
\arguments{
\item{g}{Observed genotype vector.}

\item{p}{genotype frequencies.}

\item{RRlist}{A list of RRs.}

\item{beta}{Regression coefficient.}

\item{se}{Standard error for beta.}

\item{N}{Sample size.}

\item{RR}{Relative risk.}

\item{MAF}{Minar allele frequency.}

\item{unit}{Unit to exponentiate for homozygote.}
}
\description{
These are a set of functions specifically for single nucleotide
polymorphisms (SNPs), which are biallelic markers. This is
particularly relevant to the genomewide association studies (GWAS)
using GeneChips and in line with the classic generalised single-locus
model. snpHWE is from Abecasis's website and yet to be adapted for
chromosome X.
}
\details{
\code{snpHWE} gives an exact Hardy-Weinberg Equilibrium (HWE) test and it return -1 in the case of misspecification of genotype counts.

\code{snpPAR} calculates the the population attributable risk (PAR) for a particular SNP.
Internally, it calls for an internal function PARn, given a
set of frequencies and associate relative risks (RR). Other
2x2 table statistics familiar to epidemiologists can be added when
necessary.

\code{snpPVE} provides proportion of variance explained (PVE) estimate based on the linear regression coefficient and standard error.
For logistic regression, we can have similar idea for log(OR) and log(SE(OR)).
}
\author{
Jing Hua Zhao, Shengxu Li
}
\keyword{utilities}
