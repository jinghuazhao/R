\name{gap-package}
\alias{gap-package}
\alias{gap}
\docType{package}
\title{
Genetic analysis package
}
\description{
It is designed as an integrated package for genetic data analysis
of both population and family data. Currently, it contains
functions for sample size calculations of both population-based
and family-based designs, classic twin models, 
probability of familial disease aggregation, kinship calculation, 
some statistics in linkage analysis, and association analysis 
involving one or more genetic markers including haplotype 
analysis with or without environmental covariates.
}
\details{
\tabular{ll}{
Package: \tab gap\cr
Version: \tab 1.2.2\cr
Depends: \tab R(>= 2.1.0)\cr
Suggests: \tab BradleyTerry2, MASS, Matrix, MCMCglmm, R2jags, bdsmatrix, \cr
          \tab calibrate, circlize, coda, coxme, dplyr, foreign, forestplot, \cr
          \tab gap.datasets, grid, haplo.stats, kinship2, lattice, \cr
          \tab magic, matrixStats, meta, mets, nlme, pedigree, pedigreemm, plotrix, \cr
          \tab qqman, regress, reshape, rmeta, rms, survival \cr
License: \tab GPL (>=2)\cr
URL: \tab https://jinghuazhao.github.io/R\cr
}

Index:
\tabular{ll}{
\cr
* ANALYSIS *\cr
\cr
AE3               \tab      AE model using nuclear family trios\cr
bt                \tab      Bradley-Terry model for contingency table\cr
ccsize            \tab      Power and sample size for case-cohort design\cr
cs                \tab      Credibel set\cr
fbsize            \tab      Sample size for family-based linkage and association design\cr
gc.em             \tab      Gene counting for haplotype analysis\cr
gcontrol          \tab      genomic control\cr
gcontrol2         \tab      genomic control based on p values\cr
gcp               \tab      Permutation tests using GENECOUNTING\cr
gc.lambda         \tab      Estimation of the genomic control inflation statistic (lambda)\cr
genecounting      \tab      Gene counting for haplotype analysis\cr
gif               \tab      Kinship coefficient and genetic index of familiality\cr
hap               \tab      Haplotype reconstruction\cr
hap.em            \tab      Gene counting for haplotype analysis\cr
hap.score         \tab      Score statistics for association of traits with haplotypes\cr
htr               \tab      Haplotype trend regression\cr
h2.jags           \tab      Heritability estimation based on genomic relationship matrix using JAGS\cr
hwe               \tab      Hardy-Weinberg equilibrium test for a multiallelic marker\cr
hwe.cc            \tab      A likelihood ratio test of population Hardy-Weinberg equilibrium \cr
hwe.hardy         \tab      Hardy-Weinberg equilibrium test using MCMC\cr
hwe.jags          \tab      Hardy-Weinberg equlibrium test for a multiallelic marker using JAGS\cr
invnormal         \tab      inverse Normal transformation\cr
kin.morgan        \tab      kinship matrix for simple pedigree\cr
LD22              \tab      LD statistics for two diallelic markers\cr
LDkl              \tab      LD statistics for two multiallelic markers\cr
lambda1000        \tab      A standardized estimate of the genomic inflation scaling to a study of\cr
                  \tab      1,000 cases and 1,000 controls\cr
log10p            \tab      log10(p) for a standard normal deviate\cr
logp              \tab      log(p) for a normal deviate\cr
masize            \tab      Sample size calculation for mediation analysis\cr
MCMCgrm           \tab	    Mixed modeling with genetic relationship matrices\cr
mia               \tab      multiple imputation analysis for hap\cr
mtdt              \tab      Transmission/disequilibrium test of a multiallelic marker\cr
mtdt2             \tab      Transmission/disequilibrium test of a multiallelic marker\cr
                  \tab      by Bradley-Terry model\cr
mvmeta            \tab      Multivariate meta-analysis based on generalized least squares\cr
pbsize            \tab      Power for population-based association design\cr
pbsize2           \tab      Power for case-control association design\cr
pfc               \tab      Probability of familial clustering of disease\cr
pfc.sim           \tab      Probability of familial clustering of disease\cr
pgc               \tab      Preparing weight for GENECOUNTING\cr
print.hap.score   \tab      Print a hap.score object\cr
s2k               \tab      Statistics for 2 by K table\cr
sentinels         \tab      Sentinel identification from GWAS summary statistics\cr
tscc              \tab      Power calculation for two-stage case-control design\cr

\cr
* GRAPHICS *\cr
\cr
asplot            \tab      Regional association plot\cr
ESplot            \tab      Effect-size plot\cr
circos.cis.vs.trans.plot \tab circos plot of cis/trans classification\cr
circos.cnvplot    \tab      circos plot of CNVs\cr
circos.mhtplot 	  \tab      circos Manhattan plot with gene annotation\cr
cnvplot           \tab      genomewide plot of CNVs\cr
makeRLEplot       \tab      make relative log expression plot\cr
METAL_forestplot  \tab      forest plot as R/meta's forest for METAL outputs\cr
mhtplot           \tab      Manhattan plot\cr
mhtplot2          \tab      Manhattan plot with annotations\cr
mhtplot2d         \tab      2D Manhattan plot\cr
mhtplot.trunc     \tab      truncated Manhattan plot\cr
miamiplot         \tab      Miami plot (Experimental)\cr
pedtodot          \tab      Converting pedigree(s) to dot file(s)\cr
plot.hap.score    \tab      Plot haplotype frequencies versus haplotype score statistics\cr
qqfun             \tab      Quantile-comparison plots\cr
qqunif            \tab      Q-Q plot for uniformly distributed random variable\cr

\cr
* UTITLITIES *\cr
\cr
SNP               \tab      Functions for single nucleotide polymorphisms (SNPs)\cr
BFDP              \tab      Bayesian false-discovery probability\cr
FPRP              \tab      False-positive report probability\cr
ab                \tab      Test/Power calculation for mediating effect\cr
b2r               \tab      Obtain correlation coefficients and their variance-covariances\cr
chow.test         \tab      Chow's test for heterogeneity in two regressions\cr
chr_pos_a1_a2     \tab      Form SNPID from chromosome, posistion and alleles\cr
cis.vs.trans.classification \tab a cis/trans classifier\cr
comp.score        \tab      score statistics for testing genetic linkage of quantitative trait\cr
h2                \tab      Heritability estimation according to twin correlations\cr
                  \tab      for case-control studies\cr
klem              \tab      Haplotype frequency estimation based on a genotype table\cr
                  \tab      of two multiallelic markers\cr
makeped           \tab      A function to prepare pedigrees in post-MAKEPED format\cr
metap             \tab      Meta-analysis of p values\cr
metareg           \tab      Fixed and random effects model for meta-analysis\cr
muvar             \tab      Means and variances under 1- and 2- locus (diallelic) QTL model\cr
read.ms.output    \tab      A utility function to read ms output\cr
snptest_sample    \tab      A utility to generate SNPTEST sample file\cr
whscore           \tab      Whittemore-Halpern scores for allele-sharing\cr
GRM functions     \tab      ReadGRM, ReadGRMBin, ReadGRMPLINK, ReadGRMPCA, WriteGRM, WriteGRMBin, WriteGRMSAS\cr
                  \tab      handle genomic relationship matrix involving other software\cr
heritability functions \tab h2G, VR, h2GC, h2l give point estimates as with their variances\cr
                       \tab for continuous traits and binary traits under liability threshold model and\cr
                       \tab case-control sampling\cr
\cr
}

We have incorporated functions for a wide range of problems. Nevertheless,
this largely remains as a preliminary work to be consolidated in the
near future.

}
\author{
Author: Jing Hua Zhao in collaboration with other colleagues, and with
        help from Kurt Hornik and Brian Ripley of the R core development team

Maintainer: Jing Hua Zhao <jinghuazhao@hotmail.com>
}
\references{
Zhao JH, gap: genetic analysis package. Journal of Statistical Software 2007, 23(8):1-18
}
\keyword{package}
