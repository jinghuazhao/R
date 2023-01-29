#' @title Internal functions for gap
#'
#' @param a1 Allele 1.
#' @param a2 Allele 2.
#' @param g A genotype identifier.
#' @param p1 Frequency of the risk allele.
#' @param K Prevalence of disease in the population.
#' @param loci A vector of number of alleles at all loci.
#' @param hapid Haplotype identifier.
#'
#' @name gap-internal
#' @alias{allDuplicated}
#' @alias{hg18}
#' @alias{hg19}
#' @alias{hg38}
#' @alias{lambda1000}
#' @alias{HapDesign}
#' @alias{HapFreqSE}
#' @alias{ci}
#' @alias{cov.invlogit}
#' @alias{g2a.c}
#' @alias{gc.control}
#' @alias{gcode}
#' @alias{geno.p}
#' @alias{getb1star}
#' @alias{getPTE}
#' @alias{grec2g}
#' @alias{hap.score.glm}
#' @alias{hap.score.podds}
#' @alias{hw}
#' @alias{invlogit}
#' @alias{is.miss}
#' @alias{k}
#' @alias{logit}
#' @alias{m2plem}
#' @alias{makeRLEplot}
#' @alias{micombine}
#' @alias{mr.boot}
#' @alias{plem2m}
#' @alias{revhap}
#' @alias{revhap.i}
#' @alias{ReadGRMPLINK}
#' @alias{ReadGRMPCA}
#' @alias{se.exp}
#' @alias{se.invlogit}
#' @alias{solve_skol}
#' @alias{textbox}
#' @alias{toETDT}
#' @alias{ungcode}
#' @alias{VR}
#' @alias{weighted.median}
#' @alias{WriteGRMSAS}
#' @alias{x2}
#' @alias{z}
#'
#' @section Usage:
#' ReadGRMPLINK(prefix, diag=1)
#' ReadGRMPCA(prefix)
#' revhap(loci,hapid)
#' VR(v1,vv1,v2,vv2,c12)
#' WriteGRMSAS(grmlist, outfile="gwas")
#'
#' @details
#' These functions are not so frequently called by users.
#'
#' `g2a.c` is the C version of `g2c`.
#'
#' `gc.control` is used by `gc.em()`.
#'
#' `gcode` is as `a2g`.
#'
#' `grec2g` is undocumented.
#'
#' `HapDesign` and `HapFreqSE` both accept a \code{\link[haplo.stats]{haplo.em}} object to derieve a design/dosage
#' matrix and standard error of haplotype frequency estimates. The former is appropriate for haplotype trend 
#' regression (HTR), e.g., within the generalized linear model (GLM) framework to be equivllant to a formal 
#' approach as implemented in the package haplo.stats and hap.score. However, they are expected to be compatible 
#' with objects from gc.em() \code{\link[gap]{gc.em}} and \code{\link[gap]{hap.em}}. The two functions are 
#' included as courtesy of Prof Andrea Foulkes from the useR!2008 tutorial.
#'
#' `hap.score.glm`, `hap.score.podds` are used by hap.score().
#'
#' `invlogit`, inverse logit transformation.
#'
#' `is.miss` is undocumented.
#'
#' `k` obtains 1st and 2nd order culumants for correlation coefficient.
#'
#' `m2plem` is an experimental function for PLEM format.
#'
#' `makeRLEplot` for RLE plot.
#'
#' `micombine` is used to combine imputation results.
#'
#' `plem2m` is also experimental for PLEM format.
#'
#' `ReadGRMPLINK` is a function to read PLINK PI_HAT as a genomic relationship matrix.
#'
#' `ReadGRMPCA` is a function to read .eigenval and .eigenvec files from gcta --pca.
#'
#' `revhap` recovers the allele indices for a given haplotype ID in a multiallelic system.
#'
#' `revhap.i` is similar to revhap.
#'
#' `solve.skol` is a function used by tscc.
#'
#' `toETDT` a function used to experiment with ETDT.
#'
#' `ungcode` recovers alleles from genotype(s).
#'
#' `VR` is a utility function for calculating variance of a ratio.
#'
#' `weighted.median` is a function for obtaining weighted median with interpolation.
#'
#' `WriteGRMSAS` is a utility function to write a GRM object to SAS PROCs MIXED/GLIMMIX ldata.
#'
#' `x2` is a simple chi-squared test of two proportions.
#'
#' `z` is a normal z-test of two proportions used by tscc.
#'
#' @keyword internal
