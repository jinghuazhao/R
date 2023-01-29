#' Get b and se from AF, n, and z
#'
#' The function obtains effect size and its standard error.
#'
#' @param f Allele frequency.
#' @param n Sample size.
#' @param z z-statistics.
#' @export
#' @return
#' b and se.
#' @examples
#' \dontrun{
#'   library(dplyr)
#'   # eQTLGen
#'   cis_pQTL <- merge(read.delim('eQTLGen.lz') %>%
#'               filter(GeneSymbol=="LTBR"),read.delim("eQTLGen.AF"),by="SNP") %>%
#'               mutate(data.frame(get_b_se(AlleleB_all,NrSamples,Zscore)))
#'   head(cis_pQTL,1)
#'         SNP    Pvalue SNPChr  SNPPos AssessedAllele OtherAllele Zscore
#'   rs1003563 2.308e-06     12 6424577              A           G 4.7245
#'              Gene GeneSymbol GeneChr GenePos NrCohorts NrSamples         FDR
#'   ENSG00000111321       LTBR      12 6492472        34     23991 0.006278872
#'   BonferroniP hg19_chr hg19_pos AlleleA AlleleB allA_total allAB_total
#'             1       12  6424577       A       G       2574        8483
#'   allB_total AlleleB_all          b          se
#'         7859   0.6396966 0.04490488 0.009504684
#' }

get_b_se <- function(f,n,z)
{
  v <- 2*f*(1-f)
  se <- 1 / sqrt(v*(n + z^2))
  b <- z * se
  cbind(b,se)
}
