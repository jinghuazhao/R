#' Get b and se from AF, n, and z
#'
#' The function obtains effect size and its standard error.
#'
#' @md
#' @param f Allele frequency.
#' @param n Sample size.
#' @param z z-statistics.
#' @export
#' @return b and se.
#' @examples
#' \dontrun{
#'   get_b_se(0.6396966,23991,4.7245)
#'#               b          se
#'# [1,] 0.04490488 0.009504684
#'# eQTLGen
#'   cis_pQTL <- merge(read.delim('eQTLGen.lz') %>%
#'               filter(GeneSymbol=="LTBR"),read.delim("eQTLGen.AF"),by="SNP") %>%
#'               mutate(data.frame(get_b_se(AlleleB_all,NrSamples,Zscore)))
#' }

get_b_se <- function(f,n,z)
{
  v <- 2*f*(1-f)
  se <- 1 / sqrt(v*(n + z^2))
  b <- z * se
  cbind(b,se)
}
