#' Credible set from summary statistics
#'
#' Computes a Bayesian credible set from GWAS summary statistics
#' using Wakefield-style approximate Bayes factors.
#'
#' @param tbl A data.frame containing summary statistics.
#' @param b Name of column containing effect sizes.
#' @param se Name of column containing standard errors.
#' @param log_p Optional column name containing log p-values.
#'   If supplied, z-scores are derived from p-values instead of
#'   effect sizes and standard errors.
#' @param cutoff Cumulative posterior probability threshold.
#'   Default is 0.95 (95% credible set).
#'
#' @return A subset of `tbl` containing variants in the credible set,
#' ordered by decreasing posterior probability of association.
#'
#' @details
#' Credible set is often used in fine-mapping.
#'
#' Posterior probabilities are computed from z-statistics using a
#' numerically stable log-sum-exp implementation from matrixStats.
#'
#' @export
#' @examples
#' \dontrun{
#' \preformatted{
#'   zcat ~/rds/results/private/proteomics/scallop-inf1/4E.BP1-1.tbl.gz | \
#'   awk 'NR==1 || ($1==4 && $2 >= 187158034 - 1e6 && $2 < 187158034 + 1e6)' > 4E.BP1.z
#' }
#'   tbl <- within(read.delim("4E.BP1.z"),{logp <- logp(Effect/StdErr)})
#'   z <- cs(tbl)
#'   l <- cs(tbl,log_p="logp")
#' }
#'
cs <- function(tbl, b="Effect", se="StdErr", log_p=NULL, cutoff=0.95)
# credible set based on METAL sumstats
{
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop(
      "The 'cs()' function requires the 'matrixStats' package.\n",
      "Please install it with install.packages('matrixStats').",
      call. = FALSE
    )
  }
  tbl <- within(tbl, {
           if (is.null(log_p)) z <- tbl[[b]]/tbl[[se]]
           else z <- qnorm(tbl[[log_p]]-log(2), lower.tail=FALSE, log.p=TRUE)
           z2 <- z * z / 2
           d <- matrixStats::logSumExp(z2)
           log_ppa <- z2 - d
           ppa <- exp(log_ppa)
        })
  ord <- with(tbl, order(ppa,decreasing = TRUE))
  tbl_ord <- within(tbl[ord,], {cppa <- cumsum(ppa)})
  last <- which(with(tbl_ord,cppa) >= cutoff)[1]
  tbl[ord[1:last],]
}
