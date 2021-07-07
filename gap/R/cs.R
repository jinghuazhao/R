#' Credible set
#'
#' The function implements credible set as in fine-mapping.
#'
#' @md
#' @param tbl Input data.
#' @param b Effect size.
#' @param se Standard error.
#' @param log_p if not NULL it will be used to derive z-statistic
#' @param cutoff Threshold for inclusion.
#' @export
#' @return Credible set.
#' @examples
#' \dontrun{
#' \preformatted{
#'   zcat METAL/4E.BP1-1.tbl.gz | \
#'   awk 'NR==1 || ($1==4 && $2 >= 187158034 - 1e6 && $2 < 187158034 + 1e6)' > 4E.BP1.z
#' }
#'   tbl <- within(read.delim("4E.BP1.z"),{logp <- logp(Effect/StdErr)})
#'   z <- cs(tbl)
#'   l <- cs(tbl,log_p="logp")
#' }

cs <- function(tbl, b="Effect", se="StdErr", log_p=NULL, cutoff=0.95)
# credible set based on METAL sumstats
{
  requireNamespace("matrixStats")
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
