#' Estionmation of the genomic control inflation statistic (lambda)
#'
#' @param p A vector of p values.
#' @param logscale A logical variable for `log.p`
#' @export
#' @return Estimate of inflation factor.
#' @examples
#' set.seed(12345)
#' p <- runif(100)
#' gc.lambda(p)
#' lp <- -log10(p)
#' gc.lambda(lp,logscale=TRUE)

# gc.lambda and miamiplot functions hosted at CEU by Daniel R Barnes
# A simplified version is as follows,
# obs <- median(chisq)
# exp <- qchisq(0.5, 1) # 0.4549364
# lambda <- obs/exp
# see also estlambda from GenABEL and qq.chisq from snpStats
# Note when `logscale=TRUE` one assumes a -log10(p) is assumed

gc.lambda <- function(p, logscale=FALSE) {
  p <- p[!is.na(p)]
  n <- length(p)

  if (!logscale)
  {
    obs <- qchisq(p,1,lower.tail=FALSE)
    exp <- qchisq(1:n/n,1,lower.tail=FALSE)
  } else {
    obs <- qchisq(-log(10)*p,1,lower.tail=FALSE,log.p=TRUE)
    exp <- qchisq(log(1:n/n),1,lower.tail=FALSE,log.p=TRUE)
  }

  lambda <- median(obs)/median(exp)
  return(lambda)
}
