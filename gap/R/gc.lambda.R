#' Estionmation of the genomic control inflation statistic (lambda)
#'
#' @md
#' @param p A vector of p values.
#' @export
#' @return Estimate of inflation factor.
#' @examples
#' set.seed(12345)
#' p <- runif(100)
#' lambda <- gc.lambda(p)

# gc.lambda and miamiplot functions hosted at CEU by Daniel R Barnes
# A simplified version is as follows,
# obs <- median(chisq)
# exp <- qchisq(0.5, 1) # 0.4549364
# lambda <- obs/exp
# see also estlambda from GenABEL and qq.chisq from snpStats

gc.lambda <- function(p) {
  p <- p[!is.na(p)]
  n <- length(p)

  obs <- qchisq(p,1,lower.tail=FALSE)
  exp <- qchisq(1:n/n,1,lower.tail=FALSE)

  lambda <- median(obs)/median(exp)
  return(lambda)
}
