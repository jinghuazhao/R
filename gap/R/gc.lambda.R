#' Estimation of the genomic control inflation statistic (lambda)
#'
#' @param x A real vector (p or z).
#' @param logscale A logical variable such that x as -log10(p).
#' @param z A flag to indicate x as a vector of z values.
#' @export
#' @return Estimate of inflation factor.
#' @examples
#' set.seed(12345)
#' p <- runif(100)
#' gc.lambda(p)
#' lp <- -log10(p)
#' gc.lambda(lp,logscale=TRUE)
#' z <- qnorm(p/2)
#' gc.lambda(z,z=TRUE)

# gc.lambda and miamiplot functions hosted at CEU by Daniel R Barnes
# A simplified version is as follows,
# obs <- median(chisq)
# exp <- qchisq(0.5, 1) # 0.4549364
# lambda <- obs/exp
# see also GenABEL::estlambda and snpStats::qq.chisq

gc.lambda <- function(x, logscale=FALSE, z=FALSE) {
  v <- x[!is.na(x)]
  n <- length(v)
  if (z) {
     obs <- v^2
     exp <- qchisq(log(1:n/n),1,lower.tail=FALSE,log.p=TRUE)
  } else {
    if (!logscale)
    {
      obs <- qchisq(v,1,lower.tail=FALSE)
      exp <- qchisq(1:n/n,1,lower.tail=FALSE)
    } else {
      obs <- qchisq(-log(10)*v,1,lower.tail=FALSE,log.p=TRUE)
      exp <- qchisq(log(1:n/n),1,lower.tail=FALSE,log.p=TRUE)
    }
  }

  lambda <- median(obs)/median(exp)
  return(lambda)
}
