#' log10(p) for a normal deviate z
#'
#' @param z normal deviate.
#' @export
#' @return
#' `log10(P)`
#' @examples
#' log10p(100)
#' @author James Peters

log10p <- function(z)
  log(2, base=10)+pnorm(-abs(z), lower.tail=TRUE, log.p=TRUE)/log(10)
