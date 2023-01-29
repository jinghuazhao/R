#' log(p) for a normal deviate z
#'
#' @param z normal deviate.
#' @export
#' @return
#' `log(P)`
#' @examples
#' logp(100)

logp <- function(z) log(2) + pnorm(-abs(z), lower.tail = TRUE, log.p = TRUE)
