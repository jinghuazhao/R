#' log(p) for a normal deviate z
#'
#' @md
#' @param z normal deviate.
#' @export
#' @return log10(P)
#' @examples
#' logp(100)

logp <- function(z) log(2) + pnorm(-abs(z), lower.tail = TRUE, log.p = TRUE)
