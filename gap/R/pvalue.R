#' P value for a normal deviate
#'
#' @md
#' @param z normal deviate.
#' @param decimals number of decimal places.
#' @export
#' @return P value as a string variable.
#' @examples
#' pvalue(-1.96)

pvalue <- function (z, decimals = 2)
{
    lp <- -log10p(z)
    exponent <- ceiling(lp)
    base <- 10^-(lp - exponent)
    paste0(round(base, decimals), "e", -exponent)
}
