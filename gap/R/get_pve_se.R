#' Get pve and its standard error from n, z
#'
#' This function obtains proportion of explained variance of a continuous outcome.
#' The standard error is based on \eqn{r^2} and appears \eqn{\sqrt{2}} less than that
#' from variance of the ratio.
#'
#' @md
#' @param n Sample size.
#' @param z z-statistic, i.e., b/se when they are available instead.
#' @param correction if TRUE an correction based on t-statistic is applied.
#' @export
#' @return pve.

get_pve_se <- function(n,z,correction=TRUE)
{
  pve <- ifelse(correction,z^2/(z^2+n-2),z^2/(z^2_+n))
  se <- 1/sqrt(n-1)
  cbind(pve,se)
}
