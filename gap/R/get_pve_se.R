#' Get pve and its standard error from n, z
#'
#' @param n Sample size.
#' @param z z-statistic, i.e., b/se when they are available instead.
#' @param correction if TRUE an correction based on t-statistic is applied.
#'
#' @details
#' This function obtains proportion of explained variance of a continuous outcome.
#'
#' @export
#' @return pve and its se.

get_pve_se <- function(n,z,correction=TRUE)
{
  pve <- ifelse(correction,z^2/(z^2+n-2),z^2/(z^2+n))
  se <- ifelse(correction,sqrt(2)/(n-1),sqrt(2)/(n+1))
  cbind(pve,se)
}
