#' genomic control based on p values
#'
#' @param p a vector of observed p values.
#' @param col colour for points in the Q-Q plot.
#' @param lcol colour for the diagonal line in the Q-Q plot.
#' @param ... other options for plot.
#'
#' @details
#' The function obtains 1-df \eqn{\chi^2}{chi-squared} statistics (observed) according
#' to a vector of p values, and the inflation factor (lambda) according
#' to medians of the observed and expected statistics. The latter is based
#' on the empirical distribution function (EDF) of 1-df \eqn{\chi^2}{chi-squared} statstics.
#'
#' It would be appropriate for genetic association analysis as of 1-df Armitage trend test
#' for case-control data; for 1-df additive model with continuous outcome one has to
#' consider the compatibility with p values based on z-/t- statistics.
#'
#' @export
#' @return A list containing:
#' - x the expected \eqn{\chi^2}{chi-squared} statistics.
#' - y the observed \eqn{\chi^2}{chi-squared} statistics.
#' - lambda the inflation factor.
#'
#' @references
#' Devlin B, Roeder K (1999) Genomic control for association studies. 
#' Biometrics 55:997-1004
#'
#' @examples
#' \dontrun{
#' x2 <- rchisq(100,1,.1)
#' p <- pchisq(x2,1,lower.tail=FALSE)
#' r <- gcontrol2(p)
#' print(r$lambda)
#' }
#'
#' @author Jing Hua Zhao
#' @keywords models

gcontrol2 <- function(p,col=palette()[4],lcol=palette()[2],...)
{
   p <- p[!is.na(p)]
   n <- length(p)
   x2obs <- qchisq(p,1,lower.tail=FALSE)
   x2exp <- qchisq(1:n/n,1,lower.tail=FALSE)
   lambda <- median(x2obs)/median(x2exp)
   qqunif(p,col=col,lcol=lcol,...)
   invisible(list(x=x2exp,y=x2obs,lambda=lambda))
}
