#' Q-Q plot for uniformly distributed random variables
#'
#' @param u A vector of uniformly distributed random variables.
#' @param type Distribution type: `"unif"` for uniform order statistics or `"exp"` for exponential order statistics.
#' @param logscale Logical; use log scale.
#' @param base Base of the logarithm.
#' @param col Color for points.
#' @param lcol Color for the diagonal reference line.
#' @param ci Logical; show confidence intervals.
#' @param alpha Significance level for confidence intervals.
#' @param ... Additional graphical arguments passed to `qqplot()`.
#'
#' @details
#' This function produces a Q-Q plot for a random variable following a uniform distribution, optionally on a logarithmic scale.
#'
#' For `type = "exp"`, the plot is based on exponential order statistics, which is generally more appropriate than directly log-transforming the expected uniform order statistics.
#'
#' @return
#' Invisibly returns the list produced by `qqplot()` with components:
#' \describe{
#'   \item{x}{Expected quantiles.}
#'   \item{y}{Observed quantiles.}
#' }
#'
#' @references
#' \insertRef{balakrishnan03}{gap}
#' \insertRef{casella02}{gap}
#' \insertRef{davison03}{gap}
#'
#' @seealso [`qqfun`]
#'
#' @examples
#' \dontrun{
#' u_obs <- runif(1000)
#' r <- qqunif(u_obs,pch=21,bg="blue",bty="n")
#' u_exp <- r$y
#' hits <- u_exp >= 2.30103
#' points(r$x[hits],u_exp[hits],pch=21,bg="green")
#' legend("topleft",sprintf("GC.lambda = %.4f",gc.lambda(u_obs)))
#' }
#'
#' @author Jing Hua Zhao
#' @keywords hplot distribution univar
#'
#' @export
qqunif <- function(u,type="unif",logscale=TRUE,base=10,col=palette()[4],lcol=palette()[2],ci=FALSE,alpha=0.05,...)
{
    u <- sort(u[!is.na(u)])
    if (!length(u)) stop("'u' contains no valid observations")
    if (any(u <= 0 | u >= 1)) stop("'u' must contain values strictly between 0 and 1")
    if (!is.numeric(base) || length(base) != 1 || base <= 0 || base == 1) stop("'base' must be a positive number different from 1")
    if (!type %in% c("unif","exp")) stop("invalid 'type'")
    n <- length(u)
    xlabel <- ifelse(logscale,paste("-log",base,"(Expected)",sep=""),"Expected")
    ylabel <- ifelse(logscale,paste("-log",base,"(Observed)",sep=""),"Observed")
    if (ci) crit <- qnorm(1 - alpha / 2)
    if (type == "exp") {
        n1 <- 1 / (n:1)
        n2 <- n1^2
        lambda <- 1 / log(base)
        expected <- cumsum(n1) * lambda
        observed <- -log(u,base)
        z <- qqplot(expected,observed,xlab=xlabel,ylab=ylabel,col=col,...)
        if (ci) {
            variance <- cumsum(n2) * lambda^2
            se <- sqrt(variance)
            lcl <- expected - crit * se
            ucl <- expected + crit * se
            lines(expected,lcl)
            lines(expected,ucl)
        }
    } else {
        expected <- (1:n) / (n + 1)
        if (logscale) {
            x <- -log(expected,base)
            y <- -log(u,base)
        } else {
            x <- expected
            y <- u
        }
        z <- qqplot(x,y,xlab=xlabel,ylab=ylabel,col=col,...)
        if (ci) {
            variance <- (1:n) * (n - (1:n) + 1) / ((n + 1)^2 * (n + 2))
            se <- sqrt(variance)
            lcl <- expected - crit * se
            ucl <- expected + crit * se
            lid <- lcl > 0
            uid <- ucl <= 1
            if (logscale) {
                x_lcl <- -log(expected[lid],base)
                y_lcl <- -log(lcl[lid],base)
                x_ucl <- -log(expected[uid],base)
                y_ucl <- -log(ucl[uid],base)
            } else {
                x_lcl <- expected[lid]
                y_lcl <- lcl[lid]
                x_ucl <- expected[uid]
                y_ucl <- ucl[uid]
            }
            lines(x_lcl,y_lcl)
            lines(x_ucl,y_ucl)
        }
    }
    abline(0,1,col=lcol)
    invisible(z)
}

#21-05-2021 check with ChatGPT
#09-11-2009 refine U(0,1) and exp(lambda)
#27-08-2008 first attempt of CI
