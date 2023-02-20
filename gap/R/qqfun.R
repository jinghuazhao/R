#' Quantile-comparison plots
#'
#' @param x vector of numeric values.
#' @param distribution root name of comparison distribution -- e.g., `norm` for the normal distribution; `t` for the t-distribution.
#' @param ylab label for vertical (empirical quantiles) axis.
#' @param xlab label for horizontal (comparison quantiles) axis.
#' @param main label for plot.
#' @param envelope confidence level for point-wise confidence envelope, or `FALSE` for no envelope.
#' @param labels vector of point labels for interactive point identification, or `FALSE` for no labels.
#' @param las if `0`, ticks labels are drawn parallel to the axis; set to `1` for horizontal labels (see [`graphics::par`]).
#' @param col color for points; the default is the \emph{fourth} entry in the current color palette (see [`grDevices::palette`] and [`graphics::par`]).
#' @param lcol color for lines; the default is the \emph{second} entry as above.
#' @param xlim the x limits (x1, x2) of the plot. Note that x1 > x2 is allowed and leads to a reversed axis.
#' @param ylim the y limits of the plot.
#' @param pch plotting character for points; default is `1` (a circle, see [`graphics::par`]).
#' @param bg background color of points.
#' @param cex factor for expanding the size of plotted symbols; the default is `.4.
#' @param lwd line width; default is `1` (see [`graphics::par`]). Confidence envelopes are drawn at half this line width.
#' @param line `"quartiles"` to pass a line through the quartile-pairs, or `"robust"` for a robust-regression line; the latter uses the `rlm` function in the `MASS` package. Specifying `line = "none"` suppresses the line.
#' @param \dots arguments such as \code{df} to be passed to the appropriate quantile function.
#'
#' @details
#' Plots empirical quantiles of a variable against theoretical quantiles of a comparison distribution.
#'
#' Draws theoretical quantile-comparison plots for variables and for studentized residuals from a linear model. A comparison line is drawn on the plot either through the quartiles of the two distributions, or by robust regression.
#'  
#' Any distribution for which quantile and density functions exist in R (with prefixes `q` and `d`, respectively) may be used. Studentized residuals are plotted against the appropriate t-distribution.
#'
#' This is adapted from [`car::qq.plot`] with different values for points and lines, more options, more transparent code and examples in the current setting. Another similar but sophisticated function is [`lattice::qqmath`].
#'
#' @export
#' @return
#' These functions are used only for their side effect (to make a graph).
#'
#' @references
#' \insertRef{davison03}{gap}
#'
#' \insertRef{leemis08}{gap}
#'
#' @author John Fox, Jing Hua Zhao
#' @seealso [`stats::qqnorm`], [`qqunif`], [`gcontrol2`]
#'
#' @examples
#' \dontrun{
#' p <- runif(100)
#' alpha <- 1/log(10)
#' qqfun(p,dist="unif")
#' qqfun(-log10(p),dist="exp",rate=alpha,pch=21)

#' library(car)
#' qq.plot(p,dist="unif")
#' qq.plot(-log10(p),dist="exp",rate=alpha)
#'
#' library(lattice)
#' qqmath(~ -log10(p), distribution = function(p) qexp(p,rate=alpha))
#' }
#'
#' @keywords distribution univar regression

qqfun <- function(x, distribution="norm", ylab=deparse(substitute(x)),
            xlab=paste(distribution, "quantiles"), main=NULL, las=par("las"),
            envelope=.95, labels=FALSE, col=palette()[4], lcol=palette()[2], 
            xlim=NULL, ylim=NULL, lwd=1, pch=1, bg=palette()[4], cex=.4,
            line=c("quartiles", "robust", "none"), ...)
{
    result <- NULL
    line <- match.arg(line)
    good <- !is.na(x)
    ord <- order(x[good])
    ord.x <- x[good][ord]
    qfun <- eval(parse(text=paste("q",distribution,sep="")))
    dfun <- eval(parse(text=paste("d",distribution,sep="")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- qfun(P, ...)
    plot(z, ord.x, xlab=xlab, ylab=ylab, main=main, las=las, col=col, pch=pch, cex=cex, bg=bg, xlim=xlim, ylim=ylim)
    if (line=="quartiles") {
        Qx <- quantile(ord.x, c(.25,.75))
        Qz <- qfun(c(.25,.75), ...)
        b <- (Qx[2]-Qx[1])/(Qz[2]-Qz[1])
        a <- Qx[1]-b*Qz[1]
        abline(a, b, col=lcol, lwd=lwd)
    }
    if (line=="robust") {
        for(p in c("MASS")) {
           if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
              if (!requireNamespace(p, quietly = TRUE))
              warning(paste("qqfun needs package `", p, "' to be fully functional; please install", sep=""))
           }
        }
        coef <- coefficients(MASS::rlm(ord.x~z))
        a <- coef[1]
        b <- coef[2]
        abline(a,b,col=palette()[2])
    }
    if (line != 'none' & envelope != FALSE) {
        zz <- qnorm(1-(1-envelope)/2)
        SE <- (b/dfun(z, ...))*sqrt(P*(1-P)/n)
        fit.value <- a+b*z
        upper <- fit.value+zz*SE
        lower <- fit.value-zz*SE
        lines(z, upper, lty=2, lwd=lwd/2, col=lcol)
        lines(z, lower, lty=2, lwd=lwd/2, col=lcol)
    }
    if (labels[1]==TRUE & length(labels)==1) labels<-seq(along=z)
    if (labels[1] != FALSE) {
        selected <- identify(z, ord.x, labels[good][ord])
        result <- seq(along=x)[good][ord][selected]
    }
    if (is.null(result)) invisible(result) else sort(result)
}
