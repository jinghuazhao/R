\name{qqfun}
\alias{qqfun}
\title{Quantile-comparison plots}
\description{
 Plots empirical quantiles of a variable against theoretical quantiles of a comparison distribution.
}
\usage{
qqfun(x, distribution="norm", ylab=deparse(substitute(x)),
      xlab=paste(distribution, "quantiles"), main=NULL, las=par("las"),
      envelope=.95, labels=FALSE, col=palette()[4], lcol=palette()[2],
      xlim=NULL, ylim=NULL, lwd=1, pch=1, bg=palette()[4], cex=.4,
      line=c("quartiles", "robust", "none"), ...)
}
\arguments{
  \item{x}{vector of numeric values.}
  \item{distribution}{root name of comparison distribution -- e.g., \code{norm} for the
    normal distribution; \code{t} for the t-distribution.}
  \item{ylab}{label for vertical (empirical quantiles) axis.}
  \item{xlab}{label for horizontal (comparison quantiles) axis.}
  \item{main}{label for plot.}
  \item{envelope}{confidence level for point-wise confidence envelope, or 
    \code{FALSE} for no envelope.}
  \item{labels}{vector of point labels for interactive point identification,
    or \code{FALSE} for no labels.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{col}{color for points; the default is the \emph{fourth} entry
    in the current color palette (see \code{\link[grDevices]{palette}}
    and \code{\link[graphics]{par}}).}
  \item{lcol}{color for lines; the default is the \emph{second} entry as above.}
  \item{xlim}{the x limits (x1, x2) of the plot. Note that x1 > x2 is allowed and leads to a reversed axis.}
  \item{ylim}{the y limits of the plot}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).}
  \item{bg}{background color of points}
  \item{cex}{factor for expanding the size of plotted symbols; the default is
    \code{.4}.}
  \item{lwd}{line width; default is \code{1} (see \code{\link[graphics]{par}}).
    Confidence envelopes are drawn at half this line width.}
  \item{line}{\code{"quartiles"} to pass a line through the quartile-pairs, or
    \code{"robust"} for a robust-regression line; the latter uses the \code{rlm}
    function in the \code{MASS} package. Specifying \code{line = "none"} suppresses the line.}
  \item{\dots}{arguments such as \code{df} to be passed to the appropriate quantile function.}
}
\details{
  Draws theoretical quantile-comparison plots for variables and for studentized residuals
  from a linear model. A comparison line is drawn on the plot either through the quartiles
  of the two distributions, or by robust regression. 
  
  Any distribution for which quantile and density functions exist in R (with prefixes 
  \code{q} and \code{d}, respectively) may be used. 
  Studentized residuals are plotted against the appropriate t-distribution.

  This is adapted from qq.plot of package car with different values for points and lines,
  more options, more transparent code and examples in the current setting. Another similar but
  sophisticated function is qqmath of package lattice.
}
\value{
  \code{NULL}. These functions are used only for their side effect (to make a graph).
}
\references{
  Davison, A. C. (2003)  \emph{Statistical Models.} Cambridge University Press.

  Leemis, L. M., J. T. Mcqueston (2008) \emph{Univariate distribution relationships.}
  The American Statistician 62:45-53
}
\author{John Fox, Jing Hua Zhao}

\seealso{\code{\link[stats]{qqnorm}}, \code{\link[gap]{qqunif}}, \code{\link[gap]{gcontrol2}}}

\examples{
\dontrun{
p <- runif(100)
alpha <- 1/log(10)
qqfun(p,dist="unif")
qqfun(-log10(p),dist="exp",rate=alpha,pch=21)

#library(car)
#qq.plot(p,dist="unif")
#qq.plot(-log10(p),dist="exp",rate=alpha)

#library(lattice)
#qqmath(~ -log10(p), distribution = function(p) qexp(p,rate=alpha))
}
}
\keyword{distribution}
\keyword{univar}
\keyword{regression}
