#' Conversion of chromosome names to strings
#'
#' This function converts x=1:24 to 1:22, X, Y
#'
#' @param x (alpha)numeric value indicating chromosome.
#'
#' @export
#' @return As indicated.

xy <- function(x) if (x<23) x else if (x==23) "X" else if (x==24) "Y"

#' Conversion of chrosome name from strings
#'
#' This function converts 1:22, X, Y back to 1:24.
#'
#' @param x Chromosome name in strings
#'
#' @export
#' @return As indicated.

ixy <- function(x) if (x=="X") 23 else if (x=="Y") 24 else x

#' Two-dimensional grid
#'
#' This function build 2-d grids
#'
#' @param chrlen Lengths of chromosomes; e.g., hg18, hg19 or hg38.
#' @param plot A flag for plot.
#' @param cex.labels A scaling factor for labels.
#' @param xlab X-axis title.
#' @param ylab Y-axis title.
#'
#' @export
#' @return A list with two variables.
#' \describe{
#'   \item{n}{Number of chromosomes.}
#'   \item{CM}{Cumulative lengths starting from 0.}
#' }

grid2d <- function(chrlen, plot=TRUE, cex.labels=0.6, xlab="QTL position", ylab="Gene position")
{
  CM <- cumsum(chrlen)
  n <- length(chrlen)
  xy <- xy.coords(c(0,CM), c(0,CM))
  if (plot)
  {
    par(xaxt = "n", yaxt = "n")
    plot(xy$x, xy$y, type = "n", ann = FALSE, axes = FALSE)
    par(xaxt = "s", yaxt = "s", xpd = TRUE)
    for (x in 1:n) {
        segments(CM[x],0,CM[x],CM[n],col="black")
        segments(0,CM[x],CM[n],CM[x],col="black")
        text(ifelse(x == 1, CM[x]/2, (CM[x] + CM[x-1])/2), 0, pos = 1, offset = 0.5, xy(x), cex=cex.labels)
        text(0, ifelse(x == 1, CM[x]/2, (CM[x] + CM[x-1])/2), pos = 2, offset = 0.5, xy(x), cex=cex.labels)
    }
    segments(0,0,CM[n],0)
    segments(0,0,0,CM[n])
    title(xlab=xlab,ylab=ylab,line=2)
  }
  invisible(list(n=n, CM=c(0,CM)))
}
