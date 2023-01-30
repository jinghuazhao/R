#' RLE plot
#' @noRd

makeRLEplot <- function(E, log2.data=TRUE, groups=NULL, col.group=NULL,
                        showTitle=FALSE, title="Relative log expression (RLE) plot",...){
  
  # log2 the data
  if (log2.data==TRUE){
    message("log2'ing the data")
    E <- log2(E)
  }
  
  # handle the colours
  if (is.null(groups)) mycol <- NULL
  else {
    if (length(levels(groups)) > length(col.group)) {
      warning("More groups than colours: will plot plain")
      mycol <- NULL
    } else mycol <- col.group[groups]
  }
  
  # calculate the median of log2 counts for each gene
  g.medians <- apply(E, 1, FUN = median, na.rm = TRUE)
  
  # substract the median of each gene
  E.new <- sweep(E, MARGIN=1, STATS=g.medians, FUN='-')

  # generate boxplots for all individuals
  boxplot(E.new, xaxt="n", las=2, col= mycol, ...)
  xtick <- seq(1, ncol(E.new), by=1)
  axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
  text(x=xtick, par("usr")[3], labels = colnames(data),
       srt=90, pos=1, xpd = TRUE, col = mycol, ...)
  if (showTitle) mtext(text=title, side=3, line=0.1)
  
}
