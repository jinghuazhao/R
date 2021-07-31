#' genomewide plot of CNVs
#'
#' The function generates a plot containing genomewide copy number variants (CNV) chr, start, end, freq(uencies).
#' @md
#' @param data Data to be used.
#' @export
#' @return None.
#' @examples
#' knitr::kable(cnv,caption="A CNV dataset")
#' cnvplot(cnv)

cnvplot <- function(data)
# cnvplot(cnv)
{
  d <- within(data,{chr<-replace(chr,chr=="X",23); chr<-replace(chr,chr=="Y",24)})
  pos <- vector("numeric")
  n <- length(table(with(data,chr)))
  for (x in 1:n) pos[x] <- with(subset(d,chr==paste(x)),{max(end)})
  CM <- cumsum(pos)
  par(xaxt = "n", yaxt = "n")
  xy <- xy.coords(c(0,CM), seq(1,90,by=90/(1+n)))
  plot(xy$x, xy$y, type = "n", ann = FALSE, axes = FALSE)
  colors <- rep(c("red","blue"),n)
  par(xaxt = "s", yaxt = "s", xpd = TRUE)
  xy <- function(x) if (x<23) x else if (x==23) "X" else if (x==24) "Y";
  for (x in 1:n) with(subset(d,chr==paste(x)), {
      l <- ifelse(x==1,0,CM[x-1])
      segments(l+start,freq,l+end,freq,lwd="3",col=colors[x])
      text(ifelse(x == 1, CM[x]/2, (CM[x] + CM[x-1])/2), 0, pos = 1, offset = 0.5, xy(x), cex=0.4)
  })
  segments(0,0,CM[x],0)
  axis(2,line=-0.5)
  title(xlab="Chromosome",ylab="Frequency",line=2)
}
