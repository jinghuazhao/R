#' circos plot of CNVs.
#'
#' @param data CNV data containing chromosome, start, end and freq.
#'
#' @details
#' The function plots frequency of CNVs.
#'
#' @export
#' @return None.
#' @examples
#' \dontrun{
#' circos.cnvplot(cnv)
#' }

circos.cnvplot <- function(data)
{
  for(p in c("circlize")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("circos.cnvplot needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  start <- end <- NA
  cnv <- within(subset(data,start<=end),{chr=paste0("chr",chr)})
  requireNamespace("circlize")
  circlize::circos.par(start.degree = 50, track.height = 0.3, cell.padding = c(0, 0, 0, 0))
  circlize::circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
  circlize::circos.genomicTrackPlotRegion(cnv, ylim = c(0, 50), panel.fun = function(region,value,...) {
                      color <- as.numeric(gsub("chr","",circlize::get.current.chromosome()))
                      with(cbind(region,value),circlize::circos.segments(start,freq,end,freq,col=color,lwd=1))
  })
  circlize::circos.clear()
}
