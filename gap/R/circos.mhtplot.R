#' circos Manhattan plot with gene annotation
#'
#' The function generates circos Manhattan plot with gene annotation.
#'
#' @md
#' @param data Data to be used.
#' @param glist A gene list.
#' @export
#' @return
#' None.
#' @examples
#' \dontrun{
#' require(gap.datasets)
#' glist <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3","WISP3",
#'            "PPP1R3B","RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
#' circos.mhtplot(mhtdata,glist)
#' }
#'

circos.mhtplot <- function(data, glist)
# circos.mhtplot(mhtdata,g)
{
  pos <- gene <- NULL
  for(p in c("circlize")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("circos.mhtplot needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  requireNamespace("circlize")
  d <- within(data, {
    chr <- paste0("chr",chr)
    start <- pos - 1
    end <- pos
  })[c("chr","start","end","p","gene")]
  hd <-	subset(data, gene %in% glist)[c("chr","start","end","gene")]
  hd <-	within(hd, {chr	<- paste0("chr",chr)})
  ann <- data.frame()
  for (g in glist)
  {
     m <- subset(hd, gene==g)
     n <- round(nrow(m) / 2 + 0.5)
     ann <- rbind(ann,m[n,])
  }
  circlize::circos.par(start.degree = 90, track.height = 0.4, cell.padding = c(0, 0, 0, 0))
  circlize::circos.initializeWithIdeogram(species = "hg18", track.height = 0.05, ideogram.height = 0.06)
  circlize::circos.genomicTrackPlotRegion(d[c("chr","start","end","p")], ylim = c(0, 15), 
         panel.fun = function(region, value, ...) {
           color <- as.numeric(gsub("chr", "", circlize::get.current.chromosome()))
           with(cbind(region, value), circlize::circos.genomicPoints(region, -log10(value), cex=0.3, col = color))
  })
  circlize::circos.genomicLabels(ann, labels.column = 4, side = "inside")
  circlize::circos.clear()
}
