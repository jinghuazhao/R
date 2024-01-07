#' Another circos Manhattan plot
#'
#' @param dat Data to be plotted with variables chr, pos, log10p.
#' @param labs Data on labels.
#' @param species Genome build.
#' @param ticks Tick positions.
#' @param ymax maximum value for y-axis.
#'
#' @details
#' This is adapted from work for a recent publication. It enables a y-axis to the -log10(P) for association statistics
#'
#' @export
#' @return
#' There is no return value but a plot.
#'
#' @examples
#' \dontrun{
#' require(gap.datasets)
#' library(dplyr)
#' glist <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3","WISP3",
#'            "PPP1R3B","RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
#' testdat <- mhtdata[c("chr","pos","p","gene","start","end")] %>%
#'            rename(log10p=p) %>%
#'            mutate(chr=paste0("chr",chr),log10p=-log10(log10p))
#' dat <- mutate(testdat,start=pos,end=pos) %>%
#'        select(chr,start,end,log10p)
#' labs <- subset(testdat,gene %in% glist) %>%
#'         group_by(gene,chr,start,end) %>%
#'         summarize() %>%
#'         mutate(cols="blue") %>%
#'         select(chr,start,end,gene,cols)
#' labs[2,"cols"] <- "red"
#' ticks <- 0:3*5
#' circos.mhtplot2(dat,labs,ticks=ticks,ymax=max(ticks))
#' # https://www.rapidtables.com/web/color/RGB_Color.html
#' }

circos.mhtplot2 <- function(dat,labs,species="hg18",ticks=0:3*10,ymax=30)
{
  requireNamespace("circlize")
  circlize::circos.clear()
  chrs <- with(dat,table(chr))
  chr.index <- paste0("chr",sort(as.numeric(substring(names(chrs),4))))
  circlize::circos.par("start.degree" = 90, gap.degree = c(rep(c(0.7), length(chrs)-1), 8),
                        track.margin = c(0.005, 0.005),
                        cell.padding = c(0.001, 0.01, 0.01, 0.001))
  circlize::circos.initializeWithIdeogram(plotType = NULL, species = species, chromosome.index = chr.index)
  circlize::circos.genomicLabels(labs, labels=labs[["gene"]], side = "outside",
                                 cex = 0.7, font = 4, line_lwd = 0.7, padding=0.5,
                                 connection_height = circlize::convert_height(8, "mm"),
                                 line_col = labs[["cols"]],
                                 col = labs[["cols"]])
  circlize::circos.track(ylim = c(0, 1),
               panel.fun = function(x, y) {
                 chr  = substring(circlize::CELL_META$sector.index, 4)
                 xlim = circlize::CELL_META$xlim
                 ylim = circlize::CELL_META$ylim
                 circlize::circos.rect(xlim[1], 0, xlim[2], 1, col = "white", cex = 0.2, lwd = 0.5 )
                 circlize::circos.text(mean(xlim), mean(ylim), chr, cex = 0.4, col = "black", facing = "inside", niceFacing = TRUE)
               },
               track.height = 0.03,  bg.border = NA)
  circlize::circos.track(ylim=c(0,1), track.height=0.05, bg.border=NA,
               panel.fun=function(x, y) {
                 chr=gsub("chr", "", circlize::CELL_META$sector.index)
                 xlim=circlize::CELL_META$xlim
                 ylim=circlize::CELL_META$ylim
                 circlize::circos.genomicAxis(h="top", direction="inside", labels.cex=0.25, major.at=seq(0,1e10,5e7))
               })
  circlize::circos.genomicTrackPlotRegion(dat, numeric.column = 4, panel.fun = function(region, value,  ...)
                 circlize::circos.genomicPoints(region, value, pch = 16, col = "magenta", cex = 0.3),
                                                track.height = 0.55, bg.border = 1, bg.col = "white", ylim = c(0, ymax))
  circlize::circos.yaxis(side = "left", at = ticks, labels = ticks,
              sector.index = circlize::get.all.sector.index()[1], labels.cex = 0.3, lwd = 0.3,
              tick.length = 0.5*(circlize::convert_x(1, "mm", circlize::CELL_META$sector.index,circlize::CELL_META$track.index)))
  circlize::circos.genomicText(data.frame(start=1,end=1),sector.index=circlize::get.all.sector.index()[1],
                               labels = "-log10(P)", h = "bottom", cex = 0.6, font = 2, y = 2/3 * ymax, adj = c(0.2, 1.5), facing = "clockwise")
}
