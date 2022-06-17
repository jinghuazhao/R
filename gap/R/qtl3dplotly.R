#' 3D QTL plot
#'
#' @md
#' @param d Data in qtl2d() format.
#' @param chrlen Lengths of chromosomes for specific build: hg18, hg19, hg38.
#' @param zmax Maximum value (e.g., -log10p) to truncate, above which they would be set to this value.
#' @param qtl.id QTL id.
#' @param qtl.prefix QTL prefix.
#' @param qtl.gene QTL target gene.
#' @param target.type Type of target, e.g., protein.
#' @param TSS to use TSS when TRUE.
#' @param xlab X-axis title.
#' @param ylab Y-axis title.
#' @param ... Additional arguments, e.g., to qtl2dplot().
#'
#' @export
#' @return A plotly figure.
#'
#' @examples
#' \dontrun{
#' suppressMessages(library(dplyr))
#' INF <- Sys.getenv("INF")
#' d <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE) %>%
#'      mutate(log10p=-log10p)
#' r <- qtl3dplotly(d,zmax=300)
#' htmlwidgets::saveWidget(r,file=file.path(INF,"INF1.qtl3dplotly.html"))
#' r
#' }

qtl3dplotly <- function(d, chrlen=gap::hg19, zmax=300, qtl.id="SNPid:", qtl.prefix="QTL:", qtl.gene="Gene:", target.type="Protein",
                        TSS=FALSE, xlab="QTL position", ylab="Gene position",...)
{
  n <- CM <- snpid <- pos_qtl <- pos_gene <- target_gene <- lp <- chr1 <- pos1 <- chr2 <- pos2 <- target <- gene <- value <- cistrans <- y <- NA
  t2d <- qtl2dplot(d, chrlen, TSS=TSS, plot=FALSE, ...)
  n <- with(t2d, n)
  CM <- with(t2d, CM)
  tkvals <- tktxts <- vector()
  for (x in 1:n) {
       tkvals[x] <- ifelse(x == 1, CM[x]/2, (CM[x] + CM[x-1])/2)
       tktxts[x] <- xy(x)
  }
  t2d_pos <- with(t2d, data) %>%
             dplyr::mutate(snpid=paste(qtl.id,id),pos_qtl=paste0(qtl.prefix,chr1,":",pos1),
                           pos_gene=paste0(qtl.gene,chr2,":",pos2),
                           target_gene=paste0(target.type," (gene):", target, "(", gene, ")"),
                           lp=paste("value:", value),
                           text=paste(snpid, pos_qtl, pos_gene, target_gene, lp, sep="\n")) %>%
             dplyr::mutate(z=if_else(value<=zmax,value,zmax)) %>%
             dplyr::select(x,y,z,cistrans,text)
  fig <- with(t2d_pos,
         plotly::plot_ly(t2d_pos, x = ~x, y = ~y, z = ~z, color = ~cistrans, colors = c('#BF382A', '#0C4B8E')) %>%
         plotly::add_markers(type="scatter3d", text=text) %>%
         plotly::layout(scene = list(xaxis = list(title = xlab,
                                                  tickmode = "array",
                                                  autotick = FALSE,
                                                  tick0 = 1,
                                                  dtick = 1,
                                                  ticklen = n,
                                                  tickwidth = 0,
                                                  tickfont = list (size = 10),
                                                  tickvals = tkvals,
                                                  ticktext = tktxts),
                                     yaxis = list(title = ylab,
                                                  tickmode = "array",
                                                  autotick = FALSE,
                                                  tick0 = 1,
                                                  dtick = 1,
                                                  ticklen = n,
                                                  tickwidth = 0,
                                                  tickfont = list (size = 10),
                                                  tickvals = tkvals,
                                                  ticktext = tktxts),
                                     zaxis = list(title = "-log10(p)", tickfont = list(size = 10)),
                                     aspectratio = list(x = 0.9, y = 1, z = 0.6)
                                ),
                        title = "Scatterplot of sentinels",
                        showlegend = TRUE
         )
  )
}
