#' 3D Manhattan plot
#'
#' @md
#' @param d Data in mhtplot2d() format.
#' @export
#' @return A plotly figure.
#' @examples
#'
#' INF <- Sys.getenv("INF")
#' d <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE)
#' r <- mhtplot3d(d)
#' r
#' htmlwidgets::saveWidget(r,file=file.path(INF,"INF1.latest.html"))

mhtplot3d <- function(d)
{
  n <- CM <- snpid <- pos_pqtl <- pos_prot <- prot_gene <- lp <- chr1 <- pos1 <- chr2 <- pos2 <- target <- gene <- log10p <- NA
  t2d <- mhtplot2d(d, plot=FALSE)
  n <- with(t2d, n)
  CM <- with(t2d, CM)
  tkvals <- tktxts <- vector()
  for (x in 1:n) {
       tkvals[x] <- ifelse(x == 1, CM[x]/2, (CM[x] + CM[x-1])/2)
       tktxts[x] <- xy(x)
  }
  t2d_pos <- with(t2d, data)
  t2d_pos <- t2d_pos %>% dplyr::mutate(snpid=paste("SNPid:",id),pos_pqtl=paste0("pQTL: ",chr1,":",pos1),
                                       pos_prot=paste0("Protein: ",chr2,":",pos2),
                                       prot_gene=paste0("target (gene):", target, "(", gene, ")"),
                                       lp=paste("-log10(P):", -log10p))
  fig <- plotly::plot_ly(t2d_pos, x = ~x, y = ~y, z = ~-log10p, color = ~col, colors = c('#BF382A', '#0C4B8E')) %>%
         plotly::add_markers(type="scatter3d", text=paste(snpid, pos_pqtl, pos_prot, prot_gene, lp, sep="\n")) %>%
         plotly::layout(scene = list(xaxis = list(title = "pQTL position",
                                                  tickmode = "array",
                                                  autotick = FALSE,
                                                  tick0 = 1,
                                                  dtick = 1,
                                                  ticklen = n,
                                                  tickwidth = 0,
                                                  tickfont = list (size = 10),
                                                  tickvals = tkvals,
                                                  ticktext = as.list(c(1:22,"X","Y"))
                                             ),
                                     yaxis = list(title = "Gene position",
                                                  tickmode = "array",
                                                  autotick = FALSE,
                                                  tick0 = 1,
                                                  dtick = 1,
                                                  ticklen = n,
                                                  tickwidth = 0,
                                                  tickfont = list (size = 10),
                                                  tickvals = tkvals,
                                                  ticktext = as.list(c(1:22,"X","Y"))
                                             ),
                                     zaxis = list(title = "-log10(p)", tickfont = list(size = 10)),
                                     aspectratio = list(x = 0.9, y = 1, z = 0.6)
                                ),
                        legend = list(x = 10, y = 0.5),
                        xaxis = list(domain=list(0,1)),
                        yaxis = list(domain=list(0,1)),
                        title = "Scatterplot of sentinels",
                        showlegend = TRUE
         )
}
