#' 3D pQTL plot
#'
#' @md
#' @param d Data in pqtl2d() format.
#' @param chrlen Lengths of chromosomes for specific build: hg18, hg19, hg38.
#' @param zmax Maximum -log10p to truncate, above which they would be set to this value.
#' @export
#' @return A plotly figure.
#' @examples
#'
#' \dontrun{
#' INF <- Sys.getenv("INF")
#' d <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE)
#' r <- pqtl3dplotly(d,zmax=300)
#' htmlwidgets::saveWidget(r,file=file.path(INF,"INF1.pqtl3dplotly.html"))
#' r
#' }

pqtl3dplotly <- function(d, chrlen=gap::hg19, zmax=300)
{
  n <- CM <- snpid <- pos_pqtl <- pos_prot <- prot_gene <- lp <- chr1 <- pos1 <- chr2 <- pos2 <- target <- gene <- log10p <- cistrans <- y <- NA
  t2d <- pqtl2dplot(d, chrlen, plot=FALSE)
  n <- with(t2d, n)
  CM <- with(t2d, CM)
  tkvals <- tktxts <- vector()
  for (x in 1:n) {
       tkvals[x] <- ifelse(x == 1, CM[x]/2, (CM[x] + CM[x-1])/2)
       tktxts[x] <- xy(x)
  }
  t2d_pos <- with(t2d, data) %>%
             dplyr::mutate(snpid=paste("SNPid:",id),pos_pqtl=paste0("pQTL: ",chr1,":",pos1),
                           pos_prot=paste0("Protein: ",chr2,":",pos2),
                           prot_gene=paste0("target (gene):", target, "(", gene, ")"),
                           lp=paste("-log10(P):", -log10p),
                           text=paste(snpid, pos_pqtl, pos_prot, prot_gene, lp, sep="\n")) %>%
             dplyr::mutate(z=if_else(-log10p<=zmax,-log10p,zmax)) %>%
             dplyr::select(x,y,z,cistrans,text)
  fig <- with(t2d_pos,
         plotly::plot_ly(t2d_pos, x = ~x, y = ~y, z = ~z, color = ~cistrans, colors = c('#BF382A', '#0C4B8E')) %>%
         plotly::add_markers(type="scatter3d", text=text) %>%
         plotly::layout(scene = list(xaxis = list(title = "pQTL position",
                                                  tickmode = "array",
                                                  autotick = FALSE,
                                                  tick0 = 1,
                                                  dtick = 1,
                                                  ticklen = n,
                                                  tickwidth = 0,
                                                  tickfont = list (size = 10),
                                                  tickvals = tkvals,
                                                  ticktext = tktxts),
                                     yaxis = list(title = "Gene position",
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
