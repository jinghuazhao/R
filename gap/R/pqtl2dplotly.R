#' 2D pQTL plotly
#'
#' @md
#' @param d Data in pqtl2dplot() format.
#' @param chrlen Lengths of chromosomes for specific build: hg18, hg19, hg38.
#' @export
#' @return A plotly figure.
#' @examples
#'
#' \dontrun{
#' INF <- Sys.getenv("INF")
#' d <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE)
#' r <- pqtl2dplotly(d)
#' htmlwidgets::saveWidget(r,file=file.path(INF,"INF1.pqtl2dplotly.html"))
#' r
#' }

pqtl2dplotly <- function(d, chrlen=gap::hg19)
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
             dplyr::select(x,y,cistrans,text)
  axes <- list(tickmode = "array",
               tick0 = 1,
               dtick = 1,
               ticklen = 1,
               tickwidth = 0,
               tickfont = list(family = "arial", size = 12, color = "#7f7f7f"),
               tickvals = tkvals,
               ticktext = tktxts)
  xaxis = c(title = "pQTL position", axes)
  yaxis = c(title = "Gene position", axes)
  cols <- c('#BF382A','#0C4B8E')
# cols <- setNames(cols, c("cis","trans"))
  fig <- with(t2d_pos, 
              plotly::plot_ly(x = ~x, y = ~y, color = ~cistrans, colors = cols, 
                              marker=list(size=11), mode="markers") %>%
              plotly::add_markers(type="scatter",text=text) %>%
              plotly::layout(xaxis=xaxis, yaxis=yaxis, showlegend = TRUE))
}
