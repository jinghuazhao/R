#' 2D QTL plotly
#'
#' @param d Data in qtl2dplot() format.
#' @param chrlen Lengths of chromosomes for specific build: hg18, hg19, hg38.
#' @param qtl.id QTL id.
#' @param qtl.prefix QTL prefix.
#' @param qtl.gene QTL gene.
#' @param target.type Type of target, e.g., protein.
#' @param TSS to use TSS when TRUE.
#' @param xlab X-axis title.
#' @param ylab Y-axis title.
#' @param ... Additional arguments, e.g., target, log10p, to qtl2dplot.
#'
#' @export
#' @return A plotly figure.
#'
#' @examples
#' \dontrun{
#' INF <- Sys.getenv("INF")
#' d <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE)
#' r <- qtl2dplotly(d)
#' htmlwidgets::saveWidget(r,file=file.path(INF,"INF1.qtl2dplotly.html"))
#' r
#' }

qtl2dplotly <- function(d, chrlen=gap::hg19, qtl.id="SNPid:", qtl.prefix="QTL:", qtl.gene="Gene:", target.type="Protein",
                        TSS=FALSE, xlab="QTL position", ylab="Gene position",...)
{
  n <- CM <- snpid <- pos_qtl <- pos_gene <- target_gene <- lp <- NA
  t2d_id <- id <- chr1 <- pos1 <- chr2 <- pos2 <- target <- gene <- value <- cistrans <- y <- NA
  t2d <- qtl2dplot(d, chrlen, TSS=TSS, plot=FALSE, ...)
  names(t2d$data)[names(t2d$data) == "id"] <- "t2d_id"
  n <- with(t2d, n)
  CM <- with(t2d, CM)
  tkvals <- tktxts <- vector()
  for (x in 1:n) {
       tkvals[x] <- ifelse(x == 1, CM[x]/2, (CM[x] + CM[x-1])/2)
       tktxts[x] <- xy(x)
  }
  t2d_pos <- with(t2d, data) %>%
             dplyr::mutate(snpid=paste(qtl.id,t2d_id),
                           pos_qtl=paste0(qtl.prefix,chr1,":",pos1),
                           pos_gene=paste0(qtl.gene,chr2,":",pos2),
                           target_gene=paste0(target.type, " (gene): ", target, " (", gene, ")"),
                           lp=paste("value:", value),
                           text=paste(snpid, pos_qtl, pos_gene, target_gene, lp, sep="\n")) %>%
             dplyr::select(x,y,cistrans,text)
  axes <- list(tickmode = "array",
               tick0 = 1,
               dtick = 1,
               ticklen = 1,
               tickwidth = 0,
               tickfont = list(family = "arial", size = 12, color = "#7f7f7f"),
               tickvals = tkvals,
               ticktext = tktxts)
  xaxis = c(title = xlab, axes)
  yaxis = c(title = ylab, axes)
  cols <- c('#BF382A','#0C4B8E')
# cols <- setNames(cols, c("cis","trans"))
  fig <- with(t2d_pos, 
              plotly::plot_ly(x = ~x, y = ~y, color = ~cistrans, colors = cols, 
                              marker=list(size=11), mode="markers") %>%
              plotly::add_markers(type="scatter",text=text) %>%
              plotly::layout(xaxis=xaxis, yaxis=yaxis, showlegend = TRUE))
}
