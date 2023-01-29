#' 2D QTL plot
#'
#' @param d Data to be used.
#' @param chrlen lengths of chromosomes for specific build: hg18, hg19, hg38.
#' @param snp_name variant name.
#' @param snp_pos variant position.
#' @param snp_chr variant chromosome.
#' @param gene_chr gene chromosome.
#' @param gene_start gene start position.
#' @param gene_end gene end position.
#' @param trait trait name.
#' @param gene gene name.
#' @param TSS to use TSS when TRUE.
#' @param cis cis variant when TRUE.
#' @param value A specific value to show.
#' @param plot to plot when TRUE.
#' @param cex.labels Axis label extension factor.
#' @param cex.points Data point extension factor.
#' @param xlab X-axis title.
#' @param ylab Y-axis title.
#'
#' @details
#' This function is both used as its own for a 2d plot and/or generate data for a plotly counterpart.
#'
#' @export
#' @return positional information.
#' @examples
#' \dontrun{
#' INF <- Sys.getenv("INF")
#' d <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE)
#' r <- qtl2dplot(d)
#' }

qtl2dplot <- function(d, chrlen=gap::hg19, snp_name="SNP", snp_chr="Chr", snp_pos="bp",
                      gene_chr="p.chr", gene_start="p.start", gene_end="p.end",
                      trait="p.target.short", gene="p.gene", TSS=FALSE,
                      cis="cis",value="log10p",
                      plot=TRUE,
                      cex.labels=0.6, cex.points=0.6, xlab="QTL position", ylab="Gene position")
{
  r <- grid2d(chrlen, plot=plot, cex.labels=cex.labels, xlab=xlab, ylab=ylab)
  n <- with(r, n)
  CM <- with(r, CM)
  chr1 <- d[[snp_chr]]
  chr1[chr1=="X"] <- 23
  chr1[chr1=="Y"] <- 24
  pos1 <- CM[chr1] + d[[snp_pos]]
  chr2 <- d[[gene_chr]]
  chr2[chr2=="X"] <- 23
  chr2[chr2=="Y"] <- 24
  pos <- (d[[gene_start]] + d[[gene_end]])/2
  if (TSS) pos <- d[[gene_start]]
  pos2 <- CM[chr2] + pos
  if (plot) {
     points(pos1, pos2, cex=cex.points, col=ifelse(d[[cis]],"red","blue"), pch=19)
     legend("top", legend=c("cis","trans"), box.lty=0, cex=cex.points, col=c("red","blue"),
            horiz=TRUE, inset=c(0,1), xpd=TRUE, pch=19)
  }
  return(list(n=n, CM=CM, data=data.frame(id=d[[snp_name]],
                                          chr1=chr1, pos1=d[[snp_pos]],
                                          chr2=chr2, pos2=pos,
                                          x=pos1, y=pos2, value=d[[value]],
                                          target=d[[trait]], gene=d[[gene]],
                                          cistrans=ifelse(d[[cis]],"cis","trans")
  )))
}
