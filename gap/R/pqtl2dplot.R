#' 2D pQTL plot
#'
#' @md
#' @param d Data to be used.
#' @param chrlen lengths of chromosomes for specific build: hg18, hg19, hg38.
#' @param snp_name variant name.
#' @param snp_pos variant position.
#' @param snp_chr variant chromosome.
#' @param gene_chr gene chromosome.
#' @param gene_start gene start position.
#' @param gene_end gene end position.
#' @param protein protein name.
#' @param gene gene name.
#' @param lp log10(p).
#' @param cis cis variant when TRUE.
#' @param plot to plot when TRUE.
#' @param cex extension factor.
#' @export
#' @return positional information.
#' @examples
#' \dontrun{
#' INF <- Sys.getenv("INF")
#' d <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE)
#' r <- pqtl2dplot(d)
#' }

pqtl2dplot <- function(d, chrlen=gap::hg19, snp_name="SNP", snp_chr="Chr", snp_pos="bp",
                       gene_chr="p.chr", gene_start="p.start", gene_end="p.end",
                       protein="p.target.short", gene="p.gene", lp="log10p",
                       cis="cis",
                       plot=TRUE, cex=0.6)
{
  r <- grid2d(chrlen, plot=plot)
  n <- with(r, n)
  CM <- with(r, CM)
  chr1 <- d[[snp_chr]]
  chr1[chr1=="X"] <- 23
  chr1[chr1=="Y"] <- 24
  pos1 <- CM[chr1] + d[[snp_pos]]
  chr2 <- d[[gene_chr]]
  chr2[chr2=="X"] <- 23
  chr2[chr2=="Y"] <- 24
  mid <- (d[[gene_start]] + d[[gene_end]])/2
  pos2 <- CM[chr2] + mid
  if (plot) {
     points(pos1, pos2, cex=cex, col=ifelse(d[[cis]],"red","blue"), pch=19)
     legend("top", legend=c("cis","trans"), box.lty=0, cex=cex, col=c("red","blue"),
            horiz=TRUE, inset=c(0,1), xpd=TRUE, pch=19)
  }
  return(list(n=n, CM=CM, data=data.frame(id=d[[snp_name]],
                                          chr1=chr1, pos1=d[[snp_pos]],
                                          chr2=chr2, pos2=mid, x=pos1, y=pos2,
                                          target=d[[protein]], gene=d[[gene]], log10p=d[[lp]],
                                          cistrans=ifelse(d[[cis]],"cis","trans")
  )))
}
