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
#' # A qtlClaaifier/qtl2dplot coupling example:
#' ucsc_modified <- bind_rows(ucsc,APOC,AMY,C4B,HIST,HBA)
#' pqtls <- select(merged,prot,SNP,log.P.) %>%
#'          mutate(log10p=-log.P.) %>%
#'          left_join(caprion_modified) %>%
#'          select(Gene,SNP,prot,log10p)
#' posSNP <- select(merged,SNP,Chr,bp)
#' cis.vs.trans <- qtlClassifier(pqtls,posSNP,ucsc_modified,1e6) %>%
#'                 mutate(geneChrom=as.integer(geneChrom),
#'                        cis=if_else(Type=="cis",TRUE,FALSE))
#' head(cis.vs.trans)
#'     Gene             SNP  prot log10p geneChrom geneStart  geneEnd SNPChrom    SNPPos  cis
#'    YWHAB 8:111907280_A_T 1433B   7.38        20  43530174 43535076        8 111907280 FALSE
#'      A2M 14:34808001_A_T  A2MG   7.51        12   9220421  9268445       14  34808001 FALSE
#'     APEH  1:12881809_A_G  ACPH   7.83         3  49711834 49720772        1  12881809 FALSE
#'      PGD 2:121896327_A_G  6PGD   7.79         1  10459174 10479803        2 121896327 FALSE
#' SERPINF2  17:1618262_C_T  A2AP  12.59        17   1648289  1657825       17   1618262  TRUE
#'     PGLS 19:54327869_G_T  6PGL   9.87        19  17622481 17631887       19  54327869 FALSE
#' qtl2dplot(cis.vs.trans,chrlen=gap::hg19,
#'           snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
#'           gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",
#'           trait="prot",gene="Gene",
#'           TSS=TRUE,cis="cis",plot=TRUE,cex.labels=0.6,cex.points=0.6,
#'           xlab="pQTL position",ylab="Gene position")
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
