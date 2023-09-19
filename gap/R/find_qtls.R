#' Distance-based signal identification
#'
#' @param d input data.
#' @param Chromosome chromosome.
#' @param Position position.
#' @param MarkerName RSid or SNPid.
#' @param Allele1 effect allele.
#' @param Allele2 other allele.
#' @param EAF effect allele frequency.
#' @param Effect b.
#' @param StdErr SE.
#' @param log10P -log(P).
#' @param N sample size.
#' @param radius a flanking distance.
#'
#' @details
#' This function implements an iterative merging algorithm to identify signals.
#' The setup follows output from METAL.
#'
#' @export
#' @return The function lists QTLs and meta-information.
#' @examples
#' \dontrun{
#'   f <- "ZPI_dr.p.gz"
#'   varlist=c("Chromosome","Position","MarkerName","Allele1","Allele2",
#'             "Freq1","FreqSE","MinFreq","MaxFreq",
#'             "Effect","StdErr","log10P","Direction",
#'             "HetISq","HetChiSq","HetDf","logHetP","N")
#'   d <- read.table(f,col.names=varlist,check.names=FALSE)
#'   find_qtls(d)
#' }

find_qtls <- function(d, Chromosome="Chromosome",Position="Position",
                         MarkerName="MarkerName",Allele1="Allele1",Allele2="Allele2",
                         EAF="Freq1",
                         Effect="Effect",StdErr="StdErr",log10P="log10P",
                         N="N", radius=1e6)
{
  for (q in c("dplyr","valr")) {
     if (length(grep(paste("^package:", q, "$", sep=""), search())) == 0) {
        if (!requireNamespace(q, quietly = TRUE))
        warning(paste("mhtplot.trunc needs package `", q, "' to be fully functional; please install", sep=""))
     }
  }
  chrom <- start.gene <- end.gene <- geneStart <- geneEnd <- mlog10p <- NULL
  p <- d %>%
       dplyr::transmute(chrom=paste0("chr",d[[Chromosome]]),start=d[[Position]],end=d[[Position]],
                        rsid=d[[MarkerName]],a1=d[[Allele1]],a2=d[[Allele2]],
                        EAF=d[[EAF]],b=d[[Effect]],SE=d[[StdErr]],mlog10p=-d[[log10P]],n=d[[N]])
  m <- valr::bed_merge(p,max_dist=radius)
  i <- valr::bed_intersect(p,m,suffix = c("", ".gene")) %>%
       dplyr::rename(geneStart=start.gene,geneEnd=end.gene)
  pQTLs <- i %>%
           dplyr::group_by(chrom,geneStart,geneEnd) %>%
           dplyr::slice_max(mlog10p)
}
