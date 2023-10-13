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
#' @param collapse.hla a flag to collapse signals in the HLA region.
#' @param build genome build to define the HLA region.
#'
#' @details
#' This function implements an iterative merging algorithm to identify signals.
#' The setup follows output from METAL. When collapse.hla=TRUE, a single most
#' significant signal in the HLA region is chosen. The Immunogenomics paper gives
#' hg19/GRCh37: chr6:28477797-33448354 (6p22.1-21.3),
#' hg38/GRCh38: chr6:28510020-33480577.
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
#'   qtlFinder(d)
#' }

qtlFinder <- function(d, Chromosome="Chromosome",Position="Position",
                         MarkerName="MarkerName",Allele1="Allele1",Allele2="Allele2",
                         EAF="Freq1",
                         Effect="Effect",StdErr="StdErr",log10P="log10P",
                         N="N",radius=1e6,collapse.hla=TRUE,build="hg19")
{
  for (q in c("dplyr","valr")) {
     if (length(grep(paste("^package:", q, "$", sep=""), search())) == 0) {
        if (!requireNamespace(q, quietly = TRUE))
        warning(paste("qtlFinder needs package `", q, "' to be fully functional; please install", sep=""))
     }
  }
  chrom <- start.region <- end.region <- regionStart <- regionEnd <- mlog10p <- hla_qtls <- other_qtls <- pQTLs <- NULL
  p <- data.frame(chrom=paste0("chr",d[[Chromosome]]),start=d[[Position]],end=d[[Position]],
                  rsid=d[[MarkerName]],a1=d[[Allele1]],a2=d[[Allele2]],
                  EAF=d[[EAF]],b=d[[Effect]],SE=d[[StdErr]],mlog10p=-d[[log10P]],n=d[[N]]) %>%
       valr::bed_sort()
  m <- valr::bed_merge(p,max_dist=radius)
  i <- valr::bed_intersect(p,m,suffix = c("", ".region")) %>%
       dplyr::rename(regionStart=start.region,regionEnd=end.region)
  if (nrow(i)>0)
  {
    if (!collapse.hla)
       pQTLs <- i %>%
                dplyr::group_by(chrom,regionStart,regionEnd) %>%
                dplyr::slice(which.max(mlog10p))
    else
    {
      if (build=="hg19") {startHLA <- 28477797; endHLA <- 33448354} else {startHLA <- 28510020; endHLA <- 33480577}
      hla <- with(i, chrom=="chr6" & start >= startHLA & end <= endHLA)
      hasHLA <- dplyr::filter(i,hla)
      hasOther <- dplyr::filter(i,!hla)
      if (nrow(hasHLA)>0) hla_qtls <- hasHLA %>%
                                      dplyr::group_by(chrom,regionStart,regionEnd) %>%
                                      dplyr::slice(which.max(mlog10p))
      if (nrow(hasOther)>0) other_qtls <- hasOther %>%
                                          dplyr::group_by(chrom,regionStart,regionEnd) %>%
                                          dplyr::slice(which.max(mlog10p)) %>%
                                          bind_rows(hla_qtls)
      pQTLs <- bind_rows(hla_qtls,other_qtls)
    }
  }
  pQTLs %>% group_by(chrom,regionStart,regionEnd) %>% dplyr::slice(which.max(mlog10p))
}
