#' Distance-based signal identification
#'
#' @param d Input data frame containing summary association statistics.
#' @param Chromosome Chromosome column.
#' @param Position Genomic position column (bp).
#' @param MarkerName Variant identifier (e.g. rsID or SNP ID).
#' @param Allele1 Effect allele.
#' @param Allele2 Other allele.
#' @param EAF Effect allele frequency.
#' @param Effect Effect estimate.
#' @param StdErr Standard error of the effect estimate.
#' @param log10P Column containing log10(P) values, typically as produced by
#'   METAL. Values are converted internally to -log10(P) for ranking signals.
#' @param N Sample size.
#' @param radius Maximum gap (bp) allowed between consecutive variants when
#'   defining a locus.
#' @param collapse.hla Logical flag indicating whether variants within the HLA
#'   region should be collapsed.
#' @param build Genome build used to define HLA boundaries ("hg19" or "hg38").
#'
#' @details
#' Variants are merged into non-overlapping loci using
#' `valr::bed_merge(max_dist = radius)`. Variants connected through a chain of
#' gaps less than or equal to `radius` belong to the same locus.
#'
#' The most significant variant (largest -log10(P)) within each merged locus is
#' selected as the lead signal.
#'
#' When `collapse.hla=TRUE`, signals within the HLA region are collapsed using
#' HLA boundaries:
#'
#' * hg19/GRCh37: chr6:28477797-33448354
#' * hg38/GRCh38: chr6:28510020-33480577
#'
#' @return
#' A tibble containing one lead variant per non-overlapping genomic locus.
#'
#' Loci are defined by single-linkage clustering of variants separated by at most
#' `radius` base pairs. Each locus is represented by its most significant variant
#' (largest `-log10(P)`).
#'
#' ## Output columns
#'
#' - **chrom**: chromosome
#' - **start, end**: position of the lead variant
#' - **rsid**: variant identifier
#' - **a1, a2**: effect and other allele
#' - **EAF**: effect allele frequency
#' - **b**: effect size
#' - **SE**: standard error
#' - **mlog10p**: -log10(P) used for ranking
#' - **n**: sample size
#' - **regionStart, regionEnd**: boundaries of the merged locus
#'
#' ## Notes
#'
#' Each row corresponds to a non-overlapping locus defined using
#' `valr::bed_merge(max_dist = radius)`.
#'
#' The lead variant is selected as the SNP with the highest `mlog10p`.
#'
#' @examples
#' \dontrun{
#' threshold <- log10(5*10^-8)
#' varlist <- c(
#'   "Chromosome","Position","MarkerName","Allele1","Allele2",
#'   "Freq1","FreqSE","MinFreq","MaxFreq",
#'   "Effect","StdErr","log10P","Direction",
#'   "HetISq","HetChiSq","HetDf","logHetP","N"
#' )
#' tbl <- "~/Caprion/analysis/METAL_dr/ZPI_dr-1.tbl.gz"
#' sumstats <- read.table(tbl, col.names = varlist, check.names = FALSE)
#' res1 <- qtlFinder(sumstats)
#' res1
#' subset(res1,mlog10p>=-threshold)
#' d <- subset(sumstats, log10P <= threshold)
#' res2 <- qtlFinder(d)
#' res2
#' sentinels <- "~/Caprion/analysis/METAL_dr/sentinels/ZPI_dr.p.gz"
#' d2 <- read.table(sentinels, col.names = varlist, check.names = FALSE)
#' res3 <- qtlFinder(d2)
#' res3
#' # A tibble: 2 × 14
#' # Groups:   chrom, regionStart, regionEnd [2]
#'   chrom     start       end rsid   a1    a2      EAF      b     SE mlog10p     n
#'   <chr>     <dbl>     <dbl> <chr>  <chr> <chr> <dbl>  <dbl>  <dbl>   <dbl> <int>
#' 1 chr13 113800622 113800622 13:11… t     c     0.373 -0.541 0.0306    69.0  2491
#' 2 chr14  94750486  94750486 14:94… t     c     0.989  1.52  0.136     28.5  2491
#' # 3 more variables: regionStart <dbl>, regionEnd <dbl>, .overlap <int>
#' }
#'
#' @export
#'
qtlFinder <- function(
    d,
    Chromosome = "Chromosome",
    Position = "Position",
    MarkerName = "MarkerName",
    Allele1 = "Allele1",
    Allele2 = "Allele2",
    EAF = "Freq1",
    Effect = "Effect",
    StdErr = "StdErr",
    log10P = "log10P",
    N = "N",
    radius = 1e6,
    collapse.hla = TRUE,
    build = "hg19"
)
{
  for (q in c("dplyr", "valr"))
  {
    if (length(grep(paste("^package:", q, "$", sep = ""), search())) == 0)
    {
      if (!requireNamespace(q, quietly = TRUE))
      {
        warning(
          paste(
            "qtlFinder needs package `",
            q,
            "' to be fully functional; please install",
            sep = ""
          )
        )
      }
    }
  }

  chrom <- start.region <- end.region <- NULL
  regionStart <- regionEnd <- mlog10p <- NULL

  p <- data.frame(
         chrom = paste0("chr", d[[Chromosome]]),
         start = d[[Position]],
         end   = d[[Position]],
         rsid  = d[[MarkerName]],
         a1    = d[[Allele1]],
         a2    = d[[Allele2]],
         EAF   = d[[EAF]],
         b     = d[[Effect]],
         SE    = d[[StdErr]],
         mlog10p = -d[[log10P]],
         n     = d[[N]]
       ) %>%
       valr::bed_sort()

  m <- valr::bed_merge(
         p,
         max_dist = radius
       )

  i <- valr::bed_intersect(
         p,
         m,
         suffix = c("", ".region"),
         min_overlap = 0L
       ) %>%
       dplyr::rename(
         regionStart = start.region,
         regionEnd   = end.region
       )

  if (nrow(i) == 0)
    return(i)

  if (!collapse.hla)
  {
    pQTLs <- i %>%
      dplyr::group_by(chrom, regionStart, regionEnd) %>%
      dplyr::slice(which.max(mlog10p))
  } else {

    if (build == "hg19")
    {
      startHLA <- 28477797
      endHLA   <- 33448354
    } else {
      startHLA <- 28510020
      endHLA   <- 33480577
    }

    hla <- with(
      i,
      chrom == "chr6" &
      start >= startHLA &
      end <= endHLA
    )

    hasHLA   <- dplyr::filter(i, hla)
    hasOther <- dplyr::filter(i, !hla)

    hla_qtls <- NULL
    other_qtls <- NULL

    if (nrow(hasHLA) > 0)
    {
      hla_qtls <- hasHLA %>%
        dplyr::group_by(
          chrom,
          regionStart,
          regionEnd
        ) %>%
        dplyr::slice(which.max(mlog10p))
    }

    if (nrow(hasOther) > 0)
    {
      other_qtls <- hasOther %>%
        dplyr::group_by(
          chrom,
          regionStart,
          regionEnd
        ) %>%
        dplyr::slice(which.max(mlog10p))
    }

    pQTLs <- dplyr::bind_rows(
      hla_qtls,
      other_qtls
    )
  }

  pQTLs %>%
    dplyr::group_by(
      chrom,
      regionStart,
      regionEnd
    ) %>%
    dplyr::slice(which.max(mlog10p))
}
