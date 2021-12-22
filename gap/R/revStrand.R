#' Allele on the reverse strand
#'
#' The function obtains allele on the reverse strand.
#'
#' @md
#' @param allele Allele to reverse.
#' @export
#' @return Allele on the reverse strand.
#' @examples
#' \dontrun{
#'   alleles <- c("a","c","G","t")
#'   reverse_strand(alleles)
#' }

revStrand <- function(allele)
{
  m <- sapply(allele,function(x) {
  if (x %in% LETTERS)
  { 
    Allele <- toupper(allele)
    forward <- c("A","C","G","T")
    reverse <- c("T","G","C","A")
  } else {
    Allele <- tolower(allele)
    forward <- c("a","c","g","t")
    reverse <- c("t","g","c","a")
  }
  i <- (1:4)[Allele==forward]
  reverse[i]
  })
  diag(m)
}
