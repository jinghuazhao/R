#' Cis/trans classification of pQTL signals
#'
#' Classify genetic association signals as cis or trans relative to the gene
#' encoding the target protein.
#'
#' A variant is classified as cis if it lies on the same chromosome as the gene
#' and within a specified window around the gene boundaries.
#'
#' @param hits Data frame of association signals containing at least:
#'   - identifier column (specified by `id`)
#'   - chromosome column `Chr`
#'   - position column `bp`
#'
#' @param panel Annotation data frame containing:
#'   - `uniprot` (or other id column)
#'   - `gene`, `chr`, `start`, `end`
#'
#' @param id Join key (default: `"uniprot"`).
#' @param radius Cis window size in base pairs (default: `1e6`).
#' @param verbose Logical; print summary if TRUE.
#'
#' @return A list containing:
#'
#' - **data**: annotated data frame with cis/trans classification
#'   - p.gene, p.chr, p.start, p.end
#'   - cis (logical)
#'   - cis.trans ("cis" / "trans")
#'   - cis.start, cis.end
#'
#' - **table**: gene × cis/trans contingency table
#'
#' - **total**: named vector of total cis and trans counts
#'
#' @examples
#' #----------------------------------------------------------
#' # minimal reproducible example
#' #----------------------------------------------------------
#'
#' hits <- data.frame(
#'   uniprot = c("P1", "P2", "P3", "P1"),
#'   Chr = c("1", "1", "2", "1"),
#'   bp = c(100000, 2000000, 500000, 105000)
#' )
#'
#' panel <- data.frame(
#'   uniprot = c("P1", "P2", "P3"),
#'   gene = c("G1", "G2", "G3"),
#'   chr = c("1", "1", "2"),
#'   start = c(90000, 1500000, 400000),
#'   end = c(110000, 2500000, 600000)
#' )
#'
#' ct <- cis.vs.trans.classification(hits, panel)
#'
#' ct$total
#' head(ct$data)
#' with(ct$data, table(p.gene, cis.trans))
#'
#' @export
#'
cis.vs.trans.classification <- function(hits,
                      panel,
                      id = "uniprot",
                      radius = 1e6,
                      verbose = TRUE)
{
  norm_chr <- function(x)
  {
    x <- as.character(x)
    sub("^chr", "", x, ignore.case = TRUE)
  }
  if (!id %in% names(hits))
    stop("id column not found in hits: ", id)
  if (!id %in% names(panel))
    stop("id column not found in panel: ", id)
  required_panel <- c("gene", "chr", "start", "end")
  missing_panel <- setdiff(required_panel, names(panel))
  if (length(missing_panel))
    stop("Missing panel columns: ",
         paste(missing_panel, collapse = ", "))
  panel2 <- panel[, c(
    id, "gene", "chr", "start", "end"
  ), drop = FALSE]
  names(panel2)[names(panel2) == "gene"]  <- "p.gene"
  names(panel2)[names(panel2) == "chr"]   <- "p.chr"
  names(panel2)[names(panel2) == "start"] <- "p.start"
  names(panel2)[names(panel2) == "end"]   <- "p.end"
  dat <- merge(hits, panel2, by = id, all.x = TRUE, sort = FALSE)
  dat$cis.start <- pmax(0, dat$p.start - radius)
  dat$cis.end   <- dat$p.end + radius
  dat$Chr  <- norm_chr(dat$Chr)
  dat$p.chr <- norm_chr(dat$p.chr)
  dat$cis <-
    dat$Chr == dat$p.chr &
    dat$bp >= dat$cis.start &
    dat$bp <= dat$cis.end
  dat$cis[is.na(dat$cis)] <- FALSE
  dat$cis.trans <- ifelse(dat$cis, "cis", "trans")
  tab <- with(dat, table(p.gene, cis.trans))
  total <- table(dat$cis.trans)
  if (verbose)
  {
    cat("\nCis/trans classification\n")
    cat("------------------------\n")
    cat("Variants :", nrow(dat), "\n")
    cat("Proteins :", length(unique(dat$p.gene)), "\n\n")
    print(total)
    cat("\n")
  }
  list(
    data = dat,
    table = tab,
    total = total
  )
}
