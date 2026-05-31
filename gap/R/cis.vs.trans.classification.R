#' Cis/trans classification of pQTL signals
#'
#' Classify genetic association signals as cis or trans relative to the gene
#' encoding the target protein.
#'
#' A variant is classified as cis if it lies on the same chromosome as the
#' gene and within a specified window around the gene start and end positions.
#'
#' @param hits Data frame of association signals. Must contain:
#'   - identifier column (default: `uniprot`)
#'   - chromosome column (`Chr`)
#'   - position column (`bp`)
#' @param panel Annotation data frame with columns:
#'   - `uniprot`, `gene`, `chr`, `start`, `end`
#' @param id Join key (default: `"uniprot"`).
#' @param chr Chromosome column in `hits` (default: `"Chr"`).
#' @param pos Position column in `hits` (default: `"bp"`).
#' @param radius Cis window size in base pairs (default: `1e6`).
#' @param verbose Logical; print summary if TRUE.
#'
#' @return A list containing:
#' - **data**: annotated data frame with cis/trans classification
#' - **table**: gene × cis/trans contingency table
#' - **overall**: total counts of cis and trans variants
#'
#' @examples
#' \dontrun{
#' ct <- cis.vs.trans.classification(hits, inf1)
#' ct$overall
#' head(ct$data)
#' with(ct$data, table(p.gene, cis.trans))
#' subset(ct$data, cis.trans == "cis")[1:10,
#'        c("uniprot", "p.gene", "Chr", "bp", "cis.trans")]
#' }
#'
#' @export
#'
cis.vs.trans.classification <- function(hits,
                      panel,
                      id = "uniprot",
                      chr = "Chr",
                      pos = "bp",
                      radius = 1e6,
                      verbose = TRUE)
{
  norm_chr <- function(x)
  {
    x <- as.character(x)
    sub("^chr", "", x, ignore.case = TRUE)
  }
  if (!id %in% names(hits))
    stop("Missing column in hits: ", id)
  if (!id %in% names(panel))
    stop("Missing column in panel: ", id)
  if (!chr %in% names(hits))
    stop("Missing column in hits: ", chr)
  if (!pos %in% names(hits))
    stop("Missing column in hits: ", pos)
  required_panel <- c("gene", "chr", "start", "end")
  missing_panel <- setdiff(required_panel, names(panel))
  if (length(missing_panel))
    stop("Missing panel columns: ",
         paste(missing_panel, collapse = ", "))
  panel2 <- panel[, c(
    id,
    "gene", "chr", "start", "end"
  ), drop = FALSE]
  names(panel2)[names(panel2) == "gene"]  <- "p.gene"
  names(panel2)[names(panel2) == "chr"]   <- "p.chr"
  names(panel2)[names(panel2) == "start"] <- "p.start"
  names(panel2)[names(panel2) == "end"]   <- "p.end"
  dat <- merge(hits, panel2, by = id, all.x = TRUE)
  dat$cis.start <- pmax(0, dat$p.start - radius)
  dat$cis.end   <- dat$p.end + radius
  dat[[chr]] <- norm_chr(dat[[chr]])
  dat$p.chr  <- norm_chr(dat$p.chr)
  dat$cis <-
    dat[[chr]] == dat$p.chr &
    dat[[pos]] >= dat$cis.start &
    dat[[pos]] <= dat$cis.end
  dat$cis[is.na(dat$cis)] <- FALSE
  dat$cis.trans <- ifelse(dat$cis, "cis", "trans")
  summary_table <- with(dat, table(p.gene, cis.trans))
  overall <- table(dat$cis.trans)
  if (verbose)
  {
    cat("\nCis/trans classification\n")
    cat("------------------------\n")
    cat("Variants :", nrow(dat), "\n")
    cat("Proteins :", length(unique(dat$p.gene)), "\n\n")
    print(overall)
    cat("\n")
  }
  list(
    data = dat,
    table = summary_table,
    overall = overall
  )
}
