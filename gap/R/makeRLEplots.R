#' Create a Relative Log Expression (RLE) Plot
#'
#' Generates a Relative Log Expression (RLE) plot for quality assessment of
#' gene expression data. For each feature (row), the median expression across
#' samples is subtracted, and boxplots of the resulting relative expression
#' values are displayed for each sample.
#'
#' @param E A numeric matrix or data frame with features in rows and samples
#'   in columns.
#' @param log2.data Logical; if TRUE, expression values are log2-transformed
#'   before computing RLE values.
#' @param groups Optional factor or vector defining sample groups. Must have
#'   length equal to ncol(E).
#' @param col.group Optional named vector of colours matching group levels.
#' @param showTitle Logical; if TRUE, display a plot title.
#' @param title Character string specifying the plot title.
#' @param ... Additional arguments passed to graphics::boxplot().
#'
#' @return Invisibly returns the matrix of relative log expression values.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' E <- matrix(rlnorm(987*200, meanlog=5, sdlog=1), nrow=987, ncol=200)
#' colnames(E) <- paste0("Sample_", seq_len(200))
#' rownames(E) <- paste0("Gene_", seq_len(987))
#' group <- rep(c("Control", "Treatment"), each=100)
#' col.group <- c(Control="steelblue", Treatment="tomato")
#' par(mar = c(10,4,4,2))
#' makeRLEplot(E,log2.data=TRUE,groups=group,col.group=col.group,cex=0.4,
#'   cex.axis=0.5,showTitle=TRUE,title="RLE plot (987 genes×200 samples)")
#' }
#'
#' @export
#'
#'
makeRLEplot <- function(
  E,
  log2.data = TRUE,
  groups = NULL,
  col.group = NULL,
  showTitle = FALSE,
  title = "Relative log expression (RLE) plot",
  ...
) {
  E <- as.matrix(E)
  storage.mode(E) <- "numeric"
  if (nrow(E) == 0L || ncol(E) < 2L) {
    stop("'E' must have at least one row and two columns")
  }
  if (!is.null(groups)) {
    if (length(groups) != ncol(E)) {
      stop("'groups' must have length equal to ncol(E)")
    }
    groups <- factor(groups)
  }
  if (isTRUE(log2.data)) {
    if (any(E <= 0, na.rm = TRUE)) {
      warning("Non-positive values detected; log2 may produce -Inf/NaN")
    }
    message("Applying log2 transformation")
    E <- log2(E)
  }
  mycol <- NULL
  if (!is.null(groups) && !is.null(col.group)) {
    mycol <- col.group[as.character(groups)]
    if (any(is.na(mycol))) {
      warning("Some groups have no matching colours in col.group")
    }
  }
  g.medians <- matrixStats::rowMedians(E, na.rm = TRUE)
  E.rle <- sweep(E, 1, g.medians, FUN = "-")
  dots <- list(...)
  n <- ncol(E.rle)
  cex.axis.use <- if (!is.null(dots$cex.axis)) {
    dots$cex.axis
  } else {
    if (n > 100) 0.4 else if (n > 50) 0.6 else 0.8
  }
  dots$cex.axis <- NULL
  sample.labels <- colnames(E.rle)
  if (is.null(sample.labels)) {
    sample.labels <- seq_len(ncol(E.rle))
  }
  args <- c(
    list(
      E.rle,
      cex.axis = cex.axis.use,
      las = 2,
      col = mycol,
      names = sample.labels
    ),
    dots
  )
  do.call(boxplot, args)
  if (isTRUE(showTitle)) {
    mtext(title, side = 3, line = 0.1)
  }
  invisible(E.rle)
}
