#' Truncated Manhattan plot
#'
#' @param x A data.frame.
#' @param chr Chromosome column name.
#' @param bp Base-pair position column name.
#' @param p P-value column name.
#' @param log10p Column containing -log10(P).
#' @param z Z-statistic column name.
#' @param snp SNP/variant identifier column name.
#' @param col Point colours alternating by chromosome.
#' @param chrlabs Optional chromosome labels.
#' @param suggestiveline Horizontal suggestive significance line.
#' @param genomewideline Horizontal genome-wide significance line.
#' @param highlight SNPs to highlight.
#' @param annotatelog10P Threshold for annotation.
#' @param annotateTop Annotate only top SNP per chromosome.
#' @param cex.mtext Axis title size.
#' @param cex.text SNP label size.
#' @param mtext.line Axis title line offset.
#' @param y.ax.space Y-axis tick spacing.
#' @param y.brk1 Lower truncation breakpoint.
#' @param y.brk2 Upper truncation breakpoint.
#' @param trunc.yaxis Enable truncated y-axis.
#' @param cex.axis Axis tick label size.
#' @param delta Fractional window around highlighted SNPs.
#' @param ... Additional graphical parameters passed to points().
#'
#' @details
#' Draws a Manhattan plot with optional y-axis truncation for
#' extremely significant associations commonly observed in
#' large-scale GWAS or protein GWAS analyses.
#'
#' @return
#' Invisibly returns the processed plotting data.
#'
#' Example FTO locus dataset and truncated Manhattan plot
#'
#' This example demonstrates the use of `mhtplot.trunc()` on a
#' single-chromosome FTO locus dataset with and without y-axis truncation.
#'
#' @examples
#' txt <- "
#' CHR POS SNP Z log10P
#' 16 53804965 rs10852521 -39.75000 344.8039
#' 16 53805207 rs11075985 43.88235 419.8925
#' 16 53819877 rs11075989 43.94118 421.0149
#' 16 53819893 rs11075990 -43.94118 421.0149
#' 16 53809247 rs1121980 43.76471 417.6523
#' 16 53845487 rs11642841 38.52941 324.0426
#' 16 53842908 rs12149832 42.23529 389.0756
#' 16 53800954 rs1421085 -45.76471 456.5538
#' 16 53803574 rs1558902 45.88235 458.8962
#' 16 53813367 rs17817449 -46.87500 478.8994
#' 16 53828066 rs17817964 44.29412 427.7808
#' 16 53804340 rs1861866 39.81250 345.8844
#' 16 53818460 rs3751812 44.41176 430.0480
#' 16 53822651 rs7185735 -43.94118 421.0149
#' 16 53810686 rs7193144 -44.29412 427.7808
#' 16 53821862 rs7201850 42.35294 391.2377
#' 16 53821615 rs7202116 -43.94118 421.0149
#' 16 53813450 rs8043757 -42.22222 388.8357
#' 16 53798523 rs8047395 37.76471 311.3650
#' 16 53816275 rs8050136 46.75000 476.3569
#' 16 53816752 rs8051591 -44.00000 422.1388
#' 16 53803156 rs8055197 39.81250 345.8844
#' 16 53812614 rs8057044 39.75000 344.8039
#' 16 53806280 rs9922047 -39.62500 342.6480
#' 16 53831771 rs9922619 43.37500 410.2743
#' 16 53831146 rs9922708 43.25000 407.9218
#' 16 53801549 rs9923147 43.82353 418.7716
#' 16 53819198 rs9923233 43.94118 421.0149
#' 16 53801985 rs9923544 43.82353 418.7716
#' 16 53820503 rs9926289 40.66667 360.8208
#' 16 53799905 rs9928094 -43.88235 419.8925
#' 16 53799977 rs9930333 -44.00000 422.1388
#' 16 53830452 rs9930501 -40.94118 365.6883
#' 16 53830465 rs9930506 -43.31250 409.0972
#' 16 53827179 rs9931494 -42.94118 402.1387
#' 16 53830491 rs9932754 -41.00000 366.7356
#' 16 53816838 rs9935401 46.81250 477.6273
#' 16 53819169 rs9936385 -44.23529 426.6494
#' 16 53799507 rs9937053 43.94118 421.0149
#' 16 53820527 rs9939609 44.05882 423.2642
#' 16 53800568 rs9939973 43.88235 419.8925
#' 16 53800754 rs9940128 43.82353 418.7716
#' 16 53800629 rs9940646 -43.82353 418.7716
#' 16 53825488 rs9941349 42.58824 395.5801
#' "
#' FTO <- read.table(text = txt, header = TRUE, stringsAsFactors = FALSE)
#' par(mar = c(5, 6, 2, 1))
#' mhtplot.trunc(
#'   x = FTO,
#'   chr = "CHR",
#'   bp = "POS",
#'   log10p = "log10P",
#'   snp = "SNP",
#'   cex = 1.0,
#'   cex.axis = 1.2,
#'   cex.mtext = 1.5,
#'   col = "navy"
#' )
#' title("FTO locus without truncation")
#'
#' par(mar = c(5, 6, 2, 1))
#' mhtplot.trunc(
#'   x = FTO,
#'   chr = "CHR",
#'   bp = "POS",
#'   log10p = "log10P",
#'   snp = "SNP",
#'   trunc.yaxis = TRUE,
#'   y.brk1 = 200,
#'   y.brk2 = 350,
#'   y.ax.space = 50,
#'   genomewideline = -log10(5e-8),
#'   suggestiveline = NULL,
#'   highlight = c("rs1421085", "rs1558902", "rs17817449",
#'                 "rs8050136", "rs9939609"),
#'   annotatelog10P = 420,
#'   annotateTop = FALSE,
#'   cex = 1.2,
#'   cex.text = 0.9,
#'   cex.axis = 1.2,
#'   cex.mtext = 1.5,
#'   col = "steelblue4"
#' )
#' title("FTO locus with truncated y-axis")
#' @export
#'
mhtplot.trunc <- function(
    x,
    chr = "CHR",
    bp = "BP",
    p = NULL,
    log10p = NULL,
    z = NULL,
    snp = "SNP",
    col = c("gray10", "gray60"),
    chrlabs = NULL,
    suggestiveline = -log10(1e-5),
    genomewideline = -log10(5e-8),
    highlight = NULL,
    annotatelog10P = NULL,
    annotateTop = FALSE,
    cex.mtext = 1.5,
    cex.text = 0.7,
    mtext.line = 2,
    y.ax.space = 5,
    y.brk1 = NULL,
    y.brk2 = NULL,
    trunc.yaxis = TRUE,
    cex.axis = 1.2,
    delta = 0.05,
    ...
) {
  CHR <- BP <- SNP <- BP.x <- BP.y <- log10P <- index <- pos <- NULL
  for (pkg in c("calibrate", "plotrix")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(
        "Package `", pkg,
        "` is recommended for full functionality."
      )
    }
  }
  required_cols <- c(chr, bp)
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols)) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (!is.numeric(x[[chr]])) {
    stop(
      chr,
      " must be numeric. Convert X/Y/MT chromosomes to integers first."
    )
  }
  if (!is.null(p)) {
    if (!(p %in% names(x))) {
      stop("Column ", p, " not found.")
    }
    log10P <- -log10pvalue(x[[p]])
  } else if (!is.null(log10p)) {
    if (!(log10p %in% names(x))) {
      stop("Column ", log10p, " not found.")
    }
    log10P <- as.numeric(x[[log10p]])
  } else if (!is.null(z)) {
    if (!(z %in% names(x))) {
      stop("Column ", z, " not found.")
    }
    log10P <- -log10p(as.numeric(x[[z]]))
  } else {
    stop("One of `p`, `log10p`, or `z` must be supplied.")
  }
  use_truncation <- (
    trunc.yaxis &&
    !is.null(y.brk1) &&
    !is.null(y.brk2)
  )
  if (use_truncation && y.brk2 <= y.brk1) {
    stop("y.brk2 must be larger than y.brk1")
  }
  d <- data.frame(
    CHR = x[[chr]],
    BP = as.numeric(x[[bp]]),
    log10P = log10P,
    stringsAsFactors = FALSE
  )
  if (snp %in% names(x)) {
    d$SNP <- x[[snp]]
  }
  d <- subset(
    d,
    !is.na(CHR) &
    !is.na(BP) &
    !is.na(log10P) &
    is.finite(log10P)
  )
  d <- d[order(d$CHR, d$BP), ]
  d$index <- match(d$CHR, unique(d$CHR))
  nchr <- length(unique(d$CHR))
  if (nchr == 1) {
    d$pos <- d$BP
    ticks <- median(d$pos)
    xlabel <- paste(
      "Chromosome",
      unique(d$CHR),
      "position"
    )
    labs <- unique(d$CHR)
  } else {
    chr_lengths <- tapply(d$BP, d$index, max)
    offsets <- c(
      0,
      cumsum(chr_lengths)[-length(chr_lengths)]
    )
    d$pos <- d$BP + offsets[d$index]
    ticks <- tapply(
      d$pos,
      d$index,
      function(z) (min(z) + max(z)) / 2
    )
    xlabel <- "Chromosome"
    labs <- unique(d$CHR)
  }
  if (use_truncation) {
    max.y.original <- ceiling(max(d$log10P, na.rm = TRUE))
    if (y.brk2 > max.y.original) {
      stop(
        "Upper breakpoint exceeds maximum observed -log10(P)"
      )
    }
    offset <- y.brk2 - y.brk1
    gapped <- d$log10P > y.brk1 &
              d$log10P < y.brk2
    above <- d$log10P >= y.brk2
    d$log10P[gapped] <- NA
    d$log10P[above]  <- d$log10P[above] - offset
  } else {
    offset <- 0
  }
  xmax <- ceiling(max(d$pos) * 1.03)
  xmin <- floor(min(d$pos) * 0.97)
  ymax <- ceiling(max(d$log10P, na.rm = TRUE))
  def_args <- list(
    xaxt = "n",
    yaxt = "n",
    bty = "n",
    xaxs = "i",
    las = 1,
    pch = 20,
    xlim = c(xmin, xmax),
    ylim = c(0, ymax),
    xlab = "",
    ylab = ""
  )
  dotargs <- list(...)
  do.call(
    "plot",
    c(
      list(NA),
      dotargs,
      def_args[!names(def_args) %in% names(dotargs)]
    )
  )
  mtext(
    text = xlabel,
    side = 1,
    line = mtext.line,
    cex = cex.mtext,
    font = 2
  )
  mtext(
    text = expression(-log[10](italic(P))),
    side = 2,
    line = mtext.line,
    cex = cex.mtext,
    font = 2
  )
  if (use_truncation) {
    y.tick.pos <- seq(
      0,
      ymax,
      by = y.ax.space
    )
    lower.labels <- seq(
      0,
      y.brk1,
      by = y.ax.space
    )
    upper.n <- length(y.tick.pos) - length(lower.labels)
    upper.labels <- seq(
      y.brk2,
      by = y.ax.space,
      length.out = upper.n
    )
    y.labels <- ceiling(c(lower.labels, upper.labels))
    axis(
      side = 2,
      at = y.tick.pos,
      labels = y.labels,
      las = 1,
      cex.axis = cex.axis
    )
    plotrix::axis.break(
      axis = 2,
      breakpos = y.brk1,
      style = "slash"
    )
  } else {
    axis(
      side = 2,
      las = 1,
      cex.axis = cex.axis
    )
  }
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      } else {
        warning(
          "Length of `chrlabs` does not match chromosome count."
        )
      }
    } else {
      warning("`chrlabs` must be a character vector.")
    }
  }
  if (nchr == 1) {
    axis(1, cex.axis = cex.axis)
  } else {
    axis(
      1,
      at = ticks,
      labels = labs,
      cex.axis = cex.axis
    )
  }
  chr_cols <- rep(
    col,
    length.out = length(unique(d$CHR))
  )
  for (i in seq_along(unique(d$index))) {
    idx <- d$index == i
    points(
      d$pos[idx],
      d$log10P[idx],
      col = grDevices::adjustcolor(
        chr_cols[i],
        alpha.f = 0.7
      ),
      pch = 20,
      ...
    )
  }
  if (!is.null(suggestiveline)) {
    abline(h = suggestiveline, col = "blue")
  }
  if (!is.null(genomewideline)) {
    abline(h = genomewideline, col = "red")
  }
  if (!is.null(highlight) && "SNP" %in% names(d)) {
    missing_snps <- setdiff(highlight, d$SNP)
    if (length(missing_snps)) {
      warning(
        "Some highlighted SNPs not found: ",
        paste(missing_snps, collapse = ", ")
      )
    }
    d.highlight <- subset(d, SNP %in% highlight)
    points(
      d.highlight$pos,
      d.highlight$log10P,
      col = "red",
      pch = 20,
      ...
    )
    d.column <- subset(
      merge(
        d,
        d.highlight[c("CHR", "BP")],
        by = "CHR"
      ),
      BP.x > (1 - delta) * BP.y &
      BP.x < (1 + delta) * BP.y
    )
    if (nrow(d.column) > 0) {
      points(
        d.column$pos,
        d.column$log10P,
        col = "red",
        pch = 20,
        ...
      )
    }
  }
  if (
    !is.null(annotatelog10P) &&
    "SNP" %in% names(d)
  ) {
    topHits <- subset(
      d,
      log10P >= annotatelog10P
    )
    if (!annotateTop) {
      if (!is.null(highlight)) {
        topHits <- subset(
          topHits,
          SNP %in% highlight
        )
      }
    } else {
      topHits <- do.call(
        rbind,
        lapply(
          split(topHits, topHits$CHR),
          function(z) z[which.max(z$log10P), ]
        )
      )
    }
    if (nrow(topHits) > 0) {
      calibrate::textxy(
        topHits$pos,
        topHits$log10P,
        offset = 0.625,
        pos = 3,
        labs = topHits$SNP,
        cex = cex.text,
        font = 4
      )
    }
  }
  invisible(d)
}
