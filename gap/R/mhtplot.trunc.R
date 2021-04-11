mhtplot.trunc <- function (x, chr = "CHR", bp = "BP", z = "Z", snp = "SNP",
                           col = c("gray10", "gray60"),
                           chrlabs = NULL, suggestiveline = -log10(1e-05), 
                           genomewideline = -log10(5e-08), highlight = NULL,
                           annotatelog10P = NULL, annotateTop = TRUE, cex.mtext=0.4, cex.text=0.4,
                           mtext.line = 2, cex.y = 1, y.ax.space = 5, y.brk1, y.brk2, ...)
{
  for (q in c("calibrate","plotrix")) {
     if (length(grep(paste("^package:", q, "$", sep=""), search())) == 0) {
        if (!requireNamespace(q, quietly = TRUE))
        warning(paste("mhtplot.trunc needs package `", q, "' to be fully functional; please install", sep=""))
     }
  }
  if (y.brk2 <= y.brk1) stop("y.brk2 must be larger than y.brk1")
  CHR <- BP <- Z <- SNP <- log10P <- index <- NULL
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x)))  stop(paste("Column", bp, "not found!"))
  if (!(z %in% names(x)))  stop(paste("Column", z, "not found!"))
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]]))  stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[z]]))   stop(paste(z, "column should be numeric."))
  d <- data.frame(CHR = x[[chr]], BP = x[[bp]], Z = x[[z]], log10P = -log10p(x[[z]]))
  d <- subset(d, !is.na(CHR) & !is.na(BP) & !is.na(Z))
  if (!is.null(x[[snp]])) d <- transform(d, SNP = x[[snp]])
  d <- d[order(d$CHR, d$BP), ]
  d$pos <- NA
  d$index <- NA
  ind <- 0
  for (i in unique(d$CHR)) {
    ind <- ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr <- length(unique(d$CHR))
  if (nchr == 1) {
    d$pos <- d$BP
    ticks <- floor(length(d$pos))/2 + 1
    xlabel <- paste("Chromosome", unique(d$CHR), "position")
    labs <- ticks
  } else {
    lastbase <- 0
    ticks <- NULL
    for (i in unique(d$index)) {
      if (i == 1) d[d$index == i, ]$pos = d[d$index == i, ]$BP
      else {
        lastbase = lastbase + tail(subset(d, index == i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + lastbase
      }
      ticks <- c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == i, ]$pos))/2 + 1)
    }
    xlabel <- "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax <- ceiling(max(d$pos) * 1.03)
  xmin <- floor(max(d$pos) * -0.03)
  zz <- d$log10P
  max.y <- ceiling(max(d$log10P))
  if (y.brk2 > max.y ){
    message(paste("max.y is", max.y))
    stop("User error: Upper breakpoint must be lower than maximum -log10 P-value")
  }
  zz[which(zz > y.brk1 & zz < y.brk2)] <- NA
  offset = y.brk2-y.brk1
  z[which(zz > y.brk2)] <- zz[which(zz > y.brk2)] - offset
  d$log10P <- zz
  def_args <- list(xaxt = "n", yaxt="n", bty = "n", xaxs = "i", # yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax),
                   ylim = c(0, ceiling(max(d$log10P, na.rm=TRUE))),
                   xlab = "", ylab = "")
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  mtext(text = xlabel, side = 1, line = mtext.line, cex = cex.mtext)
  mtext(text = expression(-log[10](italic(p))), side=2, line = mtext.line, cex = cex.mtext)
  y.lab.tick.pos <- seq(from = 0, by = y.ax.space, to = ceiling(max(d$log10P, na.rm = TRUE)) + y.ax.space)
  pre.brk.labs <- seq(from = 0, by = y.ax.space, to = y.brk1-y.ax.space)
  post.brk.labs <- seq(from = y.brk2, by=y.ax.space, to = max(y.lab.tick.pos))
  axis(side=2, at=c(pre.brk.labs, post.brk.labs-offset), labels=c(pre.brk.labs, post.brk.labs), cex.axis=cex.y, las=1)
# axis(side=2, at=y.lab.tick.pos, labels=c(pre.brk.labs,
#      seq(from=y.brk2, by=y.ax.space, length.out=length(y.lab.tick.pos)-length(pre.brk.labs))), cex.axis=cex.y, las=1)
  plotrix::axis.break(axis = 2, breakpos = y.brk1, style = "slash")
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
       if (length(chrlabs) == length(labs)) labs <- chrlabs
       else warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
    }
    else warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
  }
  if (nchr == 1) axis(1, ...) else axis(1, at = ticks, labels = labs, ...)
  col <- rep(col, max(d$CHR))
  if (nchr == 1) with(d, points(pos, log10P, pch = 20, col = col[1], ...))
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos, log10P, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) abline(h = suggestiveline, col = "blue")
  if (genomewideline) abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, log10P, col = "red", pch = 20, ...))
  }
  if (!is.null(annotatelog10P)) {
    topHits = subset(d, log10P >= annotatelog10P)
    if (!annotateTop) {
      with(subset(topHits,SNP %in% highlight), calibrate::textxy(pos, log10P, offset = 0.625, labs = SNP, cex = 0.45), ...)
    }
    else {
      topHits <- topHits[order(topHits$log10P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      with(topSNPs,calibrate::textxy(pos, log10P, offset = 0.625, labs = SNP, cex = cex.text),...)
    }
  }
}
