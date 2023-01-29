#' Miami plot
#'
#' @param x Input data.
#' @param chr Chromsome.
#' @param bp Position.
#' @param p P value.
#' @param pr P value of the other GWAS.
#' @param snp Marker.
#' @param col Colors.
#' @param col2 Colors.
#' @param ymax Max y.
#' @param highlight Highlight flag.
#' @param highlight.add Highlight meta-data.
#' @param pch Symbol.
#' @param cex cex.
#' @param cex.lab cex for labels.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param lcols Colors.
#' @param lwds lwd.
#' @param ltys lty.
#' @param main Main title.
#' @param ... Additional options.
#'
#' @details
#' The function allows for contrast of genomewide P values from two GWASs. It is conceptually simpler than at the first sight since it involves only one set of chromosomal positions.
#'
#' @export
#' @return None.
#' @examples
#' \dontrun{
#'   mhtdata <- within(mhtdata,{pr=p})
#'   miamiplot(mhtdata,chr="chr",bp="pos",p="p",pr="pr",snp="rsn")
#' }

miamiplot <- function (x, chr = "CHR", bp = "BP", p = "P", pr = "PR", snp = "SNP", 
    col = c("midnightblue", "chartreuse4"), col2 = c("royalblue1", 
        "seagreen1"), ymax = NULL, highlight = NULL, highlight.add = NULL, 
    pch = 19, cex = 0.75, cex.lab = 1, xlab = "Chromosome", ylab = "-log10(P) [y>0]; log10(P) [y<0]", 
    lcols = c("red", "black"), lwds = c(5, 2), ltys = c(1, 2), 
    main = "", ...) 
{
    P = index = NULL
    PR = index = NULL
    if (!(chr %in% names(x))) 
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) 
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) 
        stop(paste("Column", p, "not found!"))
    if (!(pr %in% names(x))) 
        stop(paste("Column", pr, "not found!"))
    if (!(snp %in% names(x))) 
        if (!is.numeric(x[[chr]])) 
            stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) 
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) 
        stop(paste(p, "column should be numeric."))
    if (!is.numeric(x[[pr]])) 
        stop(paste(pr, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
        PR = x[[pr]])
    if (!is.null(x[[snp]])) 
        d = transform(d, SNP = x[[snp]])
    d = subset(d[order(d$CHR, d$BP), ], (P > 0 & P <= 1 & is.numeric(P) & 
        PR > 0 & PR <= 1 & is.numeric(PR)))
    d$logp = -log10(d$P)
    d$logpr = log10(d$PR)
    d$pos = NA
    ymax = ceiling(max(d$logp) + 2)
    ymin = floor(min(d$logpr) - 2)
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        d$pos = d$BP
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position")
        labs = ticks
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index == 
                  i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
                  lastbase
            }
            ticks = c(ticks, d[d$index == i, ]$pos[floor(length(d[d$index == 
                i, ]$pos)/2) + 1])
        }
        xlabel = "Chromosome"
        labs = unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
        xlim = c(xmin, xmax), ylim = c(ymin, ymax), main = main, 
        xlab = xlab, ylab = ylab, las = 1, pch = pch, cex.lab = cex.lab, 
        ...)
    if (nchr == 1) {
        axis(1, ...)
    }
    else {
        axis(1, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    col2 = rep(col2, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logpr, cex = cex, pch = pch, ...))
        with(d, points(pos, logp, cex = cex, pch = pch, ...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos, 
                logpr, col = col2[icol], cex = cex, pch = pch, 
                ...))
            with(d[d$index == unique(d$index)[i], ], points(pos, 
                logp, col = col[icol], cex = cex, pch = pch, 
                ...))
            icol = icol + 1
        }
    }
    abline(h = -log10(5e-08), col = lcols[2], lwd = lwds[2], 
        lty = ltys[2])
    abline(h = log10(5e-08), col = lcols[2], lwd = lwds[2], lty = ltys[2])
    abline(h = 0, col = lcols[1], lwd = lwds[1], lty = ltys[1])
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight = d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col = "red3", cex = cex, 
            pch = pch, ...))
        with(d.highlight, points(pos, logpr, col = "red3", cex = cex, 
            pch = pch, ...))
    }
    if (!is.null(highlight.add)) {
        print("yessssssssssssssssss")
        if (any(!(highlight.add %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight.add = d[which(d$SNP %in% highlight.add), 
            ]
        with(d.highlight.add, points(pos, logp, col = "darkgreen", 
            cex = cex, pch = pch, ...))
        with(d.highlight.add, points(pos, logpr, col = "darkgreen", 
            cex = cex, pch = pch, ...))
    }
}
