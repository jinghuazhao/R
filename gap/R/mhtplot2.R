#' Manhattan plot with annotations
#'
#' @param data a data frame with three columns representing chromosome, position and p values.
#' @param control A control function named mht.control() with the following arguments:
#' - type a flag with value "p" or "l" indicating if points or lines are to be drawn.
#' - usepos a flag to use real chromosomal positions as composed to ordinal positions with default value FALSE.
#' - logscale a flag to indicate if p value is to be log-transformed with default value TRUE.
#' - base the base of the logarithm with default value 10.
#' - cutoffs the cut-offs where horizontal line(s) are drawn with default value NULL.
#' - colors the color for different chromosome(s), and random if unspecified with default values NULL.
#' - labels labels for the ticks on x-axis with default value NULL.
#' - srt degree to which labels are rotated with default value of 45.
#' - gap gap between chromosomes with default value NULL.
#' - cex cex for the data points.
#' - yline Margin line position.
#' - xline Margin line position.
#' @param hcontrol A control function named hmht.control() with the following arguments:
#' - data chunk of data to be highlighted with default value NULL.
#' - colors colors for annotated genes.
#' - yoffset offset above the data point showing most significant p value with default value 0.5.
#' - cex shrinkage factor for data points with default value 1.5.
#' - boxed if the label for the highlited region with default value FALSE.
#' @param ... other options in compatible with the R plot function.
#'
#' @details
#' To generate Manhattan plot with annotations. The function is generic and for instance could be used for genomewide
#' p values or any random variable that is uniformly distributed. By default, a log10-transformation is applied.
#' Note that with real chromosomal positions, it is also appropriate to plot and some but not all chromosomes.
#'
#' It is possible to specify options such as xlab, ylim and font family when the plot is requested for data in other
#' context.
#'
#' To maintain back compatibility options as in [`mhtplot`] are used. The positions of the horizontal
#' labels are now in the middle rather than at the beginning of their bands in the plot.
#'
#' @export
#' @return
#' The plot is shown on or saved to the appropriate device.
#' @examples
#' \dontrun{
#' The following example uses only chromosomes 14 and 20 of the Nat Genet paper.
#'
#' mdata <- within(hr1420,{
#'   c1<-colour==1
#'   c2<-colour==2
#'   c3<-colour==3
#'   colour[c1] <- 62
#'   colour[c2] <- 73
#'   colour[c3] <- 552
#' })
#' mdata <- mdata[,c("CHR","POS","P","gene","colour")]
#' ops <- mht.control(colors=rep(c("lightgray","gray"),11),yline=1.5,xline=2,srt=0)
#' hops <- hmht.control(data=subset(mdata,!is.na(gene)))
#' v <- "Verdana"
#' ifelse(Sys.info()['sysname']=="Windows", windowsFonts(ffamily=windowsFont(v)),
#'        ffamily <- v)
#' tiff("mh.tiff", width=.03937*189, height=.03937*189/2, units="in", res=1200,
#'      compress="lzw")
#' par(las=2, xpd=TRUE, cex.axis=1.8, cex=0.4)
#' mhtplot2(with(mdata,cbind(CHR,POS,P,colour)),ops,hops,pch=19,
#'          ylab=expression(paste(plain("-"),log[10],plain("p-value"),sep=" ")),
#'          family="ffamily")
#' axis(2,pos=2,at=seq(0,25,5),family="ffamily",cex=0.5,cex.axis=1.1)
#' dev.off()
#'
#' # To exemplify the use of chr, pos and p without gene annotation
#' # in response to query from Vallejo, Roger <Roger.Vallejo@ARS.USDA.GOV>
#' opar <- par()
#' par(cex=0.4)
#' ops <- mht.control(colors=rep(c("lightgray","lightblue"),11),srt=0,yline=2.5,xline=2)
#' mhtplot2(data.frame(mhtdata[,c("chr","pos","p")],gene=NA,color=NA),ops,xlab="",ylab="",srt=0)
#' axis(2,at=1:16)
#' title("data in mhtplot used by mhtplot2")
#' par(opar)
#' }
#' @references
#' \insertRef{denhoed13}{gap}
#' @author Jing Hua Zhao
#' @keywords hplot

mhtplot2 <- function (data, control = mht.control(), hcontrol = hmht.control(), ...)
{
    for(p in c("grid")) {
       if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
          if (!requireNamespace(p, quietly = TRUE))
          warning(paste("mhtplot2 needs package `", p, "' to be fully functional; please install", sep=""))
       }
    }
    nonmiss <- !apply(is.na(data[,1:3]), 1, any)
    data2 <- data[nonmiss, ]
    chr <- data2[, 1]
    pos <- newpos <- data2[, 2]
    p <- data2[, 3]
    rc <- data2[, 4]
    tablechr <- table(chr)
    allchr <- as.vector(tablechr)
    n.chr <- length(allchr)
    type <- control$type
    usepos <- control$usepos
    logscale <- control$logscale
    base <- control$base
    cutoffs <- control$cutoffs
    colors <- control$colors
    labels <- control$labels
    srt <- control$srt
    gap <- control$gap
    pcex <- control$cex
    yline <- control$yline
    xline <- control$xline
    colorlist <- colors()
    if (is.null(colors)) 
        colors <- sample(colorlist, n.chr)
    tablechr <- unique(chr)
    if (is.null(labels)) 
        labels <- tablechr
    if (is.null(gap)) 
        gap <- 0
    if (!is.null(hcontrol$data)) {
        hdata <- hcontrol$data
        hdata2 <- hdata[!apply(is.na(hdata), 1, any), ]
        hchr <- hdata2[, 1]
        hpos <- hnewpos <- hdata2[, 2]
        hp <- hdata2[, 3]
        hname <- hdata2[, 4]
        hcol <- hdata2[, 5]
        namecol <- table(hname,hcol)
        namecol[namecol!=0] <- 1
        hcolors <- colnames(namecol)
        gmat <- namecol%*%(1:length(hcolors))
        hyoffs <- hcontrol$yoffs
        hboxed <- hcontrol$boxed
        hcex <- hcontrol$cex
        htablechr <- unique(hchr)
        hn.chr <- length(htablechr)
        hlabels <- unique(hname)
        htablechrname <- unique(data.frame(hchr, hname))
        if (is.null(hcolors)) 
            hcolors <- rep("red", length(hlabels))
        else hcolors <- hcontrol$colors
    }
    CMindex <- cumsum(allchr)
    for (i in 1:n.chr) {
        u <- CMindex[i]
        l <- CMindex[i] - allchr[i] + 1
        chr <- l:u
        if (usepos) 
            d <- diff(pos[chr])
        else d <- rep(1, allchr[i] - 1)
        newpos[chr] <- c(gap, d)
    }
    CM <- cumsum(as.numeric(newpos))
    args <- list(...)
    if ("ylim" %in% names(args)) 
        dp <- seq(args$ylim[1], args$ylim[2], length = sum(allchr))
    else dp <- seq(min(p), max(p), length = sum(allchr))
    if (logscale) 
        y <- -log(dp, base)
    else y <- dp
    y1 <- min(y)
    par(xaxt = "n", yaxt = "n")
    xy <- xy.coords(CM, y)
    plot(xy$x, xy$y, type = "n", ann = FALSE, axes = FALSE, ...)
    axis(1)
    axis(2)
    par(xaxt = "s", yaxt = "s")
    for (i in 1:n.chr) {
        u <- CMindex[i]
        l <- CMindex[i] - allchr[i] + 1
        chr <- l:u
        cat("Plotting points ", l, "-", u, "\n")
        if (logscale)
            y <- -log(p[chr], base)
        else y <- p[chr]
        col.chr <- colors()[rc[chr]]
        col.chr[is.na(rc[chr])]<-colors[i]
        if (type == "l")
            lines(CM[chr], y, col = col.chr, cex = pcex, ...)
        else points(CM[chr], y, col = col.chr, cex = pcex, ...)
        text(ifelse(i == 1, CM[1]+(CM[u]-CM[l])/2, CM[l]+(CM[u]-CM[l])/2), y1, pos = 1, offset = 1.5, 
             labels[i], srt = srt, ...)
    }
    j <- 1
    for (i in 1:n.chr) {
        u <- CMindex[i]
        l <- CMindex[i] - allchr[i] + 1
        chr <- l:u
        if (logscale) 
            y <- -log(p[chr], base)
        else y <- p[chr]
        col.chr <- colors()[rc[chr]]
        if (type =="l"&all(is.na(rc[chr])))
          lines(CM[chr], y, col = col.chr, cex = pcex, ...)
        else points(CM[chr], y, col = col.chr, cex = pcex, ...)
        if (!is.null(hcontrol$data)) {
            chrs <- htablechrname[tablechr[i] == htablechrname[, 
                1], ]
            if (dim(chrs)[1] > 0) {
                hchrs <- as.character(chrs[, 2])
                for (k in 1:length(hchrs)) {
                  hregion <- hpos[hchr == chrs[k, 1] & hname == 
                    hchrs[k]]
                  hl <- chr[pos[chr] == hregion[1]]
                  hu <- chr[pos[chr] == hregion[length(hregion)]]
                  cat("  ... highlighting", hl, "-", hu, hchrs[k], 
                    "\n")
                  l1 <- hl - l + 1
                  l2 <- hu - l + 1
                  col.index <- as.integer(colnames(namecol)[gmat[rownames(gmat)==hchrs[k]]])
                  col.label <- colors()[col.index]
                  if (hboxed) textbox(hchrs[k], name="tbt", vp=grid::viewport(x = CM[chr][l1]/max(CM), y = (max(y[l1:l2]) +  hyoffs)/max(y)))
                  else text(CM[chr][l1], max(y[l1:l2] + hyoffs), hchrs[k],
                       col=col.label, cex = hcex, font=3, ...)
 #                 points(CM[l + (l1:l2)], y[l1:l2], col = col.label, cex = pcex, ...)
                  j <- j + 1
                }
            }
        }
    }
    if (!is.null(cutoffs)) 
        abline(h = cutoffs)
    if ("ylab" %in% names(args)) 
        mtext(args$ylab, 2, line = yline, las = 0, cex=0.5, ...)
    else mtext(ifelse(logscale, paste("-log", base, "(Observed value)", 
        sep = ""), "Observed value"), 2, line = yline, las = 0, cex=0.5, ...)
    if ("xlab" %in% names(args)) 
        xlabel <- args$xlab
    else xlabel <- ifelse(is.null(names(chr)), "Chromosome", 
        names(chr))
    mtext(xlabel, 1, line = xline, las = 0, cex=0.5, ...)
}
