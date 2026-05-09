#' Manhattan plot
#'
#' Draw a Manhattan plot for genome-wide association studies (GWAS) or any
#' genome-indexed numeric score.
#'
#' @param data A matrix or data frame with at least three columns.
#'   The first three columns are interpreted as:
#'   1. chromosome
#'   2. position
#'   3. value (typically p-value)
#'   Column names are ignored. Factors are automatically converted.
#'
#' @param control A list produced by [mht.control()] controlling plot behaviour.
#'
#' @param hcontrol A list produced by [hmht.control()] defining highlighted
#'   markers or regions and labels.
#'
#' @param ... Additional graphical arguments passed to [graphics::plot()]
#'   (e.g. `pch`, `bg`, `xlab`, `ylab`, `ylim`).
#'   If `mar` is supplied, the default margins are not modified.
#'
#' @details
#' Produces a Manhattan plot in which genomic markers are arranged along the
#' x-axis by chromosome and position.
#'
#' By default the y-axis shows `-log10(p)`, the standard GWAS significance
#' scale. Set `logscale = FALSE` in [mht.control()] to plot raw values instead.
#'
#' Rows containing missing values are removed automatically before plotting.
#'
#' ## Chromosome handling
#' Chromosomes may be numeric or character. The prefix `"chr"` is removed
#' automatically. Chromosomes are ordered numerically first, followed by
#' character chromosomes (e.g. X, Y, MT).
#'
#' ## X-axis spacing
#' * `usepos = FALSE` (default): markers are spaced evenly within chromosomes.
#' * `usepos = TRUE`: real chromosomal positions are used.
#' * `gap` inserts spacing between chromosomes when using real positions.
#'
#' ## Colours
#' Chromosome colours alternate by default. Custom colours can be supplied via
#' `colors` in [mht.control()]; colours are recycled as needed.
#'
#' ## Axis tuning
#' Axis appearance can be controlled using:
#' * `axis.cex` — tick label size
#' * `axis.lwd` — axis and tick thickness
#' * `axis.tck` — tick length and direction
#'
#' These settings are particularly useful when exporting high-resolution figures.
#'
#' ## Horizontal thresholds
#' Horizontal reference lines can be drawn using the `cutoffs` parameter.
#'
#' ## Highlighting regions
#' Highlighted regions are defined using [hmht.control()]. Matching markers are
#' recoloured, labelled, and optionally boxed. Matching is performed by
#' chromosome and position range.
#'
#' @section Mathematical background:
#'
#' A Manhattan plot visualises genome-wide association statistics by mapping
#' genetic markers to a two-dimensional coordinate system.
#'
#' \strong{Y-axis transformation}
#'
#' P-values are transformed to improve visibility of small values:
#'
#' \deqn{y = -\log_{10}(p)}
#'
#' Small p-values therefore appear as large positive values and form the
#' characteristic peaks used to identify strong associations.
#'
#' \strong{Genome linearisation}
#'
#' Markers originate from multiple chromosomes and must be mapped onto a
#' single continuous axis. For chromosome \eqn{c}:
#'
#' \deqn{x_i = pos_i + offset_c}
#'
#' where the chromosome offset is
#'
#' \deqn{offset_c = \sum_{k < c} \max(pos_k)}
#'
#' This produces a continuous genome coordinate system in which chromosomes
#' appear sequentially without overlap.
#'
#' \strong{Chromosome label placement}
#'
#' Chromosome labels are positioned at the midpoint of each chromosome:
#'
#' \deqn{mid_c = (\min(x_c) + \max(x_c)) / 2}
#'
#' This centres labels regardless of chromosome length.
#'
#' \strong{Plot limits}
#'
#' The plotting region is defined using axis tick locations so that ticks and
#' plotting limits coincide exactly.
#'
#' Conceptually, a Manhattan plot is defined by two transformations:
#'
#' \deqn{x = \text{linearised genome coordinate}}
#' \deqn{y = -\log_{10}(p)}
#'
#' All other elements (colouring, thresholds, highlighting) are visual
#' annotations applied to these transformed coordinates.
#'
#' @return invisibly returns `NULL`.
#'
#' @seealso [mht.control()], [hmht.control()]
#' @export
#'
#' @examples
#' \donttest{
#' ## -----------------------------------------------------------
#' ## 1. Minimal example
#' ## -----------------------------------------------------------
#' test <- matrix(c(
#'   1,1,4,
#'   1,1,6,
#'   1,10,3,
#'   2,1,5,
#'   2,2,6,
#'   2,4,8),
#'   byrow=TRUE, ncol=3)
#' mhtplot(test)
#' ## Raw values instead of -log10
#' mhtplot(test, mht.control(logscale = FALSE))
#'
#' ## -----------------------------------------------------------
#' ## 2. Simulated GWAS dataset
#' ## -----------------------------------------------------------
#' set.seed(1)
#' affy <- c(40220,41400,33801,32334,32056,31470,25835,27457,22864,
#'           28501,26273,24954,19188,15721,14356,15309,11281,14881,
#'           6399,12400,7125,6207)
#' n.markers <- sum(affy)
#' n.chr <- length(affy)
#' gwas <- data.frame(
#'   chr = rep(1:n.chr, affy),
#'   pos = 1:n.markers,
#'   p   = runif(n.markers)
#' )
#' mhtplot(gwas)
#'
#' ## -----------------------------------------------------------
#' ## 3. Publication-style figure
#' ## -----------------------------------------------------------
#' ops <- mht.control(
#'   axis.cex = 1.4,
#'   axis.lwd = 2,
#'   axis.tck = -0.03
#' )
#' mhtplot(gwas, control = ops, pch = 19)
#'
#' ## -----------------------------------------------------------
#' ## 4. Real positions + chromosome gaps
#' ## -----------------------------------------------------------
#' ops2 <- mht.control(usepos = TRUE, gap = 10000)
#' mhtplot(gwas, control = ops2, pch = 19)
#'
#' ## -----------------------------------------------------------
#' ## 5. Genome-wide significance thresholds
#' ## -----------------------------------------------------------
#' ops3 <- mht.control(cutoffs = c(5, 7.3))
#' mhtplot(gwas, control = ops3, pch = 19)
#'
#' ## -----------------------------------------------------------
#' ## 6. Highlight selected genes
#' ## -----------------------------------------------------------
#' hdata <- data.frame(
#'   chr  = c(1,3,5),
#'   pos  = c(10000,50000,90000),
#'   p    = c(1e-8,1e-7,1e-9),
#'   gene = c("GENE1","GENE2","GENE3")
#' )
#' hops <- hmht.control(
#'   data = hdata,
#'   colors = "red",
#'   boxed = TRUE
#' )
#' mhtplot(gwas, control = ops, hcontrol = hops, pch = 19)
#' ## A real study (data from gap.datasets package)
#' if (requireNamespace("gap.datasets", quietly = TRUE)) {
#'   data("mhtdata", package = "gap.datasets")
#'
#'   data <- with(mhtdata, cbind(chr, pos, p))
#'   glist <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3",
#'              "WISP3","PPP1R3B","RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
#'   hdata <- subset(mhtdata, gene %in% glist)[c("chr","pos","p","gene")]
#'
#'   ops  <- mht.control(colors = rep(c("lightgray","gray"),11),
#'                       labels = paste("chr",1:22,sep=""),
#'                       yline = 1.5, xline = 3)
#'   hops <- hmht.control(data = hdata, colors = "red", boxed = TRUE)
#'
#'   mhtplot(data, ops, hops, pch = 19)
#'   title("Manhattan plot with genes highlighted")
#' }
#'
#' ## -----------------------------------------------------------
#' ## 7. Export high-resolution PNG
#' ## -----------------------------------------------------------
#' png("manhattan.png", height = 3600, width = 6000, res = 600)
#' opar <- par()
#' par(cex = 0.7)
#' mhtplot(gwas, control = ops, pch = 19)
#' par(opar)
#' dev.off()
#'
#' ## -----------------------------------------------------------
#' ## 8. Miamiplot (see vignette for polished ones)
#' ## -----------------------------------------------------------
#' gwas <- within(gwas, {pr=1-p})
#' miamiplot(test,chr="chr",bp="pos",p="p",pr="pr")
#' }
#'
#' @author Jing Hua Zhao
#'
mhtplot <- function(data, control=mht.control(), hcontrol=hmht.control(), ...) {
  data2 <- data[!apply(is.na(data),1,any),]
  chr_raw <- data2[,1]
  chr_clean <- gsub("^chr", "", as.character(chr_raw), ignore.case=TRUE)
  chr_unique <- unique(chr_clean)
  chr_num <- suppressWarnings(as.numeric(chr_unique))
  num_chr  <- sort(unique(chr_num[!is.na(chr_num)]))
  char_chr <- sort(unique(chr_unique[is.na(chr_num)]))
  chr_levels <- c(as.character(num_chr), char_chr)
  data2 <- as.data.frame(data2)
  data2[,1] <- factor(chr_clean, levels = chr_levels, ordered = TRUE)
  data2 <- data2[order(data2[,1], data2[,2]), ]
  chr <- droplevels(data2[[1]])
  pos <- data2[,2]
  p <- data2[,3]
  tablechr_counts <- tabulate(chr)
  allchr <- tablechr_counts[tablechr_counts > 0]
  n.chr <- length(allchr)
  type <- control$type
  usepos <- control$usepos
  logscale <- control$logscale
  base <- control$base
  cutoffs <- control$cutoffs
  labels <- control$labels
  gap <- control$gap
  pcex <- control$cex
  labcex <- control$lab.cex
  yline <- control$yline
  xline <- control$xline
  verbose <- control$verbose
  axis.cex <- control$axis.cex
  axis.lwd <- control$axis.lwd
  axis.tck <- control$axis.tck
  full_chr_levels <- levels(chr)
  full_palette <- control$colors
  if (is.null(full_palette)) {
      full_palette <- rep(c("grey30","grey70"), length.out = length(full_chr_levels))
  } else {
      full_palette <- rep(full_palette, length.out = length(full_chr_levels))
  }
  names(full_palette) <- full_chr_levels
  chr_colors <- full_palette[levels(chr)[tablechr_counts > 0]]
  tablechr <- levels(chr)[tablechr_counts > 0]
  if(is.null(labels)) labels <- tablechr
  labels <- rep(labels, length.out = n.chr)
  if(is.null(gap)) gap <- 0
  if (!is.null(hcontrol$data))
  {
     hdata <- hcontrol$data
     hdata2 <- hdata[!apply(is.na(hdata),1,any),]
     hchr <- hdata2[,1]
     hchr <- gsub("^chr", "", as.character(hchr), ignore.case=TRUE)
     hpos <- hdata2[,2]
     hname <- hdata2[,4]
     hcolors <- hcontrol$colors
     hyoffs <- hcontrol$yoffset
     hboxed <- hcontrol$boxed
     hcex <- hcontrol$cex
     htablechr <- unique(hchr)
     hlabels <- unique(hname)
     htablechrname <- unique(data.frame(hchr, hname))
     if (is.null(hcolors)) {
         hcolors <- rep("red", length(hlabels))
     } else {
         hcolors <- rep(hcolors, length.out = length(hlabels))
     }
  }
  CMindex <- cumsum(allchr)
  CM <- numeric(length(pos))
  offset <- 0
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     idx <- l:u
     if(usepos)
        CM[idx] <- pos[idx] - min(pos[idx]) + offset
     else
        CM[idx] <- seq_along(idx) + offset
     offset <- max(CM[idx]) + gap
  }
  args <- list(...)
  y_all <- if (logscale) -log(p, base=base) else p
  xlim <- c(0, max(CM))
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  if (!"mar" %in% names(args)) {
    par(mar = c(7, 6, 3, 2))
  }
  par(xaxt="n", yaxt="n")
  ymax <- max(y_all, na.rm = TRUE)
  if (!is.null(cutoffs))
      ymax <- max(ymax, cutoffs, na.rm = TRUE)
  yticks <- pretty(c(0, ymax))
  yticks <- yticks[yticks >= 0]
  ylim <- range(yticks)
  xy <- xy.coords(CM,y_all)
  plot(xy$x,xy$y,type="n",ann=FALSE,axes=FALSE,...)
  axis(1)
  axis(2)
  par(xaxt = "s", yaxt = "s") 
  chr_mid <- numeric(n.chr)
  for(i in 1:n.chr) {
    u <- CMindex[i]
    l <- CMindex[i] - allchr[i] + 1
    chr_mid[i] <- (min(CM[l:u]) + max(CM[l:u])) / 2
  }
  axis(1, at = chr_mid, labels = labels, las = 2,
       cex.axis = axis.cex * labcex,
       lwd = axis.lwd,
       lwd.ticks = axis.lwd,
       tck = axis.tck,
       col = "black",
       col.axis = "black")
  axis(2, at = yticks, las = 1,
       cex.axis = axis.cex,
       lwd = axis.lwd,
       lwd.ticks = axis.lwd,
       tck = axis.tck,
       col = "black",
       col.axis = "black")
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     idx <- l:u
     if (logscale) y <- -log(p[idx],base=base) else y <- p[idx]
     if (verbose) message("Plotting points ",l,"-",u,"\n")
     col.chr <- chr_colors[i]
     if(type=="l") lines(CM[idx],y,col=col.chr,cex=pcex)
     else points(CM[idx],y,col=col.chr,cex=pcex)
  }
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     idx <- l:u
     if (logscale) y <- -log(p[idx],base=base) else y <- p[idx]
     col.chr <- chr_colors[i]
     if (!is.null(hcontrol$data))
     {
        chrs <- htablechrname[as.character(tablechr[i]) == as.character(htablechrname[,1]), ]
        if(dim(chrs)[1]>0) {
          hchrs <- as.character(chrs[,2])
          for(k in 1:length(hchrs))
          {
             hregion <- hpos[hchr==chrs[k,1]&hname==hchrs[k]]
             l1 <- which(pos[idx] >= min(hregion))[1]
             l2 <- tail(which(pos[idx] <= max(hregion)),1)
             if (length(l1)==0 || length(l2)==0) next
             if (verbose) message("  ... highlighting", hchrs[k], "\n")
             col.chr[l1:l2] <- hcolors[ match(hchrs[k], hlabels) ]
             hx <- CM[idx][l1]
             hy <- max(y[l1:l2]) + hyoffs
           # if (hboxed) {text(hx, hy, hchrs[k], cex=1, font=3, bg="white")} else {text(hx, hy, hchrs[k], cex=1, font=3)}
             hx <- CM[idx][l1]
             hy <- max(y[l1:l2]) + hyoffs
             lab <- hchrs[k]
             if (hboxed) {
                usr <- par("usr")
                pin <- par("pin")
                xinch <- diff(usr[1:2]) / pin[1]
                yinch <- diff(usr[3:4]) / pin[2]
                w <- strwidth(lab, units="inches") * xinch * 1.3
                h <- strheight(lab, units="inches") * yinch * 1.6
                rect(hx - w/2, hy - h/2, hx + w/2, hy + h/2,
                     col="white", border="black")
                text(hx, hy, lab, cex=hcex, font=3)
             } else {
                text(hx, hy, lab, cex=hcex, font=3)
             }
             points(CM[idx][l1:l2],y[l1:l2],col=col.chr[l1:l2],cex=pcex)
          }
        }
     }
  }
  if(!is.null(cutoffs)) segments(0, cutoffs, max(CM), cutoffs) # abline(h=cutoffs)
  if ("ylab" %in% names(args)) {
      ylabel <- args$ylab
  } else {
      ylabel <- if (logscale)
          paste0("-log", base, "(Observed value)")
      else
          "Observed value"
  }
  mtext(ylabel, 2, line=yline, las=0, cex=par("cex.lab"))
  if ("xlab" %in% names(args)) {
      xlabel <- args$xlab
  } else {
      xlabel <- "Chromosome"
  }
  mtext(xlabel, 1, line=xline, las=0, cex=par("cex.lab"))
  invisible(NULL)
}
#' Controls for Manhattan plot
#'
#' Parameter specification helper for [mhtplot()]. This function creates a list
#' of graphical and behavioural settings used when generating Manhattan plots.
#'
#' @param type Character. Either `"p"` (points) or `"l"` (lines).
#' @param usepos Logical. Use real chromosomal positions instead of ordinal
#'   marker order.
#' @param logscale Logical. If `TRUE`, values are transformed using
#'   `-log(base)(value)` before plotting.
#' @param base Numeric. Base of logarithm used when `logscale = TRUE`.
#' @param cutoffs Numeric vector of horizontal reference lines to draw.
#' @param colors Vector of chromosome colours. Recycled as needed.
#' @param labels Optional chromosome labels for the x-axis.
#' @param gap Numeric. Gap inserted between chromosomes on the x-axis.
#' @param cex Numeric. Scaling factor for plotted points.
#' @param lab.cex Numeric. Scaling factor for chromosome labels on the x-axis.
#' @param axis.cex Numeric. Scaling factor for axis tick labels.  
#'   Increase when exporting high-resolution figures.
#' @param axis.lwd Numeric. Line width for axes and tick marks.
#' @param axis.tck Numeric. Length and direction of tick marks.  
#'   Negative values draw ticks outward (recommended for publication plots).
#' @param yline Numeric. Margin line for the y-axis label.
#' @param xline Numeric. Margin line for the x-axis label.
#' @param verbose Logical. Print plotting progress messages.
#'
#' @examples
#' mht.control()
#' @return a named list of control parameters for [mhtplot()].
#'
#' @seealso [mhtplot()], [hmht.control()]
#' @export
#'
mht.control <- function(type="p", usepos=FALSE, logscale=TRUE, base=10,
                        cutoffs=NULL, colors=NULL, labels=NULL, gap=NULL,
                        cex=0.4, lab.cex=1,
                        axis.cex=1.2, axis.lwd=1.2, axis.tck=-0.02,
                        yline=3, xline=3, verbose=FALSE){
  list(type=type,usepos=usepos,logscale=logscale,base=base,
       cutoffs=cutoffs,colors=colors,labels=labels,gap=gap,
       cex=cex,lab.cex=lab.cex,
       axis.cex=axis.cex,axis.lwd=axis.lwd,axis.tck=axis.tck,
       yline=yline,xline=xline,verbose=verbose)
}

#' Controls for highlighted regions in mhtplot
#'
#' Helper function defining highlighted markers or genomic regions to be
#' emphasised in [mhtplot()].
#'
#' @param data Data frame with **four columns**:
#'   chromosome, position, value and label (e.g. gene name).
#' @param colors Character vector of colours used for highlighted regions.
#'   Colours are recycled if multiple labels are present.
#' @param yoffset Numeric vertical offset added to the highest highlighted
#'   point before placing the label.
#' @param cex Numeric scaling factor for label text.
#' @param boxed Logical. If `TRUE`, labels are drawn inside a white box with
#'   black border.
#'
#' @return A list of highlight control parameters for [mhtplot()].
#' @seealso [mhtplot()], [mht.control()]
#' @export
#'
#' @examples
#' ## Example highlight specification
#' hdata <- data.frame(
#'   chr  = c(1,3,5),
#'   pos  = c(10000,50000,90000),
#'   p    = c(1e-8,1e-7,1e-9),
#'   gene = c("GENE1","GENE2","GENE3")
#' )
#' hmht.control(data = hdata, colors = "red", boxed = TRUE)
#'
hmht.control <- function(data=NULL, colors="red", yoffset=0.25, cex=1.2, boxed=FALSE) {
  list(data=data,colors=colors,yoffset=yoffset,cex=cex,boxed=boxed)
}
