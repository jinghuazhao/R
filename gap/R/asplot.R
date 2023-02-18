#' Regional association plot
#'
#' @param locus Data frame with columns c("CHR", "POS", "NAME", "PVAL", "RSQR") containing association results.
#' @param map Genetic map, i.e, c("POS","THETA","DIST").
#' @param genes Gene annotation with columns c("START", "STOP", "STRAND", "GENE").
#' @param flanking Flanking length.
#' @param best.pval Best p value for the locus of interest.
#' @param sf scale factors for p values and recombination rates, smaller values are necessary for gene dense regions.
#' @param logpmax Maximum value for -log10(p).
#' @param pch Plotting character for the SNPs to be highlighted, e.g., 21 and 23 refer to circle and diamond.
#'
#' @details
#' This function obtains regional association plot for a particular locus, based on
#' the information about recombinatino rates, linkage disequilibria between the
#' SNP of interest and neighbouring ones, and single-point association tests p values.
#'
#' Note that the best p value is not necessarily within locus in the original design.
#'
#' @export
#' @references
#' \insertRef{saxena07}{gap}
#'
#' @examples
#' \dontrun{
#' require(gap.datasets)
#' asplot(CDKNlocus, CDKNmap, CDKNgenes)
#' title("CDKN2A/CDKN2B Region")
#' asplot(CDKNlocus, CDKNmap, CDKNgenes, best.pval=5.4e-8, sf=c(3,6))
#'
#' ## NCBI2R
#'
#' options(stringsAsFactors=FALSE)
#' p <- with(CDKNlocus,data.frame(SNP=NAME,PVAL))
#' hit <- subset(p,PVAL==min(PVAL,na.rm=TRUE))$SNP
#'
#' library(NCBI2R)
#' # LD under build 36
#' chr_pos <- GetSNPInfo(with(p,SNP))[c("chr","chrpos")]
#' l <- with(chr_pos,min(as.numeric(chrpos),na.rm=TRUE))
#' u <- with(chr_pos,max(as.numeric(chrpos),na.rm=TRUE))
#' LD <- with(chr_pos,GetLDInfo(unique(chr),l,u))
#' # We have complaints; a possibility is to get around with 
#' # https://ftp.ncbi.nlm.nih.gov/hapmap/
#' hit_LD <- subset(LD,SNPA==hit)
#' hit_LD <- within(hit_LD,{RSQR=r2})
#' info <- GetSNPInfo(p$SNP)
#' haldane <- function(x) 0.5*(1-exp(-2*x))
#' locus <- with(info, data.frame(CHR=chr,POS=chrpos,NAME=marker,
#'                     DIST=(chrpos-min(chrpos))/1000000,
#'                     THETA=haldane((chrpos-min(chrpos))/100000000)))
#' locus <- merge.data.frame(locus,hit_LD,by.x="NAME",by.y="SNPB",all=TRUE)
#' locus <- merge.data.frame(locus,p,by.x="NAME",by.y="SNP",all=TRUE)
#' locus <- subset(locus,!is.na(POS))
#' ann <- AnnotateSNPList(p$SNP)
#' genes <- with(ann,data.frame(ID=locusID,CLASS=fxn_class,PATH=pathways,
#'                              START=GeneLowPoint,STOP=GeneHighPoint,
#'                              STRAND=ori,GENE=genesymbol,BUILD=build,CYTO=cyto))
#' attach(genes)
#' ugenes <- unique(GENE)
#' ustart <- as.vector(as.table(by(START,GENE,min))[ugenes])
#' ustop <- as.vector(as.table(by(STOP,GENE,max))[ugenes])
#' ustrand <- as.vector(as.table(by(as.character(STRAND),GENE,max))[ugenes])
#' detach(genes)
#' genes <- data.frame(START=ustart,STOP=ustop,STRAND=ustrand,GENE=ugenes)
#' genes <- subset(genes,START!=0)
#' rm(l,u,ugenes,ustart,ustop,ustrand)
#' # Assume we have the latest map as in CDKNmap
#' asplot(locus,CDKNmap,genes)
#' }
#'
#' @author{Paul de Bakker, Jing Hua Zhao, Shengxu Li}
#' @keywords hplot

asplot <- function (locus, map, genes, flanking=1e3, best.pval=NULL, sf=c(4,4), logpmax=10, pch=21)
{
    NAME <- locus$NAME
    PVAL <- locus$PVAL
    snp <- NAME[PVAL==min(PVAL)]
    chr <- locus$CHR[1]
    hit <- locus[NAME==snp, ]
    if (is.null(best.pval)) best.pval <- hit$PVAL
    lb <- min(locus$POS) - flanking
    ub <- max(locus$POS) + flanking
    lu <- ub - lb
    center <- lb + lu / 2
    center.100kb <- round(center / 1e5) * 1e5
    offset.100kb <- round(lu / 5 / 1e5) * 1e5
    offset <- logpmax / sf[1]
    ylim <- logpmax + offset
    ylim4 <- ylim / 4
    yadj <- -offset + ylim / sf[2]
    keep <- subset(map, map$POS > lb & map$POS < ub)
    genes <- with(genes, subset(genes, (START > lb & START < ub) | (STOP > lb & STOP < ub)))
    par(mar = c(4, 4, 3, 4))
    xy <- xy.coords(keep$POS, yadj + (keep$THETA/60) * (3 * ylim4))
    plot(xy$x, xy$y, type = "l", col = "lightblue", lwd = 1, ylim = c(-offset, logpmax), xlab = "", ylab = "", axes = F)
    box()
    if (offset.100kb!=0)  center5 <- center.100kb + (-2:2) * offset.100kb
    else
    {
       p1 <- min(xy$x)
       p5 <- max(xy$x)
       p3 <- (p1+p5)/2
       p2 <- (p1+p3)/2
       p4 <- (p3+p5)/2
       center5 <- c(p1,p2,p3,p4,p5)
    }
    axis(1, at = center5, labels = round(center5 / 1e3, 2), las = 1)
    mtext(paste("Chromosome", chr, "position (kb)", sep = " "), side = 1, line = 2.5)
    axis(2, at = seq(0, logpmax, 2), labels = seq(0, logpmax, 2), las = 1)
    mtext("-log10(Observed p)", side = 2, at = logpmax / 2, line = 2)
    axis(4, at = yadj + (0:4) * ylim4, labels = paste((0:4)*20), las = 1, col.axis="blue")
    mtext("Recombination rate (cM/Mb)", side = 4, at = -offset + 2 * ylim4, line = 2, col="blue")
    lines(c(lb, ub), c(0, 0), lty = "dotted", lwd = 1, col = "black")
    points(hit$POS, -log10(hit$PVAL), pch = pch, cex = 2.5, bg = "red")
    text(hit$POS, -log10(hit$PVAL), labels = hit$NAME, pos = 3, offset = 1)
    if (-log10(best.pval) < logpmax)
    {
       points(hit$POS, -log10(best.pval), pch = pch, cex = 2.5, bg = "blue")
       text(hit$POS, -log10(best.pval), labels = c(paste("P=", best.pval, sep = "")), pos = 4, offset = 2)
    } else {
       points(hit$POS, logpmax, pch = pch, cex = 2.5, bg = "blue")
       text(hit$POS, logpmax, labels = c(paste("P=", best.pval, sep = "")), pos = 4, offset = 1)
    }
    RSQR <- locus$RSQR
    strong.ld <- subset(locus, NAME != snp & RSQR >= 0.8)
    moderate.ld <- subset(locus, NAME != snp & RSQR >= 0.5 & RSQR < 0.8)
    weak.ld <- subset(locus, NAME != snp & RSQR >= 0.2 & RSQR < 0.5)
    not.in.ld <- subset(locus, NAME != snp & RSQR < 0.2)
    na.ld <- subset(locus, NAME != snp & is.na(RSQR))
    colors <- c("red","orange","yellow","grey","white")
    points(strong.ld$POS, -log10(strong.ld$PVAL), pch = pch, cex = 1.25, bg = colors[1])
    points(moderate.ld$POS, -log10(moderate.ld$PVAL), pch = pch, cex = 1.25, bg = colors[2])
    points(weak.ld$POS, -log10(weak.ld$PVAL), pch = pch, cex = 1.25, bg = colors[3])
    points(not.in.ld$POS, -log10(not.in.ld$PVAL), pch = pch, cex = 1, bg = colors[4])
    points(na.ld$POS, -log10(na.ld$PVAL), pch = pch, cex = 1, bg = colors[5])
    for (i in 1:nrow(genes))
    {
        gi <- genes[i, ]
        GENE <- gi$GENE
        cat("-",GENE,"\n")
        start <- gi$START
        finish <- gi$STOP
        center <- (start + finish)/2
        lb <- min(xy$x)
        ub <- max(xy$x)
        adj <- -offset+2*(i%%4)/3
        if (!is.na(GENE))
        {  
           if (start<lb) text(lb, adj +  ylim4 / 10, labels = GENE, cex = 0.7)
           else if (finish>ub) text(ub, adj +  ylim4 / 10, labels = GENE, cex = 0.7)
           else text(center, adj +  ylim4 / 10, labels = GENE, cex = 0.7)
        }
        STRAND <- gi$STRAND
        if (STRAND == "+") arrows(start, adj, finish, adj, length = 0.05, lwd = 2, code = 2, lty = "solid", col = "darkgreen")
        else arrows(start, adj, finish, adj, length = 0.05, lwd = 2, code = 1, lty = "solid", col = "darkgreen")
    }
    ltext <- rbind("0.8-","0.5-","0.2-","0.0-","NA")
    legend(min(xy$x),logpmax,legend=ltext,title=expression(r^2),fill=colors)
}
