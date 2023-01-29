#' Sentinel identification from GWAS summary statistics
#'
#' @param p an object containing GWAS summary statistics.
#' @param pid a phenotype (e.g., protein) name in pGWAS.
#' @param st row number as in p.
#' @param debug a flag to show the actual data.
#' @param flanking the width of flanking region.
#' @param chr Chromosome name.
#' @param pos Position.
#' @param b Effect size.
#' @param se Standard error.
#' @param log_p log(P).
#' @param snp Marker name.
#' @param sep field delimitor.
#'
#' @details
#' This function accepts an object containing GWAS summary statistics for
#' signal identification as defined by flanking regions. As the associate P value
#' could be extremely small, the effect size and its standard error are used.
#'
#' A distance-based approach was consequently used and reframed as an algorithm here. It takes as input signals multiple correlated variants in 
#' particular region(s) which reach genomewide significance and output three types of sentinels in a region-based manner. For a given protein and a 
#' chromosome, the algorithm proceeds as follows:
#'
#' Algorithm sentinels
#'
#' Step 1. for a particular collection of genomewide significant variants on a chromosome, the width of the region is calculated according to the start 
#' and end chromosomal positions and if it is smaller than the flanking distance, the variant with the smallest P value is taken as sentinel (I) 
#' otherwise goes to step 2.
#'
#' Step 2. The variant at step 1 is only a candidate and a flanking region is generated. If such a region contains no variant the candidate is recorded 
#' as sentinel (II) and a new iteration starts from the variant next to the flanking region.
#'
#' Step 3.  When the flanking is possible at step 2 but the P value is still larger than the candidate at step 2, the candidate is again recorded as 
#' sentinel (III) but next iteration starts from the variant just after the variant at the end position; otherwise the variant is updated as a new 
#' candidate where the next iteration starts.
#'
#' Note Type I signals are often seen from variants in strong LD at a cis region, type II results seen when a chromosome contains two trans signals, 
#' type III results seen if there are multiple trans signals.
#'
#' Typically, input to the function are variants reaching certain level of significance and the functtion identifies minimum p value at the flanking 
#' interval; in the case of another variant in the flanking window has smaller p value it will be used instead.
#'
#' For now key variables in p are "MarkerName", "End", "Effect", "StdErr", "P.value", where "End" is as in a bed file indicating marker position,
#' and the function is set up such that row names are numbered as 1:nrow(p); see example below. When log_p is specified, log(P) is used instead, which 
#' is appropriate with output from METAL with LOGPVALUE ON. In this case, the column named log(P) in the output is actually log10(P).
#'
#' @export
#' @return The function give screen output.
#'
#' @examples
#' \dontrun{
#'  ## OPG as a positive control in our pGWAS
#'  require(gap.datasets)
#'  data(OPG)
#'  p <- reshape::rename(OPGtbl, c(Chromosome="Chrom", Position="End"))
#'  chrs <- with(p, unique(Chrom))
#'  for(chr in chrs)
#'  {
#'    ps <- subset(p[c("Chrom","End","MarkerName","Effect","StdErr")], Chrom==chr)
#'    row.names(ps) <- 1:nrow(ps)
#'    sentinels(ps, "OPG", 1)
#'  }
#'  subset(OPGrsid,MarkerName=="chr8:120081031_C_T")
#'  subset(OPGrsid,MarkerName=="chr17:26694861_A_G")
#'  ## log(P)
#'  p <- within(p, {logp <- log(P.value)})
#'  for(chr in chrs)
#'  {
#'    ps <- subset(p[c("Chrom","End","MarkerName","logp")], Chrom==chr)
#'    row.names(ps) <- 1:nrow(ps)
#'    sentinels(ps, "OPG", 1, log_p="logp")
#'  }
#'  ## to obtain variance explained
#'  tbl <- within(OPGtbl, chi2n <- (Effect/StdErr)^2/N)
#'  s <- with(tbl, aggregate(chi2n,list(prot),sum))
#'  names(s) <- c("prot", "h2")
#'  sd <- with(tbl, aggregate(chi2n,list(prot),sd))
#'  names(sd) <- c("p1","sd")
#'  m <- with(tbl, aggregate(chi2n,list(prot),length))
#'  names(m) <- c("p2","m")
#'  h2 <- cbind(s,sd,m)
#'  ord <- with(h2, order(h2))
#'  sink("h2.dat")
#'  print(h2[ord, c("prot","h2","sd","m")], row.names=FALSE)
#'  sink()
#'  png("h2.png", res=300, units="in", width=12, height=8)
#'  np <- nrow(h2)
#'  with(h2[ord,], {
#'      plot(h2, cex=0.4, pch=16, xaxt="n", xlab="protein", ylab=expression(h^2))
#'      xtick <- seq(1, np, by=1)
#'      axis(side=1, at=xtick, labels = FALSE)
#'      text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=0.5)
#'  })
#'  dev.off()
#'  write.csv(tbl,file="INF1.csv",quote=FALSE,row.names=FALSE)
#' }
#'
#' @keywords utilities

sentinels <- function(p,pid,st,debug=FALSE,flanking=1e+6,chr="Chrom",pos="End",b="Effect",se="StdErr",
                      log_p=NULL,snp="MarkerName",sep=",")
{
  nr <- nrow(p)
  u <- p[st:nr,]
  z <- within(u,{
    d <- c(0,diff(u[[pos]]))
    s <- cumsum(d)
    if (is.null(log_p)) log10p <- -log10p(u[[b]]/u[[se]])
    else log10p <- -(u[[log_p]]/log(10))
  })
  if (debug) print(z[c(chr,pos,"d","s",snp,"log10p")])
  if (tail(z[,"s"], 1) <= flanking) {
    l <- head(z[[pos]], 1)
    u <- tail(z[[pos]], 1)
    log10p1 <- with(z, max(log10p))
    x <- subset(z, log10p==log10p1)
    r1 <- row.names(x)[1]
    m <- tail(x[[pos]], 1)
    n <- tail(x[[snp]], 1)
    cat(pid, n, l, u, u-l, log10p1, r1, "I\n", sep=sep)
  } else {
    s <- subset(z, s <= flanking)
    l <- head(s[[pos]], 1)
    u <- tail(s[[pos]], 1)
    log10p1 <- with(s, max(log10p))
    x <- subset(s, log10p==log10p1)
    r1 <- tail(row.names(x), 1)
    m <- tail(x[[pos]], 1)
    n <- tail(x[[snp]], 1)
    t <- subset(z, z[[pos]] > m & z[[pos]] <= m + flanking)
    if (nrow(t)==0) {
      cat(pid, n, l, u, u-l, log10p1, r1, "II\n", sep=sep)
      message(paste0("No variants +1 MB downstream so move to next block (",pid,")"))
      r2 <- as.numeric(r1) + 1
      sentinels(p, pid, r2)
    } else {
      log10p2 <- with(t, max(log10p))
      y <- subset(t, log10p==log10p2)
      u <- tail(t[[pos]], 1)
      r2 <- as.numeric(tail(row.names(t), 1))
      if (log10p1 > log10p2) {
        cat(pid, n, l, u, u-l, log10p1, r1, "III\n", sep=sep)
        if (r2 < nr) sentinels(p, pid, r2+1)
      } else {
        r2 <- as.numeric(tail(row.names(y),1))
        if (r2 < nr) sentinels(p, pid, r2)
      }
    }
  }
}
