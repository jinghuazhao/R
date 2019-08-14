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
