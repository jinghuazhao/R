sentinels <- function(p,pid,st,debug=FALSE,flanking=1e+6,chr="Chrom",pos="End",b="Effect",se="StdErr",snp="MarkerName")
{
  Effect <- NA
  End <- NA
  StdErr <- NA
  nr <- nrow(p)
  z <- within(p[st:nr,],{
    d <- c(0,diff(End))
    s <- cumsum(d)
    logp <- -logp(Effect/StdErr)
  })
  if (debug) print(z[c(chr,pos,"d","s",snp,"logp")])
  if (tail(z[,"s"], 1) <= flanking) {
    l <- head(z[[pos]], 1)
    u <- tail(z[[pos]], 1)
    logp1 <- with(z, max(logp))
    x <- subset(z, logp==logp1)
    r1 <- row.names(x)[1]
    m <- tail(x[[pos]], 1)
    n <- tail(x[[snp]], 1)
    cat(pid, n, l, u, u-l, logp1, r1, "I\n", sep=",")
  } else {
    s <- subset(z, s <= flanking)
    l <- head(s[[pos]], 1)
    u <- tail(s[[pos]], 1)
    logp1 <- with(s, max(logp))
    x <- subset(s, logp==logp1)
    r1 <- tail(row.names(x), 1)
    m <- tail(x[[pos]], 1)
    n <- tail(x[[snp]], 1)
    t <- subset(z, z[[pos]] > m & z[[pos]] <= m + flanking)
    if (nrow(t)==0) {
      cat(pid, n, l, u, u-l, logp1, r1, "II\n", sep=",")
      message(paste0("No variants +1 MB downstream so move to next block (",pid,")"))
      r2 <- as.numeric(r1) + 1
      sentinels(p, pid, r2)
    } else {
      logp2 <- with(t, max(logp))
      y <- subset(t, logp==logp2)
      u <- tail(t[[pos]], 1)
      r2 <- as.numeric(tail(row.names(t), 1))
      if (logp1 > logp2) {
        cat(pid, n, l, u, u-l, logp1, r1, "III\n", sep=",")
        if (r2 < nr) sentinels(p, pid, r2+1)
      } else {
        r2 <- as.numeric(tail(row.names(y),1))
        if (r2 < nr) sentinels(p, pid, r2)
      }
    }
  }
}
