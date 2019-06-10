sentinels <- function(p,pid,st,debug=FALSE,flanking=1e+6)
{
  Effect <- End <- StdErr <- NA
  nr <- nrow(p)
  z <- within(p[st:nr,],{
    d <- c(0,diff(End))
    s <- cumsum(d)
    log10p <- -log10p(Effect/StdErr)
  })
  if (debug) print(z[c("Chrom","End","d","s","MarkerName","P.value")])
  if (tail(z[,"s"], 1) <= flanking) {
    l <- head(z[,"End"], 1)
    u <- tail(z[,"End"], 1)
    log10p1 <- with(z, max(log10p))
    x <- subset(z, log10p==log10p1)
    r1 <- row.names(x)[1]
    m <- tail(x[,"End"], 1)
    n <- tail(x[,"MarkerName"], 1)
    cat(pid, n, l, u, u-l, log10p1, r1, "I\n", sep=",")
  } else {
    s <- subset(z, s <= flanking)
    l <- head(s[,"End"], 1)
    u <- tail(s[,"End"], 1)
    log10p1 <- with(s, max(log10p))
    x <- subset(s, log10p==log10p1)
    r1 <- tail(row.names(x), 1)
    m <- tail(x[,"End"], 1)
    n <- tail(x[,"MarkerName"], 1)
    t <- subset(z, End > m & End < m + flanking)
    if (nrow(t)==0) {
      cat(pid, n, l, u, u-l, log10p1, r1, "II\n", sep=",")
      message(paste0("No variants +1 MB downstream so move to next block (",pid,")"))
      r2 <- as.numeric(r1) + 1
      sentinels(p, pid, r2)
    } else {
      log10p2 <- with(t, max(log10p))
      y <- subset(t, log10p==log10p2)
      u <- tail(t[,"End"], 1)
      r2 <- tail(row.names(t), 1)
      if (log10p1 > log10p2) {
        cat(pid, n, l, u, u-l, log10p1, r1, "III\n", sep=",")
        if (r2 <nr) sentinels(p, pid, as.numeric(r2)+1)
      } else {
        cat("Switching sentinel ...\n")
        r2 <- tail(row.names(y),1)
        if (r2 < nr) sentinels(p, pid, as.numeric(r2))
      }
    }
  }
}
