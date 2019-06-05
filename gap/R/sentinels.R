sentinels <- function(p,st,debug=FALSE,flanking=1e+6)
{
  Effect <- End <- StdErr <- prot <- NA
  nr <- nrow(p)
  z <- within(p[st:nr,],{
    d <- c(0,diff(End))
    s <- cumsum(d)
    log10p <- -log10p(Effect/StdErr)
  })
  if (debug) print(z[c("Chrom","End","d","s","MarkerName","P.value")])
  if (z[nrow(z),"s"] <= flanking) {
    l <- z[1, "End"]
    u <- z[nrow(z), "End"]
    log10p1 <- with(z, max(log10p))
    x <- subset(z, log10p==log10p1)
    r1 <- row.names(x)[1]
    m <- x[1,"End"]
    n <- x[1, "MarkerName"]
    cat(prot, n, l, u, u-l, log10p1, r1, "I\n", sep=",")
  } else {
    s <- subset(z, s <= flanking)
    l <- s[1, "End"]
    u <- s[nrow(s), "End"]
    log10p1 <- with(s, max(log10p))
    x <- subset(s, log10p==log10p1)
    r1 <- row.names(x)[1]
    m <- x[1, "End"]
    n <- x[1, "MarkerName"]
    t <- subset(z, End > m & End < m + flanking)
    if (nrow(t)==0) {
      # cat(prot, n, l, u, u-l, log10p1, r1, "II\n", sep=",")
      message(paste0("No variants +1 MB downstream so move to next block (",prot,")"))
      r2 <- as.numeric(r1) + 1
      sentinels(p, r2)
    } else {
      log10p2 <- with(t, max(log10p))
      y <- subset(t, log10p==log10p2)
      u <- p[nrow(t), "End"]
      r2 <- row.names(t)[nrow(t)]
      if (log10p2 < log10p1) {
        cat(prot, n, l, u, u-l, log10p1, r1, "III\n", sep=",")
        if (r2 < nr) sentinels(p, r2)
      } else {
        r2 <- as.numeric(row.names(y)[nrow(y)])
        if(r2 < nr) sentinels(p, r2)
      }
    }
  }
}