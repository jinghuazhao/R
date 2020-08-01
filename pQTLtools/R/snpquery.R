snpquery <- function(snplist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)
{
  ref_a1 <- ref_a2 <- ref_hg19_coordinates <- NULL
  s <- t <- list()
  batches <- split(snplist,ceiling(seq_along(snplist)/100))
  for(i in 1:length(batches))
  {
    cat("Block ",i,"\n")
    q <- phenoscanner::phenoscanner(snpquery=batches[[i]], catalogue=catalogue, proxies=proxies, pvalue=p, r2=r2, build=build)
    s[[i]] <- with(q,snps)
    t[[i]] <- with(q,results)
  }
  snps <- do.call(rbind,s)
  results <- do.call(rbind,t)
  results <- within(results,
  {
     a1 <- ref_a1
     a2 <- ref_a2
     swap <- ref_a1 > ref_a2
     a1[swap] <- ref_a2[swap]
     a2[swap] <- ref_a1[swap]
     ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
  })
  list(snps=snps,results=results)
}
