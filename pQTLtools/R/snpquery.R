snpquery <- function(snplist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)
{
  ref_a1 <- ref_a2 <- ref_hg19_coordinates <- NULL
  batches <- split(snplist,ceiling(seq_along(snplist)/100))
  s <- r <- vector('list',length(batches))
  for(i in 1:length(batches))
  {
    cat("Block",i,"\n")
    q <- phenoscanner::phenoscanner(snpquery=batches[[i]], catalogue=catalogue,
                                    proxies=proxies, pvalue=p, r2=r2, build=build)
    s[[i]] <- with(q,snps)
    r[[i]] <- with(q,results)
  }
  snps <- do.call("rbind",s)
  results <- within(do.call("rbind",r),
  {
     a1 <- as.character(ref_a1)
     a2 <- as.character(ref_a2)
     swap <- ref_a1 > ref_a2
     a1[swap] <- ref_a2[swap]
     a2[swap] <- ref_a1[swap]
     ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
  })
  list(snps=snps,results=results)
}
