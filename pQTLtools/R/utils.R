genequeries <- function(genelist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)
{
  ref_a1 <- ref_a2 <- ref_hg19_coordinates <- NULL
  batches <- split(genelist,ceiling(seq_along(genelist)/10))
  g <- r <- vector('list',length(batches))
  for(i in 1:length(batches))
  {
    cat("Block",i,"\n")
    q <- phenoscanner::phenoscanner(genequery=batches[[i]], catalogue=catalogue,
                                    proxies=proxies, pvalue=p, r2=r2, build=build)
    g[[i]] <- with(q,genes)
    r[[i]] <- with(q,results)
  }
  genes <- do.call("rbind",g)
  results <- within(do.call("rbind",r),
  {
     a1 <- as.character(ref_a1)
     a2 <- as.character(ref_a2)
     swap <- ref_a1 > ref_a2
     a1[swap] <- ref_a2[swap]
     a2[swap] <- ref_a1[swap]
     ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
  })
  list(genes=genes,results=results)
}

regionqueries <- function(regionlist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)
{
  ref_a1 <- ref_a2 <- ref_hg19_coordinates <- NULL
  batches <- split(regionlist,ceiling(seq_along(regionlist)/10))
  s <- r <- vector('list',length(batches))
  for(i in 1:length(batches))
  {
    cat("Block",i,"\n")
    q <- phenoscanner::phenoscanner(regionquery=batches[[i]], catalogue=catalogue,
                                    proxies=proxies, pvalue=p, r2=r2, build=build)
    s[[i]] <- with(q,regions)
    r[[i]] <- with(q,results)
  }
  regions <- do.call("rbind",s)
  results <- within(do.call("rbind",r),
  {
     a1 <- as.character(ref_a1)
     a2 <- as.character(ref_a2)
     swap <- ref_a1 > ref_a2
     a1[swap] <- ref_a2[swap]
     a2[swap] <- ref_a1[swap]
     ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
  })
  list(regions=regions,results=results)
}

snpqueries <- function(snplist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)
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
