genequeries <- function(genelist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)
{
  ref_a1 <- ref_a2 <- ref_hg19_coordinates <- ref_hg38_coordinates <- NULL
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
     ref_a1 <- as.character(a1)
     ref_a2 <- as.character(a2)
     swap <- ref_a1 > ref_a2
     a1[swap] <- ref_a2[swap]
     a2[swap] <- ref_a1[swap]
     ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
     if (build==38) ref_snpid <- paste0(ref_hg38_coordinates,"_",a1,"_",a2)
  })
  list(genes=genes,results=results)
}

regionqueries <- function(regionlist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37,wait=TRUE)
{
  ref_a1 <- ref_a2 <- ref_hg19_coordinates <- ref_hg38_coordinates <- NULL
  lrl <- strsplit(regionlist,":|-")
  chr <- as.character(lapply(lrl,"[[",1))
  start <- as.integer(lapply(lrl,"[[",2))
  end <- as.integer(lapply(lrl,"[[",3))
  gr <- GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start,end))
  print(GenomicRanges::width(gr))
  tiles <- GenomicRanges::tile(gr, width=1e+6)
  print(GenomicRanges::width(tiles))
  regionlist_ext <- with(as.data.frame(tiles),paste0(seqnames,":",start,"-",end))
  cat("Conducting queries for",length(regionlist_ext),"regions.\n")
  batches <- split(regionlist_ext,ceiling(seq_along(regionlist_ext)/10))
  s <- r <- vector('list',length(batches))
  for(i in 1:length(batches))
  {
    cat("Block",i,"\n")
    if (wait) if (i%%6==0) Sys.sleep(60*60)
    q <- phenoscanner::phenoscanner(regionquery=batches[[i]], catalogue=catalogue,
                                    proxies=proxies, pvalue=p, r2=r2, build=build)
    s[[i]] <- with(q,regions)
    r[[i]] <- with(q,results)
  }
  regions <- do.call("rbind",s)
  results <- within(do.call("rbind",r),
  {
     ref_a1 <- as.character(a1)
     ref_a2 <- as.character(a2)
     swap <- ref_a1 > ref_a2
     a1[swap] <- ref_a2[swap]
     a2[swap] <- ref_a1[swap]
     ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
     if (build==38) ref_snpid <- paste0(ref_hg38_coordinates,"_",a1,"_",a2)
  })
  list(tiles=tiles,regions=regions,results=results)
}

snpqueries <- function(snplist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)
{
  ref_a1 <- ref_a2 <- ref_hg19_coordinates <- ref_hg38_coordinates <- NULL
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
     ref_a1 <- as.character(a1)
     ref_a2 <- as.character(a2)
     swap <- ref_a1 > ref_a2
     a1[swap] <- ref_a2[swap]
     a2[swap] <- ref_a1[swap]
     ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
     if (build==38) ref_snpid <- paste0(ref_hg38_coordinates,"_",a1,"_",a2)
  })
  list(snps=snps,results=results)
}
