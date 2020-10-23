genequeries <- function(genelist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.8,build=37,wait=TRUE)
{
  a1 <- a2 <- hg19_coordinates <- hg38_coordinates <- NULL
  batches <- split(genelist,ceiling(seq_along(genelist)/10))
  g <- r <- vector('list',length(batches))
  for(i in 1:length(batches))
  {
    cat("Block",i,"\n")
    if (wait) if (i%%6==0) Sys.sleep(60*60)
    q <- phenoscanner::phenoscanner(genequery=batches[[i]], catalogue=catalogue,
                                    proxies=proxies, pvalue=p, r2=r2, build=build)
    g[[i]] <- with(q,genes)
    r[[i]] <- with(q,results)
  }
  genes <- do.call("rbind",g)
  results <- within(do.call("rbind",r),
  {
    swap_a1 <- as.character(a1)
    swap_a2 <- as.character(a2)
    swap <- swap_a1 > swap_a2
    a1[swap] <- swap_a2[swap]
    a2[swap] <- swap_a1[swap]
    if (build==37) snpid <- paste0(hg19_coordinates,"_",a1,"_",a2)
    else if (build==38) snpid <- paste0(hg38_coordinates,"_",a1,"_",a2)
    rm(swap_a1,swap_a2,swap)
  })
  list(genes=genes,results=results)
}

regionqueries <- function(regionlist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.8,build=37,wait=TRUE)
{
  a1 <- a2 <- hg19_coordinates <- hg38_coordinates <- NULL
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
    swap_a1 <- as.character(a1)
    swap_a2 <- as.character(a2)
    swap <- swap_a1 > swap_a2
    a1[swap] <- swap_a2[swap]
    a2[swap] <- swap_a1[swap]
    if (build==37) snpid <- paste0(hg19_coordinates,"_",a1,"_",a2)
    else if (build==38) snpid <- paste0(hg38_coordinates,"_",a1,"_",a2)
    rm(swap_a1,swap_a2,swap)
  })
  list(tiles=tiles,regions=regions,results=results)
}

snpqueries <- function(snplist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.8,build=37,wait=TRUE)
{
  a1 <- a2 <- hg19_coordinates <- hg38_coordinates <- NULL
  batches <- split(snplist,ceiling(seq_along(snplist)/100))
  s <- r <- vector('list',length(batches))
  for(i in 1:length(batches))
  {
    cat("Block",i,"\n")
    if (wait) if (i%%6==0) Sys.sleep(60*60)
    q <- phenoscanner::phenoscanner(snpquery=batches[[i]], catalogue=catalogue,
                                    proxies=proxies, pvalue=p, r2=r2, build=build)
    s[[i]] <- with(q,snps)
    r[[i]] <- with(q,results)
  }
  snps <- within(do.call("rbind",s),
  {
    swap_a1 <- as.character(a1)
    swap_a2 <- as.character(a2)
    swap <- swap_a1 > swap_a2
    a1[swap] <- swap_a2[swap]
    a2[swap] <- swap_a1[swap]
    if (build==37) snpid <- paste0(hg19_coordinates,"_",a1,"_",a2)
    else if (build==38) snpid <- paste0(hg38_coordinates,"_",a1,"_",a2)
    rm(swap_a1,swap_a2,swap)
  })
  results <- within(do.call("rbind",r),
  {
    swap_a1 <- as.character(a1)
    swap_a2 <- as.character(a2)
    swap <- swap_a1 > swap_a2
    a1[swap] <- swap_a2[swap]
    a2[swap] <- swap_a1[swap]
    if (build==37) snpid <- paste0(hg19_coordinates,"_",a1,"_",a2)
    else if (build==38) snpid <- paste0(hg38_coordinates,"_",a1,"_",a2)
    rm(swap_a1,swap_a2,swap)
  })
  list(snps=snps,results=results)
}

pqtlMR <- function(Ins,Ids,prefix="INF1")
{
  Ins <- data <- Ins
  Ins <- TwoSampleMR::format_data(Ins, type = "exposure", header = TRUE,
                     phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                     se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                     other_allele_col = "other_allele", pval_col = "pval")
  ao <- TwoSampleMR::available_outcomes(access_token=NULL)
  ids <- Ids
  outcome_dat <- TwoSampleMR::extract_outcome_data(snps = with(Ins,SNP), outcomes = ids)
  dat <- TwoSampleMR::harmonise_data(exposure_dat = Ins, outcome_dat = outcome_dat)
  mr_results <- mr_hetero <- mr_pleio <- mr_single <- NULL
  try(mr_results <- TwoSampleMR::mr(dat, method_list=c("mr_wald_ratio", "mr_ivw"))) # main MR analysis
  mr_hetero <- TwoSampleMR::mr_heterogeneity(dat) # heterogeneity test across instruments
  mr_pleio <- TwoSampleMR::mr_pleiotropy_test(dat) # MR-Egger intercept test
  try(mr_single <- TwoSampleMR::mr_singlesnp(dat)) #single SNP MR using Wald ratio
  options(width=200)
  mr_results <- within(mr_results,outcome <- sub(" [|]* id:ieu-a-[0-9]*| [|]* id:ukb-a-[0-9]*", "\\1", outcome, perl = TRUE))
  filename <- c("harmonise","mr","mr_hetero","mr_pleio","mr_single")
  ext <- "txt"
  result_files <- paste(prefix,filename,ext,sep=".")
  write.table(dat,file=result_files[1],sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  write.table(format(mr_results,digits=3),file=result_files[2],sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  write.table(mr_hetero,file=result_files[3],sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  write.table(mr_pleio,file=result_files[4],sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  write.table(mr_single,file=result_files[5],sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

uniprot2ids <- function(uniprotid="ACC+ID",to,query)
{
  rt <- find.package("pQTLtools")
  f <- file.path(rt ,"python","uniprot2ids.py")
  reticulate::source_python(f)
  invisible(uniprot2ids(uniprotid,to,query))
}
