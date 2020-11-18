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

swap <- function(x,y)
   eval(parse(text = paste("swap_unique_var_a <-", substitute(x), ";",
   substitute(x), "<-", substitute(y), ";",
   substitute(y), "<-swap_unique_var_a")), env=parent.frame())

pqtlMR <- function(Ins,Ids,prefix="INF1",reverse=FALSE)
{
  Ins <- data <- Ins
  Ins <- TwoSampleMR::format_data(Ins, type = "exposure", header = TRUE,
                     phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                     se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                     other_allele_col = "other_allele", pval_col = "pval")
# ao <- TwoSampleMR::available_outcomes(access_token=NULL)
  ids <- Ids
  outcome_dat <- TwoSampleMR::extract_outcome_data(snps = with(Ins,SNP), outcomes = ids)
  harmonise <- TwoSampleMR::harmonise_data(exposure_dat = Ins, outcome_dat = outcome_dat)
  if (reverse) harmonise <- within(harmonise,
  {
    swap(exposure,outcome)
    swap(id.exposure,id.outcome)
    swap(effect_allele.exposure,effect_allele.outcome)
    swap(other_allele.exposure,other_allele.outcome)
    swap(eaf.exposure,eaf.outcome)
    if(exists("samplesize.exposure") & exists("samplesize.outcome")) swap(samplesize.exposure,samplesize.outcome)
    swap(beta.exposure,beta.outcome)
    swap(se.exposure,se.outcome)
    swap(pval.exposure,pval.outcome)
  })
  result <- heterogeneity <- pleiotropy <- single <- NULL
  try(result <- TwoSampleMR::mr(harmonise, method_list=c("mr_wald_ratio", "mr_ivw"))) # main MR analysis
  heterogeneity <- TwoSampleMR::mr_heterogeneity(harmonise) # heterogeneity test across instruments
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonise) # MR-Egger intercept test
  try(single <- TwoSampleMR::mr_singlesnp(harmonise)) #single SNP MR using Wald ratio
# result <- within(result,outcome <- sub(" [|]* id:ieu-a-[0-9]*| [|]* id:ukb-a-[0-9]*", "\\1", outcome, perl = TRUE))
  ext <- ".txt"
  invisible(lapply(c("harmonise","result","heterogeneity","pleiotropy","single"), function(x) {
                   v <- lapply(x, function(x) tryCatch(get(x), error=function(e) NULL))[[1]]
                   if (!is.null(v)) write.table(format(v,digits=3),file=paste0(prefix,"-",x,ext),quote=FALSE,row.names=FALSE,sep="\t")
            })
  )
}

uniprot2ids <- function(uniprotid="ACC+ID",to,query)
{
  rt <- find.package("pQTLtools")
  f <- file.path(rt ,"python","uniprot2ids.py")
  reticulate::source_python(f)
  invisible(uniprot2ids(uniprotid,to,query))
}
