\name{regionqueries}
\alias{regionqueries}
\title{phenoscanner regionqueries in batches}
\usage{regionqueries(regionlist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37,wait=TRUE)}
\description{
R/phenoscanner only allows for certain number of items supplied. This simple function return
a large number of calls in batches as well as generating SNPIDs.
}
\arguments{
\item{regionlist}{a list of SNPs}
\item{catalogue}{"None","eQTL","mQTL","methQTL","pQTL","GWAS"}
\item{proxies}{"None", "AFR","AMR","EAS","EUR","SAS"}
\item{p}{p value threshold}
\item{r2}{r2 for LD}
\item{build}{37, 38}
\item{wait}{a flag to wait for 1hr for every 50 regions}
}
\details{
Batches are generated and queries are combined into one.
}

\value{
The returned value is a list containing tiles, regions and results:
}

\section{References}{
Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73-79.
}
\seealso{
\code{\link[phenoscanner]{phenoscanner}}
}

\examples{
\dontrun{
# single region
regionqueries("chr17:26691290-26700110")

flank <- 1e+6

# SCALLOP
INF <- Sys.getenv("INF")
INF1_merge <- merge(inf1,
                   read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE),
                   by="prot")
INF1_merge_uniprot <- with(INF1_merge,unique(uniprot))
INF1_uniprot <- setdiff(with(inf1,uniprot),"P23560")

# SomaLogic -- everything but stick to those in INF1
SomaLogic_INF1_merge <- subset(SomaLogic160410,UniProt \%in\% INF1_merge_uniprot)
SomaLogic_uniprot <- setdiff(unique(with(SomaLogic_INF1_merge,UniProt)),"P23560")
setdiff(INF1_merge_uniprot,SomaLogic_uniprot)
regions <- within(subset(INF1_merge,uniprot \%in\% SomaLogic_uniprot),
{
  chr <- Chrom
  start <- Start
  end <- End
})
singletons <- with(regions, start==end)
regions[singletons,"start"] <- regions[singletons,"start"] - flank
regions[singletons,"end"] <- regions[singletons,"start"] + flank
reset <- with(regions,start<0)
regions[reset,"start"] <- 0
r <- regionqueries(with(regions,paste0(chr,":",start,"-",end)))
save(r,file="INF1_merge.rda",compress='xz')
r2 <- with(r,
{
  region_ext <- cbind(tiles,regions)
  results_ext <- merge(region_ext,results,by="region")
  ord <- with(results_ext,order(group))
  results_ext[ord,]
})
results <- subset(r2,pmid=="29875488")
grp <- names(table(with(results,group)))
sink("INF1_merge.txt")
options(width=250)
for(g in as.numeric(grp))
{
  uniprot <- regions[g,"uniprot"]
  SNP <- regions[g,"SNP"]
  print(regions[g,])
  s <- subset(results,group==g&rsid==SNP)
  if(nrow(s)>1) print(s[c("region","group","rsid","hg19_coordinates","hgnc","beta","se","p","ref_snpid")])
}
sink()

SomaLogic_INF1 <- subset(SomaLogic160410,UniProt \%in\% INF1_uniprot)
SomaLogic_uniprot <- setdiff(unique(with(SomaLogic_INF1,UniProt)),"P23560")
setdiff(INF1_merge_uniprot,SomaLogic_uniprot)

# SomaLogic -- gene-based, also at GitHub, which may miss some from INF1
HOME <- Sys.getenv("HOME")
INTERVAL_box <- read.delim(file.path(HOME,"SomaLogic","doc","INTERVAL-box.tsv"),as.is=TRUE)
regions <- within(subset(INTERVAL_box, UniProt \%in\% INF_uniprot),
{
  start <- start-flank
  end <- end+flank
})
reset <- with(regions,start<0)
regions[reset,"start"] <- 0
boxlist <- unique(with(regions, paste0("chr",chr,":",start,"-",end)))

# SomaLogic -- sentinel-based, which may miss some from INF1
st4list <- with(subset(st4, UniProt \%in\% INF1_uniprot),
                paste0("chr",Chr,":",Region.start,"-",Region.end))
r <- regionqueries(st4list)

}
}
\note{adapted from custom codings}
\keyword{utilities}
