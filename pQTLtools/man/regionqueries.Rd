\name{regionqueries}
\alias{regionqueries}
\title{phenoscanner regionqueries in batches}
\usage{regionqueries(regionlist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)}
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

# SCALLOP
INF <- Sys.getenv("INF")
INF1_merge <- merge(read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE),
                    inf1, by="prot")
inf1_uniprot <- setdiff(with(INF1_merge,uniprot),"P23560")

# SomaLogic -- gene-based, also at GitHub
HOME <- Sys.getenv("HOME")
flank <- 1e+6
INTERVAL_box <- read.delim(file.path(HOME,"SomaLogic","doc","INTERVAL-box.tsv"),as.is=TRUE)

regions <- within(subset(INTERVAL_box, UniProt \%in\% inf1_uniprot),
{
  start <- start-flank
  end <- end+flank
})
reset <- with(regions,start<0)
regions[reset,"start"] <- 0
regionlist <- unique(with(regions, paste0("chr",chr,":",start,"-",end)))

# SomaLogic -- sentinel-based
regionlist <- with(subset(st4, UniProt \%in\% inf1_uniprot),
                   paste0("chr",Chr,":",Region.start,"-",Region.end))
r <- regionqueries(regionlist)

}
}
\note{adapted from custom codings}
\keyword{utilities}
