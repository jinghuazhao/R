\name{ps}
\alias{snpquery}
\title{phenoscanner snpqueries in batches}
\usage{snpquery(query,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)}
\description{
R/phenoscanner only allows for certain number of items supplied. This simple function return
a large number of calls by doing calls in batches as well as generating SNPIDs.
}
\arguments{
\item{query}{a list of SNPs}
\item{catalogue}{"None","eQTL","mQTL","methQTL","pQTL","GWAS"}
\item{proxies}{"None", "AFR","AMR","EAS","EUR","SAS"}
\item{p}{p value threshold}
\item{r2}{r2 for LD}
\item{build}{37, 38}
}
\details{
Batches are generated and queries are combined with do.call.
}

\value{
The returned value is a list containing snps and results:
}

\section{References}{
Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. Nature 558: 73-79.
}
\seealso{
\code{\link[phenoscanner]{phenoscanner}}
}

\examples{
\dontrun{
# single SNP
snpquery(query="rs704")
# SCALLOP/INF
INF <- Sys.getenv("INF")
rsid <- scan(paste(INF,'work','INF1.merge.snp',sep='/'),"")
r <- snpquery(query=rsid,catalogue='pQTL',p=1e-11)
INTERVAL_Olink <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
save(INTERVAL_Olink,file='INTERVAL_Olink.rda',compress='xz')
# Sun, et al. (2018) Supplementary tables
dir <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/'
file <- '41586_2018_175_MOESM4_ESM.xlsx'
xlsx <- paste0(dir,file)
t4.1 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE, 
                            cols=c(1:16,26:28), rows=c(5:1986))
t4.2 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE,
                            cols=c(17:25,29:31), rows=c(6:1986))
t4 <- cbind(t4.1,t4.2)
r <- snpquery(query=t4[,6],catalogue='pQTL',p=1e-11)
INTERVAL_SomaLogic <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
save(INTERVAL_SomaLogic,file='INTERVAL_SomaLogic.rda',compress='xz')
t5 <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE,
                          cols=c(1:19), rows=c(3:2746))
t6 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE,
                          cols=c(1:20), rows=c(3:167))
}
}
\note{adapted from custom codings}
\keyword{utilities}
