\name{snpquery}
\alias{snpquery}
\title{phenoscanner snpqueries in batches}
\usage{snpquery(snplist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)}
\description{
R/phenoscanner only allows for certain number of items supplied. This simple function return
a large number of calls in batches as well as generating SNPIDs.
}
\arguments{
\item{snplist}{a list of SNPs}
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
r <- snpquery(rsid,catalogue='pQTL',p=1e-11)
INTERVAL_Olink <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
save(INTERVAL_Olink,file='INTERVAL_Olink.rda',compress='xz')

# SomaLogic
dir <- 'https://static-content.springer.com/esm/art\%3A10.1038\%2Fs41586-018-0175-2/MediaObjects/'
file <- '41586_2018_175_MOESM4_ESM.xlsx'
xlsx <- paste0(dir,file)
t4.1 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE, 
                            cols=c(1:16,26:28), rows=c(5:1986))
t4.2 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE,
                            cols=c(17:25,29:31), rows=c(6:1986))
t4 <- cbind(t4.1,t4.2)
r <- snpquery(t4[,6],catalogue='pQTL',p=1e-11)
INTERVAL_SomaLogic <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
save(INTERVAL_SomaLogic,file='INTERVAL_SomaLogic.rda',compress='xz')
t5 <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE,
                          cols=c(1:19), rows=c(3:2746))
t6 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE,
                          cols=c(1:20), rows=c(3:167))

# Olink-SomaLogic overlap
p <- list(setdiff(inf1$uniprot,"P23560"),
          setdiff(SomaLogic160410$UniProt[!is.na(SomaLogic160410$UniProt)],"P23560"))
cnames <- c("INF1","SomaLogic")
VennDiagram::venn.diagram(x = p, category.names=cnames,
                          filename='os.png', imagetype="png", output=TRUE)
m <- merge(inf1,SomaLogic160410,by.x="uniprot",by.y="UniProt")
u <- setdiff(with(m,unique(uniprot)),"P23560")
options(width=220)
o <- subset(inf1,uniprot %in% u)
dim(o)
o
vars <- c("UniProt","chr","start","end","extGene","Target","TargetFullName")
s <- subset(SomaLogic160410[vars], UniProt %in% u)
dim(s)
us <- s[!duplicated(s),]
dim(us)
us

# --- query intersect proteins ---
intersect_list <- setdiff(intersect(t4$UniProt,inf1$uniprot),"P23560")

# SomaLogic intersect
SomaLogic_overlap <- subset(t4,UniProt \%in\% intersect_list)
r <- snpquery(SomaLogic_overlap[,6],catalogue='pQTL',p=1e-11)
SomaLogic_overlap <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
save(SomaLogic_overlap,file='INTERVAL_overlap.rda',compress='xz')

# Olink intersect
INF1_merge_rsid <- read.delim(paste(INF,"work","INF1.merge-rsid",sep="/"))
INF1_merge_rsid_uniprot <- merge(INF1_merge_rsid,inf1,by="prot")
Olink_overlap <- subset(INF1_merge_rsid_uniprot,uniprot \%in\% intersect_list)
r <- snpquery(with(Olink_overlap,MarkerName),catalogue='pQTL',p=1e-11)
Olink_overlap <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
save(Olink_overlap,file='Olink_overlap.rda',compress='xz')
}
}
\note{adapted from custom codings}
\keyword{utilities}
