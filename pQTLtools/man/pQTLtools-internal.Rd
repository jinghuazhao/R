\name{pQTLtools-internal}
\title{Internal data/functions for pQTLtools}
\alias{biomaRt}
\alias{caprion}
\alias{hg19}
\alias{hg19Tables}
\alias{inf1}
\alias{Olink_NGS}
\alias{Olink_qPCR}
\alias{SomaLogic160410}
\alias{swath_ms}
\description{
A list of protein panels: caprion, hg10Tables, Olink_NGS, Olink_qPCR (inf1 is heavily edited), SomaLogic160410, swath_ms.

In cases of incomplete information, datasets called hg19Tables based on UCSC and biomaRt from Ensembl are provided.
}
\arguments{
  \item{chr}{chromosome}
  \item{start}{Start position}
  \item{end}{End position}
  \item{gene}{Gene name}
  \item{UniProt}{UniProt ID}
}
\examples{
\dontrun{
# datasets
head(biomaRt)
head(caprion)
head(hg19)
head(hg19Tables)
head(inf1)
head(Olink_NGS)
subset(Olink_qPCR,Panel=="Inflammation")
head(SomaLogic160410)
head(swath_ms)

# INTERVAL SomaLogic results
dir <- 'https://static-content.springer.com/esm/art\%3A10.1038\%2Fs41586-018-0175-2/MediaObjects/'
file <- '41586_2018_175_MOESM4_ESM.xlsx'
xlsx <- paste0(dir,file)
st4.1 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE,
                             cols=c(1:16,26:28), rows=c(5:1986))
st4.2 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE,
                             cols=c(17:25,29:31), rows=c(6:1986))
st4 <- cbind(st4.1,st4.2)
r <- snpqueries(st4[,6],catalogue='pQTL',p=1e-11)
INTERVAL_SomaLogic <- subset(with(r,results),efo=='EFO_0004747' & pmid=='29875488')
save(INTERVAL_SomaLogic,file='INTERVAL_SomaLogic.rda',compress='xz')
st5 <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE,
                           cols=c(1:19), rows=c(3:2746))
st6 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE,
                           cols=c(1:20), rows=c(3:167))

intersect_list <- setdiff(intersect(st4$UniProt,inf1$uniprot),"P23560")

# Olink-SomaLogic overlap
p <- list(setdiff(inf1$uniprot,"P23560"),
          setdiff(SomaLogic160410$UniProt[!is.na(SomaLogic160410$UniProt)],"P23560"))
cnames <- c("INF1","SomaLogic")
VennDiagram::venn.diagram(x = p, category.names=cnames,
                          filename='os.png', imagetype="png", output=TRUE)
m <- merge(inf1,SomaLogic160410,by.x="uniprot",by.y="UniProt")
u <- setdiff(with(m,unique(uniprot)),"P23560")
options(width=220)
o <- subset(inf1,uniprot \%in\% u)
dim(o)
o
vars <- c("UniProt","chr","start","end","extGene","Target","TargetFullName")
s <- subset(SomaLogic160410[vars], UniProt \%in\% u)
dim(s)
us <- s[!duplicated(s),]
dim(us)
us
}
}
\details{This seeds for prototype for extensions}
\keyword{internal}
