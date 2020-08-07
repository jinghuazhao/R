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
\alias{st4}
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

# SomaLogic variations
SomaLogic_INF1 <- subset(SomaLogic160410,UniProt \%in\% INF1_uniprot)
SomaLogic_uniprot <- setdiff(unique(with(SomaLogic_INF1,UniProt)),"P23560")
setdiff(INF1_merge_uniprot,SomaLogic_uniprot)

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
r <- regionqueries(st4list)}
}
\details{This seeds for prototype for extensions}
\keyword{internal}
