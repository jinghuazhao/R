\name{sentinels}
\alias{sentinels}
\title{Sentinel identification from GWAS summary statistics}
\usage{sentinels(p,st,debug=FALSE,flanking=1e+6)}
\description{
This function accepts an object containing GWAS summary statistics for
signal identification as defined by flanking regions.
}
\arguments{
\item{p}{an object containing GWAS summary statistics}
\item{st}{row number as in p}
\item{debug}{a flag to show the actual data}
\item{flanking}{the width of flanking region}
}
\details{
Typically, input to the function are variants reaching certain level
of significance and the functtion identifies minimum p value at the
flanking interval; in the case of another variant in the flanking
window has smaller p value then it will be replaced.
}

\value{
The function give screen output.
}

\examples{
\dontrun{
# prot is the phenotype
# tag features data specifics, e.g., HLA
prot <- Sys.getenv("prot")
tag <- Sys.getenv("tag")
p <- read.delim(paste0("work/",prot,tag,".p"),as.is=TRUE)
chrs <- with(p,unique(Chrom))
for(chr in chrs)
{
  ps <- subset(p,Chrom==chr)
  row.names(ps) <- 1:nrow(ps)
  pp(ps, 1)
}
}
}
\keyword{utilities}
