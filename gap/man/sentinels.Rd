\name{sentinels}
\alias{sentinels}
\title{Sentinel identification from GWAS summary statistics}
\usage{sentinels(p,st,debug=FALSE,flanking=1e+6,prot="OPG")}
\description{
This function accepts an object containing GWAS summary statistics for
signal identification as defined by flanking regions.
}
\arguments{
\item{p}{an object containing GWAS summary statistics}
\item{pid}{a phenotype (e.g., protein) name in pGWAS}
\item{st}{row number as in p}
\item{debug}{a flag to show the actual data}
\item{flanking}{the width of flanking region}
}
\details{
Typically, input to the function are variants reaching certain level
of significance and the functtion identifies minimum p value at the
flanking interval; in the case of another variant in the flanking
window has smaller p value it will be used instead.

For now key variables in p are "MarkerName", "End", "Effect", "StdErr",
"P.value", where "End" is as in a bed file indicating marker position,
 and the function is set up such that row names are numbered as
1:nrow(p); see example below.
}

\value{
The function give screen output.
}

\examples{
\dontrun{
## OPG as a positive control in our pGWAS
require(gap.datasets)
data(OPG)
p <- reshape::rename(OPGtbl, c(Chromosome="Chrom"))
chrs <- with(p, unique(Chrom))
for(chr in chrs)
{
  ps <- subset(p[[c("Chrom","End","MarkerName","Effect","StdErr","P.value")]], Chrom==chr)
  row.names(ps) <- 1:nrow(ps)
  sentinels(ps, "OPG", 1)
}
subset(OPGrsid,MarkerName=="chr8:120081031_C_T")
subset(OPGrsid,MarkerName=="chr17:26694861_A_G")
}
}
\keyword{utilities}
