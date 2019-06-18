\name{sentinels}
\alias{sentinels}
\title{Sentinel identification from GWAS summary statistics}
\usage{sentinels(p,pid,st,debug=FALSE,flanking=1e+6,
                 chr="Chrom",pos="End",b="Effect",se="StdErr",snp="MarkerName")}
\description{
This function accepts an object containing GWAS summary statistics for
signal identification as defined by flanking regions. As the associate P value
could be extremely small, the effect size and its standard error are used.
}
\arguments{
\item{p}{an object containing GWAS summary statistics}
\item{pid}{a phenotype (e.g., protein) name in pGWAS}
\item{st}{row number as in p}
\item{debug}{a flag to show the actual data}
\item{flanking}{the width of flanking region}
\item{chr}{Chromosome name}
\item{pos}{Position}
\item{b}{Effect size}
\item{se}{Standard error}
\item{snp}{Marker name}
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
p <- reshape::rename(OPGtbl, c(Chromosome="Chrom", Position="End"))
chrs <- with(p, unique(Chrom))
for(chr in chrs)
{
  ps <- subset(p[c("Chrom","End","MarkerName","Effect","StdErr")], Chrom==chr)
  row.names(ps) <- 1:nrow(ps)
  sentinels(ps, "OPG", 1)
}
subset(OPGrsid,MarkerName=="chr8:120081031_C_T")
subset(OPGrsid,MarkerName=="chr17:26694861_A_G")
## to obtain variance explained
tbl <- within(OPGtbl, chi2n <- (Effect/StdErr)^2/N)
s <- with(tbl, aggregate(chi2n,list(prot),sum))
names(s) <- c("prot", "h2")
sd <- with(tbl, aggregate(chi2n,list(prot),sd))
names(sd) <- c("p1","sd")
m <- with(tbl, aggregate(chi2n,list(prot),length))
names(m) <- c("p2","m")
h2 <- cbind(s,sd,m)
sink("h2.dat")
print(h2[c("prot","h2","sd","m")], row.names=FALSE)
sink()
png("h2.png", res=300, units="in", width=12, height=8)
np <- nrow(h2)
ord <- with(h2, order(h2))
with(h2[ord,], {
    plot(h2, cex=0.4, xaxt="n", xlab="protein", ylab=expression(h^2))
    xtick <- seq(0, np, by=1)
    axis(side=1, at=xtick, labels = FALSE)
    text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=0.5)
})
dev.off()
write.csv(tbl,file="INF1.csv",quote=FALSE,row.names=FALSE)
}
}
\keyword{utilities}
