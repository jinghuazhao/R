\name{mhtplot.trunc}
\alias{mhtplot.trunc}
\title{Truncated Manhattan plot}
\usage{mhtplot.trunc(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
                   "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05),
                    genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE,
                    annotatePval = NULL, annotateTop = TRUE, cex.mtext=0.6, cex.text=0.8,
                    mtext.line=2, cex.y= 1, y.ax.space=5, y.brk1, y.brk2, ...)}
\description{
To generate truncated Manhattan plot, e.g., of genomewide significance (p values) and
a random variable that is uniformly distributed.
}
\arguments{
\item{data}{a data frame with three columns representing chromosome, position and p values}
\item{control}{A control function named mht.control() with the following arguments,
    \enumerate{
    \item x. data.frame.
    \item chr. Chromosome.
    \item bp. position.
    \item p. p value.
    \item snp. SNP.
    \item col. colours.
    \item chrlabs. chromosome labels.
    \item suggestivelnine. suggestive line.
    \item genomewideline. genomewide line.
    \item highlight. highlight.
    \item logp. log(P).
    \item annotatePval. Annotate p.
    \item annotateTop. Annotate top.
    \item cex.mtext. cex.mtext.
    \item cex.text. cex.text.
    \item mtext.line. mtext.line.
    \item cex.y. cex.y.
    \item y.ax.space. y.ax.space.
    \item y.brk1. y.brk1.
    \item y.brk2. y.brk2.
    \item ... others.
}

\value{
The plot is shown on or saved to the appropriate device.
}

\seealso{
\code{\link[gap]{mhtplot}}
}

\examples{
\dontrun{
}
# region to cut from truncation
brk.low = 120
brk.upp = 270

if (brk.upp > max(-log10(x$P)) ){
  stop("Upper breakpoint must be lower than maximum -log10 P-value")
}

source("mhtplot-trunc.R")
outFile <- "IL12B"
png(file=paste0(outFile, ".png"), width= 1500)

par(oma=c(0,0,0,0), mar=c(5,6.5,1,1))
x <- read.delim(gzfile("METAL/IL.12B-1.tbl.gz")[c("Chromosome","Position","P.value","MarkerName")],as.is=TRUE)
mhplot.trunc(x, chr="Chromosome", bp="Position", p="P.value", snp="MarkerName",
             suggestiveline =F, genomewideline = -log10(5e-10), logp = TRUE,
             cex.mtext=2, cex.text=0.7, 
             mtext.line=4, y.brk1= brk.low, y.brk2= brk.upp, cex.axis=2, cex.y=2, cex=2,
             highlight= toHighlight,
             y.ax.space=20,
             col = c("blue4", "skyblue")
             )
dev.off()
}
\author{James Peters}
\keyword{hplot}
