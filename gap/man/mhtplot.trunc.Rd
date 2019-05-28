\name{mhtplot.trunc}
\alias{mhtplot.trunc}
\title{Truncated Manhattan plot}
\usage{mhtplot.trunc(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
                   "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05),
                    genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE,
                    annotatePval = NULL, annotateTop = TRUE, cex.mtext=0.6, cex.text=0.8,
                    mtext.line=2, cex.y= 1, y.ax.space=5, y.brk1, y.brk2, ...)}
\description{
To generate truncated Manhattan plot, e.g., of genomewide significance (P values) or a
random variable that is uniformly distributed.
}
\arguments{
\item{x}{A data.frame}
\item{chr}{Chromosome}
\item{bp}{Position}
\item{p}{P value}
\item{snp}{SNP}
\item{col}{Colours}
\item{chrlabs}{Chromosome labels}
\item{suggestiveline}{Suggestive line}
\item{genomewideline}{Genomewide line}
\item{highlight}{Highlight}
\item{logp}{log(P)}
\item{annotatePval}{Annotate P}
\item{annotateTop}{Annotate top}
\item{cex.mtext}{cex.mtext}
\item{cex.text}{cex.text}
\item{mtext.line}{mtext.line}
\item{cex.y}{cex.y}
\item{y.ax.space}{y.ax.space}
\item{y.brk1}{y.brk1}
\item{y.brk2}{y.brk2}
\item{...}{other options}
}

\value{
The plot is shown on or saved to the appropriate device.
}

\seealso{
\code{\link[gap]{mhtplot}}
}

\examples{
\dontrun{
#
require(gap.datasets)
mhtplot.trunc(mhtdata,chr = "chr", bp = "pos", p = "p", snp = "rsn", y.brk1=10, y.brk2=12)
#
z <- gzfile("METAL/IL.12B-1.tbl.gz")
x <- read.delim(z, as.is=TRUE)[c("Chromosome","Position","P.value","MarkerName")]
png("IL.12B.png", res=300, units="in", width=9, height=6)
par(oma=c(0,0,0,0), mar=c(5,6.5,1,1))
mhtplot.trunc(x, chr="Chromosome", bp="Position", p="P.value", snp="MarkerName",
             suggestiveline=FALSE, genomewideline=-log10(5e-8), logp = TRUE,
             cex.mtext=2, cex.text=0.7, 
             mtext.line=4, y.brk1=120, y.brk2=270, cex.axis=2, cex.y=2, cex=2,
             y.ax.space=20,
             col = c("blue4", "skyblue")
)
dev.off()
}
}
\author{James Peters}
\keyword{hplot}
