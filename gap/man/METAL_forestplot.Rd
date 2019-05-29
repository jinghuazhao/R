\name{METAL_forestplot}
\alias{METAL_forestplot}
\title{forest plot as R/meta's forest for METAL outputs}
\usage{
METAL_forestplot(tbl,all,rsid,pdf="INF1.fp.pdf",package="meta",...)
}
\description{
This functions takes a meta-data from METAL (tbl) and data from contributing studies (all)
for forest plot. It also takes a SNPID-rsid mapping (rsid) as contributing studies often
involve discrepancies in rsid so it is appropriate to use SNPID, i.e., chr:pos_A1_A2 (A1<=A2).
}
\arguments{
  \item{tbl}{Meta-anslysis summary statistics}
  \item{all}{statistics from all contributing studies}
  \item{rsid}{SNPID-rsid mapping file}
  \item{pdf}{name of a PDF file}
  \item{package}{style of plot as in meta, rmeta or forestplot}
  \item{...}{options to use for the pdf device}
}

\value{
It will generate a forest plot specified by pdf for direction-adjusted effect sizes.
}

\section{References}{
Scharzer G. (2007). meta: An R package for meta-analysis. R News, 7:40-5, https://cran.r-project.org/doc/Rnews/Rnews_2007-3.pdf, 
https://CRAN.R-project.org/package=meta.

Willer CJ, Li Y, Abecasis GR. (2010). METAL: fast and efficient meta-analysis of genomewideassociation scans. Bioinformations. 26:2190-1,
https://github.com/statgen/METAL, https://genome.sph.umich.edu/wiki/METAL.
}
\seealso{
\code{\link[gap]{METAL_forestplot}}
}
\examples{
\dontrun{
require(gap.datasets)
data(OPG)
METAL_forestplot(OPGtbl,OPGall,OPGrsid,width=8.75,height=5)
}
}
\author{Jing Hua Zhao}
\keyword{hplot}
\keyword{distribution}
