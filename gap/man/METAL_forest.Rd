\name{METAL_forest}
\alias{METAL_forest}
\title{forest plot as R/meta's forest for METAL outputs}
\usage{
METAL_forest(tbl,all,rsid,pdf)
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
}

\value{
It will generate a forest plot specified by pdf for direction-adjusted effect sizes.
}

\section{References}{
METAL paper and software.

R/meta package
}
\seealso{
\code{\link[gap]{METAL_forestplot}}
}
\examples{
\dontrun{
# examples will be added here.
}
}
\author{Jing Hua Zhao}
\keyword{hplot}
\keyword{distribution}
