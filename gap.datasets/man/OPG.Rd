\name{OPG}
\docType{data}
\alias{OPG}
\alias{OPGtbl}
\alias{OPGall}
\alias{OPGrsid}
\title{An example data for forest plot using METAL output}
\description{
This example contains METAL outputs (OPGtbl) as with association statistics from contributing studies (OPGall).
It is appropriate to use chr:pos_A1_A2 (A1<=A2) (SNPID) rather than reference id (rsid) due to its 
variability -- therefore a SNPID-rsid mapping file (OPGrsid) is also provided.
}
\usage{data(OPG)}
\format{Three data frames}
\examples{
data(OPG)
head(OPGtbl)
head(OPGall)
head(OPGrsid)
}
\source{SCALLOP consortium}
\references{
The SCALLOP paper.
}
\keyword{datasets}
