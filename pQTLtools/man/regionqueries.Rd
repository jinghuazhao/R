\name{regionqueries}
\alias{regionqueries}
\title{phenoscanner regionqueries in batches}
\usage{regionqueries(regionlist,catalogue="pQTL",proxies="EUR",p=5e-8,r2=0.7,build=37)}
\description{
R/phenoscanner only allows for certain number of items supplied. This simple function return
a large number of calls in batches as well as generating SNPIDs.
}
\arguments{
\item{regionlist}{a list of SNPs}
\item{catalogue}{"None","eQTL","mQTL","methQTL","pQTL","GWAS"}
\item{proxies}{"None", "AFR","AMR","EAS","EUR","SAS"}
\item{p}{p value threshold}
\item{r2}{r2 for LD}
\item{build}{37, 38}
}
\details{
Batches are generated and queries are combined into one.
}

\value{
The returned value is a list containing regions and results:
}

\section{References}{
Sun BB, et al. (2018). Genomic atlas of the human plasma proteome. *Nature* 558: 73-79.
}
\seealso{
\code{\link[phenoscanner]{phenoscanner}}
}

\examples{
\dontrun{
# single region
regionqueries("chr17:26691290-26700110")
}
}
\note{adapted from custom codings}
\keyword{utilities}
