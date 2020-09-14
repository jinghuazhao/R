\name{pqtlMR}
\alias{pqtlMR}
\title{pQTL MR analysis}
\usage{pqtlMR(Ins,Ids,prefix="INF1")}
\description{
Given a list of pQTLs as instrumental variables (IVs), this function
performs MR analysis according to epigraphdb.
}
\arguments{
\item{Ins}{a data.frame containing a list of pQTLs and their effect sizes}
\item{Ids}{a vector containing a list of MR-Base IDs}
\item{prefix}{a prefix for output files}
}
\details{
This function is based on TwoSampleMR.
}

\value{
No value is returned but several files.
}

\section{References}{
Zheng, J. et al.(2020). Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases. *Nature Genetics*.
}

\examples{
\dontrun{
prefix <- Sys.getenv("prefix")
ivs <- within(read.table(paste0(prefix,".ins"),as.is=TRUE,header=TRUE),{pval=10^pval})
ids <- scan("ieu-id.txt",what="")
pQTLtools::pqtlMR(ivs,ids)
pqtlMR("INF1.ins","ieu-id.txt")
}
}
\note{adapted from script by Jie Zheng}
\keyword{utilities}
