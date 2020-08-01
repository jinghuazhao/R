\name{pQTLtools-internal}
\title{Internal data/functions for pQTLtools}
\alias{biomaRt}
\alias{caprion}
\alias{hg19}
\alias{hg19Tables}
\alias{inf1}
\alias{Olink_NGS}
\alias{Olink_qPCR}
\alias{SomaLogic160410}
\alias{swath_ms}
\description{
A list of protein panels: caprion, hg10Tables, Olink_NGS, Olink_qPCR (inf1 is heavily edited), SomaLogic160410, swath_ms.

In cases of incomplete information, datasets called hg19Tables based on UCSC and biomaRt from Ensembl are provided.
}
\arguments{
  \item{chr}{chromosome}
  \item{start}{Start position}
  \item{end}{End position}
  \item{gene}{Gene name}
  \item{UniProt}{UniProt ID}
}
\examples{
\dontrun{
head(biomaRt)
head(caprion)
head(hg19)
head(hg19Tables)
head(inf1)
head(Olink_NGS)
subset(Olink_qPCR,Panel=="Inflammation")
head(SomaLogic160410)
head(swath_ms)
}
}
\details{This seeds for prototype for extensions}
\keyword{internal}
