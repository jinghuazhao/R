\name{pQTLtools-internal}
\title{Internal data/functions for pQTLtools}
\alias{caprion}
\alias{inf1}
\alias{Olink_NGS}
\alias{hgTables}
\alias{CVD_II}
\alias{CVD_III}
\alias{Cardiometabolic}
\alias{Cell_Regulation}
\alias{Development}
\alias{Immune_Response}
\alias{Immuno_Oncology}
\alias{Inflammation}
\alias{Metabolism}
\alias{Neurology}
\alias{Oncology_II}
\alias{Organ_Damage}
\alias{SomaLogic160410}
\alias{swath_ms}
\usage{
head(inf1)
}
\description{
A list of protein panels: caprion, inf1, Olink_NGS, Olink_qPCR (CVD_II, ...), SomaLogic160410, swath_ms.

In cases of incomplete information, a dataset called hgTables based on UCSC is provided.
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
head(caprion)
head(inf1)
head(Olink_NGS)
head(CVD_II)
head(SomaLogic160410)
head(swath_ms)
}
}
\details{This seeds for prototype for extensions}
\keyword{internal}
