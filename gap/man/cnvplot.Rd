% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnvplot.R
\name{cnvplot}
\alias{cnvplot}
\title{genomewide plot of CNVs}
\usage{
cnvplot(data)
}
\arguments{
\item{data}{Data to be used.}
}
\value{
The cis/trans classification.
}
\description{
The function generates a plot containing genomewide copy number variants (CNV) chr, start, end, freq(uencies).
}
\examples{
knitr::kable(cnv,caption="A CNV dataset")
cnvplot(cnv)
circos.cnvplot(cnv)
}
