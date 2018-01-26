\name{CDKN}
\docType{data}
\alias{CDKN}
\alias{CDKNgenes}
\alias{CDKNmap}
\alias{CDKNlocus}
\title{An example data for regional association plot}
\description{
These data are adapted from the DGI study on CDKN2A/CDKN2B region.
}
\usage{data(CDKN)}
\format{
There are three data objects in the dataset: \code{CDKNgenes}, the gene list from the Chromosome 9
according to UCSC browser (\url{http://genome.ucsc.edu/}); \code{CDKNmap}, the genetic map as from 
the HapMap website (\url{http://www.hapmap.org/downloads/recombination/2006-10_rel21_phaseI+II/rates/});
\code{CDKNlocus}, the results from the association analysis of the locus based on DGI data.
}
\source{
  The data were obtained from the Harvard-MIT Broad Institute
  (see \url{http://www.broad.mit.edu/diabetes/})
}
\references{
 Diabetes Genetics Initiative of Broad Institute of Havard and MIT, Lund University and Novartis Institute
 for BioMedical Research. \emph{Whole-genome association analysis identifies novel loci for type 2 diabetes
 and triglyceride levels} Science 2007;316(5829):1331-6
}
\examples{
data(CDKN)
head(CDKNlocus)
}
\keyword{datasets}
