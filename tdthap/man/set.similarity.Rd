\name{set.similarity}
\alias{set.similarity}
\title{
Set similarity function for Geary-Moran type TDT tests
}
\usage{
set.similarity(nloci=1, spacing=rep(1, nloci + 1), focus=1, power=1)
}
\arguments{
\item{nloci}{
The number of loci.
}
\item{spacing}{
A numeric array of length (nloci+1) giving marker spacings and "off-end" 
distances.
}
\item{focus}{
An integer in the range 1:nloci indicating the "focus" for the similarity 
function.
}
\item{power}{
The power to which the shared haplotype length is raised.
}}
\value{
A list of the values loaded.
}
\section{Side Effects}{
Sets constants accessed by tdt.quad() when calculating Geary-Moran type 
statistics.
}
\description{
Haplotypes are only similar if they are IBS at the focal locus. The extent of 
the similar region to each side is determined by stepping outwards until 
the haplotypes are no longer IBS, the region being assumed to end midway 
between the last IBS locus and the first non-IBS locus. If the haplotypes are 
IBS at the last locus, half the "off-end" distance is scored. The similarity
is defined as the total length of this shared region raised to some power. 
}
\references{
Clayton, D. and Jones, H. (1999) Transmission/disequilibrium tests for extended marker 
haplotypes.  Am.J.Hum.Gen., 65:1161-1169.
}
\seealso{
\code{\link{tdt.quad}}, \code{\link{get.similarity}}
}
\examples{
\dontrun{
# To do a Geary_Moran test on a 10 marker haplotype
	gaps <- c(0, 50, 60, 80, 20, 30, 50, 40, 50, 100, 0)
	set.similarity(nloci=10, spacing=gaps, power=0.5)
	test <- tdt.quad(hap.use, funct=T)
}}
\keyword{htest}
% Converted by Sd2Rd version 0.2-a4.
