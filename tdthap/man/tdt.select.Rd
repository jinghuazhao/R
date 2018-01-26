\name{tdt.select}
\alias{tdt.select}
\title{
Select informative transmissions of sub-haplotypes for the TDT test
}
\usage{
tdt.select(hap.data, markers=1:((ncol(hap.data) - 4)/4), complete=TRUE)
}
\arguments{
\item{hap.data}{
The input dataframe. This will usually have been created by
\code{\link{hap.transmit}}. 
}
\item{markers}{
An integer array indicating which loci make up the relevant part of the 
haplotype.
}
\item{complete}{
If TRUE, only "complete" haplotypes are used (ie no zero's will be included).
}}
\value{
A list of two arrays of class "factor". The first (trans) contains transmitted
haplotypes and the second (untrans) contains untransmitted haplotypes. 
Rownames identify the transmission in terms of pedigree id, offspring id, 
father's id, mother's id, and whether it is a paternal transmission ("f") or 
a maternal transmission ("m"). 
}
\description{
This function is just a data handling intermediary between
\code{\link{hap.transmit}},  
which computes haplotypes, and \code{\link{tdt.quad}} and
\code{\link{tdt.rr}} which do TDT tests.
}
\references{
Clayton, D. and Jones, H. (1999) Transmission/disequilibrium tests for extended marker 
haplotypes.  Am.J.Hum.Gen., 65:1161-1169.
}
\seealso{
\code{\link{hap.transmit}}, \code{\link{tdt.rr}}, \code{\link{tdt.quad}}
}
\examples{
\dontrun{
# Select the sub-haplotype made up from the first two markers and print 
# tables of frequencies of transmitted and untransmitted haplotypes


	hap.use <- tdt.select(haps, markers=1:2)
	table(hap.use$trans)
	table(hap.use$untrans)
}}
\keyword{htest}
% Converted by Sd2Rd version 0.2-a4.
