\name{hap.transmit}
\alias{hap.transmit}
\title{
Build parental haplotypes in nuclear families
}
\usage{
hap.transmit(pedfile, markers=1:((ncol(pedfile) - 6)/2),
             multiple.cases=0, use.affected=TRUE)
}
\arguments{
\item{pedfile}{
The input dataframe. The first six columns contain the pedigree id, the member
id, the two parental id's, the sex, and the affectation status. Subsequent 
fields are in pairs and represent alleles at marker loci. All variables must 
take integer values, with zero being taken as "missing".
}
\item{markers}{
Integer array indicating markers to be used and their order. 
}
\item{multiple.cases}{
The action to be taken if multiple affected offspring in any pedigree. Options
are (0) include all, (1) include all, but whole family is duplicated and only
one offspring is treated as affected in each repeated family, and (2) use only
first affected offspring.
}
\item{use.affected}{
If TRUE, data from affected offspring is used when imputing any missing 
parental data. Otherwise it is ignored.
}}
\value{
A dataframe with one row for each affected offspring. The first four columns 
identify the offspring by pedigree id, member id, and parental id's. The next 
block of columns hold the transmitted paternal haplotype. Following blocks 
contain the untransmitted paternal haplotype and maternal transmitted and 
untransmitted haplotypes. 
}
\description{
This version only computes parental haplotypes in so far as they can be 
derived with complete certainty. Any locus with is uncertain in the final 
haplotype is coded as zero.  
}
\references{
Clayton, D. and Jones, H. (1999) Transmission/disequilibrium tests for extended marker
haplotypes. Am. J. Hum. Gen., 65:1161-1169.
}
\seealso{
\code{\link{tdt.select}}, \code{\link{tdt.rr}}, \code{\link{tdt.quad}}
}
\examples{
\dontrun{
# Read a pedfile (which includes the variable names in the top line) 
# and build haplotypes using the markers which appear third, second, and 
# first in the pedfile. 

  filespec <- system.file("tests/test.ped", package="tdthap")
  ped <- read.table(filespec)

  haps <- hap.transmit(ped, markers=c(3,2,1))

}}
\keyword{htest}
% Converted by Sd2Rd version 0.2-a4.
