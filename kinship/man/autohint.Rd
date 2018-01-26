\name{autohint}
\alias{autohint}
\title{
  Align a pedigree to print well
}
\description{
A pedigree structure can contain a \code{hints} matrix
which helps to reorder the pedigree (e.g. left-to-right order of
children within family) so as to plot with minimal distortion.
This routine is called by the \code{pedigree} function
to create an intial hints matrix.
}
\usage{
autohint(ped)
}
\arguments{
\item{ped}{
a pedigree structure
}
}
\value{
a two column hints matrix
}
\details{
This routine would not normally be called by a user. 
It moves children within families, so that marriages are on the "edge"
of a set children, closest to the spouse.  
For pedigrees that have only a single connection between two families 
this simple-minded approach works surprisingly well.
For more complex structures either hand-tuning of the hints matrix, or
use of the \code{besthint} routine will usually
be required.
}
\seealso{
\code{\link{pedigree}, \link{besthint}}
}
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 1.21.
