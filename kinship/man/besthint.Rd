\name{besthint}
\alias{besthint}
\title{
  Create a hints matrix for a pedigree.
}
\description{
A pedigree structure can contain a \code{hints} matrix
which helps to reorder the pedigree (e.g. left-to-right order of
children within family) so as to plot with minimal distortion.
This routine trys out a large number of configurations, finding the
best by brute force.
}
\usage{
besthint(ped, wt=c(1000, 10, 1), tolerance=0)
}
\arguments{
\item{ped}{
a pedigree object
}
\item{wt}{
relative weights for three types of "distortion" in a plotted pedigree.
The final score for a pedigree is the weighted sum of these; the lowest
score is considered the best.
The three components are 1: the number of dotted lines, connecting two
instances of the same person; 2: the lengths of those dotted lines; and
3: the horizontal offsets between parent/child pairs.
}
\item{tolerance}{
the threshold for acceptance.  
If any of the orderings that are attempted have a score that is less than
or equal to this value, the routine ceases searching for a better one.
}
}
\value{
a hints matrix
}
\details{
Assume that a pedigree has k founding couples, i.e., husband-wife pairs for
which neither has a parent in the pedigree.
The routine tries all k!/2 possible left to right orderings of the founders
(in random order), uses the \code{autohint} function
to optimize the order of children within each family, and computes a
score.
The hints matrix for the first pedigree to match the tolerance
level is returned, or that for the best score found if none match the
tolerance.
}
\seealso{
\code{\link{pedigree}, \link{plot.pedigree}, \link{autohint}}
}
\examples{
\dontrun{
# Find a good plot, only trying to avoid dotted connectors
myped$hints <- besthint(myped, wt=c(1000,100,0))
}
}
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 1.21.
