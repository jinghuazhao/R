\name{align.pedigree}
\alias{align.pedigree}
\alias{alignped1}
\alias{alignped2}
\alias{alignped3}
\alias{alignped4}
\title{
  Generate plotting information for a pedigree
}
\description{
Given a pedigree, this function creates helper matrices that descibe
the layout of a plot of the pedigree.
}
\usage{
align.pedigree(ped, packed=T, hints=ped$hints, width=6, align=T)
}
\arguments{
\item{ped}{
a pedigree object
}
\item{packed}{
should the pedigree be compressed, i.e., to allow diagonal lines connecting
parents to children in order to have a smaller overall width for the plot.
}
\item{hints}{
two column hints matrix.  
The first column determines the relative order of subjects within a sibship,
as well as the relative order of processing for the founder couples.
(For this latter, the female founders are ordered as though they were sisters).
The second column contains spouse information, e.g., if hints[2,6] = 17,
then subject number 17 of the pedigree is a spouse of number 2, and is
preferentially plotted to the right of number 2.  Negative numbers plot
the spouse preferentially to the left.
}
\item{width}{
for a packed output, the minimum width
}
\item{align}{
should iterations of the `springs' algorithm be used to improve the plotted
output.  If True, a default number of iterations is used.  If numeric, this
specifies the number of iterations.
}
}
\value{
a structure with components
\item{n}{
a vector giving the number of subjects on each horizonal level of the
plot
}
\item{nid}{
a matrix with one row for each level, giving the numeric id of each subject
plotted.  (An value of 17 means the 17th subject in the pedigree).
}
\item{pos}{
a matrix giving the horizontal position of each plot point
}
\item{fam}{
a matrix giving the family id of each plot point.  A value of "3" would mean
that the two subjects in positions 3 and 4, in the row above, are this
subject's parents.
}
\item{spouse}{
a matrix with values 1= subject plotted to the immediate right is a spouse,
2= subject plotted to the immediate right is an inbred spouse, 0 = not a 
spouse
}
\item{twins}{
optional matrix which will only be present if the pedigree contains twins.
It has values 1= sibling to the right is a monozygotic twin, 2= sibling to
the right is a dizygotic twin, 3= sibling to the right is a twin of unknown
zygosity, 0 = not a twin
}
}
\details{
This is an internal routine, used almost exclusively by
\code{plot.pedigree}.
The subservient functions \code{alignped1}, 
\code{alignped2}, 
\code{alignped3}, 
and \code{alignped4} contain the bulk of the computation.
}
\seealso{
\code{\link{plot.pedigree}}
}
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 1.21.
