\name{pedigree}
\alias{pedigree}
\title{
  Create pedigree structure 
}
\description{
Create pedigree structure in format needed for plotting function.
}
\usage{
pedigree(id, dadid, momid, sex, affected, status, relations)
}
\arguments{
\item{id}{
Identification variable for individual
}
\item{dadid}{
Identification variable for father
}
\item{momid}{
Identification variable for mother
}
\item{sex}{
Gender of individual noted in
`id'. Character("male","female","unknown",
"terminated")
or numeric (1="male", 2="female", 3="unknown", 4="terminated") allowed.
}
\item{affected}{
One variable, or a matrix, indicating affection status.  Assumed that
1="unaffected", 2="affected", NA or 0 = "unknown". 
}
\item{status}{
Status (0="censored", 1="dead")
}
\item{relations}{
A matrix with 3 columns (id1, id2, code) specifying special relationship
between pairs of individuals. Codes: 1=Monozygotic
twin,  2=Dizygotic twin, 3=Twin of unknown zygosity, 4=Spouse and no children in pedigree  
}
}
\value{
An object of class pedigree.
}
\seealso{
\code{\link{plot.pedigree}}
}
\examples{
\dontrun{
ptemp <- pedigree(id=d10$upn, dadid=d10$dadid,momid=d10$momid, 
                  sex=d10$sex, affected=d10$affect)
plot(ptemp)
}
}
\keyword{dplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
