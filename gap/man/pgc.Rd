\name{pgc}
\alias{pgc}
\title{Preparing weight for GENECOUNTING}
\usage{pgc(data,handle.miss=1,is.genotype=0,with.id=0)}
\description{
This function is a R port of the GENECOUNTING/PREPARE program which takes
an array of genotyep data and collapses individuals with the same multilocus
genotype. This function can also be used to prepare for the genotype table in testing
Hardy-Weinberg equilibrium.
}

\arguments{
  \item{data}{the multilocus genotype data for a set of individuals}
  \item{handle.miss}{a flag to indicate if missing data is kept, 0 = no, 1 = yes}
  \item{is.genotype}{a flag to indicate if the data is already in the form of genotype identifiers}
  \item{with.id}{a flag to indicate if the unique multilocus genotype identifier is generated}
}

\value{
The returned value is a list containing:

\item{cdata}{the collapsed genotype data}
\item{wt}{the frequency weight}
\item{obscom}{the observed number of combinations or genotypes}
\item{idsave}{optional, available only if with.id = 1}
}

\section{References}{
Zhao JH, Sham PC (2003). Generic number system and haplotype analysis. Comp Prog Meth Biomed 70:1-9
}

\seealso{
\code{\link[gap]{genecounting}},\code{\link[gap]{hwe.hardy}}
}

\examples{
\dontrun{

data(hla)
x <- hla[,3:8]

# do not handle missing data
y<-pgc(x,handle.miss=0,with.id=1)
hla.gc<-genecounting(y$cdata,y$wt,handle.miss=0)

# handle missing but with multilocus genotype identifier
pgc(x,handle.miss=1,with.id=1)

# handle missing data with no identifier
pgc(x,handle.miss=1,with.id=0)
}
}

\author{Jing Hua Zhao}
\note{Built on pgc.c}
\keyword{utilities}
