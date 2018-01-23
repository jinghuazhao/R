\name{mia}
\alias{mia}
\title{multiple imputation analysis for hap}
\usage{mia(hapfile,assfile,miafile,so,ns,mi,allsnps,sas)}
\description{
This command reads outputs from hap session that uses multiple imputations, i.e. -mi\# option. To
simplify matters it assumes -ss option is specified together with -mi option there.

This is a very naive version of MIANALYZE, but can produce results for PROC MIANALYZE of SAS
}
\arguments{
\item{hapfile}{hap haplotype output file name}
\item{assfile}{hap assignment output file name}
\item{miafile}{mia output file name}
\item{so}{to generate results according to subject order}
\item{ns}{do not sort in subject order}
\item{mi}{number of multiple imputations used in hap}
\item{allsnps}{all loci are SNPs}
\item{sas}{produce SAS data step program}
}
\details{
It simply extracts outputs from hap

}

\value{
The returned value is a list containing:

}

\section{References}{

Zhao JH and W Qian (2003) Association analysis of unrelated individuals
using polymorphic genetic markers. RSS 2003, Hassalt, Belgium

Clayton DG (2001) SNPHAP. http://www-gene.cimr.cam.ac.uk/clayton/software.


}
\seealso{
\code{\link[gap]{hap}}
}

\examples{
\dontrun{
# 4 SNP example, to generate hap.out and assign.out alone
data(fsnps)
hap(id=fsnps[,1],data=fsnps[,3:10],nloci=4)

# to generate results of imputations
control <- hap.control(ss=1,mi=5)
hap(id=fsnps[,1],data=fsnps[,3:10],nloci=4,control=control)

# to extract information from the second run above
mia(so=1,ns=1,mi=5)
file.show("mia.out")

## commands to check out where the output files are as follows:
## Windows
# system("command.com")
## Unix
# system("csh")
}
}
\note{adapted from hap, in fact cline.c and cline.h are not used}
\keyword{utilities}
