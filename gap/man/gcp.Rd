\name{gcp}
\alias{gcp}
\alias{pm}
\alias{pmplus}
\title{Permutation tests using GENECOUNTING}
\usage{gcp(y, cc, g, handle.miss=1, miss.val=0, n.sim=0,
       locus.label=NULL, quietly=FALSE)}
\description{
This function is a R port of the GENECOUNTING/PERMUTE program which
generates EHPLUS-type statistics including z-tests for individual haplotypes
}

\arguments{
 \item{y}{A column of 0/1 indicating cases and controls}
 \item{cc}{analysis indicator, 0 = marker-marker, 1 = case-control}
 \item{g}{the multilocus genotype data}
 \item{handle.miss}{a flag with value 1 indicating missing data are allowed}
 \item{miss.val}{missing value}
 \item{n.sim}{the number of permutations}
 \item{locus.label}{label of each locus}
 \item{quietly}{a flag if TRUE will suppress the screen output}
}

\value{
The returned value is a list containing (p.sim and ph when n.sim > 0):

 \item{x2obs}{the observed chi-squared statistic}
 \item{pobs}{the associated p value}
 \item{zobs}{the observed z value for individual haplotypes}
 \item{p.sim}{simulated p value for the global chi-squared statistic}
 \item{ph}{simulated p values for individual haplotypes}

}

\section{References}{
Zhao JH, Curtis D, Sham PC (2000). Model-free analysis and permutation 
tests for allelic associations. Human Heredity 50(2): 133-139

Zhao JH (2004). 2LD, GENECOUNTING and HAP: Computer programs for linkage 
disequilibrium analysis. Bioinformatics 20: 1325-1326

Zhao JH, Qian WD Association analysis of unrelated individuals using
polymorphic genetic markers -- methods, implementation and application, Royal
Statistical Society 2003, Hassallt-Diepenbeek, Belgium.
}

\seealso{
\code{\link[gap]{genecounting}}
}

\examples{
\dontrun{

data(fsnps)
y<-fsnps$y
cc<-1
g<-fsnps[,3:10]

gcp(y,cc,g,miss.val="Z",n.sim=5)
hap.score(y,g,method="hap",miss.val="Z")
}
}
\author{Jing Hua Zhao}
\note{Built on gcp.c}
\keyword{models}
\keyword{htest}
